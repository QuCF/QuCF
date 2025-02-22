#include "../third_party/helper_cuda.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <filesystem>
#include <math.h>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <time.h>
#include <limits>
#include <chrono>
#include <stdarg.h>
#include <memory>
#include <numeric>
#include <unistd.h>
#include <cuda_runtime_api.h>
#include "H5Cpp.h"

using namespace std;

struct FUNC_
{
    string  sel;
    int32_t parity; // 0 - even, 1 - odd;
    double  coef_norm;
};

// used for the computation of angles
const vector<FUNC_> avail_functions_ = { 
   //   sel,       parity, coef_norm
   {"inversion",        1,     0.125}, // id = 0
   {"gaussian",         0,     0.980}, // id = 1
   {"gaussian-arcsin",  0,     0.980}, // id = 2
   {"arcsin",           1,     0.980}, // id = 3
   {"xgaussian",        1,     0.980}, // id = 4
   {"sqrt_inv_arcsin2", 0,     0.980}, // id = 5
   {"vH",               1,     1.000}, // id = 6
   {"ase-init",         0,     1.000}, // id = 7
};

struct FUNC_DATA_
{
    int32_t id;      // ID of the function: specifies the function type: inversion, gaussian, exponential etc.
    int32_t parity;
    double param;
    double coef_norm;
};
__constant__ FUNC_DATA_ function_d_;

void calculate_coefficients(
    const uint32_t& Nd, const uint32_t& N_coefs_avail, const FUNC_DATA_& function_h,
    double*& coefs_real, double*& coefs_imag, uint32_t& N_coefs
);
void construct_polynomial(
    const double* coefs_real, const uint32_t& N_coefs,
    const FUNC_DATA_& function_h, const uint32_t& Nx_half, 
    double*& x, double*& pol, double*& orig_func, double& err_res
);
void save_coefs(
    const FUNC_& pars_func,
    const double& param, 
    const double& err_res, 
    const uint32_t& N_coefs, 
    double* coefs_real, 
    double* coefs_imag,
    const uint32_t& Nx_half, 
    const double* x, const double* pol, const double* orig_func,
    const std::string& work_path,
    const std::string& filename_out
);

__global__ 
void calc_coefs_odd(
    uint32_t Nd, uint32_t N_coefs_device, double *coefs_real, double *coefs_imag
);
__global__ 
void calc_coefs_even(
    uint32_t Nd, uint32_t N_coefs_device, double *coefs_real, double *coefs_imag
);

bool compare_strings(const string& line1, const string& line2);
void get_current_date_time(string& line_date_time);



void remove_hdf5_ext(std::string& str) 
{
    const std::string extension = ".hdf5";
    if (
        str.size() >= extension.size() && 
        str.rfind(extension) == (str.size() - extension.size())
    )
        str = str.substr(0, str.size() - extension.size());
}


/**
 * The resulting .hdf5 file is stored in the launch folder.
 * Assume that 
 *      $qucf_pol = [QuCF]/build_angles/approx_polyn
 * To launch:
 * > cd [launch-folder]
 * > $qucf_pol -sel_function [sel-function] 
 *             -param [parameter-value] 
 *             -Nd [number-of-coefficients-in-Fourier] 
 *             -work_path [where-output-will-be-saved]  (optional)
*/
int main(int argc, char *argv[])
{
    int nDevices;
    uint32_t Nd;         // initial number of coefficients in the polynomial;
    double param;        // function main parameter: kappa, time etc.
    string sel_function; // ID of the function to approximate;
    FUNC_DATA_ function_h;
    uint32_t id_arg;
    string work_path = "./";
    string filename_out = "";

    cout << "--- Fourier approach ---" << endl;

    // --- INPUT parameters ---
    if(argc < 7)
    {
        cout << "Error: some input parameters are missing." << endl;
        return -1;
    }

    id_arg = 1;
    while(id_arg < (argc - 1))
    {
        if(compare_strings(argv[id_arg], "-sel_function"))
        {
            id_arg += 1;
            sel_function = string (argv[id_arg]);
        }
        if(compare_strings(argv[id_arg], "-param"))
        {
            id_arg += 1;
            param = stod(string (argv[id_arg]));
        }
        if(compare_strings(argv[id_arg], "-Nd"))
        {
            id_arg += 1;
            Nd = stoi(string (argv[id_arg]));
        }
        if(compare_strings(argv[id_arg], "-work_path"))
        {
            id_arg += 1;
            work_path = string (argv[id_arg]);
        }
        if(compare_strings(argv[id_arg], "-filename_out"))
        {
            id_arg += 1;
            filename_out = string (argv[id_arg]);
            remove_hdf5_ext(filename_out);
        }
        ++id_arg;
    }

    function_h.id = -1;
    for(int32_t ii = 0; ii < avail_functions_.size(); ii++)
        if(compare_strings(avail_functions_[ii].sel, sel_function))
            function_h.id = ii;
    if(function_h.id == -1)
    {
        cout << "Error: the function with ID = {" << sel_function << "} is not defined." << endl;
        return -1;
    }
    function_h.parity    = avail_functions_[function_h.id].parity;
    function_h.param     = param;
    function_h.coef_norm = avail_functions_[function_h.id].coef_norm;

    cudaGetDeviceCount(&nDevices);
    if(nDevices == 0)
    {
        cout << "Error: GPU devices that support CUDA are not found." << endl;
        return -1;
    }

    string current_path = filesystem::current_path();

    cout << "Function to approximate: " << avail_functions_[function_h.id].sel << "\n";
    cout << "Its parity: \t" << function_h.parity << "\n";
    cout << "Parameter: \t" << param << "\n";
    cout << "Rescaling parameter: \t" << function_h.coef_norm << "\n";
    cout << "Nd: \t" << Nd << "\n";
    cout << "N of avail. GPU devices: " << nDevices << endl;
    cout << "work directory: " << current_path << "/" << work_path << endl;

    if(!filename_out.empty())
        cout << "output filename: " << filename_out << endl;

    // assume that all devices have the same properties:
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    double coef_mem = 0.9;
    uint32_t N_coefs_avail = coef_mem * prop.totalGlobalMem / (2*sizeof(double));
    cout << "\nAvail. GPU mem for coefs. (MB): " << 
        coef_mem * prop.totalGlobalMem / (1024 * 1024.) << "\n";

    for(int id_device = 0; id_device < nDevices; id_device++)
    {
        cudaSetDevice(id_device);
        cudaMemcpyToSymbol(function_d_, &function_h, sizeof(FUNC_DATA_));
        cudaDeviceSynchronize();
        checkCudaErrors(cudaGetLastError());
    }

    // calculate the polynomial coefficients:
    double* coefs_real;
    double* coefs_imag;
    uint32_t N_coefs;
    calculate_coefficients(Nd, N_coefs_avail, function_h, coefs_real, coefs_imag, N_coefs);
    
    // construct the resulting polynomial and estimate the error:
    uint32_t Nx_half = 4001; 
    double *x, *pol, *orig_fun;
    double err_res;
    construct_polynomial(coefs_real, N_coefs, function_h, Nx_half, x, pol, orig_fun, err_res);
  
    // --- Save the coefficients to the .hdf5 file ---
    save_coefs(
        avail_functions_[function_h.id],
        param, err_res, N_coefs, 
        coefs_real, coefs_imag, 
        Nx_half, x, 
        pol, orig_fun,
        work_path,
        filename_out
    );
    return 0;
}


/**
 * If modified, change also F_CALC_HOST.
*/
__device__ __forceinline__
double F_CALC(const double& x) 
{
    auto& dd = function_d_;

    // inversion function:
    if(dd.id == 0)
    {
        double kappa = dd.param;
        return (dd.coef_norm/kappa) * (1. - exp(-pow(5*kappa*x, 2))) / x;
    }
    // Gaussian(x):
    if(dd.id == 1)
    {
        double mu = dd.param;
        return dd.coef_norm * exp(-x*x/(2*mu*mu));
    }
    // Gaussian(arcsin(x)):
    if(dd.id == 2)
    {
        double mu = dd.param;
        return dd.coef_norm * exp(-asin(x)*asin(x)/(2*mu*mu));
    }
    // arcsin(x):
    if(dd.id == 3)
    {
        double a = dd.param;
        return dd.coef_norm * asin(a*x);
    }
    // x*Gaussian(x):
    if(dd.id == 4)
    {
        double mu = dd.param;
        return dd.coef_norm * x * exp(-x*x/(2*mu*mu));
    }
    // sqrt(1 + par^2 * arcsin(x)^2)^{-1}:
    if(dd.id == 5)
    {
        double x_max = dd.param;
        double asin2 = asin(x) * asin(x);
        return dd.coef_norm / sqrt(1 + x_max*x_max * asin2);
    }
    // par * asin(x) * exp(pred_coef*asin(x)^2):
    if(dd.id == 6)
    {
        double ampl = dd.param;
        double yy  = asin(x);
        double yy2 = yy*yy;
        double pred_coef = -7.999999938811e+00;
        return dd.coef_norm * ampl * yy * exp(pred_coef * yy2);
    }
    // par * exp(ww*asin(x)^2):
    if(dd.id == 7)
    {
        // double ww = -1.249984330652e+03;
        double ww = -5.065995676665e+02;
        double yy  = asin(x);
        double asin2 = yy*yy;
        return dd.coef_norm * dd.param * exp(ww * asin2);
    }
    return 0; // function is missing;
}


/**
 * If modified, change also F_CALC.
*/
double F_CALC_HOST(const double& x, const FUNC_DATA_& dd)  
{
    // inversion function:
    if(dd.id == 0)
    {
        double kappa = dd.param;
        return (dd.coef_norm/kappa) * (1. - exp(-pow(5*kappa*x, 2))) / x;
    }
    // Gaussian(x):
    if(dd.id == 1)
    {
        double mu = dd.param;
        return dd.coef_norm * exp(-x*x/(2*mu*mu));
    }
    // Gaussian(arcsin(x)):
    if(dd.id == 2)
    {
        double mu = dd.param;
        return dd.coef_norm * exp(-asin(x)*asin(x)/(2*mu*mu));
    }
    // arcsin(x):
    if(dd.id == 3)
    {
        double a = dd.param;
        return dd.coef_norm * asin(a*x);
    }
    // x*Gaussian(x):
    if(dd.id == 4)
    {
        double mu = dd.param;
        return dd.coef_norm * x * exp(-x*x/(2*mu*mu));
    }
    // sqrt(1 + par^2 * arcsin(x)^2)^{-1}:
    if(dd.id == 5)
    {
        double x_max = dd.param;
        double asin2 = asin(x) * asin(x);
        return dd.coef_norm / sqrt(1 + x_max*x_max * asin2);
    }
    // par * asin(x) * exp(pred_coef*asin(x)^2):
    if(dd.id == 6)
    {
        double ampl = dd.param;
        double yy  = asin(x);
        double yy2 = yy*yy;
        double pred_coef = -7.999999938811e+00;
        return dd.coef_norm * ampl * yy * exp(pred_coef * yy2);
    }
    // par * exp(ww*asin(x)^2):
    if(dd.id == 7)
    {
        // double ww = -1.249984330652e+03;
        double ww = -5.065995676665e+02;
        double yy  = asin(x);
        double asin2 = yy*yy;
        return dd.coef_norm * dd.param * exp(ww * asin2);
    }
    return 0; // function is missing;
}


void calculate_coefficients(
    const uint32_t& Nd, const uint32_t& N_coefs_avail, const FUNC_DATA_& function_h,
    double*& coefs_real, double*& coefs_imag, uint32_t& N_coefs
){
    N_coefs = int(Nd/2.);

    // ATTENTION: the case where 
    //  (the number of polynomial coefficients) > (number of coefficients that can be stored on a single GPU) 
    // is not debugged:
    uint32_t N_coefs_device = (N_coefs > N_coefs_avail) ? N_coefs_avail: N_coefs;
    uint32_t N_iter         = int(N_coefs/N_coefs_avail) + 1;
    double* coefs_real_d;
    double* coefs_imag_d;

    // ATTENTION: N_coefs maybe too large, in this case it will be necessary 
    // to create new arrays coefs_real and coefs_imag at each iteration:
    coefs_real = new double[N_coefs]; 
    coefs_imag = new double[N_coefs];

    cout << "\n---\n";
    cout << "Required mem. for coefs. (MB): " << 
        2*sizeof(double) * N_coefs / (1024 * 1024.) << endl;

    cudaMalloc((void**) &(coefs_real_d), N_coefs_device * sizeof(double));
    cudaMalloc((void**) &(coefs_imag_d), N_coefs_device * sizeof(double));
    uint32_t N_coefs_device_init = N_coefs_device;
    for(auto counter_iter = 0; counter_iter < N_iter; counter_iter++)
    {
        if(counter_iter == N_iter - 1)
            N_coefs_device = N_coefs - N_coefs_device * counter_iter;

        cudaMemset(coefs_real_d, 0.0, N_coefs_device * sizeof(double));
        cudaMemset(coefs_imag_d, 0.0, N_coefs_device * sizeof(double));

        uint32_t N_threads_per_block = (N_coefs_device < 1024)?  N_coefs_device: 1024;
        uint32_t N_cuda_block = 
            (N_coefs_device <= N_threads_per_block) ? 
                1:
                int(N_coefs_device / N_threads_per_block) + 1;
        
        cout << "\n\t *** Iteration: " << counter_iter        << " ***\n";
        cout << "\tN-coefs per device: " << N_coefs_device << "\n";
        cout << "\tN-blocks: "  << N_cuda_block        << "\n";
        cout << "\tN-threads: " << N_threads_per_block << "\n";

        // CALCULATE the polynomial coefficients for the ODD function:
        if(function_h.parity == 1)
            calc_coefs_odd<<<N_cuda_block, N_threads_per_block>>>(
                Nd, N_coefs_device, coefs_real_d, coefs_imag_d
            );
        // CALCULATE the polynomial coefficients for the EVEN function:
        if(function_h.parity == 0)
            calc_coefs_even<<<N_cuda_block, N_threads_per_block>>>(
                Nd, N_coefs_device, coefs_real_d, coefs_imag_d
            );

        cudaMemcpy(
            coefs_real + N_coefs_device_init * counter_iter, 
            coefs_real_d,
            N_coefs_device * sizeof(double),
            cudaMemcpyDeviceToHost
        );
        cudaMemcpy(
            coefs_imag + N_coefs_device_init * counter_iter, 
            coefs_imag_d,
            N_coefs_device * sizeof(double),
            cudaMemcpyDeviceToHost
        );

        cudaDeviceSynchronize();
        checkCudaErrors(cudaGetLastError());
    }
}


void construct_polynomial(
    const double* coefs_real, const uint32_t& N_coefs,
    const FUNC_DATA_& function_h, const uint32_t& Nx_half, 
    double*& x, double*& pol, double*& orig_func, double& err_res
){
    // ATTENTION: here we assume that all coefficients can be saved on a single GPU.

    cout << "\nConstruction of the polynomial...\n";
    x         = new double[2*Nx_half];
    pol       = new double[2*Nx_half];
    orig_func = new double[2*Nx_half];

    for(auto ii = 0; ii < (2*Nx_half); ii++)
        x[ii] = cos((2*ii + 1)*M_PI / (4*Nx_half));

    if(
        function_h.id == 2 || function_h.id == 3 || function_h.id == 5 || 
        function_h.id == 6 || function_h.id == 7
    ) 
        for(auto ii = 0; ii < (2*Nx_half); ii++)
            x[ii] = sin(x[ii]);

    // // x-grid:
    // if(function_h.id == 0) // for the inversion function
    // {
    //     double inv_param = 1./function_h.param;
    //     double dx = (1. - inv_param) / (Nx_half - 1.);
    //     for(auto ii = 0; ii < Nx_half; ii++)
    //         x[ii] = -1 + dx * ii;
    //     for(auto ii = 0; ii < Nx_half; ii++)
    //         x[Nx_half + ii] = inv_param + dx * ii;
    // }
    // else 
    // {
    //     // double dx = 2. / (2*Nx_half-1);
    //     // for(auto ii = 0; ii < (2*Nx_half); ii++)
    //     //     x[ii] = -1 + dx * ii;

    //     for(auto ii = 0; ii < (2*Nx_half); ii++)
    //         x[ii] = cos((2*ii + 1)*M_PI / (2*Nx_half));

    //     // for functions depending on arcsin:
    //     if(function_h.id == 2 || function_h.id == 3 || function_h.id == 5) 
    //         for(auto ii = 0; ii < (2*Nx_half); ii++)
    //             x[ii] = sin(x[ii]);
    // }

    // Construct the polynomial (using only real coefficients) for the ODD function:
    if(function_h.parity == 1)
        for(auto id_x = 0; id_x < (2*Nx_half); id_x++)
        {
            double pol_one = 0;
            double x_one = x[id_x];
            for(auto id_coef = 1; id_coef < N_coefs+1; id_coef++)
                pol_one += coefs_real[id_coef - 1] * cos((2*id_coef-1)*acos(x_one));
            pol[id_x] = pol_one;
        }

    // Construct the polynomial (using only real coefficients) for the EVEN function:
    if(function_h.parity == 0)
        for(auto id_x = 0; id_x < (2*Nx_half); id_x++)
        {
            double pol_one = 0;
            double x_one = x[id_x];
            for(auto id_coef = 0; id_coef < N_coefs; id_coef++)
                pol_one += coefs_real[id_coef] * cos((2*id_coef)*acos(x_one));
            pol[id_x] = pol_one;
        }

    // Construct the original function:
    for(auto id_x = 0; id_x < (2*Nx_half); id_x++)
        orig_func[id_x] = F_CALC_HOST(x[id_x], function_h);

    // Estimate the absolute error of the polynomial:
    err_res = 0;
    for(auto id_x = 0; id_x < (2*Nx_half); id_x++)
    {
        double err1 = abs(pol[id_x] - orig_func[id_x]);
        if(err1 > err_res)
            err_res = err1;
    }
    cout << "Resulting approximation error: " << std::scientific << setprecision(3) << err_res << endl;
}


__global__ void calc_coefs_odd(
    uint32_t Nd, uint32_t N_coefs_device, double *coefs_real, double *coefs_imag 
){
    // printf("id-thread = %u, parity = %d\n", threadIdx.x, function_d_.parity);
    auto idx = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t N_quad = 2 * Nd; 

    if(idx < N_coefs_device)
    {
        uint32_t ii = 2 * idx + 1;
        double coef_in_front = pow(-1, ii) / N_quad;
        double sum_temp_real = 0;
        double sum_temp_imag = 0;
        double temp;
        double th;
        for(uint32_t kk = 0; kk < 2*N_quad; kk++)
        {
            th = M_PI * kk / N_quad;
            temp = F_CALC(-cos(th));
            sum_temp_real += temp * cos(ii * th);
            sum_temp_imag += temp * sin(ii * th);
        }
        coefs_real[idx] = coef_in_front * sum_temp_real;
        coefs_imag[idx] = coef_in_front * sum_temp_imag;
    }
}


__global__ void calc_coefs_even(
    uint32_t Nd, uint32_t N_coefs_device, double *coefs_real, double *coefs_imag 
){
    // printf("id-thread = %u, parity = %d\n", threadIdx.x, function_d_.parity);
    auto idx = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t N_quad = 2 * Nd; 

    if(idx < N_coefs_device)
    {
        uint32_t ii = 2 * idx;
        double coef_in_front = pow(-1, ii) / N_quad;
        if(ii == 0) coef_in_front = coef_in_front/2.;
        double sum_temp_real = 0;
        double sum_temp_imag = 0;
        double temp;
        double th;
        for(uint32_t kk = 0; kk < 2*N_quad; kk++)
        {
            th = M_PI * kk / N_quad;
            temp = F_CALC(-cos(th));
            sum_temp_real += temp * cos(ii * th);
            sum_temp_imag += temp * sin(ii * th);
        }
        coefs_real[idx] = coef_in_front * sum_temp_real;
        coefs_imag[idx] = coef_in_front * sum_temp_imag;
    }
}


void save_coefs(
    const FUNC_& pars_func,
    const double& param, 
    const double& err_res, 
    const uint32_t& N_coefs, 
    double* coefs_real, 
    double* coefs_imag,
    const uint32_t& Nx_half, 
    const double* x, const double* pol, const double* orig_func,
    const std::string& work_path,
    const std::string& filename_out
){
    std::stringstream sstr;

    // // --- filename based on eps ---
    // sstr << work_path << "/" << pars_func.sel << "_" << to_string(param) << "_" << round(-log10(err_res)) << ".hdf5";

    // --- filename based on N_coefs ---
    sstr << work_path << "/";
    if(filename_out.empty())
        sstr << pars_func.sel << "_param_" << 
            std::setprecision(3) << param << "_Nc_" << N_coefs << ".hdf5";
    else
        sstr << filename_out << ".hdf5";

    string filename_hdf5 = sstr.str(); 
    H5::H5File* f_ = new H5::H5File(filename_hdf5, H5F_ACC_TRUNC);
    H5::Group grp_basic(f_->createGroup("basic"));
    H5::Group grp_coefs(f_->createGroup("coefs"));
    H5::Group grp_functions(f_->createGroup("functions"));

    // description of the data:
    string descr = pars_func.sel;
    H5::StrType dtype_descr(H5::PredType::C_S1, descr.size()+1);
    H5::DataSet dataset_descr = grp_basic.createDataSet(
        "descr", 
        dtype_descr, 
        H5::DataSpace(H5S_SCALAR)
    );
    dataset_descr.write(descr, dtype_descr);

    // save the date of the simulation:
    string str_date_time;
    get_current_date_time(str_date_time);

    H5::StrType dtype_str_time(H5::PredType::C_S1, str_date_time.size()+1);
    H5::DataSet dataset_str_time = grp_basic.createDataSet(
        "date-of-simulation", 
        dtype_str_time, 
        H5::DataSpace(H5S_SCALAR)
    );
    dataset_str_time.write(str_date_time, dtype_str_time);


    // save the function parity:
    H5::DataSet dataset_parity = grp_basic.createDataSet(
        "parity", 
        H5::PredType::NATIVE_INT32, 
        H5::DataSpace(H5S_SCALAR)
    );
    dataset_parity.write((int*) &pars_func.parity, H5::PredType::NATIVE_INT32);

    // save the function parameter:
    H5::DataSet dataset_param = grp_basic.createDataSet(
        "param", 
        H5::PredType::NATIVE_DOUBLE, 
        H5::DataSpace(H5S_SCALAR)
    );
    dataset_param.write((int*) &param, H5::PredType::NATIVE_DOUBLE);

    // save the approximation error:
    H5::DataSet dataset_err = grp_basic.createDataSet(
        "eps", 
        H5::PredType::NATIVE_DOUBLE, 
        H5::DataSpace(H5S_SCALAR)
    );
    dataset_err.write((int*) &err_res, H5::PredType::NATIVE_DOUBLE);

    // save the rescaling factor of the function:
    H5::DataSet dataset_factor_norm = grp_basic.createDataSet(
        "coef_norm", 
        H5::PredType::NATIVE_DOUBLE, 
        H5::DataSpace(H5S_SCALAR)
    );
    dataset_factor_norm.write((int*) &(pars_func.coef_norm), H5::PredType::NATIVE_DOUBLE);

    // save the coefficients:
    hsize_t dims_coefs[] = {N_coefs};
    H5::DataSpace dspace_coefs(1, dims_coefs);
    H5::DataSet dataset_real = grp_coefs.createDataSet(
        "real", 
        H5::PredType::NATIVE_DOUBLE, 
        dspace_coefs
    );
    dataset_real.write(coefs_real, H5::PredType::NATIVE_DOUBLE);

    H5::DataSet dataset_imag = grp_coefs.createDataSet(
        "imag", 
        H5::PredType::NATIVE_DOUBLE, 
        dspace_coefs
    );
    dataset_imag.write(coefs_imag, H5::PredType::NATIVE_DOUBLE);

    // save the x-grid:
    hsize_t dims_x[] = {2*Nx_half};
    H5::DataSpace dspace_x(1, dims_x);
    H5::DataSet dataset_x = grp_functions.createDataSet(
        "x", 
        H5::PredType::NATIVE_DOUBLE, 
        dspace_x
    );
    dataset_x.write(x, H5::PredType::NATIVE_DOUBLE);

    // save the original function:
    hsize_t dims_orig[] = {2*Nx_half};
    H5::DataSpace dspace_orig(1, dims_orig);
    H5::DataSet dataset_orig = grp_functions.createDataSet(
        "orig", 
        H5::PredType::NATIVE_DOUBLE, 
        dspace_orig
    );
    dataset_orig.write(orig_func, H5::PredType::NATIVE_DOUBLE);

    // save the constructed polynomial:
    hsize_t dims_pol[] = {2*Nx_half};
    H5::DataSpace dspace_pol(1, dims_pol);
    H5::DataSet dataset_pol = grp_functions.createDataSet(
        "pol", 
        H5::PredType::NATIVE_DOUBLE, 
        dspace_pol
    );
    dataset_pol.write(pol, H5::PredType::NATIVE_DOUBLE);

    delete f_;

    cout << "Data are saved to the file: " << filename_hdf5 << endl;
}


bool compare_strings(const string& line1, const string& line2)
{
    string new_line1(line1), new_line2(line2);
    transform(new_line1.begin(), new_line1.end(), new_line1.begin(), ::tolower);
    transform(new_line2.begin(), new_line2.end(), new_line2.begin(), ::tolower);

    if(new_line1.compare(new_line2) == 0)
        return true;
    else
        return false;
}

void get_current_date_time(string& line_date_time)
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer, sizeof(buffer), "%m-%d-%Y %H:%M:%S", timeinfo);
    line_date_time = string(buffer);
}







