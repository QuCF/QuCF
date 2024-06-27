#ifndef COMPUTEANGLES_H
#define COMPUTEANGLES_H

#include "QLib.h"

struct DATA_ANGLES_
{
    std::string function_type;
    short parity; // 0 - even, 1 - odd, otherwise - without definite parity
    double coef_norm;
    double abs_error;
    double parameter;
    double stopping_criterion;

    // the maximum number of iterations in the L-BFGS algorithm:
    uint32_t max_n_iters;

    // the number of pairs (s,y) to remember to represent a Hessian in the L-BFGS:
    uint32_t N_pairs;

    // the coefficients in the polynomial approximation:
    std::shared_ptr<double[]> coefs;

    // The number of coefficients in the polynomial approximation.
    // Remark, in Dong-21, the number of angles is equal to (d+1), and Nc = \tilde(d).
    // Due to the assumed symmetry, the number of angles (d+1) is:
    // 2*Nc   if d is odd;
    // 2*Nc-1 if d is even;
    uint32_t Nc;
};


/**
 * Keeping \p N_pairs pairs (s,y) from previous L-BFGS iterations.
*/
struct LBGFS_pairs_
{
    // the number of pairs (s,y) from previous steps;
    uint32_t N_pairs; 

    // the number of the QSVT angles:
    uint32_t N_phis;

    // difference between previous and 
    // new estimation of the unknowns (i.e., QSVT angles phis);
    // s[id_pair, id_phi]
    std::shared_ptr<YMATH::VectorD_[]> s; 

    // gradients of the cost function from previous iterations;
    // y[id_pair. id_phi]
    std::shared_ptr<YMATH::VectorD_[]> y; 

    // The scalar coefficients 1./(s^T * y);
    // rho[id_pair]
    std::shared_ptr<double[]> rho;

    // the current number of saved pairs:
    uint32_t N_current;
    
    // the position of last saved pair:
    uint32_t id_last_saved;


    LBGFS_pairs_(YCU max_N_pairs, YCU in_N_phis)
    {
        using namespace std;

        N_pairs = max_N_pairs;
        N_phis  = in_N_phis;

        s = shared_ptr<YMATH::VectorD_[]>(new YMATH::VectorD_[N_pairs]);
        y = shared_ptr<YMATH::VectorD_[]>(new YMATH::VectorD_[N_pairs]);
        for(uint32_t ii = 0; ii < N_pairs; ii++)
        {
            s[ii].init(N_phis);
            y[ii].init(N_phis);
        }

        rho = shared_ptr<double[]>(new double[N_pairs]);

        N_current = 0;
        id_last_saved = -1;
    }


    /**
     * From the recently saved to the latest.
     *  [id_last_saved, ..., 0, N_current -1, N_current - 2,...]
    */
    inline uint32_t get_id_pair_first(YCU ii)
    {
        int temp = id_last_saved - ii;
        return temp >= 0 ? temp: N_current + temp;
    }
    /**
     * From the latest saved to the newest.
    */
    inline uint32_t get_id_pair_second(YCU ii)
    {  
        uint32_t temp = id_last_saved + ii + 1;
        return (temp < N_current) ? temp: temp - N_current; 
    }


    void first_loop(
        YMATH::VectorD_& d, 
        std::shared_ptr<double[]>& alpha,
        const YMATH::VectorD_& mean_grad_phi
    ){
        using namespace std;
        uint32_t id_pair;

        d = YMATH::VectorD_(mean_grad_phi);
        for(uint32_t ii = 0; ii < N_current; ii++)
        {
            id_pair   = get_id_pair_first(ii);
            alpha[ii] = rho[id_pair] * s[id_pair].dot(d);
            d.plus_mult_by(y[id_pair], -alpha[ii]);
        }
        return;
    }


    void second_loop(
        YMATH::VectorD_& d, 
        const std::shared_ptr<const double[]>& alpha
    ){
        using namespace std;
        uint32_t id_pair;
        double beta;
        for(uint32_t ii = 0; ii < N_current; ii++)
        {
            id_pair = get_id_pair_second(ii);
            beta = rho[id_pair] * y[id_pair].dot(d);
            d.plus_mult_by(s[id_pair], alpha[N_current - 1 - ii] - beta);
        }
    }
};



class ComputeAngles_{

protected:
    std::string project_name_;
    std::string work_directory_;
    std::string filename_coefs_;
    DATA_ANGLES_ dd_;

    // samples where the approximation is tested;
    std::shared_ptr<double []> x_ = nullptr; 

    // target polynomial:
    std::shared_ptr<double[]> pol_ = nullptr;

    // sought-after angles (a half of them):
    std::shared_ptr<double []> phis_ = nullptr; 

    // components of the cost functio:
    std::shared_ptr<double[]> cost_fx_;

    // gradient_phi of the cost function  (ix, i_phi):
    double** grad_cost_; 

    // name of the output file:
    std::string output_name_;

    /**
     * Output \p project_name_.hdf5 file in \p work_directory_ directory.
    */
    YMIX::H5File hfo_; 


public:
    ComputeAngles_(YCS pname) :
        project_name_(pname), work_directory_(std::filesystem::current_path())
    {
        using namespace std::complex_literals;

        output_name_ = "";

        dd_.max_n_iters = 5e4;
        dd_.N_pairs = 200;

        read_input_data();
        create_output_hdf5();

        // --- Save some input parameters into the output .hdf5 file ---
        std::cout << "\nSaving some input parameters into the output .hdf5 file..." << std::endl;
        hfo_.open_w();
        hfo_.add_scalar(dd_.function_type, "function-type",      "basic");
        hfo_.add_scalar(dd_.parity,        "function-parity",    "basic");
        hfo_.add_scalar(dd_.parameter,     "function-parameter", "basic");
        hfo_.add_scalar(dd_.abs_error,     "abs-error",          "basic");
        hfo_.add_scalar(dd_.coef_norm,     "factor-norm",        "basic");
        hfo_.close();

        // Create samples (using roots of Chebyschev polynomials - Chebyschev nodes):
        std::cout << "\nInitialising x-grid and the grid of sougth-after angles..." << std::endl;
        x_ = std::shared_ptr<double[]>(new double[dd_.Nc]);
        for(uint32_t ii = 0; ii < dd_.Nc; ii++)
            x_[ii] = cos((2*ii + 1)*M_PI / (4 * dd_.Nc));

        // for(uint32_t ii = 0; ii < dd_.Nc; ii++)
        //     x_[ii] = sin(x_[ii]);

        // initialize sought-after angles:
        phis_ = std::shared_ptr<double[]>(new double[dd_.Nc]);
        for(uint32_t ii = 0; ii < dd_.Nc; ii++)
            phis_[ii] = 0.;

        // compute the polynomial approximation of the target function:
        compute_target_polynomial();

        // initialize the gradient (ix, i_phi):
        grad_cost_ = new double*[dd_.Nc];
        for(uint32_t ix = 0; ix < dd_.Nc; ix++)
            grad_cost_[ix] = new double[dd_.Nc];
        for(uint32_t ix = 0; ix < dd_.Nc; ix++)
            for(uint32_t i_phi = 0; i_phi < dd_.Nc; i_phi++)
                grad_cost_[ix][i_phi] = 0.0;

        YMIX::YTimer timer;
        timer.StartPrint("\nPerforming the L-BFGS algorithm...\n", 0, false);
        launch_LGFGS();
        timer.StopPrint(0, false);

        hfo_.open_w();
        hfo_.add_scalar(timer.get_dur_s(), "sim-time-seconds", "basic");
        hfo_.close();

        // --- Construct a full set of angles taking into account their symmetry. ---
        // See Eq.(24) in Dong-21.
        std::shared_ptr<double[]> phis_full; 
        uint32_t Nc = dd_.Nc;
        uint32_t Nc_full;
        std::cout << "\nForming the final full set of QSVT angles..." << std::endl;
        if(dd_.parity == 1) // odd d in notations of Dong-21
        {
            Nc_full = 2*dd_.Nc;
            phis_full = std::shared_ptr<double[]>(new double[Nc_full]);
            for(uint32_t ii = 0; ii < Nc; ii++)
                phis_full[ii] = phis_[Nc - 1 - ii];
            for(uint32_t ii = 0; ii < Nc; ii++)
                phis_full[Nc+ii] = phis_[ii];
        }
        else if(dd_.parity == 0) // even d in notations of Dong-21
        {
            Nc_full = 2*dd_.Nc - 1;
            phis_full = std::shared_ptr<double[]>(new double[Nc_full]);
            for(uint32_t ii = 0; ii < Nc; ii++)
                phis_full[ii] = phis_[Nc - 1 - ii];
            for(uint32_t ii = 1; ii < Nc; ii++)
                phis_full[Nc+ii-1] = phis_[ii];
        }
        else
        {
            // for the function without denifite parity;
        }

        // Take into account initial conditions:
        phis_full[0]         += M_PI_4;
        phis_full[Nc_full-1] += M_PI_4;

        // --- Reconstruct the polynomial using the angles ---
        double* rec_pol = new double[dd_.Nc];
        for(uint32_t ix = 0; ix < dd_.Nc; ix++)
            rec_pol[ix] = compute_sequence_of_rotations_FULL(phis_full, Nc_full, x_[ix]);

        std::cout << "Saving the polynomial reconstructed using QSVT angels..." << std::endl;
        hfo_.open_w();
        hfo_.add_array(rec_pol, dd_.Nc, "pol-angles", "results");
        hfo_.close();
        delete [] rec_pol;
        
        // Modify the angles according to Eq.(15) in Dong-21:
        phis_full[0]         += M_PI_4;
        phis_full[Nc_full-1] += M_PI_4;
        for(uint32_t ii = 1; ii < Nc_full-1; ii++)
            phis_full[ii] += M_PI_2;

        // for(uint32_t ii = 0; ii < Nc_full; ii++)
        //     phis_full[ii] += M_PI_2;

        // Save the angles:
        std::cout << "Saving the full set of QSVT angles..." << std::endl;
        hfo_.open_w();
        hfo_.add_array(phis_full.get(), Nc_full, "phis", "results");
        hfo_.close();
    }

    ~ComputeAngles_()
    {
        for(uint32_t ix = 0; ix < dd_.Nc; ix++)
            delete [] grad_cost_[ix];
        delete grad_cost_;

    }


protected:

    void launch_LGFGS()
    {
        using namespace std;
        using namespace std::complex_literals;
        uint32_t Nx = dd_.Nc;
        uint32_t N_phis = dd_.Nc;
        uint32_t counter_iter;

        LBGFS_pairs_ pairs(dd_.N_pairs, N_phis);

        // used for intermediate computations during line-search procedure;
        double cost_mean;
        double cost_mean_new;

        // used to decide whether stop iterations or not; 
        double cost_max = 1.0;
        double err2 = dd_.stopping_criterion * dd_.stopping_criterion;

        YMATH::VectorD_ mean_grad_phi; // of size N_phis
        YMATH::VectorD_ mean_grad_phi_new; // of size N_phis
        YMATH::VectorD_ d; // of size N_phis
        shared_ptr<double[]> alpha       = shared_ptr<double[]>(new double[pairs.N_pairs]);
        shared_ptr<double[]> phis_new    = shared_ptr<double[]>(new double[N_phis]);
        shared_ptr<double[]> cost_fx_new = shared_ptr<double[]>(new double[Nx]);

        // compute initial cost function and its derivative w.r.t. initial angles:
        cost_fx_ = shared_ptr<double[]>(new double[Nx]);
        compute_cost_function(phis_, Nx, cost_fx_); 
        compute_grad_of_cost_function(phis_);

        // find the mean values (across x_) of the cost function and its derivative:
        YMATH::compute_mean(cost_fx_.get(), Nx, cost_mean);
        YMATH::compute_mean(grad_cost_, Nx, N_phis, mean_grad_phi);

        // --- L-BFGS technique ---
        // book Sun-06, Algorithms 5.7.1 and 5.7.2;
        // the condition to stop the iterations is defined by Theorem 4 in Dong-21; 
        counter_iter = 0;
        while(
            cost_max > err2 &&
            counter_iter < dd_.max_n_iters
        ){
            counter_iter++;

            // --- Two-loop recursion: compute d = - H*g --- 
            // --- First loop ---
            pairs.first_loop(d, alpha, mean_grad_phi);

            // multiply by first approximation of the Hessian (Step 3 in Algorithm 5.7.1 in Sun-06):
            d.mult_by(0.5);
            if(dd_.parity == 0)
                d.el(0) *= 2.;

            // --- Second loop ---
            pairs.second_loop(d, alpha);

            // --- Line search: !!! Test signs !!!---
            double w = 0.5;
            double rr = 1.e-3;
            double min_step = 1e-5;

            double step = 1./w; // alpha in backtracking line search Algorithm 2.5.2 in Sun-06.
            double diff = -1.; // f(x) - f(x + alpha*d)
            double comp = 0.; // rr * alpha * gd

            double gd = mean_grad_phi.dot(d); // sign ?

            while(diff < comp && step >= min_step)
            {
                step *= w;
                for(uint32_t ii = 0; ii < N_phis; ii++)
                    phis_new[ii] = phis_[ii] - step * d.v(ii);
                compute_cost_function(phis_new, Nx, cost_fx_new); 
                YMATH::compute_mean(cost_fx_new.get(), Nx, cost_mean_new);
                diff = cost_mean - cost_mean_new;
                comp = rr*step*gd;
            }

            // --- The rest of the algorithm ---
            YMIX::copy_array(phis_new, N_phis, phis_);
            cost_mean = cost_mean_new;
            YMATH::find_max(cost_fx_new, Nx, cost_max);
            compute_grad_of_cost_function(phis_);
            YMATH::compute_mean(grad_cost_, Nx, N_phis, mean_grad_phi_new);

            pairs.N_current = min(pairs.N_current+1, pairs.N_pairs);
            pairs.id_last_saved = (pairs.id_last_saved+1)%pairs.N_pairs; 

            pairs.y[pairs.id_last_saved].set_diff(mean_grad_phi_new, mean_grad_phi);
            pairs.s[pairs.id_last_saved].copy_mult_by(d, -step);
            pairs.rho[pairs.id_last_saved] = 
                1./(pairs.s[pairs.id_last_saved].dot(pairs.y[pairs.id_last_saved])); 
            mean_grad_phi = YMATH::VectorD_(mean_grad_phi_new);

            printf("iter = %d, cost-max = %0.3e, step = %0.3e\n", counter_iter, cost_max, step);
        }

        cout << "Done." << endl;
    }


    void compute_grad_of_cost_function(
        const std::shared_ptr<const double[]>& phis
    ){
        using namespace std;
        using namespace std::complex_literals;

        double x1, diff;
        complex<double> ix2;
        shared_ptr<YMATH::Matrix2_> Wx;
        shared_ptr<YMATH::Matrix2_> U;
        shared_ptr<YMATH::Matrix2_> first_half;
        shared_ptr<YMATH::Matrix2_> grad_U;
        uint32_t Nc = dd_.Nc;

        shared_ptr<YMATH::Matrix2_> G0 = make_shared<YMATH::Matrix2_>(
            exp(1i*M_PI_4), 0., 0., exp(-1i*M_PI_4)
        );

        shared_ptr<complex<double>[]> ephi = 
            shared_ptr<complex<double>[]>(new complex<double>[Nc]);
        for(uint32_t ii = 0; ii < Nc; ii++)
            ephi[ii] = exp(1i*phis[ii]);

        complex<double> temp_vec[2];
        auto vec_left = [&ephi, &temp_vec](YCU j){ 
            temp_vec[0] = ephi[j-1];
            temp_vec[1] = conj(temp_vec[0]);
            return temp_vec;
        };
        auto vec_right = [&ephi, &Nc, &temp_vec](YCU j){ 
            temp_vec[0] = ephi[Nc-j-1];
            temp_vec[1] = conj(temp_vec[0]);
            return temp_vec;
        };

        complex<double> vec_i[2];
        vec_i[0] = 1i;
        vec_i[1] = conj(vec_i[0]);

        shared_ptr<YMATH::Matrix2_[]> left_sequ = shared_ptr<YMATH::Matrix2_[]>(
            new YMATH::Matrix2_[Nc]
        );
        shared_ptr<YMATH::Matrix2_[]> right_sequ = shared_ptr<YMATH::Matrix2_[]>(
            new YMATH::Matrix2_[Nc]
        );

        left_sequ[0].set_I();
        right_sequ[0].copy_elements_from( 
            YMATH::Matrix2_(ephi[Nc-1], 0.0, 0.0, conj(ephi[Nc-1])).dot(G0)
        );

        for(uint32_t ix = 0; ix < Nc; ix++)
        {
            x1  = x_[ix];
            ix2 = 1i*sqrt(1-x1*x1);
            Wx  = make_shared<YMATH::Matrix2_>(x1, ix2, ix2, x1);
            for(uint32_t i_phi = 1; i_phi < Nc; i_phi++)
            {
                left_sequ[i_phi]  = *(left_sequ[i_phi-1].mult_per_el_row(vec_left(i_phi))->dot(Wx));
                right_sequ[i_phi] = *(Wx->mult_per_el_col(vec_right(i_phi))->dot(right_sequ[i_phi-1]));
            }

            U = compute_sequence_of_rotations(phis, x1);
            diff = U->get_v(0,0).real() - pol_[ix];

            if(dd_.parity == 1) // for odd function;
                first_half = right_sequ[Nc-1].transpose()->dot(Wx);
            else  // for even polynomial;
                first_half = right_sequ[Nc-2].transpose()->dot(Wx);

            for(uint32_t i_phi = 0; i_phi < Nc; i_phi++)
            {
                grad_U = first_half->dot(left_sequ[i_phi])->
                    mult_per_el_row(vec_i)->dot(right_sequ[Nc-1-i_phi]);
                grad_cost_[ix][i_phi] = 2.*real(grad_U->get_v(0,0)) * diff;
            }
            if(dd_.parity == 0)
                grad_cost_[ix][0] /= 2.;
        }
    }


    /**
     * @param phis - half of the QSVT angles, the second half is reconstructed using the assumed symmetry of angles;
    */
    std::shared_ptr<YMATH::Matrix2_> compute_sequence_of_rotations(
        const std::shared_ptr<const double[]>& phis,
        const double& x1
    ){
        // see Eqs. (13) and (24) on page 7 in Dong-21
        using namespace std::complex_literals;
        using namespace std;

        auto ix2 = 1i*sqrt(1-x1*x1);

        shared_ptr<YMATH::Matrix2_> Wx = make_shared<YMATH::Matrix2_>(
            x1, ix2, ix2, x1
        );

        // Describes the initial condition on the QSVT angles:
        // (pi/4, 0.0, 0.0, ..., 0.0, 0.0, pi/4).
        shared_ptr<YMATH::Matrix2_> G0 = make_shared<YMATH::Matrix2_>(
            exp(1i*M_PI_4), 0., 0., exp(-1i*M_PI_4)
        );

        shared_ptr<complex<double>[]> ephi = 
            shared_ptr<complex<double>[]>(new complex<double>[dd_.Nc]);
        for(uint32_t ii = 0; ii < dd_.Nc; ii++)
            ephi[ii] = exp(1i*phis[ii]);

        auto A_ephi = [&ephi](YCU j){ 
            return make_shared<YMATH::Matrix2_>(
                ephi[j], 0.0, 0.0, conj(ephi[j])
            );
        };
        auto WxA_ephi = [&ephi, &x1, &ix2](YCU j){ 
            auto temp = ephi[j];
            auto ephi_c = conj(temp);
            return make_shared<YMATH::Matrix2_>(
                x1  * temp, ix2 * ephi_c, 
                ix2 * temp,  x1 * ephi_c
            );
        };

        shared_ptr<YMATH::Matrix2_> U = make_shared<YMATH::Matrix2_>();
        if(dd_.parity == 1) // odd function
        {
            U->copy_elements_from(A_ephi(0));
            for(uint32_t ii = 1; ii < dd_.Nc; ii++)
                U = U->dot(Wx)->dot(A_ephi(ii));
            U = U->dot(G0);

            auto res_U = U->transpose()->dot(Wx)->dot(U);
            return res_U;
        }
        else if(dd_.parity == 0) // even function
        {
            U->set_I();
            for(uint32_t ii = 1; ii < dd_.Nc; ii++)
                U = U->dot(WxA_ephi(ii));
                // U = U->dot(Wx)->dot(WxA_ephi(ii));
            U = U->dot(G0);

            auto res_U = U->transpose()->dot(A_ephi(0))->dot(U);
            return res_U;
        }
        else
        {
            // if the function does not have definite parity
        }
        return nullptr;
    }


    /**
     * @param phis - all QSVT angles;
    */
    double compute_sequence_of_rotations_FULL(
        const std::shared_ptr<const double[]>& phis_full,
        YCU N_phis,
        const double& x1
    ){
        // see Eqs. (13) and (24) on page 7 in Dong-21
        using namespace std::complex_literals;
        using namespace std;

        auto ix2 = 1i*sqrt(1-x1*x1);
        shared_ptr<YMATH::Matrix2_> Wx = make_shared<YMATH::Matrix2_>(
            x1, ix2, ix2, x1
        );

        complex<double>* ephi = new complex<double>[N_phis];
        for(uint32_t ii = 0; ii < N_phis; ii++)
            ephi[ii] = exp(1i*phis_full[ii]);

        auto A_ephi = [&ephi](YCU j){ 
            return make_shared<YMATH::Matrix2_>(
                ephi[j], 0.0, 0.0, conj(ephi[j])
            );
        };
        
        shared_ptr<YMATH::Matrix2_> U = make_shared<YMATH::Matrix2_>();
        U->copy_elements_from(A_ephi(0));
        for(uint32_t ii = 1; ii < N_phis; ii++)
            U = U->dot(Wx)->dot(A_ephi(ii));

        delete [] ephi;    
        return U->get_v(0,0).real();
    }


    // Reconstruct the polynomial approximation of the target function 
    //  using the given coefficients.
    void compute_target_polynomial()
    {
        double cx1;

        std::cout << "Computing the target polynomial using the given coefficients..." << std::endl;

        pol_ = std::shared_ptr<double[]>(new double[dd_.Nc]);
        for(uint32_t ix = 0; ix < dd_.Nc; ix++)
        {
            cx1 = acos(x_[ix]);
            pol_[ix] = 0.;
            if(dd_.parity == 0) // even function
            {
                for(uint32_t ic = 0; ic < dd_.Nc; ic++)
                    pol_[ix] += dd_.coefs[ic] * cos((2*ic)*cx1);
            }
            else if(dd_.parity == 1) // odd function
            {
                for(uint32_t ic = 1; ic <= dd_.Nc; ic++)
                    pol_[ix] += dd_.coefs[ic-1] * cos((2*ic-1)*cx1);
            }
            else
            {
                // if the function does not have definite parity
            }
        }

        // Save the polynomial and the x-grid:
        std::cout << "Saving the polynomial and the x-grid..." << std::endl;
        hfo_.open_w();
        hfo_.add_array(x_.get(),   dd_.Nc, "x",   "results");
        hfo_.add_array(pol_.get(), dd_.Nc, "pol-coefs", "results");
        hfo_.close();
    }


    void create_output_hdf5()
    {
        using namespace std;

        printf("Creating output %s.hdf5 file in the directory:\n", project_name_.c_str());
        printf("%s\n", work_directory_.c_str());

        string fname;
        if(YMIX::compare_strings(output_name_, "")) 
            fname = project_name_;
        else
            fname = output_name_;

        string filename_hdf5 = work_directory_ +"/" + fname + ".hdf5";
        hfo_.create(filename_hdf5);
        hfo_.add_group("basic");
        hfo_.add_group("results");

        // date of simulation:
        string str_date_time;
        YMIX::get_current_date_time(str_date_time);
        hfo_.add_scalar(str_date_time, "date-of-simulation", "basic");

        // save the project name and the path to input files:
        hfo_.add_scalar(project_name_,   "project-name",     "basic");
        hfo_.add_scalar(work_directory_, "work-directory",   "basic");
        hfo_.close();   

        printf("Done\n");
    }


    void read_input_data()
    {
        using namespace std;
        string data;

        printf("\nReading input parameters from the %s.ca file in the directory:\n", project_name_.c_str());
        printf("%s\n", work_directory_.c_str());
        string filename_ca = work_directory_ +"/" + project_name_ + ".ca";

        // remove comments:
        YMIX::read_input_file(data, filename_ca);

        // read the data:
        string word;
        istringstream istr(data);

        while(istr >> word)
        {
            if(YMIX::compare_strings(word, "filename_coefs"))
            {
                istr >> filename_coefs_;

                printf("\n\tReading coefficients from the file: ./%s\n", filename_coefs_.c_str());
                read_coefficients();
            }

            if(YMIX::compare_strings(word, "stopping_criterion"))
            {
                istr >> dd_.stopping_criterion;
            }


            if(YMIX::compare_strings(word, "output_name"))
            {
                cout << "here";
                istr >> output_name_;
            }
        }
        
        // --- Print input parameters ---
        cout << "function type: \t" << dd_.function_type << endl;
        cout << "parity: \t" << dd_.parity << endl;
        cout << "parameter: \t" << dd_.parameter << endl;
        printf("absolute error of polynomial approximation: %0.3e\n", dd_.abs_error);
        cout << "coefficient of normalization: \t" << dd_.coef_norm << endl;
        cout << "number of coefficients: \t" << dd_.Nc << endl;
        cout << "number of L-BFGS pairs: \t" << dd_.N_pairs << endl;
        printf("stopping criterion square: %0.3e\n", 
            dd_.stopping_criterion * dd_.stopping_criterion
        );
        cout << "Output file name: " << output_name_ << ".hdf" << endl;
        cout << "Done." << endl;
    }


    void read_coefficients()
    {
        using namespace std::complex_literals;
        using namespace std;

        string date_sim;
        vector<double> coefs_real;
        // vector<double> coefs_imag;

        string filename_coefs = work_directory_ + "/" + filename_coefs_ + ".hdf5";

        YMIX::H5File ff;
        ff.set_name(filename_coefs);
        ff.open_r();
        ff.read_scalar(date_sim,          "date-of-simulation", "basic");
        ff.read_scalar(dd_.function_type, "descr",              "basic");
        ff.read_scalar(dd_.parameter,     "param",              "basic");
        ff.read_scalar(dd_.coef_norm,     "coef_norm",          "basic");
        ff.read_scalar(dd_.abs_error,     "eps",                "basic");
        ff.read_scalar(dd_.parity,        "parity",             "basic");
        ff.read_vector(coefs_real, "real", "coefs");
        // ff.read_vector(coefs_imag, "imag", "coefs");
        ff.close();

        dd_.Nc = coefs_real.size();
        dd_.coefs = shared_ptr<double[]>(new double[dd_.Nc]);
        for(uint32_t ii = 0; ii < dd_.Nc; ii++)
            // dd_.coefs[ii] = coefs_real[ii] + 1i * coefs_imag[ii];
            dd_.coefs[ii] = coefs_real[ii];

        cout << "\tThe coefficients were computed on: " << date_sim << endl;
        cout << "\tDone." << endl;
    }


    void compute_cost_function(
        const std::shared_ptr<const double[]>& phis,
        YCU Nx,
        std::shared_ptr<double[]>& cost_fx
    ){
        using namespace std;
        double x1;
        double temp;

        shared_ptr<YMATH::Matrix2_> U;
        for(uint32_t ix = 0; ix < Nx; ix++)
        {
            x1 = x_[ix];
            U = compute_sequence_of_rotations(phis, x1);
            temp = U->get_v(0,0).real() - pol_[ix];
            cost_fx[ix] = 0.5 * temp * temp;
        }
    }

};
#endif