#ifndef QLIB_H
#define QLIB_H



#include "../submodules/QuEST/QuEST/include/QuEST.h"
// #include "QuEST.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>

#include <filesystem>

#include <math.h>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>

#include <time.h>
#include <stdlib.h> 
#include <stdio.h> 

#include <limits>
#include <chrono>
#include <stdarg.h>
#include <memory>
#include <numeric>
#include <unistd.h>

#include <cuda_runtime_api.h>
#include "H5Cpp.h"

#include <complex>

#define YCUDA true

// ------------------------------------------
// --- Type synonyms --- 
// ------------------------------------------
#define YCB const bool&

#define YCsh  const short&
#define YCVsh const std::vector<short>&
#define YVsh  std::vector<short>&
#define YVshv std::vector<short>

#define YCU  const uint32_t&
#define YCVU const std::vector<uint32_t>&
#define YVU  std::vector<uint32_t>&
#define YVUv std::vector<uint32_t>

#define YCUL const uint64_t&

#define YCI  const int&
#define YCVI const std::vector<int>&
#define YVI  std::vector<int>&
#define YVIv std::vector<int>

#define YCVT const std::vector<T>&

#define YCS  const std::string&
#define YS   std::string&
#define YCVS const std::vector<std::string>&
#define YVSv std::vector<std::string>

#define YCQR const qreal&
#define YCVQ const std::vector<qreal>&
#define YVQ  std::vector<qreal>&
#define YVQv std::vector<qreal>

#define YCVCo const std::vector<Complex>&
#define YVCo  std::vector<Complex>&
#define YVCov  std::vector<Complex>

#define YPQC YSQ 
#define YISS std::istringstream& 

#define YQCP QCircuit*
#define YSQ  std::shared_ptr<QCircuit>
#define YSG  std::shared_ptr<Gate__>
#define YSB  std::shared_ptr<Box__>
#define YSM  std::shared_ptr<YMATH::YMatrix>
#define YSCQ  std::shared_ptr<const QCircuit>
#define YCCQ const std::shared_ptr<const QCircuit>
#define YCCB const std::shared_ptr<const Box__>
#define YCCM const std::shared_ptr<const YMATH::YMatrix>

#define YCVVI const std::vector<std::vector<int>>&
#define YVVI std::vector<std::vector<int>>&
#define YVVIv std::vector<std::vector<int>>

#define YMBo std::make_shared<Box__>

#define YGV YGlobalVariables

#define LOG_INDENT "    "
#define WHITESPACE " \n\r\t\f\v"

#define YSIZE_CHAR_ARRAY 100000

#define TEX_BEGIN_CIRCUIT    "\\begin{quantikz}[row sep={0.5cm,between origins},column sep={0.3cm}]\n"
#define TEX_VERTICAL_GAP     "\\vspace{20pt}"
#define TEX_CIRCUIT_LENGTH 10

#define COMMENT_ORACLE    "//"
#define FORMAT_ORACLE     ".oracle"s 
#define FORMAT_PROFILE    ".condR_profile"s
#define FORMAT_CIRCUIT    ".circuit"s
#define FORMAT_TEX        ".tex"s
#define FORMAT_ANGLES     ".angles"s
#define FORMAT_LOG        ".log"s
#define FORMAT_QSP        ".qsvt"s
#define FORMAT_HDF5       ".hdf5"
// #define FORMAT_INIT       ".init_state"s
#define FORMAT_INIT       "_INIT_STATE.hdf5"s
#define FORMAT_RANDOM     ".random"s
#define ENDING_FORMAT_OUTPUT "_OUTPUT.hdf5"s
#define ENDING_FORMAT_RESTART "_RESTART.hdf5"s

#define ZERO_ERROR 1e-14
#define PARAMETER_ACCURACY 12


enum SEL_INIT_STATE_PREP {use_init_vector, use_init_oracle};


// ------------------------------------------
// --- Structure with QSVT parameters --- 
// ------------------------------------------
struct QSVT_pars
{
    std::string name;
    std::string filename_angles;

    std::string type;

    YVQv angles_phis_even;
    YVQv angles_phis_odd;
    YVQv angles_phis_arbitrary;

    qreal eps_qsvt;
    qreal f_par;
    qreal rescaling_factor;

    int parity;

    uint32_t n_repeat;
};


// ----------------------------------------------------------
// --- Structure with parameters for a compression gadget --- 
// ----------------------------------------------------------
struct GADGET_pars
{
    std::string name;
    std::string type;
    int N_mult;
    YVIv counter_qubits;
};


// ----------------------------------------------------------
// --- QuCF complex data --- 
// ----------------------------------------------------------
struct QuCF_complex_data
{
    std::map<std::string, QSVT_pars>   qsvt;
    std::map<std::string, GADGET_pars> gadgets;

    // void check_name(YCS name_structure)
    // {
    //     // if(
    //     //     qsvt.find(name_structure) != qsvt.end() or
    //     //     gadgets.find(name_structure) != gadgets.end()
    //     // ) 
    //     //     throw std::string("In QuCF_complex_data: the name " + name_structure + " has alredy been used.");
    // }

    // True: this data has already been read;
    bool check_name(YCS name_structure)
    {
        if(
            qsvt.find(name_structure) != qsvt.end() or
            gadgets.find(name_structure) != gadgets.end()
        ) 
            return true;
        else
            return false;
    }
};


// ------------------------------------------
// --- Control qubits --- 
// ------------------------------------------
struct ControlQubits
{
    // indices of the qubits by which UNIT state an operator is controlled;
    YVIv ids_unit; 

    // indices of the qubits by which ZERO state an operator is controlled;
    YVIv ids_zero; 

    /**
     * indices of the qubits that are used as conditional qubits for an operator in a compression gadget;
     * dependence on these qubits takes place only if 
     *      the corresponding operator is called within a compression gadget;
     */
    YVIv ids_gadget; 
};


// ------------------------------------------
// --- Structure with an initial state --- 
// ------------------------------------------
struct INIT_STATE__
{
    bool flag_defined = false;
    long long b_ampl;   // position of the first element in a state vector to set
    long long n_ampls; // number of elements to set
    YVQv ampl_vec_real;
    YVQv ampl_vec_imag;
};

// ------------------------------------------
// --- Structure with Global variables --- 
// ------------------------------------------
struct YGlobalVariables{
    static int tex_circuit_length;
    static const std::string reg_whole_circuit;
    static const qreal inv_sqrt2;
    static const ComplexMatrix2 mX;
    static const ComplexMatrix2 mY;
    static const ComplexMatrix2 mZ;
    static const ComplexMatrix2 mH;

    static ComplexMatrix2 mRx(YCQR a);
    static ComplexMatrix2 mRy(YCQR a);
    static ComplexMatrix2 mRz(YCQR a);
    static ComplexMatrix2 mRc(YCQR az, YCQR ay);
    static ComplexMatrix2 mPhase(YCQR a);
    static ComplexMatrix2 mPhaseZero(YCQR a);
};

// ------------------------------------------
// --- Math functions --- 
// ------------------------------------------ 
namespace YMATH{
    /** Check if the variable \p x is zero.*/
    bool is_zero(YCQR x);

    /** Convert the integer \p x to an array of bits \p binaryNum :
    * 2 -> {1, 0}; 
    * 6 -> {1, 1, 0} etc.
    * @param[out] binaryNum resulting array of bits. 
    * The array should be initialized with a necessary size in advance.
    */
    void intToBinary(int x, YVsh binaryNum);

    /**
     * @brief Convert an array of bits to a corresponding integer.
     * @param bb array of binary. bb[0] - the most significant bit.
     * @return integer = bb[0]*2**(n - 1) + ... + bb[n-2]*2**1 + bb[n-1]
     */
    long long int binaryToInt(const std::vector<short>& bb);

    /**
     * @brief Two-dimensional matrix of a qreal (which is usually double) type.
     * It is an object that provides some basic manipulations with a two-dimensional 
     * raw pointer.
     */
    class YMatrix{
        public:
            /**
             * @brief Create an empty matrix object whithout reserving any memory.
             */
            YMatrix();

            /**
             * @brief Create a zero matrix (\p Nrows, \p Ncols).
             * @param Nrows is a number of rows in a new matrix.
             * @param Ncols is a number of columns.
             * @param flag_zero if true, set elements to zero.
             */
            YMatrix(YCU Nrows, YCU Ncols, YCB flag_zero=true);

            /**
             * @brief Copy a matrix \p M. 
             */
            YMatrix(YCCM M);
            ~YMatrix();

            /**
             * @brief Create an empty matrix.
             */
            void create(YCU Nrows, YCU Ncols);

            /**
             * @brief Create a zero matrix.
             */
            void create_zero_matrix(YCU Nrows, YCU Ncols);

            /**
             * @brief Create an identity matrix.
             */
            void create_identity_matrix(YCU Nrows, YCU Ncols);

            void create_x10_matrix(YCU Nrows, YCU Ncols);

            /**
             * @brief Gives a raw pointer of the matrix.
             */
            qreal** get_pointer();

            /**
             * @brief Create and return 1-D pointer to the matrix.
             * Matrix elements in the resulting 1D-vector are arranged in the row-major format.
             */
            qreal* get_1d_pointer();

            /**
             * @brief Get a pointer to a transposed matrix.
             * Matrix elements in the resulting 1D-vector are arranged in the column-major format.
             */
            qreal* get_1d_transposed_pointer();

            inline int get_nr(){return nr_;}
            inline int get_nc(){return nc_;}

            void set_squared_from_transposed_1d_matrix(int N, qreal* M);

            inline
            qreal& operator()(YCU id_row, YCU id_col)
            {
                if (id_row >= nr_)
                {
                    std::cerr << "\nError: id-row = " << id_row << ", while n-rows = " << nr_ << std::endl;
                    exit(-1);
                }
                if (id_col >= nc_)
                {
                    std::cerr << "\nError: id-column = " << id_col << ", while n-columns = " << nc_ << std::endl;
                    exit(-1);
                }
                return a_[id_row][id_col];
            }

            inline
            qreal operator()(YCU id_row, YCU id_col) const
            {
                if (id_row >= nr_)
                {
                    std::cerr << "\nError: id-row = " << id_row << ", while n-rows = " << nr_ << std::endl;
                    exit(-1);
                }
                if (id_col >= nc_)
                {
                    std::cerr << "\nError: id-column = " << id_col << ", while n-columns = " << nc_ << std::endl;
                    exit(-1);
                }
                return a_[id_row][id_col];
            }

            /**
             * @brief Print the matrix with precision = \p prec.
             * @param prec precision of matrix elements.
             * @param flag_scientific if true, then print in the scientific notation.
             * @param wc width of every column.
             */
            void print(int prec=3, bool flag_scientific=false, int wc=2);


            /**
             * Save a sparse version of the matrix.
            */
            void form_sparse_format(qreal* nz_values, int* nz_columns, int* nz_rows)
            {
                // count the number of nonzero values:
                int Nnz = 0;
                for(uint64_t ir = 0; ir < nr_; ir++)
                    for(uint64_t ic = 0; ic < nc_; ic++)
                        if(abs(a_[ir][ic]) > ZERO_ERROR)
                            ++Nnz;

                nz_values  = new qreal[Nnz];
                nz_columns = new int[Nnz];
                nz_rows    = new int[nr_+1];

                uint64_t counter = 0;
                for(int ir = 0; ir < nr_; ir++)
                {
                    nz_rows[ir] = counter;
                    for(int ic = 0; ic < nc_; ic++)
                        if(abs(a_[ir][ic]) > ZERO_ERROR)
                        {
                            nz_values[counter]  = a_[ir][ic];
                            nz_columns[counter] = ic;
                            ++counter;
                        }
                }
                nz_rows[nr_] = Nnz;
            }

        protected:
            /**
             * @brief Reserve memory for a new matrix of a known size.
             * The function does not check whether the matrix has been already initialized.
             */
            void create_new();

            /**
             * @brief Free the memory occupied by the matrix.
             */
            void clear();

            /**
             * @brief Set elements to zeros.
             * The function does not check where the matrix has been already initialized.
             */
            void set_zeros();

        private:
            unsigned nr_, nc_;
            qreal** a_ = nullptr;
            std::shared_ptr<qreal[]> a_1d_ = nullptr;
            std::shared_ptr<qreal[]> a_1d_transposed_ = nullptr;
    };



    // ------------------------------------------------------------
    // | e00_, e01_|
    // | e10_, e11_|
    // ------------------------------------------------------------
    class Matrix2_
    {
    protected:
        std::complex<double> e00_, e01_, e10_, e11_;

    public:
        Matrix2_(){}

        Matrix2_(
            const std::shared_ptr<const Matrix2_>& A
        ){
            e00_ = A->e00_;
            e01_ = A->e01_;
            e10_ = A->e10_;
            e11_ = A->e11_;
        }


        Matrix2_(
            const std::complex<double>& e00, 
            const std::complex<double>& e01, 
            const std::complex<double>& e10, 
            const std::complex<double>& e11
        ){
            e00_ = e00;
            e01_ = e01;
            e10_ = e10;
            e11_ = e11;
        }

        void set_I()
        {
            e00_ = 1.0;
            e01_ = 0.0;
            e10_ = 0.0;
            e11_ = 1.0;
        }


        void set_elements(
            const std::complex<double>& e00, 
            const std::complex<double>& e01, 
            const std::complex<double>& e10, 
            const std::complex<double>& e11
        ){
            e00_ = e00;
            e01_ = e01;
            e10_ = e10;
            e11_ = e11;
        }


        void copy_elements_from(
            const std::shared_ptr<const Matrix2_>& A
        ){
            e00_ = A->e00_;
            e01_ = A->e01_;
            e10_ = A->e10_;
            e11_ = A->e11_;
        }


        // Scalar multiplication.
        std::shared_ptr<Matrix2_> dot(
            const std::shared_ptr<const Matrix2_>& B
        ){
            std::shared_ptr<Matrix2_> C = std::make_shared<Matrix2_>();
            C->e00_ = e00_ * B->e00_ + e01_ * B->e10_;
            C->e01_ = e00_ * B->e01_ + e01_ * B->e11_;
            C->e10_ = e10_ * B->e00_ + e11_ * B->e10_;
            C->e11_ = e10_ * B->e01_ + e11_ * B->e11_;
            return C;
        }

        std::shared_ptr<Matrix2_> dot(const Matrix2_& B)
        {
            std::shared_ptr<Matrix2_> C = std::make_shared<Matrix2_>();
            C->e00_ = e00_ * B.e00_ + e01_ * B.e10_;
            C->e01_ = e00_ * B.e01_ + e01_ * B.e11_;
            C->e10_ = e10_ * B.e00_ + e11_ * B.e10_;
            C->e11_ = e10_ * B.e01_ + e11_ * B.e11_;
            return C;
        }


        std::shared_ptr<Matrix2_> mult_per_el_row(std::complex<double>* v)
        {
            std::shared_ptr<Matrix2_> C = std::make_shared<Matrix2_>();
            C->e00_ = e00_ * v[0];
            C->e01_ = e01_ * v[1];
            C->e10_ = e10_ * v[0];
            C->e11_ = e11_ * v[1];
            return C;
        }

        std::shared_ptr<Matrix2_> mult_per_el_col(std::complex<double>* v)
        {
            std::shared_ptr<Matrix2_> C = std::make_shared<Matrix2_>();
            C->e00_ = e00_ * v[0];
            C->e01_ = e01_ * v[0];
            C->e10_ = e10_ * v[1];
            C->e11_ = e11_ * v[1];
            return C;
        }


        // Return the Hermitian adjoint matrix.
        std::shared_ptr<Matrix2_> adjoint_H()
        {
            std::shared_ptr<Matrix2_> adjoint = std::make_shared<Matrix2_>();
            adjoint->e00_ = std::conj(e00_);
            adjoint->e01_ = std::conj(e10_);
            adjoint->e10_ = std::conj(e01_);
            adjoint->e11_ = std::conj(e11_);
            return adjoint;
        }


        std::shared_ptr<Matrix2_> transpose()
        {
            std::shared_ptr<Matrix2_> adjoint = std::make_shared<Matrix2_>();
            adjoint->e00_ = e00_;
            adjoint->e01_ = e10_;
            adjoint->e10_ = e01_;
            adjoint->e11_ = e11_;
            return adjoint;
        }


        std::complex<double> get_v(YCU ir, YCU ic)
        {
            if(ir == 0)
                if(ic == 0) return e00_;
                if(ic == 1) return e01_;
            if(ir == 1)
                if(ic == 0) return e10_;
                if(ic == 1) return e11_;
            return NAN;
        }

        void print()
        {
            using namespace std;

            cout << "\t" << e00_ << "\t" << e01_ << endl;
            cout << "\t" << e10_ << "\t" << e11_ << endl;
        }
    };


    ComplexMatrix2 inv_matrix2(const ComplexMatrix2& a);


    class VectorD_
    {
    protected:
        uint32_t N_;
        std::shared_ptr<double[]> arr_ = nullptr;

    public:
        VectorD_(){}

        VectorD_(YCU N){ init(N); }

        VectorD_(const VectorD_& v_src)
        {
            N_ = v_src.N_;
            arr_ = std::shared_ptr<double[]>(new double[N_]);
            for(uint32_t ii = 0; ii < N_; ii++)
                arr_[ii] = v_src.arr_[ii];
        }

        inline void init(YCU N)
        {
            N_ = N;
            arr_ = std::shared_ptr<double[]>(new double[N_]);
        }


        inline std::shared_ptr<double[]>& get_array(){ return arr_;}

        /**
         * this->arr_ *= v 
        */
        inline void mult_by(const double& v)
        {
            for(uint32_t ii = 0; ii < N_; ii++)
                arr_[ii] *= 0.5;
        }

        /**
         * this->arr_ += v * vec.arr_ 
        */
        inline void plus_mult_by(const VectorD_& vec, const double& v)
        {
            for(uint32_t ii = 0; ii < N_; ii++)
                arr_[ii] += v * vec.arr_[ii];
        }

        inline double& el(YCU i){ return arr_[i]; }
        inline double v(YCU i) const { return arr_[i]; }

        inline double dot(const VectorD_& vec) const
        {
            double res = 0;
            for(uint32_t ii = 0; ii < N_; ii++)
                res += arr_[ii] * vec.arr_[ii];
            return res;
        }

        void set_diff(const VectorD_& vec1, const VectorD_& vec2)
        {
            for(auto ii = 0; ii < N_; ii++)
                arr_[ii] = vec1.arr_[ii] - vec2.arr_[ii];
        }

        void copy_mult_by(const VectorD_& vec, const double& v)
        {
            for(auto ii = 0; ii < N_; ii++)
                arr_[ii] = v * vec.arr_[ii];
        }

    };

    /**
     * @brief Create a vector with integers in the interval [start, end).
     */
    YVIv get_range(YCI start, YCI end);

    /**
     * @brief Whether a string \p str can be converted to a number.
     * @param str is a string to check;
     */
    bool is_number(YCS str);

    /**
     * Compute a mean value of the array \p arr which has \p N elements;
    */
    void compute_mean(const double* arr, YCU N, double& res_mean);

    /**
     * Compute a mean value along the first dimension 
     * of the two-dimensional array \p arr with 
     * \p N1 and \p N2 values along the first and second dimensions, respectively.
    */
    void compute_mean(double** arr, YCU N1, YCU N2, std::shared_ptr<double[]>& res_mean);
    void compute_mean(double** arr, YCU N1, YCU N2, YMATH::VectorD_& res_mean);

    void find_max(const std::shared_ptr<const double[]> arr, YCU N, double& res_max);
}

// ------------------------------------------
// --- Mix functions --- 
// ------------------------------------------
namespace YMIX{
    void print_log(
        YCS line, 
        YCI n_indent = 0, 
        const bool& flag_only_file=false, 
        const bool& flag_new_line=true
    );
    void print_log_flush(YCS line, YCI n_indent=0, YCB flag_file = true, YCB flag_flux = true);
    void print_log_err(YCS line);

    /** Timer*/
    struct YTimer{
        public:           
            void Start(){
                start_ = std::chrono::steady_clock::now();
            }
            void Stop(){
                end_ = std::chrono::steady_clock::now();
            }
            void StartPrint(YCS mess, YCI n_indent=0, YCB flag_file = true, YCB flag_flux = true)
            {
                Start();
                print_log_flush(mess, n_indent, flag_file, flag_flux);
            }
            void StopPrint(YCI n_indent=0, YCB flag_file = true)
            {
                Stop();

                std::ostringstream oss;
                oss << std::scientific << std::setprecision(3) 
                    << get_dur_s() << " s" << std::endl;
                print_log_flush(oss.str(), n_indent, flag_file);
            }
            double get_dur(){
                std::chrono::duration<double> dur_seconds = end_ - start_;
                return 1000.*dur_seconds.count(); // in ms
            }
            double get_dur_s(){
                std::chrono::duration<double> dur_seconds = end_ - start_;
                return dur_seconds.count(); // in seconds
            }
            std::string get_dur_str_ms(){
                std::ostringstream ostr;
                ostr << std::scientific << std::setprecision(3) << get_dur() << " ms";
                return ostr.str();
            }
            std::string get_dur_str_s(){
                std::ostringstream ostr;
                ostr << std::scientific << std::setprecision(3) << get_dur()/1000. << " s";
                return ostr.str();
            }
        protected:
            std::chrono::time_point<std::chrono::steady_clock> start_;
            std::chrono::time_point<std::chrono::steady_clock> end_;
    };


    struct StateVectorOut
    {
        std::list<std::vector<short>> states; // output states.
        std::vector<Complex> ampls; // output state amplitudes.
        uint32_t n_low_prior_qubits; // a number of first low-priority qubits to calculate.

        std::vector<uint32_t> organize_state; // indicates how to print each state:
            // e.g. if n_low_prior_qubits = 4, and organize_state = [2,1,1], then
            // print  |q3 q2> |q1> |q0>;

        std::vector<short> state_to_choose; // array with bits that a state must have to be chosen:
            // Elements of this vector must be 0 or 1. 
            // -1 means that it can be either 0 or 1.
            // The first elements in the array correspond to the high-priority qubits.
        bool flag_str; // if true, then form a string representation of the statevector.
        std::string str_wv; // a string representation of the statevector.
        uint32_t prec; // precision of amplitudes to print in str_wv;

        StateVectorOut()
        {
            n_low_prior_qubits = 0;
            flag_str = true;
            prec = 3;
            state_to_choose = std::vector<short>{};
        }
    };


    /** Form a string representation of the statevector. */
    void getStrWavefunction(StateVectorOut& out);


    struct File
    {
        std::ofstream of; // stream connected to a file;
        
        /** Open a file.
         * @param[in] fileName name of the file;
         * @param[in] flagNew  is it a new file (if new, old data will be deleted);
         * */
        File(YCS fileName, bool flagNew=false)
        {
            if(fileName.empty())
            {
                std::cerr << "Error: File name is not set." << std::endl;
                exit(-1);
            }

            if(flagNew)
                of.open(fileName);
            else
                of.open(fileName, std::ios::app);
            if(!of)
            {
                std::cerr << "Error: It's impossible to open the file:\n" << fileName << std::endl;
                exit(-1);
            }
        }

        ~File()
        {
            if(of.is_open())
            {
                of.clear();
                of.close();
            }
        }

        template<class T>
        std::ostream& operator<<(const T& x)
        {
            of << x;
            return of;
        }
    };


    struct LogFile : File
    {
        static std::string name_global_;
        LogFile(bool flagNew=false) : File(name_global_, flagNew){}
    };

    /**
     * @brief Remove a comment from a line.
     * @param line line from where to remove a comment;
     * @return a new line;
     */
    std::string remove_comment(YCS line);

    // true if line1 == line2
    bool compare_strings(YCS line1, YCS line2);

    // true if line1 == line2 and in lines
    bool compare_strings(YCS line1, YCS line2, YCVS lines);

    // true if line1 in lines
    bool compare_strings(YCS line1, YCVS lines);


    // replace substrings in a string:
    void replace_substrings(std::string& init_line, YCS susstr, YCS new_substr);


    template<class T>
    std::vector<T> conc_vectors(YCVT v1, YCVT v2)
    {
        std::vector<T> vv(v1);
        copy(v2.begin(), v2.end(), back_inserter(vv));

        return vv;
    }
    template<class T>
    std::vector<T> conc_vectors(YCVT v1, YCVT v2, YCVT v3)
    {
        std::vector<T> vv(v1);
        copy(v2.begin(), v2.end(), back_inserter(vv));
        copy(v3.begin(), v3.end(), back_inserter(vv));
        return vv;
    }
    template<class T>
    std::vector<T> conc_vectors(YCVT v1, YCVT v2, YCVT v3, YCVT v4)
    {
        std::vector<T> vv(v1);
        copy(v2.begin(), v2.end(), back_inserter(vv));
        copy(v3.begin(), v3.end(), back_inserter(vv));
        copy(v4.begin(), v4.end(), back_inserter(vv));
        return vv;
    }

    void print(const ComplexMatrix2& a, YCI prec = 3);

    template<class T>
    void print(const std::vector<T>& a, YCI prec = 3, const bool& flag_scientific = false)
    {
        for(auto& x: a)
            if(flag_scientific)
                std::cout << std::scientific << std::setprecision(prec) <<  x << " ";
            else
                std::cout << std::setprecision(prec) << x << " ";
        std::cout << std::endl;
    }

    template<class T>
    void print(
        std::ofstream& of, const std::vector<T>& a, 
        YCI prec = 3, const bool& flag_scientific = false
    ){
        for(auto& x: a)
            if(flag_scientific)
                of << std::scientific << std::setprecision(prec) <<  x << " ";
            else
                of << std::setprecision(prec) <<  x << " ";
        of << "\n";
    }

    void print(std::vector<int> a);
    void print(short* a, const unsigned long& nr, const unsigned long& nc);

    std::string get_line(std::vector<int> a);

    std::string ltrim(YCS s);
    std::string rtrim(YCS s);
    std::string trim(YCS s);

    void insert_indent(YS line_original, YCS line_indent);

    bool is_present(YCVS v, YCS e);
    bool is_present(YCVI v, YCI e);

    void get_array_from_list(
        const std::list<std::vector<short>>& v, 
        short* array_1d, 
        const unsigned long& nr, 
        const unsigned long& nc
    );

    void get_current_date_time(YS line_date_time);
    
    struct H5File
    {
        /**
         * @brief Create and open an .hdf5 file with a name \p fname. 
         */
        void create(YCS fname);
        void close();

        inline void set_name(YCS fname){ name_ = fname; }

        /**
         * @brief Open an .hdf5 file with a name \p fname only to read it. 
         */
        void open_r();

        /**
         * @brief Open an .hdf5 file with a name \p fname to write-read it. 
         */
        void open_w();

        /**
         * @brief Add a group (folder) with a name \p gname to an already opened file.
         */
        void add_group(YCS gname);

        bool is_exist(YCS grp_name)
        {
            return H5Lexists(f_->getId(), grp_name.c_str(), H5P_DEFAULT) > 0;
        }

        /**
         * @brief Add a dataset with a name \p dname, where a scalar \p v is to be written.
         * The dataset is put to a group \p gname.
         */
        template<class T>
        void add_scalar(const T& v, YCS dname, YCS gname)
        {
            if(!flag_opened) 
                throw "HDF5 File " + name_ + 
                    " is not opened to add a dataset " + dname + " to a group " + gname;
            // add_group(gname);

            H5::Group grp(f_->openGroup(gname));
            write(v, dname, grp);
        }

        template<class T>
        void add_vector(const std::vector<T>& v, YCS dname, YCS gname)
        {
            if(!flag_opened) 
                throw "HDF5 File " + name_ + 
                    " is not opened to add a dataset " + dname + " to a group " + gname;
            // add_group(gname);

            H5::Group grp(f_->openGroup(gname));
            write(v, dname, grp);
        }

        template<class T>
        void add_array(const T* v, YCUL N, YCS dname, YCS gname)
        {
            if(!flag_opened) 
                throw "HDF5 File " + name_ + 
                    " is not opened to add a dataset " + dname + " to a group " + gname;
            // add_group(gname);

            H5::Group grp(f_->openGroup(gname));
            write(v, N, dname, grp);
        }

        template<class T>
        void add_matrix(const std::list<std::vector<T>>& v, YCS dname, YCS gname)
        {
            if(!flag_opened) 
                throw "HDF5 File " + name_ + 
                    " is not opened to add a dataset " + dname + " to a group " + gname;
            // add_group(gname);

            T* array_1d;
            unsigned long nr, nc;

            nr = v.size();
            nc = v.back().size();
            array_1d = new T[nr*nc];

            YMIX::get_array_from_list(v, array_1d, nr, nc);

            H5::Group grp(f_->openGroup(gname));
            write(array_1d, nc, nr, dname, grp);

            delete [] array_1d;
        }

        template<class T>
        void read_scalar(T& v, YCS dname, YCS gname)
        {
            if(!flag_opened) 
                throw "HDF5 File " + name_ + 
                    " is not opened to read a dataset " + dname + " from a group " + gname;    
            H5::Group grp(f_->openGroup(gname));
            read(v, dname, grp);
        }

        template<class T>
        void read_vector(std::vector<T>& v, YCS dname, YCS gname)
        {
            if(!flag_opened) 
                throw "HDF5 File " + name_ + 
                    " is not opened to add a dataset " + dname + " to a group " + gname;
            H5::Group grp(f_->openGroup(gname));
            read(v, dname, grp);
        }


        protected:
            inline void write(YCS v, YCS dname, H5::Group& grp)
            {
                auto dspace = H5::DataSpace(H5S_SCALAR);
                H5::StrType dtype(H5::PredType::C_S1, v.size()+1);
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write(v, dtype);
            }
            inline void write(YCI v, YCS dname, H5::Group& grp)
            {
                auto dspace = H5::DataSpace(H5S_SCALAR);
                auto dtype = H5::PredType::NATIVE_INT;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write((int*) &v, dtype);
            }
            inline void write(YCU v, YCS dname, H5::Group& grp)
            {
                auto dspace = H5::DataSpace(H5S_SCALAR);
                auto dtype = H5::PredType::NATIVE_UINT;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write((unsigned*) &v, dtype);
            }
            inline void write(const long unsigned int& v, YCS dname, H5::Group& grp)
            {
                auto dspace = H5::DataSpace(H5S_SCALAR);
                auto dtype = H5::PredType::NATIVE_ULONG;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write((long unsigned int*) &v, dtype);
            }
            inline void write(const double& v, YCS dname, H5::Group& grp)
            {
                auto dspace = H5::DataSpace(H5S_SCALAR);
                auto dtype = H5::PredType::NATIVE_DOUBLE;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write((int*) &v, dtype);
            }

            inline void write(YCVU v, YCS dname, H5::Group& grp)
            {
                hsize_t dims[] = {v.size()};
                H5::DataSpace dspace(1, dims);
                auto dtype = H5::PredType::NATIVE_UINT;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write(&v[0], dtype);
            }
            inline void write(YCVI v, YCS dname, H5::Group& grp)
            {
                hsize_t dims[] = {v.size()};
                H5::DataSpace dspace(1, dims);
                auto dtype = H5::PredType::NATIVE_INT;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write(&v[0], dtype);
            }

            inline void write(YCVsh v, YCS dname, H5::Group& grp)
            {
                hsize_t dims[] = {v.size()};
                H5::DataSpace dspace(1, dims);
                auto dtype = H5::PredType::NATIVE_SHORT;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write(&v[0], dtype);
            }

            inline void write(const std::vector<double>& v, YCS dname, H5::Group& grp)
            {
                hsize_t dims[] = {v.size()};
                H5::DataSpace dspace(1, dims);
                auto dtype = H5::PredType::NATIVE_DOUBLE;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write(&v[0], dtype);
            }
            inline void write(const std::vector<Complex>& v, YCS dname, H5::Group& grp)
            {
                hsize_t dims[] = {v.size()};
                H5::DataSpace dspace(1, dims);

                hid_t dtype = H5Tcreate(H5T_COMPOUND, sizeof(Complex));
                H5Tinsert (dtype, "real", HOFFSET(Complex,real), H5T_NATIVE_DOUBLE);
                H5Tinsert (dtype, "imag", HOFFSET(Complex,imag), H5T_NATIVE_DOUBLE);

                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write(&v[0], dtype);
            }

            inline void write(const double* v, YCUL N, YCS dname, H5::Group& grp)
            {
                hsize_t dims[] = {N};
                H5::DataSpace dspace(1, dims);
                auto dtype = H5::PredType::NATIVE_DOUBLE;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write(v, dtype);
            }

            inline void write(short* v, const unsigned long& nr, const unsigned long& nc, YCS dname, H5::Group& grp)
            {
                hsize_t dims[] = {nr, nc};
                H5::DataSpace dspace(2, dims);
                auto dtype = H5::PredType::NATIVE_SHORT;
                H5::DataSet dataset = grp.createDataSet(dname, dtype, dspace);
                dataset.write(v, dtype);
            }


            template<class T>
            inline void read(T& v, YCS dname, H5::Group& grp)
            {
                if(H5Lexists(grp.getId(), dname.c_str(), H5P_DEFAULT))
                {
                    H5::DataSet dataset = grp.openDataSet(dname);
                    H5::DataType dtype = dataset.getDataType();
                    dataset.read(&v, dtype);
                } 
                else
                {
                    std::cout << "\n>>> HDF5 WARNING READING: the dataset " << dname 
                        << " does not exist in the group " << grp.getObjName() << ".\n" << std::endl;
                } 
            }
            inline void read(YS v, YCS dname, H5::Group& grp)
            {
                H5::DataSet dataset = grp.openDataSet(dname);
                H5::DataType dtype = dataset.getDataType();
                v="";
                dataset.read(v, dtype);
            }
            template<class T>
            inline void read(std::vector<T>& v, YCS dname, H5::Group& grp)
            {
                H5::DataSet dataset = grp.openDataSet(dname);

                H5::DataSpace dataspace = dataset.getSpace();
                int rank = dataspace.getSimpleExtentNdims();
                hsize_t dims_out[rank];
                int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

                unsigned long long N = 1;
                for(unsigned i_dim = 0; i_dim < rank; i_dim++)
                    N *= dims_out[i_dim];
                v = std::vector<T>(N);

                H5::DataType dtype = dataset.getDataType();

                dataset.read(&v[0], dtype, dataspace, dataspace);
            }

        protected:
            std::string name_;
            bool flag_opened;
            H5::H5File* f_;
    };


    void read_init_state(YCS fname, YVQ v_real, YVQ v_imag);

    void read_input_file(YS data, YCS file_name);

    /**
     * Copy \p N elements from the array \p source_arr to the array \p res_arr.
    */
    void copy_array(
        const std::shared_ptr<const double[]>& source_arr, 
        YCU N, 
        std::shared_ptr<double[]>& res_arr
    );
}


#endif