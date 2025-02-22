#include "../include/QLib.h"
using namespace std;

std::string YMIX::LogFile::name_global_ = "output.log";


int YGV::tex_circuit_length = TEX_CIRCUIT_LENGTH;
const string YGV::reg_whole_circuit = "the_whole_circuit";
const qreal YGV::inv_sqrt2 = 1./sqrt(2.);

const ComplexMatrix2 YGV::mX = 
    {
        .real = {
            {0., 1.}, 
            {1., 0.}
        },
        .imag = {{0., 0.}, {0., 0.}}
    };

const ComplexMatrix2 YGV::mY = 
    {
        .real = {
            {0., 0.}, 
            {0., 0.}
        },
        .imag = {
            {0., -1.}, 
            {1.,  0.}
        }
    };

const ComplexMatrix2 YGV::mZ = 
    {
        .real = {
            {1.,  0.}, 
            {0., -1.}
        },
        .imag = {{0., 0.}, {0., 0.}}
    };

const ComplexMatrix2 YGV::mH = 
    {
        .real = {
            {YGV::inv_sqrt2,  YGV::inv_sqrt2}, 
            {YGV::inv_sqrt2, -YGV::inv_sqrt2}
        },
        .imag = {{0., 0.}, {0., 0.}}
    };

ComplexMatrix2 YGV::mRx(YCQR a)
{
    qreal a2 = a/2.;
    ComplexMatrix2 res = {
        .real = {
            {cos(a2),      0.},
            {     0., cos(a2)}
        },
        .imag = {
            {      0., -sin(a2)}, 
            {-sin(a2),       0.}
        }
    };
    return res;
}

ComplexMatrix2 YGV::mRy(YCQR a)
{
    qreal a2 = a/2.;
    ComplexMatrix2 res = {
        .real = {
            {cos(a2), -sin(a2)},
            {sin(a2),  cos(a2)}
        },
        .imag = {{0., 0.}, {0., 0.}}
    };
    return res;
}

ComplexMatrix2 YGV::mRz(YCQR a)
{
    qreal a2 = a/2.;
    ComplexMatrix2 res = {
        .real = {
            {cos(a2),      0.},
            {0.,      cos(a2)}
        },
        .imag = {
            {-sin(a2),      0.}, 
            {0.,       sin(a2)}
        }
    };
    return res;
}

ComplexMatrix2 YGV::mRc(YCQR az, YCQR ay)
{
    qreal az2 = az/2.;
    qreal ay2 = ay/2.;

    // Ry(ay) * Rz(az)
    ComplexMatrix2 res = {
        .real = {
            {cos(az2)*cos(ay2), -cos(az2)*sin(ay2)},
            {cos(az2)*sin(ay2),  cos(az2)*cos(ay2)}
        },
        .imag = {
            {-sin(az2)*cos(ay2),  -sin(az2)*sin(ay2)}, 
            {-sin(az2)*sin(ay2),  sin(az2)*cos(ay2)}
        }
    };
    return res;
}

ComplexMatrix2 YGV::mPhase(YCQR a)
{
    ComplexMatrix2 res = {
        .real = {
            {1., 0.},
            {0., cos(a)}
        },
        .imag = {
            {0., 0.}, 
            {0., sin(a)}
        }
    };
    return res;
}


ComplexMatrix2 YGV::mPhaseZero(YCQR a)
{
    ComplexMatrix2 res = {
        .real = {
            {cos(a), 0.},
            {0.,     1.}
        },
        .imag = {
            {sin(a), 0.}, 
            {0.,     0.}
        }
    };
    return res;
}


bool YMATH::is_zero(YCQR x)
{
    // if(abs(x) < std::numeric_limits<double>::min())
    if(abs(x) < ZERO_ERROR)
        return true;
    else
        return false;
}

void YMATH::intToBinary(int x, std::vector<short>& binaryNum)
{
    int i = 0;
    while (x > 0) {
        binaryNum[i] = x % 2;
        x = x / 2;
        i++;
    }
    std::reverse(binaryNum.begin(), binaryNum.end());
}

long long int YMATH::binaryToInt(const std::vector<short>& bb)
{
    int ii = 0;
    int count_b = bb.size();
    for(auto& b: bb)
    {
        --count_b;
        ii += b * int(pow(2, count_b));
    }
    return ii;
}

YMATH::YMatrix::YMatrix()
{
    nr_ = 0;
    nc_ = 0;
}
YMATH::YMatrix::YMatrix(YCU Nrows, YCU Ncols, YCB flag_zero)
    :nr_(Nrows), nc_(Ncols)
{
    create_new();
    if(flag_zero) set_zeros();
}
YMATH::YMatrix::YMatrix(YCCM oo)
{
    // copy the matrix oo to this object:
    nr_ = oo->nr_;
    nc_ = oo->nc_;
    create_new();

    for(int i = 0; i < nr_; ++i)
        for(int k = 0; k < nc_; ++k)
            a_[i][k] = oo->a_[i][k];
}

YMATH::YMatrix::~YMatrix()
{
    clear();
}

void YMATH::YMatrix::create_new()
{
    a_ = new qreal*[nr_];
    for(int i = 0; i < nr_; ++i)
        a_[i] = new qreal[nc_];
}

void YMATH::YMatrix::clear()
{
    if(nc_ != 0 || nr_ != 0)
    {
        // cout << "delete a matrix" << endl;
        for(unsigned i = 0; i < nr_; ++i) 
            delete [] a_[i];
        delete [] a_;

        nc_ = 0;
        nr_ = 0;

        a_1d_.reset();
    }
}

void YMATH::YMatrix::set_zeros()
{
    for(int i = 0; i < nr_; ++i)
        for(int k = 0; k < nc_; ++k)
            a_[i][k] = 0.0;
}

void YMATH::YMatrix::create(YCU Nrows, YCU Ncols)
{
    clear();
    nc_ = Ncols;
    nr_ = Nrows;
    create_new();
}

void YMATH::YMatrix::set_squared_from_transposed_1d_matrix(int N, qreal* M)
{
    clear();
    nc_ = N;
    nr_ = N;
    create_new();
    
    for(int i = 0; i < nr_; ++i)
        for(int k = 0; k < nc_; ++k)
            a_[i][k] = M[k*N + i];
}

void YMATH::YMatrix::create_zero_matrix(YCU Nrows, YCU Ncols)
{
    clear();
    nc_ = Ncols;
    nr_ = Nrows;
    create_new();

    set_zeros();
}

void YMATH::YMatrix::create_identity_matrix(YCU Nrows, YCU Ncols)
{
    clear();
    nc_ = Ncols;
    nr_ = Nrows;
    create_new();

    for(int i = 0; i < nr_; ++i)
        for(int k = 0; k < nc_; ++k)
            if(i == k) a_[i][k] = 1.0;
            else       a_[i][k] = 0.0;
}

void YMATH::YMatrix::create_x10_matrix(YCU Nrows, YCU Ncols)
{
    clear();
    nc_ = Ncols;
    nr_ = Nrows;
    create_new();

    unsigned coef_10 = 0;
    for(unsigned i = 0; i < nr_; ++i)
    {
        coef_10 = 10*i;
        for(unsigned k = 0; k < nc_; ++k)
            a_[i][k] = coef_10 + k;
    }
}

qreal** YMATH::YMatrix::get_pointer()
{
    return a_;
}

qreal* YMATH::YMatrix::get_1d_pointer()
{
    if(!a_1d_ && nc_ > 0 && nr_ > 0)
    {
        a_1d_ = shared_ptr<qreal[]>(new qreal[nc_*nr_]);
        for(unsigned ir = 0; ir < nr_; ir++)
            for(unsigned ic = 0; ic < nc_; ic++)
                a_1d_[ir*nc_ + ic] = a_[ir][ic];
    }
    return a_1d_.get();
}

qreal* YMATH::YMatrix::get_1d_transposed_pointer()
{
    if(!a_1d_transposed_ && nc_ > 0 && nr_ > 0)
    {
        a_1d_transposed_ = shared_ptr<qreal[]>(new qreal[nc_*nr_]);
        for(unsigned ic = 0; ic < nc_; ic++)
            for(unsigned ir = 0; ir < nr_; ir++)
                a_1d_transposed_[ic*nr_ + ir] = a_[ir][ic];
    }
    return a_1d_transposed_.get();
}

void YMATH::YMatrix::print(int prec, bool flag_scientific, int wc)
{
    for(unsigned i = 0; i < nr_; ++i)
    {
        for(unsigned k = 0; k < nc_; ++k)
            if(flag_scientific)
                std::cout << std::scientific << std::setprecision(prec) << std::setw(prec+wc) << a_[i][k] << " ";
            else
                std::cout << std::fixed << std::setprecision(prec) << std::setw(prec+wc+1) << a_[i][k] << " ";
        std::cout << "\n";
    }
    std::cout << std::endl;
}

ComplexMatrix2 YMATH::inv_matrix2(const ComplexMatrix2& a)
{
    ComplexMatrix2 res;
    for(unsigned i = 0; i < 2; ++i)
        for(unsigned k = 0; k < 2; ++k)
        {
            res.real[i][k] =   a.real[k][i];
            res.imag[i][k] = - a.imag[k][i];
        }
    return res;
}

YVIv YMATH::get_range(YCI start, YCI end)
{
    int n = (end - start);
    YVIv res(n);
    std::iota(res.begin(), res.end(), start);
    return res;
}

bool YMATH::is_number(YCS str)
{
    for (char const &c : str) {
        if (std::isdigit(c) == 0) return false;
    }
    return true;
}


void YMATH::compute_mean(const double* arr, YCU N, double& res_mean)
{
    res_mean = 0.0;
    for(uint32_t ii = 0; ii < N; ii++)
        res_mean += arr[ii];
    res_mean /= N;
}

void YMATH::compute_mean(
    double** arr, 
    YCU N1, YCU N2, 
    std::shared_ptr<double[]>& res_mean
){
    res_mean = std::shared_ptr<double[]>(new double[N2]);
    for(uint32_t i2 = 0; i2 < N2; i2++)
    {
        res_mean[i2] = 0.0;
        for(uint32_t i1 = 0; i1 < N1; i1++)
            res_mean[i2] += arr[i1][i2];
        res_mean[i2] /= N2;
    }
}


void YMATH::compute_mean(double** arr, YCU N1, YCU N2, YMATH::VectorD_& v_mean)
{
    v_mean.init(N2);
    auto res_mean = v_mean.get_array();
    for(uint32_t i2 = 0; i2 < N2; i2++)
    {
        res_mean[i2] = 0.0;
        for(uint32_t i1 = 0; i1 < N1; i1++)
            res_mean[i2] += arr[i1][i2];
        res_mean[i2] /= N2;
    }
}


void YMATH::find_max(const std::shared_ptr<const double[]> arr, YCU N, double& res_max)
{
    double temp;
    res_max = arr[0];
    for(uint32_t ii = 1; ii < N; ii++)
    {
        temp = arr[ii];
        if(temp > res_max)
            res_max = temp;
    }
}




void YMIX::getStrWavefunction(StateVectorOut& out)
{ 
    out.str_wv = "";
    if(!out.flag_str) return;

    unsigned n = out.states.front().size();
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(out.prec);

    unsigned count_i = 0;
    std::string str_state;
    Complex aa;
    qreal ar, ai;
    unsigned w_str  = 9 + out.prec;
    for(auto& one_state:out.states){
        str_state = "|";
        if(out.organize_state.empty())
            for(auto& one_qubit_state:one_state)
                str_state += std::to_string(one_qubit_state);
        else
        {
            unsigned count_org = 0;
            unsigned prev_sum = 0;
            for(unsigned jj = 0; jj < n; ++jj){
                str_state += std::to_string(one_state[jj]);
                if(jj == prev_sum + out.organize_state[count_org]-1){
                    prev_sum += out.organize_state[count_org];
                    ++count_org;
                    if(count_org < out.organize_state.size())
                        str_state += ">|";
                }
            }
        }
        str_state += ">";

        oss.str(std::string());

        aa = out.ampls.at(count_i);
        ar = aa.real; 
        ai = aa.imag;

        oss << std::setw(w_str) << ar;
        oss << std::setw(w_str) << ai << "j";    
        out.str_wv += oss.str() + "   " + str_state + "\n";

        ++count_i;
    }
}


void YMIX::print(const ComplexMatrix2& a, YCI prec)
{
    for(unsigned i = 0; i < 2; ++i)
    {
        for(unsigned k = 0; k < 2; ++k)
        {
            cout << scientific << setprecision(prec) << 
                a.real[i][k] << " " << a.imag[i][k] << "j    ";
        }
        cout << endl;
    }
        
}

string YMIX::remove_comment(YCS line)
{
    string new_line = line.substr(0, line.find(COMMENT_ORACLE, 0));
    return new_line;
}

bool YMIX::compare_strings(YCS line1, YCS line2)
{
    string new_line1(line1), new_line2(line2);
    transform(new_line1.begin(), new_line1.end(), new_line1.begin(), ::tolower);
    transform(new_line2.begin(), new_line2.end(), new_line2.begin(), ::tolower);

    if(new_line1.compare(new_line2) == 0)
        return true;
    else
        return false;
}

bool YMIX::compare_strings(YCS line1, YCS line2, YCVS lines)
{
    if(compare_strings(line1, line2))
        for(const auto& line_one: lines)
            if(compare_strings(line1, line_one))
                return true;
    return false;
}

bool YMIX::compare_strings(YCS line1, YCVS lines)
{
    for(const auto& line_one: lines)
        if(compare_strings(line1, line_one))
            return true;
    return false;
}


void YMIX::replace_substrings(std::string& init_line, YCS susstr, YCS new_substr)
{
    if(new_substr.empty())
        return;

    size_t start_pos = 0;
    while((start_pos = init_line.find(susstr, start_pos)) != std::string::npos) {
        init_line.replace(start_pos, susstr.length(), new_substr);
        start_pos += new_substr.length(); // Handles case where 'to' is a substring of 'from'
    }
}



void YMIX::print(std::vector<int> a)
{
    for(auto const& a1: a) std::cout << a1 << "  ";
    std::cout << std::endl;
}

string YMIX::get_line(std::vector<int> a)
{
    ostringstream inf;
    for(auto const& a1: a) inf << a1 << "  ";
    return inf.str();
}

void YMIX::print_log(
    YCS line, 
    YCI n_indent, 
    const bool& flag_only_file,
    const bool& flag_new_line
)
{
    string line_print = line;

    if(n_indent>0)
    {
        string str_indent = "";
        for(unsigned i = 0; i < n_indent; i++) str_indent += LOG_INDENT;
        insert_indent(line_print, str_indent);
    }
    
    YMIX::LogFile cf;
    cf << line_print; 
    if(flag_new_line) cf << "" << endl;

    if(!flag_only_file)
    {
        cout << line_print;
        if(flag_new_line) cout << "" << endl;
    }
}
void YMIX::print_log_flush(YCS line, YCI n_indent, YCB flag_file, YCB flag_flux)
{
    string line_print = line;
    if(n_indent>0)
    {
        string str_indent = "";
        for(unsigned i = 0; i < n_indent; i++) str_indent += LOG_INDENT;
        insert_indent(line_print, str_indent);
    }
    
    if(flag_file)
    {
        YMIX::LogFile cf;
        if(flag_flux)
            cf << line_print << flush;
        else
            cf << line_print << endl;   
    }
    if(flag_flux)
        cout << line_print << flush;
    else
        cout << line_print << endl;
}
void YMIX::print_log_err(YCS line)
{
    throw line;
}

string YMIX::ltrim(YCS s)
{
    size_t start = s.find_first_not_of(WHITESPACE);
    return (start == string::npos) ? "" : s.substr(start);
}

string YMIX::rtrim(YCS s)
{
    size_t end = s.find_last_not_of(WHITESPACE);
    return (end == string::npos) ? "" : s.substr(0, end + 1);
}

string YMIX::trim(YCS s) {
    return rtrim(ltrim(s));
}

void YMIX::insert_indent(YS line_original, YCS line_indent)
{
    stringstream sstr(line_original);
    string one_line;

    getline(sstr, one_line);
    string line_res = line_indent + one_line;

    while(getline(sstr, one_line)) line_res += "\n" + line_indent + one_line;
    line_original = line_res;
}

bool YMIX::is_present(YCVS v, YCS e)
{
    if(find(v.begin(), v.end(), e) == v.end())
        return false;
    return true;
}
bool YMIX::is_present(YCVI v, YCI e)
{
    if(find(v.begin(), v.end(), e) == v.end())
        return false;
    return true;
}

void YMIX::get_current_date_time(YS line_date_time)
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer, sizeof(buffer), "%m-%d-%Y %H:%M:%S", timeinfo);
    line_date_time = string(buffer);
}

void YMIX::H5File::create(YCS fname)
{

    if(!fname.empty()) set_name(fname);
    f_ = new H5::H5File(name_.c_str(), H5F_ACC_TRUNC);
    flag_opened = true;
}

void YMIX::H5File::open_r()
{
    f_ = new H5::H5File(name_.c_str(), H5F_ACC_RDONLY);
    flag_opened = true;
}

void YMIX::H5File::open_w()
{
    f_ = new H5::H5File(name_.c_str(), H5F_ACC_RDWR);
    flag_opened = true;
}

void YMIX::H5File::close()
{
    delete f_;
    flag_opened = false;
}

void YMIX::H5File::add_group(YCS gname)
{
    // if(!flag_opened) throw "HDF5 File " + name_ + " is not opened. One cannot add a group " + gname;

    // if(find(grp_names_.begin(), grp_names_.end(), gname) == grp_names_.end())
    // {
    //     grp_names_.push_back(gname);
    //     H5::Group grp(f_->createGroup(gname));
    // }
    H5::Group grp(f_->createGroup(gname));
}


void YMIX::get_array_from_list(
    const std::list<std::vector<short>>& v, 
    short* array_1d, 
    const unsigned long& nr, 
    const unsigned long& nc
){
    // !!! assume that every vector has the same amount of elements !!!
    unsigned long count_r = 0;
    unsigned long count_c = 0;
    for(auto const& one_row: v)
    {
        count_c = 0;
        for(auto const& el: one_row)
        {
            array_1d[count_c*nr + count_r] = el; // inversed row <-> column ordering
            ++count_c;
        }
        ++count_r;
    }
}


void YMIX::print(short* a, const unsigned long& nr, const unsigned long& nc)
{
    unsigned long long count = -1;
    for(unsigned long ir = 0; ir < nr; ir++)
    {
        for(unsigned long ic = 0; ic < nc; ic++)
        {
            ++count;
            cout << setw(6) << a[count] << " ";
        }
        cout << "\n";
    }
}


void YMIX::read_init_state(YCS fname, YVQ v_real, YVQ v_imag)
{
    YMIX::H5File ff;
    ff.set_name(fname);
    ff.open_r();
    ff.read_vector(v_real, "real", "init_state");
    ff.read_vector(v_imag, "imag", "init_state");
}



void YMIX::read_input_file(YS data, YCS file_name)
{   
    string file_name_res = file_name;
    ifstream ff(file_name_res);
    if(!ff.is_open()) throw "Error: the file [" + file_name_res + "] does not exist.";
    data = string((istreambuf_iterator<char>(ff)), istreambuf_iterator<char>());
    ff.close();

    // clean the buffer from empty lines and comments:
    istringstream istr(data);
    string data_clr = "";
    string line;
    while(getline(istr, line))
    {
        line = YMIX::remove_comment(line);
        line = YMIX::trim(line);
        if(line.find_first_not_of(' ') == string::npos)
            continue;
        data_clr += line + "\n";
    }
    // std::transform(data_clr.begin(), data_clr.end(), data_clr.begin(), ::tolower);
    data = data_clr;
}

void YMIX::copy_array(
    const std::shared_ptr<const double[]>& source_arr, 
    YCU N, 
    std::shared_ptr<double[]>& res_arr
){
    res_arr = std::shared_ptr<double[]>(new double[N]);
    for(uint32_t ii = 0; ii < N; ii++)
        res_arr[ii] = source_arr[ii];
}