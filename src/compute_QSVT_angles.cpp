
// #define NO_QUEST

#include "../include/compute_angles_class.h"

using namespace std;


/**
 * REMARKS:
 * To compute QSP angles I can use the same approach but without assuming symmetry of angles,
 *  also I need to change the form of Wx.
 *  What about parity?
*/



/**
 * Launch:
 *      ./compute_angles name_project
 * where
 *  > [name_project].ca - text file wit parameters for computing angles
 * Data are saved in the directory from where the code has been launched
*/
int main(int argc, char *argv[])
{
    uint64_t N_angles; // defined by the number of coefficients in the sequence of polynomials;

    // --- Interpret input parameters ---
    uint32_t id_arg;

    id_arg = 1;
    string pname(argv[id_arg]); // project name
    cout << "Project name: "   << pname << endl;

    try
    {
        ComputeAngles_ oo = ComputeAngles_(pname);


        // // --- read the [name_project].ca file ---

        // // choose the function for which the QSVT or QSP angles should be computed;
        // // choose the parameters of the computation: stopping criteria, file with coefficients for a polynomial
        
        // //
    

        // // set initial angles:
        // shared_ptr<double []> angles = shared_ptr<double[]>(new double[N_angles]);

        // // compute U

        // // compute the cost function (will need for a line search)
        // void compute_cost_function();

        // // compute gradient of the cost function
        // void compute_grad_of_cost_function();


        // // compute the target function
        // void compute_target_function();




    }
    catch(YCS e)
    {
        std::cerr << "\n" << e << endl;
        return -1;
    }
    catch(const std::exception& e)
    {
        std::cerr << "General error:\n" << e.what() << '\n';
        return -1;
    }
    return 0;
}


