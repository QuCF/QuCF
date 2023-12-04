
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
 * Launch from a work directory where [anme_project].ca file is placed:
 *      ./compute_angles name_project
 * where
 *  > [name_project].ca - text file with parameters for computing angles
 * Data are saved in the directory from where the code has been launched.
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
        // Compute angles for the Eq.(12) in Dong-21 or Eq.(13) in Novikau-23:
        ComputeAngles_ oo = ComputeAngles_(pname);

        // L-BFGS
        // test odd and even functions;
        // for testing use the same method as in matlab
        //      and test in a circuit using a simple matrix;
        // test for large kappa;
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


