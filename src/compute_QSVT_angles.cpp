
// #define NO_QUEST

#include "../include/compute_angles_class.h"

#ifdef _WIN32
    #include <iostream>
    #include <windows.h>
    #include <thread>
    #include <chrono>
#endif

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
 * Data are saved in the same directory from where the code has been launched.
*/
int main(int argc, char *argv[])
{
    // --- Interpret input parameters ---
    uint32_t id_arg;

    id_arg = 1;
    string pname(argv[id_arg]); // project name
    cout << "********************************************************" << endl;
    cout << "********************************************************" << endl;
    cout << "Project name: "   << pname << endl;

    // // --- for debuging !!! ---
    // #ifdef _WIN32
    //     std::cout << "Waiting for a debugger..." << std::endl;
    //     while (!IsDebuggerPresent()) {
    //         std::this_thread::sleep_for(std::chrono::milliseconds(100));
    //     }
    //     std::cout << "Debugger attached! Continuing..." << std::endl;
    // #endif

    try
    {
        // Compute angles for the Eq.(12) in Dong-21 or Eq.(13) in Novikau-23:
        ComputeAngles_ oo = ComputeAngles_(pname);
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


