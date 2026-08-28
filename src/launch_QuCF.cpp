#include "../include/QuCF.h"

#ifdef _WIN32
    #include <iostream>
    #include <windows.h>
    #include <thread>
    #include <chrono>
#endif

/**
 * @brief To launch the oracletool.
 * @param argv QuCF [project_name] [path_to_input_files] [flag_debug]
 */
int main(int argc, char *argv[])
{
    uint32_t id_arg;
    QuESTEnv env = createQuESTEnv();
    if(argc < 2)
    {
        std::cerr << "Project name is missing" << std::endl;
        exit(-1);
    }
    if(argc < 3)
    {
        std::cerr << "Path to input files is missing" << std::endl;
        exit(-1);
    }


    // --- project name ---
    id_arg = 1;
    std::string pname(argv[id_arg]);
    
    // --- path to input files ---
    id_arg += 1;
    std::string path_input(argv[id_arg]);

    // --- path to input files ---
    std::string debug_flag = "0";
    if(argc >= 4)
    {
        id_arg += 1;
        debug_flag = std::string(argv[id_arg]);
    }
    
    // --- for debuging !!! ---
    if(YMIX::compare_strings(debug_flag, "1"))
    {
        #ifdef _WIN32
            std::cout << "Waiting for a debugger..." << std::endl;
            while (!IsDebuggerPresent()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }
            std::cout << "Debugger attached! Continuing..." << std::endl;
        #endif
    }

    // file of the log file:
    YMIX::LogFile::name_global_ = path_input + "/" + pname + ".clog";
    YMIX::LogFile cf(true);

    std::cout << "\n\n********************************************************************************" << std::endl;
    std::cout << "********************************************************************************" << std::endl;
    std::cout << "Project name: "        << pname << std::endl;
    std::cout << "Path to input files: " << path_input << "/" << std::endl;

    try
    {
        QuCF__ oo = QuCF__(
            env, 
            pname, 
            path_input
        );
        oo.launch();
    }
    catch(YCS e)
    {
        std::cerr << "\n" << e << std::endl;
        YMIX::print_log("\n" + e);
        destroyQuESTEnv(env);
        return -1;
    }
    catch(const std::exception& e)
    {
        std::cerr << "General error:\n" << e.what() << '\n';
        YMIX::print_log("General error:\n" + std::string(e.what()));
        destroyQuESTEnv(env);
        return -1;
    }
    destroyQuESTEnv(env);
    return 0;
}