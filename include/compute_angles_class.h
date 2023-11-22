#ifndef COMPUTEANGLES_H
#define COMPUTEANGLES_H

#include "QLib.h"


class ComputeAngles_{

protected:
    std::string project_name_;
    std::string work_directory_;

    /**
     * Output \p project_name_.hdf5 file in \p work_directory_ directory.
    */
    YMIX::H5File hfo_; 


public:
    ComputeAngles_(YCS pname) :
        project_name_(pname), work_directory_(std::filesystem::current_path())
    {
        create_output_hdf5();


    }


protected:
    void create_output_hdf5()
    {
        using namespace std;

        printf("Creating output %s.hdf5 file in the directory:\n", project_name_.c_str());
        printf("%s\n", work_directory_.c_str());
        string filename_hdf5 = work_directory_ + project_name_ + ".hdf5";

        hfo_.create(filename_hdf5);
        hfo_.add_group("basic");

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

    }


    void compute_cost_function()
    {

    }


    void compute_grad_of_cost_function()
    {

    }


    void compute_target_function()
    {

    }




};
#endif