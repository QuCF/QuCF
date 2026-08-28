# author: Ivan Novikau
cmake_minimum_required(VERSION 3.17)

# name of the executable file:
set(OUTPUT_EXE "compute_angles")

# DEBUG or RELEASE:
set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Release)

# Set the project:
project(APP LANGUAGES CXX VERSION 0.0)


# output some information:
message(STATUS "CMake version: ${CMAKE_VERSION}")
message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")
# message(STATUS "CUDA Toolkit include: ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
# message(STATUS "CUDA COMPILER ID: ${CMAKE_CUDA_COMPILER_ID}")

# set directories, where the executable file will be stored: 
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "./")
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   "./")

# --- CUDA ---
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
else()
    message(STATUS "CUDA compiler not found.")
endif()

# --- HDF5 ---
find_package(HDF5 REQUIRED  COMPONENTS CXX HL)

if(TARGET HDF5::HDF5)
    get_target_property(HDF5_TYPE HDF5::HDF5 TYPE)
    message(STATUS "HDF5 Target Type: ${HDF5_TYPE}")
    
    # Try to extract the imported location paths
    get_target_property(HDF5_LOC HDF5::HDF5 IMPORTED_LOCATION_DEBUG)
    if(NOT HDF5_LOC)
        get_target_property(HDF5_LOC HDF5::HDF5 IMPORTED_LOCATION)
    endif()
    message(STATUS "HDF5 Main Library Location: ${HDF5_LOC}")
endif()


# create the target:
add_executable(${OUTPUT_EXE} ../src/compute_QSVT_angles.cpp ../src/QLib.cpp)

# set C++ compiler version:
set_property(TARGET ${OUTPUT_EXE} PROPERTY CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set CUDA architecture:
set_property(TARGET ${OUTPUT_EXE} PROPERTY CUDA_ARCHITECTURES 75)

#link HDF5 library
target_link_libraries(${OUTPUT_EXE} PRIVATE HDF5::HDF5)
target_compile_definitions(${OUTPUT_EXE} PRIVATE H5_BUILT_AS_DYNAMIC_LIB)

# To supress warning about std::string 
target_compile_options(${OUTPUT_EXE} PRIVATE 
    # MSVC specific warnings
    $<$<COMPILE_LANGUAGE:CXX>:/wd4251>
    
    # NVCC specific warnings
    $<$<COMPILE_LANGUAGE:CUDA>:-diag-suppress=1394>
)

# for correct debuging on device:
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # target_compile_options(${OUTPUT_EXE} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-G>)
    target_compile_options(${OUTPUT_EXE} PRIVATE "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CUDA>>:-G;-src-in-ptx>")
endif()