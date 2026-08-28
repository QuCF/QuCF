# author: Ivan Novikau
cmake_minimum_required(VERSION 3.17)

# name of the executable file:
set(OUTPUT_EXE "approx_polyn")

# # 1. SET RUNTIME POLICY FIRST
# set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")

# DEBUG or RELEASE:
set(CMAKE_BUILD_TYPE Debug)

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

# create the target:
add_executable(${OUTPUT_EXE} ../src/ApprPolynFourier.cu)

# set C++ compiler version:
set_property(TARGET ${OUTPUT_EXE} PROPERTY CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set CUDA architecture:
set_property(TARGET ${OUTPUT_EXE} PROPERTY CUDA_ARCHITECTURES 75)

#link HDF5 library
target_link_libraries(${OUTPUT_EXE} PRIVATE HDF5::HDF5)

# to avoid warning abour DLL link
target_link_options(${OUTPUT_EXE} PRIVATE 
    $<$<CONFIG:Debug>:/NODEFAULTLIB:LIBCMT>
    $<$<CONFIG:Release>:/NODEFAULTLIB:LIBCMT>
)

# To supress warning about std::string 
target_compile_options(${OUTPUT_EXE} PRIVATE 
    $<$<COMPILE_LANGUAGE:CXX>:/wd4251> # MSVC specific warnings
    $<$<COMPILE_LANGUAGE:CUDA>:-diag-suppress=1394> # NVCC specific warnings
)

# for correct debuging on device:
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # target_compile_options(${OUTPUT_EXE} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-G>)
    target_compile_options(${OUTPUT_EXE} PRIVATE "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CUDA>>:-G;-src-in-ptx>")
endif()