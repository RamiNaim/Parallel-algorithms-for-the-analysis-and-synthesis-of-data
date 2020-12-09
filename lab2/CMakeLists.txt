cmake_minimum_required(VERSION 3.17)
project(lab2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_COMPILER "/usr/bin/mpicc")
set(CMAKE_CXX_COMPILER "/usr/bin/mpic++")

include_directories(/usr/include/mpi/)
add_executable(lab2 main.cpp include/utils.h mpi_algo1.cpp)