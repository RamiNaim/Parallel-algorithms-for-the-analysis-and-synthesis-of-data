cmake_minimum_required(VERSION 3.17)
project(lab2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_COMPILER "/Users/rami/opt/usr/local/bin/mpicc")
set(CMAKE_CXX_COMPILER "/Users/rami/opt/usr/local/bin/mpic++")

include_directories(/Users/rami/opt/usr/local/include)
add_executable(lab2 main.cpp include/utils.h)
add_executable(mpi_algo1 mpi_algo1.cpp)
add_executable(mpi_algo2 mpi_algo2.cpp)