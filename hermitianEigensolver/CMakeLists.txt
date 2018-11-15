cmake_minimum_required(VERSION 3.10)
project(hermitianEigensolver)

include_directories("/usr/local/include/eigen3")

set(CMAKE_CXX_STANDARD 14)

add_executable(hermitianEigensolver main.cpp)
