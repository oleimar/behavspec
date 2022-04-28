cmake_minimum_required(VERSION 3.10)

# The Speci program runs simulations of learning in a social network
# Copyright (C) 2022  Olof Leimar
# See Readme.md for copyright notice

# set the project name
project(Speci VERSION 1.0 LANGUAGES CXX C)

# set c++ standard
set(CMAKE_CXX_STANDARD 17)

# we need HighFive to handle hdf5 files
find_package(HighFive REQUIRED)

find_package(OpenMP)

# add the executable
add_executable(Speci Speci.cpp SpeciCode.cpp hdf5code.cpp Utils.cpp)

# link target
target_link_libraries(Speci PUBLIC HighFive OpenMP::OpenMP_CXX)
