
# behavspec: C++ code for simulation of behavioural specialisation through learning


## Overview

This repository contains C++ code and example data.
The executable programs, each called `Speci` and built from the code variants, will run learning simulations of a number of social groups in which individuals are placed and interact in a social network.
The programs were used to produce the simulation results for the paper "Behavioural specialisation and learning in social networks" by Olof Leimar, Sasha R. X. Dall, Alasdair I. Houston, and John M. McNamara.


## System requirements

The programs have been compiled and run on a Linux server with Ubuntu 20.04 LTS.
The C++ compiler was g++ version 9.4.0, provided by Ubuntu, with compiler flags for c++17, and `cmake` (<https://cmake.org/>) was used to build the programs.
The programs can be run multithreaded using OpenMP, which speeds up execution times, but for the simulations reported in the paper, times are fairly short (a few minutes) also for single-threaded use.
Most likely the instructions below will work for many Linux distributions.
For single-threaded use, the programs have also been compiled and run on macOS, using the Apple supplied Clang version of g++.

The programs read input parameters from TOML files (<https://github.com/toml-lang/toml>), using the open source `cpptoml.h` header file (<https://github.com/skystrife/cpptoml>), included in this repository.

The programs use the Boost graph library to represent small-world networks, which is part of the open source Boost libraries (<https://www.boost.org/>).
The programs also store the simulated groups in HDF5 files (<https://www.hdfgroup.org/>), which is an open source binary file format.
The programs use the open source HighFive library (<https://github.com/BlueBrain/HighFive>) to read and write to such files.
These pieces of software need to be installed in order for `cmake` to successfully build the programs.


## Installation guide

Install the repository from Github to a local computer.
There is a directory `behavspec` containing this file and a license file, and subdirectories `PS`, `CS`, and `HD` with source code for producer-scrounger, caller-satellite, and Hawk-Dove simulations.
In each of these there is a subdirectory `Data` where input data and data files containing simulated groups are kept, and a subdirectory `build` used by `cmake` for files generated during building, including the executable `Speci`.


## Building the programs

The CMake build system is used.
The programs in the different project folders, `PS`, `CS`, and `HD`, are built separately.
If it does not exist, create a build subdirectory in the desired project folder (`mkdir build`) and make it the current directory (`cd build`).
If desired, for a build from scratch, delete any previous content (`rm -rf *`).
Run CMake from the build directory. For a release build:
```
cmake -D CMAKE_BUILD_TYPE=Release ../
```
and for a debug build replace Release with Debug.
If this succeeds, i.e. if the `CMakeLists.txt` file in the project folder is processed without problems, build the program:
```
make
```
This should produce an executable in the `build` directory.


## Running

Make the Data directory current.
Assuming that the executable is called `Speci` and with an input file called `Run_g99_deg04Q1a1prw00.toml`, corresponding to one of the cases in Figure 2, run the program as
```
../build/Speci Run_g99_deg04Q1a1prw00.toml
```

The programs are capable of letting some individual traits evolve over a number of generations, although this was not used in any of the simulations for the paper.
While running the program, a counter of the number of generations is displayed, but it only reaches 1 for the simulations for the paper.

## Description of the simulations

There is an input file, for instance `Run_g99_deg04Q1a1prw00.toml`, for each case, which typically simulates 500 groups of `N = 99` individuals each, optionally inputting the individuals from, e.g., the HDF5 file `Run_g99_deg04Q1a1prw00.h5` and outputting to the same file.
Without an existing `Run_g99_deg04Q1a1prw00.h5` data file, the program can start by constructing individuals with genotypes from the allelic values given by `all0`in the input file.
To make this happen, use `read_from_file = false` in the input file.
This is how the simulations for the paper were run.

### The cases in Figure 2 and Figure 3 of the paper

The input files for the procedure-scrounger cases in Figure 2 are found in the `PS/Data` directory and are as follows:
1. `Run_g99_deg04Q1a1prw00.toml`
2. `Run_g99_deg04Q1a1prw01.toml`
3. `Run_g99_deg16Q1a1prw00.toml`
4. `Run_g99_deg32Q1a1prw00.toml`
5. `Run_g99_deg98Q1a1prw00.toml`
6. `Run_g99_deg04Q1a0prw00.toml`
7. `Run_g99_deg98Q1a0prw00.toml`

The input files for the caller-satellite cases in Figure 3 are found in the `CS/Data` directory and are as follows:
1. `Run_g99_deg02Q1a1prw00.toml`
2. `Run_g99_deg02Q1a1prw01.toml`
3. `Run_g99_deg08Q1a1prw00.toml`
4. `Run_g99_deg98Q1a1prw00.toml`
5. `Run_g99_deg08Q1a1prw01.toml`
6. `Run_g99_deg02Q1a0prw00.toml`
7. `Run_g99_deg98Q1a0prw00.toml`

The input files for the Hawk-Dove cases in Figure 3 are found in the `HD/Data` directory and are as follows:
1. `Run_g99_deg02Q1a1prw00.toml`
2. `Run_g99_deg02Q1a1prw01.toml`
3. `Run_g99_deg08Q1a1prw00.toml`
4. `Run_g99_deg98Q1a1prw00.toml`
5. `Run_g99_deg08Q1a1prw01.toml`
6. `Run_g99_deg02Q1a0prw00.toml`
7. `Run_g99_deg98Q1a0prw00.toml`


## License

The `Speci` program runs simulations of  learning in a social network.

Copyright (C) 2022  Olof Leimar

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

