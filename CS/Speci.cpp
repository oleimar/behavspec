#include "SpeciCode.hpp"
#include <iostream>

// The Speci program runs simulations of learning in a social network
// Copyright (C) 2022  Olof Leimar
// See Readme.md for copyright notice

int main(int argc, char* argv[])
{
    // Open input file and read indata
    SpeciInpData sid(argv[1]);
    if (!sid.OK) {
        std::cout << "Input failed!" << "\n";
        return -1;
    }
    // Run the iteration
    Speci speci(sid);
    speci.Run();
    return 0;
}
