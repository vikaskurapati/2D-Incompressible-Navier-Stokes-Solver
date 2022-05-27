#include <iostream>
#include <string>

#include "Case.hpp"

int main(int argn, char **args) {
    // taking an extra argument in case the user wants to store the vtk files as a separate group without
    // overwriting older data sets
    if (argn > 2) {
        std::string file_name{args[1]};
        int rank{std::stoi(args[2])};
        Case problem(file_name, argn, args);
        problem.simulate(rank);
    }
    // in case the user doesn't provide any extra argument apart from the input file path
    else if (argn > 1) {
        std::string file_name{args[1]};
        Case problem(file_name, argn, args);
        problem.simulate();
    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
}
