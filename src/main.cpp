#include <iostream>
#include <string>
#include <mpi.h>

#include "Case.hpp"

int main(int argn, char **args) {
    // taking an extra argument in case the user wants to store the vtk files as a separate group without
    // overwriting older data sets

    int rank, size;
    MPI_Init(&argn, &args);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argn > 2) {
        std::string file_name{args[1]};
        int problem_rank{std::stoi(args[2])};
        Case problem(file_name, argn, args, problem_rank, rank);
        if(size == 1)
        {   
            problem.simulate_serial(problem_rank);
        }
        else
        {
            problem.simulate_serial(problem_rank);
        }
    }
    // in case the user doesn't provide any extra argument apart from the input file path
    else if (argn > 1) {
        std::string file_name{args[1]};
        Case problem(file_name, argn, args);
        if(size == 1){
        problem.simulate_serial();
        }
        else{
            problem.simulate_serial();
        } 
    }else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
    MPI_Finalize();
    return 0;
}
