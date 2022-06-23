#pragma once
#include "Datastructures.hpp"
#include "Domain.hpp"

/**
 * @brief Main Class which encapsulates the communication part of the
 * Problem
 *
 */

class Communication {
  public:
    /**
     * @brief Method to initialize MPI for the entire problem
     *
     * @param argn the problem arguments
     * @param args Extra arugments provided by the user
     * @param rank The process rank
     * @param size Total number of processors for the program
     */
    static void init_parallel(int *argn, char **args, int &rank, int &size);
    /**
     * @brief Method to finalize MPI and exit after the program is finished
     *
     */
    static void finalize();
    /**
     * @brief Method to abort the program due to bad user inputs
     *
     */
    static void abort();
    /**
     * @brief MPI method to find minimum value across dt across all processors
     *
     * @param dt dt value across different processors
     * @return double returns the minimum value of the dt across all processors and broadcasts to all processors
     */
    static double reduce_min(double dt);
    /**
     * @brief MPI method to find sum of a value across all processors
     *
     * @param res residual value to be summed across all processors
     * @return double returns the summed value of the residual across all processors and broadcasts to all processors
     */
    static double reduce_sum(double res);
    /**
     * @brief MPI method to communicate field matrixes across the boundary of different processors
     *
     * @param matrix the field matrix whose values need to be communicated
     * @param domain domain details of the processor
     * @param incoming_rank the rank of the processor
     * @param iproc total number of processors in the x-direction
     */
    static void communicate(Matrix<double> &matrix, const Domain &domain, int incoming_rank, int iproc);
};