#pragma once
#include "Datastructures.hpp"
#include "Domain.hpp"

class Communication
{
    public:

    static void init_parallel(int* argn, char** args, int& rank, int& size);
    static void finalize();
    static double reduce_min(double dt);
    static double reduce_sum(double res);
    static void communicate(Matrix<double>& matrix, const Domain& domain, int incoming_rank, int iproc);
};