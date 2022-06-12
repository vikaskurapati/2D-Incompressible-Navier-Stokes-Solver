#include "Communication.hpp"
#include <mpi.h>

void Communication::init_parallel(int* argn, char** args, int& rank, int& size)
{
    MPI_Init(argn , &args);
    MPI_Comm_rank( MPI_COMM_WORLD , &rank);
    MPI_Comm_size( MPI_COMM_WORLD , &size);
}

double Communication::reduce_min(double dt) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double reduced_dt;

    MPI_Reduce(&dt, &reduced_dt, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Bcast(&reduced_dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return reduced_dt;
}

void Communication::finalize()
{
    MPI_Finalize();
}