#include "Communication.hpp"
#include <mpi.h>

void Communication::init_parallel(int *argn, char **args, int &rank, int &size) {
    MPI_Init(argn, &args);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
}

double Communication::reduce_min(double dt) {
    double reduced_dt;

    MPI_Reduce(&dt, &reduced_dt, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Bcast(&reduced_dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return reduced_dt;
}

void Communication::finalize() { MPI_Finalize(); }

double Communication::reduce_sum(double res) {
    double reduced_res;

    // replace with weighted average later

    MPI_Reduce(&res, &reduced_res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&reduced_res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return reduced_res;
}

void Communication::communicate(Matrix<double> &matrix, const Domain &domain, int incoming_rank, int iproc) {

    std::vector<double> sender;

    //send left to right, receive right to left

    if (domain.neighbour_ranks[0] != -1) {

        sender = matrix.get_col(domain.size_x); // skipping one level to get the first fluid level without buffer
        std::vector<double> receiver(sender.size());
        MPI_Send(&sender , sender.size() , MPI_DOUBLE , incoming_rank + 1 , 1000 , MPI_COMM_WORLD);
        MPI_Status status;

        MPI_Recv( &receiver , sender.size() , MPI_DOUBLE , incoming_rank + 1 , 1001 , MPI_COMM_WORLD , &status);
        // MPI_Sendrecv(&sender, sender.size(), MPI_DOUBLE, domain.neighbour_ranks[0], 1000, &receiver, sender.size(),
        //              MPI_DOUBLE, incoming_rank, 1000, MPI_COMM_WORLD, &status);

        matrix.set_col(receiver, domain.size_x + 1);
    }
    //receive left to right, send right to left

    if (domain.neighbour_ranks[1] != -1) {

        sender = matrix.get_col(0); // skipping one level to get the first fluid level without buffer
        std::vector<double> receiver(sender.size());
        MPI_Status status;

        MPI_Recv(&receiver, sender.size(), MPI_DOUBLE, incoming_rank-1, 1000, MPI_COMM_WORLD, &status);
        matrix.set_col(receiver, 0);

        MPI_Send(&sender, sender.size(), MPI_DOUBLE, incoming_rank - 1, 1001, MPI_COMM_WORLD);
        // MPI_Sendrecv(&sender, sender.size(), MPI_DOUBLE, domain.neighbour_ranks[1], 1001, &receiver, sender.size(),
        //              MPI_DOUBLE, incoming_rank, 1001, MPI_COMM_WORLD, &status);
    }

    if (domain.neighbour_ranks[2] != -1) {

        sender = matrix.get_row(1);
        std::vector<double> receiver(sender.size());

        MPI_Send(&sender, sender.size(), MPI_DOUBLE, incoming_rank - iproc, 1002, MPI_COMM_WORLD);
        MPI_Status status;
        MPI_Recv(&receiver, sender.size(), MPI_DOUBLE, incoming_rank - iproc, 1003, MPI_COMM_WORLD, &status);
        // MPI_Sendrecv(&sender, sender.size(), MPI_DOUBLE, domain.neighbour_ranks[2], 1002, &receiver, sender.size(),
        //              MPI_DOUBLE, incoming_rank, 1002, MPI_COMM_WORLD, &status);

        matrix.set_row(receiver, 0);
    }

    if (domain.neighbour_ranks[3] != -1) {

        sender = matrix.get_row(domain.size_y);
        std::vector<double> receiver(sender.size());
        MPI_Status status;

        MPI_Recv(&receiver, sender.size(), MPI_DOUBLE, incoming_rank + iproc, 1002, MPI_COMM_WORLD, &status);
        matrix.set_row(receiver, domain.size_y + 1);
        MPI_Send(&sender, sender.size(), MPI_DOUBLE, incoming_rank + iproc, 1003, MPI_COMM_WORLD);
        // MPI_Sendrecv(&sender, sender.size(), MPI_DOUBLE, domain.neighbour_ranks[1], 1003, &receiver, sender.size(),
        //              MPI_DOUBLE, incoming_rank, 1003, MPI_COMM_WORLD, &status);
    }
}