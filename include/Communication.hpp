#pragma once


class Communication
{
    public:

    static void init_parallel(int* argn, char** args, int& rank, int& size);
    static void finalize();
    static double reduce_min(double dt);

};


// namespace Communication {
// void init_parallel(int argn, char** args, int& rank, int& size);
// void finalize();
// double reduce_min(double dt);
// }