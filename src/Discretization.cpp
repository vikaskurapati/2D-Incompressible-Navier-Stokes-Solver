#include "Discretization.hpp"

#include <cmath>

double Discretization::_dx = 0.0;
double Discretization::_dy = 0.0;
double Discretization::_gamma = 0.0;

Discretization::Discretization(double dx, double dy, double gamma) {
    _dx = dx;
    _dy = dy;
    _gamma = gamma;
}

double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) 
{
    double du2x = 0.0;
    double duvy = 0.0;
    du2x = (((U(i,j)+U(i+1,j))/2.0)*((U(i,j)+U(i+1,j))/2.0) -((U(i-1,j)+U(i,j))/2.0)*((U(i-1,j)+U(i,j))/2.0))/_dx;
    duvy = (((((U(i,j)+U(i,j+1))/2.0))*((V(i,j)+V(i+1,j))/2.0)) - ((((U(i,j-1)+U(i,j))/2.0))*((V(i,j-1)+V(i+1,j-1))/2.0)))/_dy;
    return du2x+duvy;
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) 
{
    double duvx = 0.0;
    double dv2y = 0.0;

    dv2y = (((V(i,j)+V(i,j+1))/2.0)*((V(i,j)+V(i,j+1))/2.0) -((V(i,j-1)+V(i,j))/2.0)*((V(i,j-1)+V(i,j))/2.0))/_dy;
    duvx = (((((U(i,j)+U(i+1,j))/2.0))*((V(i,j)+V(i,j+1))/2.0)) - ((((U(i-1,j)+U(i,j))/2.0))*((V(i-1,j)+V(i-1,j+1))/2.0)))/_dx;
    return dv2y+duvx;
}

//double Discretization::diffusion(const Matrix<double> &A, int i, int j) {}

double Discretization::laplacian(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                    (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

// double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {}
