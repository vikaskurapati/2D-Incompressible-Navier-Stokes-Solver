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
    double du2dx = 0.0;
    double duvdy = 0.0;

    du2dx = (0.25/_dx)*((U(i,j)+U(i+1,j))*(U(i,j)+U(i+1,j)) - (U(i-1,j)+U(i,j))*(U(i-1,j)+U(i,j)));
    du2dx = du2dx + (0.25*_gamma/_dx)*(abs(U(i,j)+U(i+1,j))*(U(i,j)+U(i+1,j)) - abs(U(i-1,j)+U(i,j))*(U(i-1,j)+U(i,j)));

    duvdy = (0.25/_dy)*((V(i,j)+V(i+1,j))*(U(i,j)+U(i,j+1)) - (V(i,j-1)+V(i+1,j-1))*(U(i,j-1)+U(i,j)));
    duvdy = duvdy + (0.25*_gamma/_dy)*(abs(V(i,j)+V(i+1,j))*(U(i,j)-U(i,j+1)) - abs(V(i,j-1)+V(i+1,j-1))*(U(i,j-1)-U(i,j)));

    return du2dx + duvdy;
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) 
{
    double duvdx = 0.0;
    double dv2dy = 0.0;

    duvdx = (0.25/_dx)*((U(i,j)+U(i,j+1))*(V(i,j)+V(i+1,j)) - ((U(i-1,j)+U(i-1,j+1))*(V(i-1,j)+V(i,j))));
    duvdx = duvdx + (0.25*_gamma/_dx)*(abs(U(i,j)+U(i,j+1))*(V(i,j)-V(i+1,j)) - (abs(U(i-1,j)+U(i-1,j+1))*(V(i-1,j)-V(i,j))));

    dv2dy = (0.25/_dy)*((V(i,j)+V(i,j+1))*(V(i,j)+V(i,j+1)) - (V(i,j-1)+V(i,j))*(V(i,j-1)+V(i,j)));
    dv2dy = dv2dy + (0.25*_gamma/_dy)*(abs(V(i,j)+V(i,j+1))*(V(i,j)-V(i,j+1)) - abs(V(i,j-1)+V(i,j))*(V(i,j-1)-V(i,j)));

    return duvdx + dv2dy;
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
