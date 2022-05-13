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

    du2dx = (1/_dx)*(interpolate(U,i,j,i+1,j)*interpolate(U,i,j,i+1,j) - interpolate(U,i-1,j,i,j)*interpolate(U,i-1,j,i,j));
    du2dx = du2dx + (_gamma/_dx)*(std::abs(interpolate(U,i,j,i+1,j))*(U(i,j)+U(i+1,j)*0.5) - std::abs(interpolate(U,i-1,j,i,j))*(U(i-1,j)+U(i,j))*0.5);

    duvdy = (1/_dy)*(interpolate(V,i,j,i+1,j)*interpolate(U,i,j,i,j+1) - interpolate(V,i,j-1,i+1,j-1)*interpolate(U,i,j-1,i,j));
    duvdy = duvdy + (_gamma/_dy)*(std::abs(interpolate(V,i,j,i+1,j))*(U(i,j)-U(i,j+1)*0.5) - std::abs(interpolate(V,i,j-1,i+1,j-1))*(U(i,j-1)-U(i,j))*0.5);

    return du2dx + duvdy;
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) 
{
    double duvdx = 0.0;
    double dv2dy = 0.0;

    duvdx = (1/_dx)*(interpolate(U,i,j,i,j+1)*interpolate(V,i,j,i+1,j) - interpolate(U,i-1,j,i-1,j+1)*interpolate(V,i-1,j,i,j));
    duvdx = duvdx + (_gamma/_dx)*(std::abs(interpolate(U,i,j,i,j+1))*(V(i,j)-V(i+1,j))*0.5 - std::abs(interpolate(U,i-1,j,i-1,j+1))*(V(i-1,j)-V(i,j))*0.5);

    dv2dy = (1/_dy)*(interpolate(V,i,j,i,j+1)*interpolate(V,i,j,i,j+1) - interpolate(V,i,j-1,i,j)*interpolate(V,i,j-1,i,j));
    dv2dy = dv2dy + (_gamma/_dy)*(std::abs(interpolate(V,i,j,i,j+1))*(V(i,j)-V(i,j+1))*0.5 - std::abs(interpolate(V,i,j-1,i,j))*(V(i,j-1)-V(i,j))*0.5);

    return duvdx + dv2dy;
}

double Discretization::convection_t(const Matrix<double> &T, const Matrix<double> &U, const Matrix<double> &V, int i, int j) 
{
    double dutdx = 0.0;
    double dvtdy = 0.0;

    dutdx = (1/_dx)*(interpolate(T,i,j,i+1,j)*U(i,j) - U(i-1,j)*interpolate(T,i-1,j,i,j));
    dutdx = dutdx + (_gamma/_dx)*(std::abs(U(i,j))*(T(i,j)-T(i+1,j))*0.5 - std::abs(U(i-1,j))*(T(i-1,j)-T(i,j))*0.5);

    dvtdy = (1/_dy)*(interpolate(T,i,j,i,j-1)*V(i,j) - V(i,j-1)*interpolate(T,i,j-1,i,j));
    dvtdy = dvtdy + (_gamma/_dy)*(std::abs(V(i,j))*(T(i,j)-T(i,j+1))*0.5 - std::abs(V(i,j-1))*(T(i,j-1)-T(i,j))*0.5);

    return dutdx + dvtdy;
}

double Discretization::diffusion(const Matrix<double> &A, int i, int j) 
{
    return (A(i + 1, j) - 2.0 * A(i, j) + A(i - 1, j)) / (_dx * _dx) +
                    (A(i, j + 1) - 2.0 * A(i, j) + A(i, j - 1)) / (_dy * _dy);
}

double Discretization::laplacian(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                    (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) 
{
    return 0.5*(A(i,j) + A(i_offset,j_offset));
}
