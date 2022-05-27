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

double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double du2dx = 0.0;
    double duvdy = 0.0;
    double result;

    du2dx = (U(i, j) + U(i + 1, j)) * (U(i, j) + U(i + 1, j));
    du2dx -= (U(i - 1, j) + U(i, j)) * (U(i - 1, j) + U(i, j));

    du2dx += _gamma * std::abs(U(i, j) + U(i + 1, j)) * (U(i, j) - U(i + 1, j));
    du2dx -= _gamma * std::abs(U(i - 1, j) + U(i, j)) * ((U(i - 1, j) - U(i, j)));

    du2dx = (0.25 / _dx) * du2dx;

    duvdy = (V(i, j) + V(i + 1, j)) * (U(i, j) + U(i, j + 1));
    duvdy -= (V(i, j - 1) + V(i + 1, j - 1)) * (U(i, j - 1) + U(i, j));

    duvdy += _gamma * std::abs(V(i, j) + V(i + 1, j)) * (U(i, j) - U(i, j + 1));
    duvdy -= _gamma * std::abs(V(i, j - 1) + V(i + 1, j - 1)) * (U(i, j - 1) - U(i, j));

    duvdy = (0.25 / _dy) * duvdy;

    result = du2dx + duvdy;

    return result;
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double duvdx = 0.0;
    double dv2dy = 0.0;
    double result;

    duvdx = (U(i, j) + U(i, j + 1)) * (V(i, j) + V(i + 1, j));
    duvdx -= (U(i - 1, j) + U(i - 1, j + 1)) * (V(i - 1, j) + V(i, j));

    duvdx += _gamma * std::abs(U(i, j) + U(i, j + 1)) * (V(i, j) - V(i + 1, j));
    duvdx -= _gamma * std::abs(U(i - 1, j) + U(i - 1, j + 1)) * (V(i - 1, j) - V(i, j));

    duvdx = (0.25 / _dx) * duvdx;

    dv2dy = (V(i, j) + V(i, j + 1)) * (V(i, j) + V(i, j + 1));
    dv2dy -= (V(i, j - 1) + V(i, j)) * (V(i, j - 1) + V(i, j));

    dv2dy += _gamma * std::abs(V(i, j) + V(i, j + 1)) * (V(i, j) - V(i, j + 1));
    dv2dy -= _gamma * std::abs(V(i, j - 1) + V(i, j)) * (V(i, j - 1) - V(i, j));

    dv2dy = (0.25 / _dy) * dv2dy;

    result = duvdx + dv2dy;
    return result;
}

double Discretization::convection_t(const Matrix<double> &T, const Matrix<double> &U, const Matrix<double> &V, int i,
                                    int j) {
    double dutdx = 0.0;
    double dvtdy = 0.0;
    double result;

    dutdx = (T(i, j) + T(i + 1, j)) * U(i, j);
    dutdx -= (T(i - 1, j) + T(i, j)) * U(i - 1, j);

    dutdx += _gamma * std::abs(U(i, j)) * (T(i, j) - T(i + 1, j));
    dutdx -= _gamma * std::abs(U(i - 1, j)) * (T(i - 1, j) - T(i, j));

    dutdx = (0.5 / _dx) * dutdx;

    dvtdy = (T(i, j) + T(i, j + 1)) * V(i, j);
    dvtdy -= (T(i, j - 1) + T(i, j)) * V(i, j - 1);

    dvtdy += _gamma * std::abs(V(i, j)) * (T(i, j) - T(i, j + 1));
    dvtdy -= _gamma * std::abs(V(i, j - 1)) * (T(i, j - 1) - T(i, j));

    dvtdy = (0.5 / _dy) * dvtdy;

    result = dutdx + dvtdy;

    return result;
}

double Discretization::diffusion(const Matrix<double> &A, int i, int j) {
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

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {
    return 0.5 * (A(i, j) + A(i_offset, j_offset));
}
