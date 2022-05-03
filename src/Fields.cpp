#include "Fields.hpp"

#include <algorithm>
#include <iostream>
#include <cmath>

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI)
    : _nu(nu), _dt(dt), _tau(tau) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

void Fields::calculate_fluxes(Grid &grid)
{
    double d2ux = 0.0;
    double d2uy = 0.0;
    double du2x = 0.0;
    double duvy = 0.0;
    double gx = 0.0;
    double gy = 0.0;
    double d2vx = 0.0;
    double d2vy = 0.0;
    double duvx = 0.0;
    double dv2y = 0.0;

    double dx = grid.dx();
    double dy = grid.dy();

    int i, j;

    for (auto currentCell : grid.fluid_cells()) 
    {
        i = currentCell->i();
        j = currentCell->j();
        d2ux = (_U(i+1,j)-2*_U(i,j)+_U(i-1,j))/(dx*dx);
        d2uy = (_U(i,j+1)-2*_U(i,j)+_U(i,j-1))/(dy*dy);
        du2x = (((_U(i,j)+_U(i+1,j))/2.0)*((_U(i,j)+_U(i+1,j))/2.0) -((_U(i-1,j)+_U(i,j))/2.0)*((_U(i-1,j)+_U(i,j))/2.0))/dx;
        duvy = (((((_U(i,j)+_U(i,j+1))/2.0))*((_V(i,j)+_V(i+1,j))/2.0)) - ((((_U(i,j-1)+_U(i,j))/2.0))*((_V(i,j-1)+_V(i+1,j-1))/2.0)))/dy;
        d2vx = (_V(i+1,j)-2*_V(i,j)+_V(i-1,j))/(dx*dx);
        d2vy = (_V(i,j+1)-2*_V(i,j)+_V(i,j-1))/(dy*dy);
        dv2y = (((_V(i,j)+_V(i,j+1))/2.0)*((_V(i,j)+_V(i,j+1))/2.0) -((_V(i,j-1)+_V(i,j))/2.0)*((_V(i,j-1)+_V(i,j))/2.0))/dy;
        duvx = (((((_U(i,j)+_U(i+1,j))/2.0))*((_V(i,j)+_V(i,j+1))/2.0)) - ((((_U(i-1,j)+_U(i,j))/2.0))*((_V(i-1,j)+_V(i-1,j+1))/2.0)))/dx;

        _F(i, j) = _U(i,j) + _dt*(_nu*(d2ux+d2uy) - du2x - duvy + gx);
        _G(i, j) = _V(i,j) + _dt*(_nu*(d2vx+d2vy) - duvx - dv2y + gy);
    }
    
}

void Fields::calculate_rs(Grid &grid) 
{
    double dx = grid.dx();
    double dy = grid.dy();
    int i, j;
    for (auto currentCell : grid.fluid_cells()) 
    {
        i = currentCell->i();
        j = currentCell->j();
        _RS(i,j) = (((_F(i,j)-_F(i-1,j))/dx)+((_G(i,j)-_G(i,j-1))/dy))/_dt;
    }
}

void Fields::calculate_velocities(Grid &grid) 
{
    double dx = grid.dx();
    double dy = grid.dy();
    int i, j;
    for (auto currentCell : grid.fluid_cells()) 
    {
        i = currentCell->i();
        j = currentCell->j();
        _U(i,j) = _F(i,j) - (_dt/dx)*(_P(i+1,j)-_P(i,j));
        _V(i,j) = _G(i,j) - (_dt/dy)*(_P(i,j+1)-_P(i,j));
    }

}

double Fields::calculate_dt(Grid &grid) { 
    double dt1, dt2, dt3;
    double dx = grid.dx();
    double dy = grid.dy();
    double umax = 0.001;
    double vmax = 0.001;
    int i, j;

    for (auto currentCell: grid.fluid_cells())
    {
        i = currentCell->i();
        j = currentCell->j();
        if (umax < std::abs(_U(i,j)))
        {
            umax = std::abs(_U(i,j));
        }
        if (vmax < std::abs(_V(i,j)))
        {
            vmax = std::abs(_V(i,j));
        }
    }

    dt1 = 0.5*(dx*dx*dy*dy)/((dx*dx + dy*dy)*_nu);
    dt2 = dx/umax;
    dt3 = dy/vmax;
    _dt = _tau*std::min({dt1, dt2, dt3});
    return _dt; 
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }