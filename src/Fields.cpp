#include "Fields.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

Fields::Fields(Grid &grid, double nu, double dt, double tau, double alpha, double beta, std::string energy_eq, int imax,
               int jmax, double UI, double VI, double PI, double TI, double gx, double gy, int process_rank, int size)
    : _nu(nu), _dt(dt), _tau(tau), _alpha(alpha), _beta(beta), _energy_eq(energy_eq), _gx(gx), _gy(gy) {

    
    _process_rank = process_rank;
    _size = size;

    _U = Matrix<double>(imax + 2, jmax + 2);
    _V = Matrix<double>(imax + 2, jmax + 2);
    _P = Matrix<double>(imax + 2, jmax + 2);
    _T = Matrix<double>(imax + 2, jmax + 2);
    _T_new = Matrix<double>(imax + 2, jmax + 2);

    int i, j;

    for (const auto &cell : grid.fluid_cells()) {
        i = cell->i();
        j = cell->j();
        _U(i, j) = UI;
        _V(i, j) = VI;
        _P(i, j) = PI;
        _T(i, j) = TI;
        _T_new(i, j) = TI;
    }
    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

void Fields::calculate_fluxes(Grid &grid) {
    int i, j;

    if (_energy_eq == "off") {
        for (const auto &currentCell : grid.fluid_cells()) {
            i = currentCell->i();
            j = currentCell->j();

            if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
                _F(i, j) = _U(i, j) + _dt * (_nu * (Discretization::laplacian(_U, i, j)) -
                                             Discretization::convection_u(_U, _V, i, j) + _gx);
                _G(i, j) = _V(i, j) + _dt * (_nu * (Discretization::laplacian(_V, i, j)) -
                                             Discretization::convection_v(_U, _V, i, j) + _gy);
            }
        }
    } else if (_energy_eq == "on") {
        for (const auto &currentCell : grid.fluid_cells()) {
            i = currentCell->i();
            j = currentCell->j();

            if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
                _F(i, j) = _U(i, j) + _dt * (_nu * (Discretization::laplacian(_U, i, j)) -
                                             Discretization::convection_u(_U, _V, i, j));
                _F(i, j) -= 0.5 * _dt * _beta * (_T(i, j) + _T(i + 1, j)) * _gx;
                _G(i, j) = _V(i, j) + _dt * (_nu * (Discretization::laplacian(_V, i, j)) -
                                             Discretization::convection_v(_U, _V, i, j));
                _G(i, j) -= 0.5 * _dt * _beta * (_T(i, j) + _T(i, j + 1)) * _gy;
            }
        }
    } else {
        std::cout << "Something went wrong with energy equation on and off\nPlease check\n";
        MPI_Finalize();
        exit(0);
    }

    for (const auto &boundary :
         {grid.fixed_wall_cells(), grid.hot_fixed_wall_cells(), grid.cold_fixed_wall_cells(),
          grid.adiabatic_fixed_wall_cells(), grid.inflow_cells(), grid.outflow_cells(), grid.moving_wall_cells()}) {
        for (const auto &currentCell : boundary) {
            i = currentCell->i();
            j = currentCell->j();
            if (currentCell->is_border(border_position::TOP)) {
                if (currentCell->is_border(border_position::RIGHT)) {
                    _F(i, j) = 0.0;
                    _G(i, j) = 0.0;
                } else if (currentCell->is_border(border_position::LEFT)) {
                    _F(i - 1, j) = 0.0;
                    _G(i, j) = 0.0;
                } else if (currentCell->is_border(border_position::BOTTOM)) {
                    _G(i, j) = 0.0;
                    _G(i, j - 1) = 0.0;
                } else {
                    _G(i, j) = _V(i, j);
                }
            } else if (currentCell->is_border(border_position::BOTTOM)) {
                if (currentCell->is_border(border_position::RIGHT)) {
                    _F(i, j) = 0.0;
                    _G(i, j - 1) = 0.0;
                } else if (currentCell->is_border(border_position::LEFT)) {
                    _F(i - 1, j) = 0.0;
                    _G(i, j - 1) = 0.0;
                } else {
                    _G(i, j - 1) = _V(i, j - 1);
                }
            } else if (currentCell->is_border(border_position::RIGHT)) {
                if (currentCell->is_border(border_position::LEFT)) {
                    _F(i, j) = 0.0;
                    _F(i - 1, j) = 0.0;
                } else {
                    _F(i, j) = _U(i, j);
                }
            } else if (currentCell->is_border(border_position::LEFT)) {
                _F(i - 1, j) = _U(i - 1, j);
            }
        }
    }
}

void Fields::calculate_rs(Grid &grid) {
    double dx = grid.dx();
    double dy = grid.dy();
    int i, j;
    for (const auto &currentCell : grid.fluid_cells()) {
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            i = currentCell->i();
            j = currentCell->j();
            _RS(i, j) = (((_F(i, j) - _F(i - 1, j)) / dx) + ((_G(i, j) - _G(i, j - 1)) / dy)) / _dt;
        }
    }
}

void Fields::calculate_velocities(Grid &grid) {
    double dx = grid.dx();
    double dy = grid.dy();
    int i, j;
    for (const auto &currentCell : grid.fluid_cells()) {
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {

            i = currentCell->i();
            j = currentCell->j();
            _U(i, j) = _F(i, j) - (_dt / dx) * (_P(i + 1, j) - _P(i, j));
            _V(i, j) = _G(i, j) - (_dt / dy) * (_P(i, j + 1) - _P(i, j));
        }
    }
}
double Fields::calculate_dt(Grid &grid) {
    double dt1, dt2, dt3, dt4;
    double dx = grid.dx();
    double dy = grid.dy();
    double umax = 0.001;
    double vmax = 0.001;
    int i, j;

    for (const auto &currentCell : grid.fluid_cells()) {
        i = currentCell->i();
        j = currentCell->j();
        if (umax < std::abs(_U(i, j))) {
            umax = std::abs(_U(i, j));
        }
        if (vmax < std::abs(_V(i, j))) {
            vmax = std::abs(_V(i, j));
        }
    }

    dt1 = 0.5 * (dx * dx * dy * dy) / ((dx * dx + dy * dy) * _nu);
    dt2 = dx / umax;
    dt3 = dy / vmax;
    if (_energy_eq == "on") {
        dt4 = 0.5 * (dx * dx * dy * dy) / ((dx * dx + dy * dy) * _alpha);
        _dt = _tau * std::min({dt1, dt2, dt3, dt4});
    } else {
        _dt = _tau * std::min({dt1, dt2, dt3});
    }
    return _dt;
}

void Fields::calculate_temperatures(Grid &grid) {
    double dx = grid.dx();
    double dy = grid.dy();
    int i, j;
    for (const auto &currentCell : grid.fluid_cells()) {
        i = currentCell->i();
        j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            _T_new(i, j) = _T(i, j) + _dt * (_alpha * Discretization::laplacian(_T, i, j) -
                                             Discretization::convection_t(_T, _U, _V, i, j));
        }
    }
    _T = _T_new;
}

double &Fields::t(int i, int j) { return _T(i, j); }
double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }
