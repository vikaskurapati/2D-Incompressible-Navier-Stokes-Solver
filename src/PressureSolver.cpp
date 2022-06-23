#include "PressureSolver.hpp"

#include <cmath>
#include <iostream>

double JACOBI::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    double dx = grid.dx();
    double dy = grid.dy();

    double coeff = 1 / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
    int i, j;

    auto p_old = field.p_matrix();
    for (auto currentCell : grid.fluid_cells()) {
        i = currentCell->i();
        j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            field.p(i, j) = coeff * (Discretization::sor_helper(p_old, i, j) - field.rs(i, j));
        }
    }

    double res = 0.0;
    double rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        i = currentCell->i();
        j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
            rloc += (val * val);
        }
    }

    return rloc;
}

SOR::SOR(double omega) : _omega(omega) {}

double SOR::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {

    double dx = grid.dx();
    double dy = grid.dy();

    double coeff = _omega / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = _omega * h^2 / 4.0, if dx == dy == h

    int i, j;

    for (auto currentCell : grid.fluid_cells()) {
        i = currentCell->i();
        j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            field.p(i, j) = (1.0 - _omega) * field.p(i, j) +
                            coeff * (Discretization::sor_helper(field.p_matrix(), i, j) - field.rs(i, j));
        }
    }

    double res = 0.0;
    double rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        i = currentCell->i();
        j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
            rloc += (val * val);
        }
    }

    return rloc;
}

WEIGHTED_JACOBI::WEIGHTED_JACOBI(double omega) : _omega(omega) {}

double WEIGHTED_JACOBI::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {

    double dx = grid.dx();
    double dy = grid.dy();

    double coeff = _omega / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = _omega * h^2 / 4.0, if dx == dy == h

    int i, j;
    auto p_old = field.p_matrix();
    for (auto currentCell : grid.fluid_cells()) {
        i = currentCell->i();
        j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            field.p(i, j) =
                (1.0 - _omega) * p_old(i, j) + coeff * (Discretization::sor_helper(p_old, i, j) - field.rs(i, j));
        }
    }

    double res = 0.0;
    double rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        i = currentCell->i();
        j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
            rloc += (val * val);
        }
    }

    return rloc;
}

double GAUSS_SEIDEL::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    double dx = grid.dx();
    double dy = grid.dy();

    double coeff = 1.0 / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = h^2 / 4.0, if dx == dy == h

    int i, j;

    for (auto currentCell : grid.fluid_cells()) {
        i = currentCell->i();
        j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            field.p(i, j) = coeff * (Discretization::sor_helper(field.p_matrix(), i, j) - field.rs(i, j));
        }
    }

    double res = 0.0;
    double rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        i = currentCell->i();
        j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
            rloc += (val * val);
        }
    }

    return rloc;
}

RICHARDSON::RICHARDSON(double omega) : _omega(omega) {
    if (omega == 0) {
        std::cout << "Invalid Omega for this scheme chosen. Setting omega to be 1.0" << std::endl;
        _omega = 1.0;
    }
}

double RICHARDSON::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    double dx = grid.dx();
    double dy = grid.dy();

    double coeff = 2.0 * (1 / (dx * dx) + 1 / (dy * dy));

    int i, j;

    auto p_old = field.p_matrix();

    for (auto currentCell : grid.fluid_cells()) {
        i = currentCell->i();
        j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            field.p(i, j) =
                p_old(i, j) + _omega * (field.rs(i, j) - Discretization::sor_helper(p_old, i, j) - coeff * p_old(i, j));
        }
    }

    double res = 0.0;
    double rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        i = currentCell->i();
        j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
            rloc += (val * val);
        }
    }

    return rloc;
}