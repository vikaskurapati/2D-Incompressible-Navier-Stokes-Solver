#include "PressureSolver.hpp"

#include <cmath>
#include <iostream>

double Jacobi::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
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

WeightedJacobi::WeightedJacobi(double omega) : _omega(omega) {}

double WeightedJacobi::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {

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

double GaussSeidel::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
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

Richardson::Richardson(double omega) : _omega(omega) {
    if (omega == 0) {
        std::cout << "Invalid Omega for this scheme chosen. Setting omega to be 1.0" << std::endl;
        _omega = 1.0;
    }
}

double Richardson::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    double dx = grid.dx();
    double dy = grid.dy();

    int i, j;

    auto p_old = field.p_matrix();

    for (auto currentCell : grid.fluid_cells()) {
        i = currentCell->i();
        j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            field.p(i, j) = p_old(i, j) + _omega * (field.rs(i, j) - Discretization::laplacian(field.p_matrix(), i, j));
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

ConjugateGradient::ConjugateGradient(Fields &field) {

    int imax = field.p_matrix().imax();
    int jmax = field.p_matrix().jmax();

    d = Matrix<double>(imax, jmax, 0.0);
    residual = Matrix<double>(imax, jmax, 0.0);
}

double ConjugateGradient::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {

    auto pressure = field.p_matrix();
    auto rhs = field.rs_matrix();

    int imax, jmax;

    imax = field.p_matrix().imax();
    jmax = field.p_matrix().jmax();

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            residual(i, j) = rhs(i, j) - Discretization::laplacian(pressure, i, j);
        }
    }

    if (iter == 0) {
        d = residual;
    }

    Matrix<double> q(imax, jmax, 0.0);

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            q(i, j) =
                Discretization::laplacian(d, i, j); // the q = Ad for this system would be laplacian of the d vector
        }
    }

    double alpha_num = 0;
    double alpha_den = 0;
    double alpha = 0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            alpha_num += residual(i, j) * residual(i, j); // num of alpha is delta_new = r^T * r
            alpha_den += d(i, j) * q(i, j);               // den of alpha is d^T * q
        }
    }

    alpha = alpha_num / alpha_den;

    double beta_den = 0;
    double beta_num = 0;
    double beta = 0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {

            field.p(i, j) = pressure(i, j) + alpha * d(i, j);
            beta_den += residual(i, j) * residual(i, j);
            if (iter == 50) {
                iter = 0;
                residual(i, j) = rhs(i, j) - Discretization::laplacian(pressure, i, j);
            } else {
                residual(i, j) -= alpha * q(i, j);
            }
            beta_num += residual(i, j) * residual(i, j);
        }
    }

    beta = beta_num / beta_den;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            d(i, j) = residual(i, j) + beta * d(i, j);
        }
    }

    iter++;

    double res = 0.0;
    double rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();
        if (i != 0 && j != 0 && i != grid.domain().size_x + 1 && j != grid.domain().size_y + 1) {
            double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
            rloc += (val * val);
        }
    }

    return rloc;
}
