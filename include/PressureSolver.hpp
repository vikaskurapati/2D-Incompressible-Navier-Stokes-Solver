#pragma once

#include "Boundary.hpp"
#include "Discretization.hpp"
#include "Fields.hpp"
#include "Grid.hpp"
#include <utility>
/**
 * @brief Abstract class for pressure Poisson equation solver
 *
 */
class PressureSolver {
  public:
    PressureSolver() = default;
    virtual ~PressureSolver() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) = 0;
};

/**
 * @brief Successive Over-Relaxation algorithm for solution of pressure Poisson
 * equation
 *
 */
class SOR : public PressureSolver {
  public:
    SOR() = default;

    /**
     * @brief Constructor of SOR solver
     *
     * @param[in] relaxation factor
     */
    SOR(double omega);

    virtual ~SOR() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:
    double _omega;
};

/**
 * @brief Jacobi iteration to solve the Pressure Poisson Equation
 *
 */

class JACOBI : public PressureSolver {
  public:
    JACOBI() = default;

    /**
     * @brief Constructor of JACOBI solver
     *
     * @param[in] relaxation factor
     */

    virtual ~JACOBI() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);
};

/**
 * @brief Weighted Jacobi iterations to solve the Pressure Poisson Equation
 *
 */
class WEIGHTED_JACOBI : public PressureSolver {
  public:
    WEIGHTED_JACOBI() = default;

    /**
     *@brief Constructor of WEIGHTED_JACOBI
     *
     *@param[in] relaxation factor
     */
    WEIGHTED_JACOBI(double omega);

    virtual ~WEIGHTED_JACOBI() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */

    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:
    double _omega;
};

/**
 * @brief Gauss Seidel iteration to solve the Pressure Poisson Equation
 *
 */
class GAUSS_SEIDEL : public PressureSolver {
  public:
    GAUSS_SEIDEL() = default;

    virtual ~GAUSS_SEIDEL() = default;

    /**
     * @brief Solver the pressure equation on given field using Gauss Seidel iterations
     *
     * @param field to be used
     * @param grid to be used
     * @param boundaries used
     * @return double the MSE residual value
     */

    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);
};

class RICHARDSON : public PressureSolver {
  public:
    RICHARDSON() = default;

    /**
     * @brief Construct a new RICHARDSON object
     *
     * @param omega
     */
    RICHARDSON(double omega);
    virtual ~RICHARDSON() = default;

    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:
    double _omega{1.0};
};

class ConjugateGradient : public PressureSolver {
  public:
    ConjugateGradient() = default;
    virtual ~ConjugateGradient() = default;

    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:
    int iter = 0;
}

Matrix<double>
residual(Matrix<double> pressure, Matrix<double> rhs, double dx, double dy) {
    Matrix<double> residual = Matrix<double>(p.imax(), p.jmax(), 0.0);

    for (int i = 1; i < p.imax() - 2; ++i) {
        /* code */
        for (int j = 1; j < p.jmax(); ++j) {
            /* code */
            double delta = Discretization::laplacian(pressure, i, j);
            residual(i, j) = rhs(i, j) - delta;
        }
    }

    return residual;
}