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

class Jacobi : public PressureSolver {
  public:
    Jacobi() = default;

    /**
     * @brief Constructor of Jacobi solver
     *
     * @param[in] relaxation factor
     */

    virtual ~Jacobi() = default;

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
class WeightedJacobi : public PressureSolver {
  public:
    WeightedJacobi() = default;

    /**
     *@brief Constructor of WeightedJacobi
     *
     *@param[in] relaxation factor
     */
    WeightedJacobi(double omega);

    virtual ~WeightedJacobi() = default;

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
class GaussSeidel : public PressureSolver {
  public:
    GaussSeidel() = default;

    virtual ~GaussSeidel() = default;

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

class Richardson : public PressureSolver {
  public:
    Richardson() = default;

    /**
     * @brief Construct a new Richardson object
     *
     * @param omega
     */
    Richardson(double omega);
    virtual ~Richardson() = default;

    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:
    double _omega{1.0};
};

class ConjugateGradient : public PressureSolver {
  public:
    ConjugateGradient() = default;
    /**
     * @brief Construct a new Conjugate Gradient object
     *
     * @param field to be used for calculations
     */
    ConjugateGradient(Fields &field);
    virtual ~ConjugateGradient() = default;

    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:
    int iter = 0;
    Matrix<double> d;
    Matrix<double> residual;
};
