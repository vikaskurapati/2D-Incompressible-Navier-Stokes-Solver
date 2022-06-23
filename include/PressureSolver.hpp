#pragma once

#include "Boundary.hpp"
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

virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);
};