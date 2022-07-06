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
    /**
     * @brief Solver the pressure equation on given field using Richardson iterations
     *
     * @param field to be used
     * @param grid to be used
     * @param boundaries used
     * @return double the MSE residual value
     */
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
    /**
     * @brief Solver the pressure equation on given field using Conjugate Gradient iterations
     *
     * @param field to be used
     * @param grid to be used
     * @param boundaries used
     * @return double the MSE residual value
     */

    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:
    int iter = 0;
    Matrix<double> d;
    Matrix<double> residual;
};

class MultiGridVCycle : public PressureSolver {
  public:
    MultiGridVCycle() = default;

    /**
     * @brief Construct a new Multi Grid V Cycle object
     *
     * @param smoothing_pre_recur number of smoothing iterations at the beginning of the solver
     * @param smoothing_post_recur number of smoothing iteratoins after the multigrid step
     */

    MultiGridVCycle(int smoothing_pre_recur, int smoothing_post_recur);

    virtual ~MultiGridVCycle() = default;
    /**
     * @brief Solver the pressure equation on given field using Multigrid V cycle
     *
     * @param field to be used
     * @param grid to be used
     * @param boundaries used
     * @return double the MSE residual value
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:
    int _smoothing_pre_recur, _smoothing_post_recur, _max_multi_grid_level;

    /**
     * @brief Recursive call to solve the multigrid problem
     *
     * @param field to be used
     * @param p pressure matrix to be prolongated or restricted
     * @param rs right hand side to be prolongated or restricted
     * @param current_level which level the multigrid calculations are
     * @param dx x-stepsize of the problem
     * @param dy y-stepsize of the problem
     * @return Matrix<double> p matrix for the next level
     */
    Matrix<double> recursiveMultiGridVCycle(Fields &field, Matrix<double> p, Matrix<double> rs, int current_level,
                                            double dx, double dy);
    /**
     * @brief Smoother function for the Multi grid scheme using Jacobi iterations (Ref Sci comp 2 for more details on
     * why Jacobi)
     *
     * @param error the error matrix on which smoothing is performed
     * @param residual residual matrix used for smoothing
     * @param iter number of smoothing iterations
     * @param dx x-stepsize of the problem
     * @param dy y-stepsize of the problem
     * @return Matrix<double> smoothed error matrix
     */
    Matrix<double> smoother(Matrix<double> error, Matrix<double> residual, int iter, double dx, double dy);
    /**
     * @brief Method to calculate the residual based on the laplacian operator
     *
     * @param p vector we're solving for
     * @param rs right hand side of the equation
     * @param dx x-stepsize for the problem
     * @param dy y-stepsize for the problem
     * @return Matrix<double> residual in the current stage
     */
    Matrix<double> residual(Matrix<double> p, Matrix<double> rs, double dx, double dy);
    /**
     * @brief Restrictor function to calculate the error on the coarser grid
     *
     * @param fine finer error values to coarsen
     * @return Matrix<double> coarse values of the error on the coarse grid
     */
    Matrix<double> restrictor(Matrix<double> fine);
    /**
     * @brief Prolongator function to calculate the error on the finer grid
     *
     * @param coarse the coarser error values to make them fine
     * @return Matrix<double> Fine values of the error on the fine grid
     */
    Matrix<double> prolongator(Matrix<double> coarse);
};
