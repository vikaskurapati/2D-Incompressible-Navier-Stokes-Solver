#pragma once

#include "Datastructures.hpp"
#include "Discretization.hpp"
#include "Grid.hpp"

/**
 * @brief Class of container and modifier for the physical fields
 *
 */
class Fields {
  public:
    Fields() = default;

    /**
     * @brief Constructor for the fields
     *
     * @param[in] kinematic viscosity
     * @param[in] initial timestep size
     * @param[in] adaptive timestep coefficient
     * @param[in] number of cells in x direction
     * @param[in] number of cells in y direction
     * @param[in] initial x-velocity
     * @param[in] initial y-velocity
     * @param[in] initial pressure
     *
     */
    Fields(double _nu, double _dt, double _tau, int imax, int jmax, double UI, double VI, double PI);

    /**
     * @brief Calculates the convective and diffusive fluxes in x and y
     * direction based on explicit discretization of the momentum equations
     *
     * @param[in] grid in which the fluxes are calculated
     * 
     * $F_(i, j) = U_(i,j) + dt*(\nu*(\nabla(U,i,j)) - convection_u + g_x)$
     * $G_(i, j) = V_(i,j) + dt*(\nu*(\nabla(V,i,j)) - convection_v + g_y)$
     * 
     * Where convection_u and convection_v are defined in Discretization
     * 
     * $convection_u = \fraction{\partial u^2}{\partial x} + \fraction{\partial uv}{\partial y}$
     * $convection_v = \fraction{\partial uv}{\partial x} + \fraction{\partial v^2}{\partial y}$
     * 
     * Also, boundary conditions for fluxes were also updated
     * 
     * $F_(0,j) = U_(0,j) for j = 1...j_max$
     * $F_(imax, j) = U_(imax, j) for j = 1...j_max$
     * $G_(i,0) = V_(i,0) for i = 1...i_max$
     * $G_(i, jmax) = V_(i, jmax) for i = 1...i_max$
     * 
     * Refer Equations 9, 10 and 17 in Worksheet
     *
     */
    void calculate_fluxes(Grid &grid);

    /**
     * @brief Right hand side calculations using the fluxes for the pressure
     * Poisson equation
     *
     * @param[in] grid in which the calculations are done
     * 
     * $RS_(i,j) = (((F_(i,j)-F_(i-1,j))/dx)+((G_(i,j)-G_(i,j-1))/dy))/dt$
     *
     */
    void calculate_rs(Grid &grid);

    /**
     * @brief Velocity calculation using pressure values
     *
     * @param[in] grid in which the calculations are done
     * 
     * $U_(i,j) = F_(i,j) - (dt/dx)*(P_(i+1,j)-P_(i,j))$
     * $V_(i,j) = G_(i,j) - (dt/dy)*(P_(i,j+1)-P_(i,j))$
     *
     */
    void calculate_velocities(Grid &grid);

    /**
     * @brief Adaptive step size calculation using x-velocity condition,
     * y-velocity condition and CFL condition
     *
     * @param[in] grid in which the calculations are done
     *
     */
    double calculate_dt(Grid &grid);

    /// x-velocity index based access and modify
    double &u(int i, int j);

    /// y-velocity index based access and modify
    double &v(int i, int j);

    /// pressure index based access and modify
    double &p(int i, int j);

    /// RHS index based access and modify
    double &rs(int i, int j);

    /// x-momentum flux index based access and modify
    double &f(int i, int j);

    /// y-momentum flux index based access and modify
    double &g(int i, int j);

    /// get timestep size
    double dt() const;

    /// pressure matrix access and modify
    Matrix<double> &p_matrix();

  private:
    /// x-velocity matrix
    Matrix<double> _U;
    /// y-velocity matrix
    Matrix<double> _V;
    /// pressure matrix
    Matrix<double> _P;
    /// x-momentum flux matrix
    Matrix<double> _F;
    /// y-momentum flux matrix
    Matrix<double> _G;
    /// right hand side matrix
    Matrix<double> _RS;

    /// kinematic viscosity
    double _nu;
    /// gravitional accelearation in x direction
    double _gx{0.0};
    /// gravitional accelearation in y direction
    double _gy{0.0};
    /// timestep size
    double _dt;
    /// adaptive timestep coefficient
    double _tau;
};