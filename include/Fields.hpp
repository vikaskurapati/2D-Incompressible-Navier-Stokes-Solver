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
    Fields(Grid &grid, double _nu, double _dt, double _tau, double _alpha, double _beta, std::string energy_eq,
           int imax, int jmax, double UI, double VI, double PI, double TI, double gx, double gy);

    /**
     * @brief Calculates the convective and diffusive fluxes in x and y
     * direction based on explicit discretization of the momentum equations
     *
     * @param[in] grid in which the fluxes are calculated
     *
     * \f$F_(i, j) = U_(i,j) + dt*(\nu*(\nabla(U,i,j)) - convection_u + g_x)\f$
     * \f$G_(i, j) = V_(i,j) + dt*(\nu*(\nabla(V,i,j)) - convection_v + g_y)\f$
     *
     * Where convection_u and convection_v are defined in Discretization
     *
     * $convection_u = \fraction{\partial u^2}{\partial x} + \fraction{\partial uv}{\partial y}$
     * $convection_v = \fraction{\partial uv}{\partial x} + \fraction{\partial v^2}{\partial y}$
     *
     * Also, boundary conditions for fluxes were also updated
     *
     * \f$F_(0,j) = U_(0,j) for j = 1...j_max\f$
     * \f$F_(imax, j) = U_(imax, j) for j = 1...j_max\f$
     * \f$G_(i,0) = V_(i,0) for i = 1...i_max\f$
     * \f$G_(i, jmax) = V_(i, jmax) for i = 1...i_max\f$
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
     * \f$RS_(i,j) = (((F_(i,j)-F_(i-1,j))/dx)+((G_(i,j)-G_(i,j-1))/dy))/dt\f$
     *
     */
    void calculate_rs(Grid &grid);

    /**
     * @brief Velocity calculation using pressure values
     *
     * @param[in] grid in which the calculations are done
     *
     * \f$U_(i,j) = F_(i,j) - (dt/dx)*(P_(i+1,j)-P_(i,j))\f$
     * \f$V_(i,j) = G_(i,j) - (dt/dy)*(P_(i,j+1)-P_(i,j))\f$
     *
     */
    void calculate_velocities(Grid &grid);

    /**
     * @brief Temperature calculation using previous and velocities values
     *
     * @param[in] grid in which the calculations are done
     *
     * Need to Write
     *
     */
    void calculate_temperatures(Grid &grid);

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

    /// Temperature index based access and modify
    double &t(int i, int j);

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

    /// getting energy equation status on or off
    bool get_energy_eq() { return _energy_eq == "on"; };

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
    /// Temperature matrix
    Matrix<double> _T;
    /// Newly calculated Temperature matrix
    Matrix<double> _T_new;

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
    /// thermal diffusivity
    double _alpha;
    /// coefficient of thermal expansion
    double _beta;
    // status if the energy_equations should be turned on or not
    std::string _energy_eq;
};