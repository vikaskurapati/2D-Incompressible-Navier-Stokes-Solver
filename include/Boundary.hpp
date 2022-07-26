#pragma once

#include <vector>

#include "Cell.hpp"
#include "Fields.hpp"
/**
 * @brief Abstact of boundary conditions.
 *
 * This class patches the physical values to the given field.
 */
class Boundary {
  public:
    /**
     * @brief Main method to patch the boundary conditons to given field and
     * grid
     *
     * @param[in] Field to be applied
     */
    virtual void apply(Fields &field) = 0;
    virtual ~Boundary() = default;
    virtual void apply_pressures(Fields &field) = 0;
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure with zero gradient
 */
class FixedWallBoundary : public Boundary {
  public:
    FixedWallBoundary(std::vector<Cell *> cells);
    FixedWallBoundary(std::vector<Cell *> cells, double wall_temperature);
    virtual ~FixedWallBoundary() = default;
    virtual void apply(Fields &field);
    virtual void apply_pressures(Fields &field);

  private:
    std::vector<Cell *> _cells;
    double _wall_temperature;
};

/**
 * @brief Adiabatic wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure with zero gradient,
 * Neumann for temperatures with zero gradient
 */
class AdiabaticWallBoundary : public Boundary {
  public:
    AdiabaticWallBoundary(std::vector<Cell *> cells);
    virtual ~AdiabaticWallBoundary() = default;
    virtual void apply(Fields &field);
    virtual void apply_pressures(Fields &field);

  private:
    std::vector<Cell *> _cells;
};

/**
 * @brief Moving wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure with zero gradient
 */
class MovingWallBoundary : public Boundary {
  public:
    MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity);
    MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity, double wall_temperature);
    virtual ~MovingWallBoundary() = default;
    virtual void apply(Fields &field);
    virtual void apply_pressures(Fields &field);

  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _wall_velocity;
    double _wall_temperature;
};

/**
 * @brief Inflow wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity at the inlet,
 * Neumann for pressure with zero gradient
 */

class InFlow : public Boundary {
  public:
    InFlow(std::vector<Cell *> cells, std::map<int, double> inlet_velocity);
    InFlow(std::vector<Cell *> cells, std::map<int, double> inlet_velocity, double wall_temperature);
    virtual ~InFlow() = default;
    virtual void apply(Fields &field);
    virtual void apply_pressures(Fields &field);

  private:
    std::vector<Cell *> _cells;
    std::map<int, double> _inlet_velocity;
    double _wall_temperature;
};

/**
 * @brief Outflow wall boundary condition for the outer boundaries of the domain.
 * Neumann for velocities with zero gradient and Dirichilet for pressure with zero gradient
 */

class OutFlow : public Boundary {
  public:
    OutFlow(std::vector<Cell *> cells, double outlet_pressure);
    OutFlow(std::vector<Cell *> cells, double wall_temperature, double outlet_pressure);
    virtual ~OutFlow() = default;
    virtual void apply(Fields &field);
    virtual void apply_pressures(Fields &field);

  private:
    std::vector<Cell *> _cells;
    double _wall_temperature;
    double _outlet_pressure;
};
