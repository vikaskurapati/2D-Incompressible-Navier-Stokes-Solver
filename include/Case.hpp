#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Boundary.hpp"
#include "Discretization.hpp"
#include "Domain.hpp"
#include "Fields.hpp"
#include "Grid.hpp"
#include "PressureSolver.hpp"

/**
 * @brief Class to hold and orchestrate the simulation flow.
 *
 */
class Case {
  public:
    /**
     * @brief Parallel constructor for the Case.
     *
     * Reads input file, creates Fields, Grid, Boundary, Solver and sets
     * Discretization parameters Creates output directory
     *
     * @param[in] Input file name
     */
    Case(std::string file_name, int argn, char **args, int process_rank = 0, int size = 1, int my_rank = 0);

    /**
     * @brief Main function to simulate the flow until the end time(Serial Implementation).
     *
     * Calculates the fluxes
     * Calculates the right hand side
     * Solves pressure
     * Calculates velocities
     * Outputs the solution files
     */
    void simulate(int my_rank = 1);

    /**
     * @brief Destroy the Case object
     *
     */
    ~Case();

  private:
    /// Plain case name without paths
    std::string _case_name;
    /// Output directiory name
    std::string _dict_name;
    // Dat File name
    std::string _datfile_name;
    /// Geometry file name
    std::string _geom_name{"NONE"};
    /// Relative input file path
    std::string _prefix;
    /// Energy equation should be consider or not
    std::string _energy_eq = "off";

    /// Simulation time
    double _t_end;
    /// Solution file outputting frequency
    double _output_freq;

    int _iproc{1};
    int _jproc{1};
    int _process_rank{0};
    int _size{1};

    std::string _solver_type;

    Fields _field;
    Grid _grid;
    Domain domain;
    Discretization _discretization;
    std::unique_ptr<PressureSolver> _pressure_solver;
    std::vector<std::unique_ptr<Boundary>> _boundaries;

    /// Solver convergence tolerance
    double _tolerance;

    /// Maximum number of iterations for the solver
    int _max_iter;

    /**
     * @brief Creating file names from given input data file
     *
     * Extracts path of the case file and creates code-readable file names
     * for outputting directory and geometry file.
     *
     * @param[in] input data file
     */
    void set_file_names(std::string file_name);

    /**
     * @brief Solution file outputter
     *
     * Outputs the solution files in .vtk format. Ghost cells are excluded.
     * Pressure is cell variable while velocity is point variable while being
     * interpolated to the cell faces
     *
     * @param[in] Timestep of the solution
     */
    void output_vtk(int t, int my_rank = 0);
    void output_log(std::string dat_file_name, double nu, double UI, double VI, double PI, double GX, double GY,
                    double xlength, double ylength, double dt, double imax, double jmax, double gamma, double omg,
                    double tau, double itermax, double eps, double TI, double alpha, double beta, double num_walls,
                    double Tc, double Th, int my_rank);
    std::ofstream simulation_log_file(int my_rank = 0);
    void build_domain(Domain &domain, int imax_domain, int jmax_domain);
};
