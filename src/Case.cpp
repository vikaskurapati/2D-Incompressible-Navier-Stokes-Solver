
#include "Case.hpp"
#include "Communication.hpp"
#include "Enums.hpp"

#include <algorithm>
#ifdef GCC_VERSION_9_OR_HIGHER
#include <filesystem>
#else
#include <experimental/filesystem>
#endif
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <typeinfo>
#include <vector>

#ifdef GCC_VERSION_9_OR_HIGHER
namespace filesystem = std::filesystem;
#else
namespace filesystem = std::experimental::filesystem;
#endif

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>

Case::Case(std::string file_name, int argn, char **args, int process_rank, int size, int my_rank) {
    // Read input parameters
    const int MAX_LINE_LENGTH = 1024;
    std::ifstream file(file_name);
    double nu;      /* viscosity   */
    double UI;      /* velocity x-direction */
    double VI;      /* velocity y-direction */
    double PI;      /* pressure */
    double GX{0.0}; /* gravitation x-direction */
    double GY{0.0}; /* gravitation y-direction */
    double xlength; /* length of the domain x-dir.*/
    double ylength; /* length of the domain y-dir.*/
    double dt;      /* time step */
    int imax;       /* number of cells x-direction*/
    int jmax;       /* number of cells y-direction*/
    double gamma;   /* uppwind differencing factor*/
    double omg;     /* relaxation factor */
    double tau;     /* safety factor for time step*/
    int itermax;    /* max. number of iterations for pressure per time step */
    double eps;     /* accuracy bound for pressure*/
    // Worksheet 2 Additions
    double TI;     /* initial Temperature */
    double alpha;  /* thermal diffusivity */
    double beta;   /* coefficient of thermal expansion */
    int num_walls; /* number of walls */
    double Tc;     /*Cold wall temperature */
    double Th;     /*Hot wall temperature */
    double T3;     /*3rd wall temperature */
    double T4;     /*4rd wall temperature */
    double T5;     /*5rd wall temperature */

    if (file.is_open()) {

        std::string var;
        while (!file.eof() && file.good()) {
            file >> var;
            if (var[0] == '#') { /* ignore comment line*/
                file.ignore(MAX_LINE_LENGTH, '\n');
            } else {
                if (var == "xlength") file >> xlength;
                if (var == "ylength") file >> ylength;
                if (var == "nu") file >> nu;
                if (var == "t_end") file >> _t_end;
                if (var == "dt") file >> dt;
                if (var == "omg") file >> omg;
                if (var == "eps") file >> eps;
                if (var == "tau") file >> tau;
                if (var == "gamma") file >> gamma;
                if (var == "dt_value") file >> _output_freq;
                if (var == "UI") file >> UI;
                if (var == "VI") file >> VI;
                if (var == "GX") file >> GX;
                if (var == "GY") file >> GY;
                if (var == "PI") file >> PI;
                if (var == "itermax") file >> itermax;
                if (var == "imax") file >> imax;
                if (var == "jmax") file >> jmax;
                // Worksheet 2 Additions
                if (var == "geo_file") file >> _geom_name;
                if (var == "TI") file >> TI;
                if (var == "alpha") file >> alpha;
                if (var == "beta") file >> beta;
                if (var == "energy_eq") file >> _energy_eq;
                if (var == "num_walls") file >> num_walls;
                if (var == "wall_temp_3") file >> T3;
                if (var == "wall_temp_4") file >> T4;
                if (var == "wall_temp_5") file >> T5;
                // Worksheet 3 Additions
                if (var == "iproc") file >> _iproc;
                if (var == "jproc") file >> _jproc;
            }
        }
    }

    file.close();

    _process_rank = process_rank;
    _size = size;

    if (_iproc * _jproc != _size) {
        MPI_Finalize();
        if (_process_rank == 0) {
            std::cerr << "iproc*jproc!=size, please check" << std::endl;
        }
        exit(0);
    }

    if (imax % _iproc != 0) {
        MPI_Finalize();
        if (_process_rank == 0) {
            std::cerr << "imax isnt divisible by iproc, please check" << std::endl;
        }
        exit(0);
    }
    if (jmax % _jproc != 0) {
        MPI_Finalize();
        if (_process_rank == 0) {
            std::cerr << "jmax isnt divisible by jproc, please check" << std::endl;
        }
        exit(0);
    }

    _datfile_name = file_name;

    std::map<int, double> wall_vel;
    if (_geom_name.compare("NONE") == 0) {
        wall_vel.insert(std::pair<int, double>(LidDrivenCavity::moving_wall_id, LidDrivenCavity::wall_velocity));
    }
    // Set file names for geometry file and output directory
    set_file_names(file_name);

    // Build up the domain
    Domain domain;
    domain.dx = xlength / static_cast<double>(imax);
    domain.dy = ylength / static_cast<double>(jmax);
    domain.domain_size_x = imax;
    domain.domain_size_y = jmax;

    build_domain(domain, imax, jmax);
    _grid = Grid(_geom_name, domain, _process_rank, _size);

    // Assigning hot and cold temperatues accordingly
    // Come back and check this once parse geometry file is written
    if (_grid.get_hot_fixed_wall_id() == 3) {
        Th = T3;
    }
    if (_grid.get_hot_fixed_wall_id() == 4) {
        Th = T4;
    }
    if (_grid.get_hot_fixed_wall_id() == 5) {
        Th = T5;
    }
    if (_grid.get_cold_fixed_wall_id() == 3) {
        Tc = T3;
    }
    if (_grid.get_cold_fixed_wall_id() == 4) {
        Tc = T4;
    }
    if (_grid.get_cold_fixed_wall_id() == 5) {
        Tc = T5;
    }

    _field = Fields(_grid, nu, dt, tau, alpha, beta, _energy_eq, _grid.domain().size_x, _grid.domain().size_y, UI, VI,
                    PI, TI, GX, GY, _process_rank, _size);

    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _pressure_solver = std::make_unique<SOR>(omg);
    _max_iter = itermax;
    _tolerance = eps;
    // Construct boundaries
    if (not _grid.moving_wall_cells().empty()) {
        _boundaries.push_back(
            std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), LidDrivenCavity::wall_velocity));
    }

    if (not _grid.fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells()));
    }

    if (not _grid.hot_fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.hot_fixed_wall_cells(), Th));
    }

    if (not _grid.cold_fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.cold_fixed_wall_cells(), Tc));
    }

    if (not _grid.adiabatic_fixed_wall_cells().empty()) {
        _boundaries.push_back(std::make_unique<AdiabaticWallBoundary>(_grid.adiabatic_fixed_wall_cells()));
    }

    if (not _grid.inflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<InFlow>(_grid.inflow_cells(),
                                                       std::map<int, double>{{PlaneShearFlow::inflow_wall_id, UI}}));
    }

    if (not _grid.outflow_cells().empty()) {
        _boundaries.push_back(std::make_unique<OutFlow>(_grid.outflow_cells(), 0.0));
    }

    if (_process_rank == 0) {
        output_log(_datfile_name, nu, UI, VI, PI, GX, GY, xlength, ylength, dt, imax, jmax, gamma, omg, tau, itermax,
                   eps, TI, alpha, beta, num_walls, Tc, Th, my_rank);
    }
}

void Case::set_file_names(std::string file_name) {
    std::string temp_dir;
    bool case_name_flag = true;
    bool prefix_flag = false;

    for (int i = file_name.size() - 1; i > -1; --i) {
        if (file_name[i] == '/') {
            case_name_flag = false;
            prefix_flag = true;
        }
        if (case_name_flag) {
            _case_name.push_back(file_name[i]);
        }
        if (prefix_flag) {
            _prefix.push_back(file_name[i]);
        }
    }

    for (int i = file_name.size() - _case_name.size() - 1; i > -1; --i) {
        temp_dir.push_back(file_name[i]);
    }

    std::reverse(_case_name.begin(), _case_name.end());
    std::reverse(_prefix.begin(), _prefix.end());
    std::reverse(temp_dir.begin(), temp_dir.end());

    _case_name.erase(_case_name.size() - 4);
    _dict_name = temp_dir;
    _dict_name.append(_case_name);
    _dict_name.append("_Output");

    if (_geom_name.compare("NONE") != 0) {
        _geom_name = _prefix + _geom_name;
    }

    // Create output directory
    filesystem::path folder(_dict_name);
    try {
        filesystem::create_directory(folder);
    } catch (const std::exception &e) {
        std::cerr << "Output directory could not be created." << std::endl;
        std::cerr << "Make sure that you have write permissions to the "
                     "corresponding location"
                  << std::endl;
    }
}

/**
 * This function is the main simulation loop. In the simulation loop, following steps are required(Serial)
 * - Calculate and apply boundary conditions for all the boundaries in _boundaries container
 *   using apply() member function of Boundary class -> done
 * - Calculate fluxes (F and G) using calculate_fluxes() member function of Fields class. -> Hope Surya done it
 * correctly Flux consists of diffusion and convection part, which are located in Discretization class
 * - Calculate right-hand-side of PPE using calculate_rs() member function of Fields class  -> Hope Surya done it
 * correctly
 * - Iterate the pressure poisson equation until the residual becomes smaller than the desired tolerance
 *   or the maximum number of the iterations are performed using solve() member function of PressureSolver class
 * - Calculate the velocities u and v using calculate_velocities() member function of Fields class
 * - Calculate the maximal timestep size for the next iteration using calculate_dt() member function of Fields class
 * - Write vtk files using output_vtk() function
 */
void Case::simulate(int my_rank) {

    double t = 0.0;
    double dt = _field.dt();
    int timestep = 0;
    int output_counter = 0;
    double t_end = _t_end;
    double err = 100;
    // starting simulation
    int iter_count = 0;

    Case::output_vtk(timestep, my_rank);
    std::ofstream output;
    int progress, last_progress;

    if (_process_rank == 0) {
        std::string str = _dict_name + "_run_log_" + std::to_string(my_rank) + ".log";

        output.open(str, std::ios_base::app);
        output << "\n\nIteration Log:\n";
        std::cout << "Simulation is Running!\nPlease Refer to " << _dict_name << "_run_log_" << my_rank
                  << ".log for Simulation log!\n";
        // Simulation Progress
        if (Case::_energy_eq == "on") {
            std::cout << "\nEnergy Equation is On\n";
        }
    }

    while (t < t_end) {
        err = 100.0;
        iter_count = 0;
        dt = _field.calculate_dt(_grid);
        dt = Communication::reduce_min(dt);
        for (size_t i = 0; i < _boundaries.size(); i++) {
            _boundaries[i]->apply(_field);
        }
        if (Case::_energy_eq == "on") {
            _field.calculate_temperatures(_grid);
            //communicate temperature
        }
        _field.calculate_fluxes(_grid);
        _field.calculate_rs(_grid);
        while (err > _tolerance && iter_count < _max_iter) {
            for (const auto &boundary : _boundaries) {
                boundary->apply_pressures(_field);
            }
            err = _pressure_solver->solve(_field, _grid, _boundaries);
            iter_count += 1;
        }
        _field.calculate_velocities(_grid);
        t += dt;
        timestep += 1;
        // break;
        if (_process_rank == 0) {
            output << std::setprecision(4) << std::fixed;
        }
        if (t - output_counter * _output_freq >= 0) {
            Case::output_vtk(timestep, my_rank);
            if (_process_rank == 0) {
                output << "Time: " << t << " Residual: " << err << " PPE Iterations: " << iter_count << std::endl;
                if (iter_count == _max_iter || std::isnan(err)) {
                    std::cout
                        << "The PPE Solver didn't converge for Time = " << t
                        << " Please check the log file and increase max iterations or other parameters for convergence"
                        << "\n";
                }
            }
            output_counter += 1;
        }
        // Printing Simulation Progress
        progress = t / t_end * 100;
        if (progress % 10 == 0 && progress != last_progress) {
            std::cout << "[";
            for (int i = 0; i < progress / 10; i++) {
                std::cout << "===";
            }
            if (progress == 100)
                std::cout << "]";
            else
                std::cout << ">";
            std::cout << " %" << progress << std::endl;
            last_progress = progress;
        }
    }

    if (t_end !=
        (output_counter - 1) * _output_freq) // Recording at t_end if the output frequency is not a multiple of t_end
    {
        Case::output_vtk(timestep, my_rank);
        output << "Time Step: " << timestep << " Residue: " << err << " PPE Iterations: " << iter_count << std::endl;
        output_counter += 1;
    }
    std::cout << "Simulation has ended\n";
    output.close();
}

void Case::output_vtk(int timestep, int my_rank) {
    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    int inflow_id = _grid.get_inflow_wall_id();
    int outflow_id = _grid.get_outflow_wall_id();

    double dx = _grid.dx();
    double dy = _grid.dy();

    double x = _grid.domain().imin * dx;
    double y = _grid.domain().jmin * dy;

    { y += dy; }
    { x += dx; }

    double z = 0;
    for (int col = 0; col < _grid.domain().size_y + 1; col++) {
        x = _grid.domain().imin * dx;
        { x += dx; }
        for (int row = 0; row < _grid.domain().size_x + 1; row++) {
            points->InsertNextPoint(x, y, z);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(_grid.domain().size_x + 1, _grid.domain().size_y + 1, 1);
    structuredGrid->SetPoints(points);

    std::vector<vtkIdType> obstacle_wall_cells;

    // Pressure Array
    vtkDoubleArray *Pressure = vtkDoubleArray::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Velocity Array
    vtkDoubleArray *Velocity = vtkDoubleArray::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    for (int i = 1; i <= _grid.imax(); ++i) {
        for (int j = 1; j <= _grid.jmax(); ++j) {
            if (_grid.cell(i, j).wall_id() != 0) {
                obstacle_wall_cells.push_back(i - 1 + (j - 1) * _grid.imax());
            }
        }
    }

    for (int i = 0; i < obstacle_wall_cells.size(); ++i) {
        structuredGrid->BlankCell(obstacle_wall_cells.at(i));
    }

    // Print pressure and temperature from bottom to top
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            double pressure = _field.p(i, j);
            Pressure->InsertNextTuple(&pressure);
        }
    }

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Print Velocity from bottom to top
    for (int j = 0; j < _grid.domain().size_y + 1; j++) {
        for (int i = 0; i < _grid.domain().size_x + 1; i++) {
            vel[0] = (_field.u(i, j) + _field.u(i, j + 1)) * 0.5;
            vel[1] = (_field.v(i, j) + _field.v(i + 1, j)) * 0.5;
            Velocity->InsertNextTuple(vel);
        }
    }

    if (_energy_eq == "on") {
        // Temperatrue Array
        vtkDoubleArray *Temperature = vtkDoubleArray::New();
        Temperature->SetName("temperature");
        Temperature->SetNumberOfComponents(1);

        // Print pressure and temperature from bottom to top
        for (int j = 1; j < _grid.domain().size_y + 1; j++) {
            for (int i = 1; i < _grid.domain().size_x + 1; i++) {
                double temperature = _field.t(i, j);
                Temperature->InsertNextTuple(&temperature);
            }
        }
        // Add Temperature to Structured Grid
        structuredGrid->GetCellData()->AddArray(Temperature);
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname = _dict_name + '/' + _case_name + "_" + std::to_string(my_rank) + "_" +
                             std::to_string(_process_rank) + "." + std::to_string(timestep) +
                             ".vtk"; // my_rank is the user's input and _process_rank is the process rank

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain) {

    int I;
    int J;
    int imin, imax, jmin, jmax, size_x, size_y;

    if (_process_rank == 0) {
        for (int i = 1; i < _size; ++i) {
            I = i % _iproc + 1;
            J = i / _iproc + 1;
            imin = (I - 1) * imax_domain / _iproc;
            imax = I * imax_domain / _iproc + 2;
            jmin = (J - 1) * jmax_domain / _jproc;
            jmax = J * jmax_domain / _jproc + 2;
            size_x = imax_domain / _iproc;
            size_y = jmax_domain / _jproc;
            MPI_Send(&imin, 1, MPI_INT, i, 999, MPI_COMM_WORLD);
            MPI_Send(&imax, 1, MPI_INT, i, 998, MPI_COMM_WORLD);
            MPI_Send(&jmin, 1, MPI_INT, i, 997, MPI_COMM_WORLD);
            MPI_Send(&jmax, 1, MPI_INT, i, 996, MPI_COMM_WORLD);
            MPI_Send(&size_x, 1, MPI_INT, i, 995, MPI_COMM_WORLD);
            MPI_Send(&size_y, 1, MPI_INT, i, 994, MPI_COMM_WORLD);
        }
        I = _process_rank % _iproc + 1;
        J = _process_rank / _iproc + 1;
        domain.imin = (I - 1) * imax_domain / _iproc;
        domain.imax = I * imax_domain / _iproc + 2;
        domain.jmin = (J - 1) * jmax_domain / _jproc;
        domain.jmax = J * jmax_domain / _jproc + 2;
        domain.size_x = imax_domain / _iproc;
        domain.size_y = jmax_domain / _jproc;
    }

    else {
        MPI_Status status;
        MPI_Recv(&domain.imin, 1, MPI_INT, 0, 999, MPI_COMM_WORLD, &status);
        MPI_Recv(&domain.imax, 1, MPI_INT, 0, 998, MPI_COMM_WORLD, &status);
        MPI_Recv(&domain.jmin, 1, MPI_INT, 0, 997, MPI_COMM_WORLD, &status);
        MPI_Recv(&domain.jmax, 1, MPI_INT, 0, 996, MPI_COMM_WORLD, &status);
        MPI_Recv(&domain.size_x, 1, MPI_INT, 0, 995, MPI_COMM_WORLD, &status);
        MPI_Recv(&domain.size_y, 1, MPI_INT, 0, 994, MPI_COMM_WORLD, &status);
    }

    if (_process_rank + 1 < _size && (_process_rank + 1) % _iproc != 0) {
        domain.neighbour_ranks[1] = _process_rank + 1;
    }
    if (_process_rank - 1 >= 0 && (_process_rank) % _iproc != 0) {
        domain.neighbour_ranks[0] = _process_rank - 1;
    }
    if (_process_rank + _iproc < _size) {
        domain.neighbour_ranks[3] = _process_rank + _iproc;
    }
    if (_process_rank - _iproc >= 0) {
        domain.neighbour_ranks[2] = _process_rank - _iproc;
    }
}

void Case::output_log(std::string dat_file_name, double nu, double UI, double VI, double PI, double GX, double GY,
                      double xlength, double ylength, double dt, double imax, double jmax, double gamma, double omg,
                      double tau, double itermax, double eps, double TI, double alpha, double beta, double num_walls,
                      double Tc, double Th, int my_rank) {

    const int MAX_LINE_LENGTH = 1024;
    std::string str = _dict_name + "_run_log_" + std::to_string(my_rank) + ".log";
    std::stringstream stream;
    // stream<<std::fixed<<std::setprecision(2)<<_pressure_solver->return_omega();
    // str += stream.str();
    std::ofstream output(str);

    output << "Log File for : " << dat_file_name << "\n";
    output << "Simulation Parameters:\n";
    output << "xlength : " << xlength << "\n";
    output << "ylength : " << ylength << "\n";
    output << "nu : " << nu << "\n";
    output << "t_end : " << _t_end << "\n";
    output << "dt : " << dt << "\n";
    output << "omg : " << omg << "\n";
    output << "eps : " << eps << "\n";
    output << "tau : " << tau << "\n";
    output << "gamma : " << gamma << "\n";
    output << "dt_value : " << _output_freq << "\n";
    output << "UI : " << UI << "\n";
    output << "VI : " << VI << "\n";
    output << "GX : " << GX << "\n";
    output << "GY : " << GY << "\n";
    output << "PI : " << PI << "\n";
    output << "itermax : " << itermax << "\n";
    output << "Energy Equation : " << _energy_eq << "\n";
    if (_energy_eq == "on") {
        output << "Temp Initial : " << TI << "\n";
        output << "alpha : " << alpha << "\n";
        output << "beta : " << beta << "\n";
        output << "No of Temperature walls: " << num_walls << "\n";
        output << "Cold Wall Temperature : " << Tc << "\n";
        output << "Hot Wall Temperature : " << Th << "\n";
    }

    output.close();
}