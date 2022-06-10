#include "Grid.hpp"
#include "Enums.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

Grid::Grid(std::string geom_name, Domain &domain) {

    MPI_Comm_rank(MPI_COMM_WORLD, &_process_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_size);

    _domain = domain;

    _cells = Matrix<Cell>(_domain.size_x + 2, _domain.size_y + 2);

    if (geom_name.compare("NONE")) {
        std::vector<std::vector<int>> geometry_data(_domain.domain_size_x + 2,
                                                    std::vector<int>(_domain.domain_size_y + 2, 0));
        parse_geometry_file(geom_name, geometry_data);
        assign_cell_types(geometry_data, geom_name);
    }

    else {
        build_lid_driven_cavity(geom_name);
    }
}

void Grid::build_lid_driven_cavity(std::string geom_name) {
    std::vector<std::vector<int>> geometry_data(_domain.size_x + 2, std::vector<int>(_domain.size_y + 2, 0));

    for (int i = 0; i < _domain.size_x + 2; ++i) {
        for (int j = 0; j < _domain.size_y + 2; ++j) {
            // Bottom, left and right walls: no-slip
            if ((i == 0 && _domain.imin == 0) || (j == 0 && _domain.jmin == 0) ||
                (i == _domain.size_x + 1 && _domain.imax == _domain.domain_size_x + 2)) {
                geometry_data.at(i).at(j) = LidDrivenCavity::fixed_wall_id;
            }
            // Top wall: moving wall
            else if (j == _domain.size_y + 1 && _domain.jmax == _domain.domain_size_y + 2) {
                geometry_data.at(i).at(j) = LidDrivenCavity::moving_wall_id;
            }
        }
    }

    assign_cell_types(geometry_data, geom_name);
}

void Grid::assign_cell_types(std::vector<std::vector<int>> &geometry_data, std::string geom_name) {

    int i = 0;
    int j = 0;

    if (geom_name.compare("NONE") == 0) {
        for (int j_geom = 0; j_geom < _domain.size_y+2; ++j_geom) {
            { i = 0; }
            for (int i_geom = 0; i_geom < _domain.size_x+2; ++i_geom) {
                if (geometry_data.at(i_geom).at(j_geom) == LidDrivenCavity::fixed_wall_id) {
                    _cells(i, j) = Cell(i, j, cell_type::FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                    _fixed_wall_cells.push_back(&_cells(i, j));
                } else if (geometry_data.at(i_geom).at(j_geom) == LidDrivenCavity::moving_wall_id) {
                    _cells(i, j) = Cell(i, j, cell_type::MOVING_WALL, geometry_data.at(i_geom).at(j_geom));
                    _moving_wall_cells.push_back(&_cells(i, j));
                } else {
                    _cells(i, j) = Cell(i, j, cell_type::FLUID);
                    _fluid_cells.push_back(&_cells(i, j));
                }
                ++i;
            }
            ++j;
        }
    } else {
    for (int j_geom = 0; j_geom < _domain.jmax - _domain.jmin; ++j_geom) {
        { i = 0; }
        for (int i_geom = 0; i_geom < _domain.imax - _domain.imin; ++i_geom) {
            if (geometry_data.at(i_geom).at(j_geom) == 0) {
                _cells(i, j) = Cell(i, j, cell_type::FLUID);
                _fluid_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == inflow_wall_id) {
                _cells(i, j) = Cell(i, j, cell_type::INFLOW, geometry_data.at(i_geom).at(j_geom));
                _inflow_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == outflow_wall_id) {
                _cells(i, j) = Cell(i, j, cell_type::OUTFLOW, geometry_data.at(i_geom).at(j_geom));
                _outflow_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == fixed_wall_id) {
                _cells(i, j) = Cell(i, j, cell_type::FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                _fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == hot_fixed_wall_id) {
                _cells(i, j) = Cell(i, j, cell_type::HOT_FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                _hot_fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == cold_fixed_wall_id) {
                _cells(i, j) = Cell(i, j, cell_type::COLD_FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                _cold_fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == adiabatic_fixed_wall_id) {
                _cells(i, j) = Cell(i, j, cell_type::ADIABATIC_FIXED_WALL, geometry_data.at(i_geom).at(j_geom));
                _adiabatic_fixed_wall_cells.push_back(&_cells(i, j));
            } else if (geometry_data.at(i_geom).at(j_geom) == moving_wall_id) {
                _cells(i, j) = Cell(i, j, cell_type::MOVING_WALL, geometry_data.at(i_geom).at(j_geom));
                _moving_wall_cells.push_back(&_cells(i, j));
            }
            ++i;
        }
        ++j;
    }
    }
    // Corner cell neighbour assigment
    // Bottom-Left Corner
    i = 0;
    j = 0;
    _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
    _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
    if (_cells(i, j).type() != cell_type::FLUID) {
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::TOP);
        }
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
    }
    // Top-Left Corner
    i = 0;
    j = _domain.size_y + 1;
    _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
    _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
    if (_cells(i, j).type() != cell_type::FLUID) {
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::BOTTOM);
        }
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
    }

    // Top-Right Corner
    i = _domain.size_x + 1;
    j = _domain.size_y + 1;
    _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
    _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
    if (_cells(i, j).type() != cell_type::FLUID) {
        if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::BOTTOM);
        }
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::LEFT);
        }
    }

    // Bottom-Right Corner
    i = _domain.size_x + 1;
    j = 0;
    _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
    _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
    if (_cells(i, j).type() != cell_type::FLUID) {
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::TOP);
        }
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::LEFT);
        }
    }

    // Bottom cells
    j = 0;
    for (int i = 1; i < _domain.size_x + 1; ++i) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::RIGHT);
        }
        if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::LEFT);
        }
        if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
            _cells(i, j).add_border(border_position::TOP);
        }
    }

    // Top Cells
    j = _domain.size_y + 1;

    for (int i = 1; i < _domain.size_x + 1; ++i) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        if (_cells(i, j).type() != cell_type::FLUID) {

            if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::RIGHT);
            }
            if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::LEFT);
            }
            if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::BOTTOM);
            }
        }
    }

    // Left Cells
    i = 0;
    for (int j = 1; j < _domain.size_y + 1; ++j) {
        _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).type() != cell_type::FLUID) {
            if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::RIGHT);
            }
            if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::BOTTOM);
            }
            if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::TOP);
            }
        }
    }

    // Right Cells
    i = _domain.size_x + 1;
    for (int j = 1; j < _domain.size_y + 1; ++j) {
        _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
        _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);
        _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
        if (_cells(i, j).type() != cell_type::FLUID) {
            if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::LEFT);
            }
            if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::BOTTOM);
            }
            if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
                _cells(i, j).add_border(border_position::TOP);
            }
        }
    }

    // Inner cells
    for (int i = 1; i < _domain.size_x + 1; ++i) {
        for (int j = 1; j < _domain.size_y + 1; ++j) {
            _cells(i, j).set_neighbour(&_cells(i + 1, j), border_position::RIGHT);
            _cells(i, j).set_neighbour(&_cells(i - 1, j), border_position::LEFT);
            _cells(i, j).set_neighbour(&_cells(i, j + 1), border_position::TOP);
            _cells(i, j).set_neighbour(&_cells(i, j - 1), border_position::BOTTOM);

            if (_cells(i, j).type() != cell_type::FLUID) {
                if (_cells(i, j).neighbour(border_position::LEFT)->type() == cell_type::FLUID) {
                    _cells(i, j).add_border(border_position::LEFT);
                }
                if (_cells(i, j).neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
                    _cells(i, j).add_border(border_position::RIGHT);
                }
                if (_cells(i, j).neighbour(border_position::BOTTOM)->type() == cell_type::FLUID) {
                    _cells(i, j).add_border(border_position::BOTTOM);
                }
                if (_cells(i, j).neighbour(border_position::TOP)->type() == cell_type::FLUID) {
                    _cells(i, j).add_border(border_position::TOP);
                }
            }
        }
    }
    if (_process_rank == 0) {
        std::cout << std::endl;
        for (int j = _domain.size_y + 1; j >= 0; --j) {
            for (int i = 0; i < _domain.size_x + 2; ++i) {
                std::cout << int(_cells(i, j).type()) << " ";
            }
            std::cout << std::endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Grid::parse_geometry_file(std::string filedoc, std::vector<std::vector<int>> &geometry_data) {

    int numcols, numrows, depth;
    char wall_id;

    std::ifstream infile(filedoc);
    std::ifstream legend_lines(filedoc);
    std::stringstream ss;
    std::string inputLine = "";

    // First line : version
    getline(infile, inputLine);
    if (inputLine.compare("P2") != 0) {
        std::cerr << "First line of the PGM file should be P2" << std::endl;
    }

    // Second line : comment
    getline(infile, inputLine);

    // Continue with a stringstream
    ss << infile.rdbuf();
    // Third line : size
    ss >> numrows >> numcols;
    // Fourth line : depth
    ss >> depth;

    // Following lines : data
    for (int col = numcols - 1; col > -1; --col) {
        for (int row = 0; row < numrows; ++row) {
            ss >> geometry_data[row][col];
        }
    }

    // To read the last few lines to get the id for the wall which are defined in the pgm file
    size_t found;
    while (getline(legend_lines, inputLine)) {
        found = inputLine.find("Inflow");
        if (found != std::string::npos) {
            inflow_wall_id = inputLine[2] - '0';
        }
        found = inputLine.find("Outflow");
        if (found != std::string::npos) {
            outflow_wall_id = inputLine[2] - '0';
        }
        found = inputLine.find("Wall/Obstacle") && !inputLine.find("Wall/Obstacle (hot)") &&
                !inputLine.find("Wall/Obstacle (cold)") && !inputLine.find("Wall/Obstacle (adiabatic)");
        bool found1 = (inputLine.find("Wall/Obstacle") != std::string::npos);
        found1 = found1 && (inputLine.find("Wall/Obstacle (hot)") == std::string::npos);
        found1 = found1 && (inputLine.find("Wall/Obstacle (cold)") == std::string::npos);
        found1 = found1 && (inputLine.find("Wall/Obstacle (adiabatic") == std::string::npos);
        if (found1) {
            fixed_wall_id = inputLine[2] - '0';
        }
        found = inputLine.find("Wall/Obstacle (hot)");
        if (found != std::string::npos) {
            hot_fixed_wall_id = inputLine[2] - '0';
        }
        found = inputLine.find("Wall/Obstacle (cold)");
        if (found != std::string::npos) {
            cold_fixed_wall_id = inputLine[2] - '0';
        }
        found = inputLine.find("Wall/Obstacle (adiabatic)");
        if (found != std::string::npos) {
            adiabatic_fixed_wall_id = inputLine[2] - '0';
        }
        found = inputLine.find("MovingWall");
        if (found != std::string::npos) {
            moving_wall_id = inputLine[2] - '0';
        }
    }

    infile.close();

    // Checking for boundary cells if they have more fluid cells and throw error to the user
    int count;
    for (int col = 1; col < numcols - 1; ++col) {
        for (int row = 1; row < numrows - 1; ++row) {
            if (geometry_data[row][col] != 0) {
                count = 0;
                if (geometry_data[row - 1][col] == 0) {
                    count++;
                }
                if (geometry_data[row][col - 1] == 0) {
                    count++;
                }
                if (geometry_data[row + 1][col] == 0) {
                    count++;
                }
                if (geometry_data[row][col + 1] == 0) {
                    count++;
                }
                if (count > 2) {
                    std::cerr
                        << "Error: Given PGM file has a boundary with more fluid neighbors\nPlease check the file\n";
                    exit(0);
                }
            }
        }
    }

    // Checking for outer boundary cells if they have fluid cells and throw error to the use
    for (int col = 1; col < numcols - 1; ++col) {
        if (geometry_data[0][col] == 0 or geometry_data[numrows - 1][col] == 0) {
            std::cerr << "Error: Given PGM file has a outer boundary cell with fluid wall id\nPlease check the file\n";
            exit(0);
        }
    }
    for (int row = 1; row < numrows - 1; ++row) {
        if (geometry_data[row][0] == 0 or geometry_data[row][numcols - 1] == 0) {
            std::cerr << "Error: Given PGM file has a outer boundary cell with fluid wall id\nPlease check the file\n";
            exit(0);
        }
    }
}

int Grid::imax() const { return _domain.size_x; }
int Grid::jmax() const { return _domain.size_y; }

int Grid::imaxb() const { return _domain.size_x + 2; }
int Grid::jmaxb() const { return _domain.size_y + 2; }

Cell Grid::cell(int i, int j) const { return _cells(i, j); }

double Grid::dx() const { return _domain.dx; }

double Grid::dy() const { return _domain.dy; }

const Domain &Grid::domain() const { return _domain; }

const std::vector<Cell *> &Grid::fluid_cells() const { return _fluid_cells; }

const std::vector<Cell *> &Grid::fixed_wall_cells() const { return _fixed_wall_cells; }

const std::vector<Cell *> &Grid::moving_wall_cells() const { return _moving_wall_cells; }

const std::vector<Cell *> &Grid::inflow_cells() const { return _inflow_cells; }

const std::vector<Cell *> &Grid::outflow_cells() const { return _outflow_cells; }

const std::vector<Cell *> &Grid::hot_fixed_wall_cells() const { return _hot_fixed_wall_cells; }

const std::vector<Cell *> &Grid::cold_fixed_wall_cells() const { return _cold_fixed_wall_cells; }

const std::vector<Cell *> &Grid::adiabatic_fixed_wall_cells() const { return _adiabatic_fixed_wall_cells; }