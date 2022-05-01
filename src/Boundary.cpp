#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {
    int imax = field.p_matrix().imax();
    int jmax = field.p_matrix().jmax();

    for (int j = 0; j < jmax; j++)
    {
        field.u(0,j) = 0.0;
        field.u(imax-1,j) = 0.0;
        field.u(imax-2,j) = 0.0;
        field.v(0,j) = -field.v(1,j);
        field.v(imax-1, j) = -field.v(imax-2,j);
    }

    for (int i = 0; i < imax; i++)
    {
        field.v(i,0) = 0.0;
        field.u(i, 0) = -field.u(i,1);
    }
}


MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {
    int imax = field.p_matrix().imax();
    int jmax = field.p_matrix().jmax();
    for (int i = 0; i < imax; i++)
    {
        field.u(i, jmax-1) = -field.u(i, jmax-2);
        field.v(i, jmax-2) = 0.0;
        field.v(i, jmax-3) = 0.0;
    }
}
