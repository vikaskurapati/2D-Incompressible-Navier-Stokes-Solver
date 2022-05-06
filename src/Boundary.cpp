#include "Boundary.hpp"
#include <cmath>
#include <iostream>
#include "Grid.hpp"
FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {

    unsigned int i,j;
    for(auto& cell : _cells){
        i = cell->i();
        j = cell->j();
        if(cell->is_border(border_position::RIGHT)){
            field.u(i,j) = 0;
            field.v(i,j) = -field.v(i+1,j);
            field.p(i,j) = field.p(i+1,j);
        }
        if(cell->is_border(border_position::LEFT)){
            field.u(i-1,j) = 0;
            field.v(i,j) = -field.v(i-1,j);
            field.p(i,j) = field.p(i-1,j);
        }
        else if(cell->is_border(border_position::TOP)){
            field.u(i,j) = -field.u(i,j+1);
            field.v(i,j) = 0;
            field.p(i,j) = field.p(i,j+1);
        }
    }
}


MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {
    unsigned int i,j;
    for(auto& cell : _cells){
        if(cell->is_border(border_position::BOTTOM)){
            i = cell->i();
            j = cell->j();
            field.u(i,j) = 2.0*LidDrivenCavity::wall_velocity - field.u(i,j-1);
            field.v(i,j-1) = 0;
            field.p(i,j) = field.p(i,j-1);
        }
    }
}