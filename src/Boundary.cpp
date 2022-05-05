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
            field.u(i,j) = 0.0;
            field.v(i,j) = - field.v(i+1,j);
            field.p(i,j) = field.p(i+1,j);
            field.f(i,j) = field.u(i,j);
        }
        if(cell->is_border(border_position::LEFT)){
            field.u(i-1,j) = 0.0;
            field.v(i,j) = - field.v(i-1,j);
            field.p(i,j) = field.p(i-1,j);
            field.f(i-1,j) = 0.0;
        }
        if(cell->is_border(border_position::TOP)){
            field.u(i,j) = -field.u(i,j+1);
            field.v(i,j) = 0.0;
            field.p(i,j) = field.p(i,j+1);
            field.g(i,j) = field.v(i,j);
        }

        //std:: cout << cell->is_border(border_position::BOTTOM) << " ";
    }
    // for (int j = 0; j < jmax; j++)
    // {
        
    //     field.u(0,j) = 0.0;
    //     field.u(imax-1,j) = 0.0;
    //     field.u(imax-2,j) = 0.0;
    //     field.v(0,j) = -field.v(1,j);
    //     field.v(imax-1, j) = -field.v(imax-2,j);
    //     field.p(0, j) = field.p(1, j);
    //     field.p(imax-1, j) = field.p(imax-2, j);
    //     field.f(0,j) = 0.0;
    //     field.f(imax-1, j) = field.u(imax-1, j);
    // }

    // for (int i = 0; i < imax; i++)
    // {
    //     field.v(i,0) = 0.0;
    //     field.u(i, 0) = -field.u(i,1);
    //     field.p(i,0) = field.p(i,1);
    //     field.g(i, 0) = field.v(i, 0);
    // }
}


MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {
    unsigned int i,j;
    //int jmax = field.p_matrix().jmax();
    for(auto& cell : _cells){
        
        if(cell->is_border(border_position::BOTTOM)){
            i = cell->i();
            j = cell->j();
            field.u(i,j) = 2.0 - field.u(i,j-1);
            field.v(i,j) = 0.0;
            field.p(i,j) = field.p(i,j-1);
            field.g(i,j) = field.v(i,j);
        }
            
    }
    //std:: cout << std::endl;
}