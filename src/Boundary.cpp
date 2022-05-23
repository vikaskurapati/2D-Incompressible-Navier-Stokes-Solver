#include "Boundary.hpp"
#include <cmath>
#include <iostream>
#include "Grid.hpp"
FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, double wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {
    int i,j;
    for(auto& cell : _cells){
        i = cell->i();
        j = cell->j();
        if(cell->is_border(border_position::TOP))
        {
            if (cell->is_border(border_position::RIGHT))
            {
                field.u(i,j) = 0.0;
                field.v(i,j) = 0.0;
                field.u(i-1,j) = -field.u(i-1,j+1);
                field.v(i,j-1) = -field.v(i+1,j-1);
                field.p(i,j) = 0.5*(field.p(i,j+1) + field.p(i+1, j));
                if(field.get_energy_eq()){
                    field.t(i,j) = 0.5*(4.0*_wall_temperature - field.t(i,j+1) - field.t(i+1,j));
                }
            }
            else if(cell->is_border(border_position::LEFT))
            {
                field.u(i-1,j) = 0.0;
                field.v(i,j) = 0.0;
                field.u(i,j) = -field.u(i,j+1);
                field.v(i,j-1) = -field.v(i-1,j-1);
                field.p(i,j) = 0.5*(field.p(i,j+1)+field.p(i-1,j));
                if(field.get_energy_eq()){
                    field.t(i,j) = 0.5*(4.0*_wall_temperature - field.t(i,j+1) - field.t(i-1,j));
                }
            }
            else if(cell->is_border(border_position::BOTTOM))
            {
                field.u(i,j) = -field.u(i,j-1);
                field.v(i,j) = 0.0;
                field.v(i,j-1) = 0.0;
                field.u(i-1,j) = 0.0;
                field.p(i,j) = 0.5*(field.p(i,j+1)+field.p(i,j-1));
                if(field.get_energy_eq()){
                    field.t(i,j) = 0.5*(4.0*_wall_temperature - field.t(i,j+1) - field.t(i,j-1));
                }
            }
            else{
                field.u(i,j) = -field.u(i,j+1);
                field.v(i,j) = 0.0;
                field.p(i,j) = field.p(i,j+1);
                if(field.get_energy_eq()){
                    field.t(i,j) = 2.0*_wall_temperature - field.t(i,j+1);
                }
            }
        }
        else if (cell->is_border(border_position::BOTTOM))
        {
            if(cell->is_border(border_position::RIGHT))
            {
                field.u(i,j) = 0.0;
                field.v(i,j-1) = 0.0;
                field.u(i-1,j) = -field.u(i-1,j-1);
                field.v(i,j) = -field.v(i+1, j);
                field.p(i,j) = 0.5*(field.p(i+1,j) + field.p(i,j-1));
                if(field.get_energy_eq()){
                    field.t(i,j) = 0.5*(4.0*_wall_temperature - field.t(i+1,j) - field.t(i,j-1));
                }
            }
            else if(cell->is_border(border_position::LEFT))
            {
                field.u(i-1,j) = 0.0;
                field.v(i,j-1) = 0.0;
                field.u(i,j) = -field.u(i,j-1);
                field.v(i,j) = -field.v(i-1,j);
                field.p(i,j) = 0.5*(field.p(i,j-1) + field.p(i-1,j));
                if(field.get_energy_eq()){
                    field.t(i,j) = 0.5*(4.0*_wall_temperature - field.t(i,j-1) - field.t(i-1,j));
                }
             }

            else{
            field.u(i,j) = -field.u(i,j-1);
            field.v(i, j-1) = 0;
            field.p(i,j) = field.p(i,j-1);
            if(field.get_energy_eq()){
                field.t(i,j) = 2.0*_wall_temperature - field.t(i,j-1);
            }
            }
        }
        else if(cell->is_border(border_position::RIGHT))
        {
            if(cell->is_border(border_position::LEFT))
            {
                field.u(i,j) = 0.0;
                field.u(i-1,j) = 0.0;
                field.v(i,j-1) = 0.5*(field.v(i+1,j-1)+field.v(i-1,j-1));
                field.v(i,j) = -0.5*(field.v(i+1,j) + field.v(i-1,j));
                field.p(i,j) = 0.5*(field.p(i+1,j) + field.p(i-1,j));
                if(field.get_energy_eq()){
                    field.t(i,j) = 0.5*(4.0*_wall_temperature - field.t(i+1,j) - field.t(i-1,j));
                }
            }
            else{
            field.u(i,j) = 0;
            field.v(i,j) = -field.v(i+1,j);
            field.p(i,j) = field.p(i+1,j);
            if(field.get_energy_eq()){
                field.t(i,j) = 2.0*_wall_temperature - field.t(i+1,j);
            }
            }
        }
        else if(cell->is_border(border_position::LEFT))
        {
            field.u(i-1,j) = 0;
            field.v(i,j) = -field.v(i-1,j);
            field.p(i,j) = field.p(i-1,j);
            if(field.get_energy_eq()){
                field.t(i,j) = 2.0*_wall_temperature - field.t(i-1,j);
            }
        }
    }
}

void FixedWallBoundary::apply_pressures(Fields &field)
{
    int i, j;
    for (const auto& cell: _cells)
    {
        i = cell->i();
        j = cell->j();
        if(cell->is_border(border_position::TOP))
        {
            if (cell->is_border(border_position::RIGHT))
            {
                field.p(i,j) = 0.5*(field.p(i,j+1) + field.p(i+1, j));
            }
            else if(cell->is_border(border_position::LEFT))
            {
                field.p(i,j) = 0.5*(field.p(i,j+1)+field.p(i-1,j));
            }
            else if(cell->is_border(border_position::BOTTOM))
            {
                field.p(i,j) = 0.5*(field.p(i,j+1)+field.p(i,j-1));
            }
            else{
                field.p(i,j) = field.p(i,j+1);
            }
        }
        else if (cell->is_border(border_position::BOTTOM))
        {
            if(cell->is_border(border_position::RIGHT))
            {
                field.p(i,j) = 0.5*(field.p(i+1,j) + field.p(i,j-1));
            }
            else if(cell->is_border(border_position::LEFT))
            {
                field.p(i,j) = 0.5*(field.p(i,j-1) + field.p(i-1,j));
             }

            else{
            field.p(i,j) = field.p(i,j-1);
            }
        }
        else if(cell->is_border(border_position::RIGHT))
        {
            if(cell->is_border(border_position::LEFT))
            {
                field.p(i,j) = 0.5*(field.p(i+1,j) + field.p(i-1,j));
            }
            else{
            field.p(i,j) = field.p(i+1,j);
            }
        }
        else if(cell->is_border(border_position::LEFT))
        {
            field.p(i,j) = field.p(i-1,j);
        }
    }
}

AdiabaticWallBoundary::AdiabaticWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

void AdiabaticWallBoundary::apply(Fields &field) {
    int i,j;
    for(auto& cell : _cells){
        i = cell->i();
        j = cell->j();
        // std::cout << "Adiabatic";
        if(cell->is_border(border_position::TOP))
        {
            if (cell->is_border(border_position::RIGHT))
            {
                field.u(i,j) = 0.0;
                field.v(i,j) = 0.0;
                field.u(i-1,j) = -field.u(i-1,j+1);
                field.v(i,j-1) = -field.v(i+1,j-1);
                field.p(i,j) = 0.5*(field.p(i,j+1) + field.p(i+1, j));
                field.t(i,j) = 0.5*(field.t(i,j+1) + field.t(i+1, j));
            }
            else if(cell->is_border(border_position::LEFT))
            {
                field.u(i-1,j) = 0.0;
                field.v(i,j) = 0.0;
                field.u(i,j) = -field.u(i,j+1);
                field.v(i,j-1) = -field.v(i-1,j-1);
                field.p(i,j) = 0.5*(field.p(i,j+1)+field.p(i-1,j));
                field.t(i,j) = 0.5*(field.t(i,j+1)+field.t(i-1,j));
            }
            else if(cell->is_border(border_position::BOTTOM))
            {
                field.u(i,j) = -field.u(i,j-1);
                field.v(i,j) = 0.0;
                field.v(i,j-1) = 0.0;
                field.u(i-1,j) = 0.0;
                field.p(i,j) = 0.5*(field.p(i,j+1)+field.p(i,j-1));
                field.t(i,j) = 0.5*(field.t(i,j+1)+field.t(i,j-1));
            }
            else{
                field.u(i,j) = -field.u(i,j+1);
                field.v(i,j) = 0.0;
                field.p(i,j) = field.p(i,j+1);
                field.t(i,j) = field.t(i,j+1);
            }
        }
        else if (cell->is_border(border_position::BOTTOM))
        {
            if(cell->is_border(border_position::RIGHT))
            {
                field.u(i,j) = 0.0;
                field.v(i,j-1) = 0.0;
                field.u(i-1,j) = -field.u(i-1,j-1);
                field.v(i,j) = -field.v(i+1, j);
                field.p(i,j) = 0.5*(field.p(i+1,j) + field.p(i,j-1));
                field.t(i,j) = 0.5*(field.t(i+1,j) + field.t(i,j-1));
            }
            else if(cell->is_border(border_position::LEFT))
            {
                field.u(i-1,j) = 0.0;
                field.v(i,j-1) = 0.0;
                field.u(i,j) = -field.u(i,j-1);
                field.v(i,j) = -field.v(i-1,j);
                field.p(i,j) = 0.5*(field.p(i,j-1) + field.p(i-1,j));
                field.t(i,j) = 0.5*(field.t(i,j-1) + field.t(i-1,j));
             }

            else{
            field.u(i,j) = -field.u(i,j-1);
            field.v(i, j-1) = 0;
            field.p(i,j) = field.p(i,j-1);
            field.t(i,j) = field.t(i,j-1);
            }
        }
        else if(cell->is_border(border_position::RIGHT))
        {
            if(cell->is_border(border_position::LEFT))
            {
                field.u(i,j) = 0.0;
                field.u(i-1,j) = 0.0;
                field.v(i,j-1) = 0.5*(field.v(i+1,j-1)+field.v(i-1,j-1));
                field.v(i,j) = -0.5*(field.v(i+1,j) + field.v(i-1,j));
                field.p(i,j) = 0.5*(field.p(i+1,j) + field.p(i-1,j));
                field.t(i,j) = 0.5*(field.t(i+1,j) + field.t(i-1,j));
            }
            else{
            field.u(i,j) = 0;
            field.v(i,j) = -field.v(i+1,j);
            field.p(i,j) = field.p(i+1,j);
            field.t(i,j) = field.t(i+1,j);
            }
        }
        else if(cell->is_border(border_position::LEFT))
        {
            field.u(i-1,j) = 0;
            field.v(i,j) = -field.v(i-1,j);
            field.p(i,j) = field.p(i-1,j);
            field.t(i,j) = field.t(i-1,j);
        }
    }
}
void AdiabaticWallBoundary::apply_pressures(Fields &field)
{
    int i, j;
    for (const auto& cell: _cells)
    {
        i = cell->i();
        j = cell->j();
        if(cell->is_border(border_position::TOP))
        {
            if (cell->is_border(border_position::RIGHT))
            {
                field.p(i,j) = 0.5*(field.p(i,j+1) + field.p(i+1, j));
            }
            else if(cell->is_border(border_position::LEFT))
            {
                field.p(i,j) = 0.5*(field.p(i,j+1)+field.p(i-1,j));
            }
            else if(cell->is_border(border_position::BOTTOM))
            {
                field.p(i,j) = 0.5*(field.p(i,j+1)+field.p(i,j-1));
            }
            else{
                field.p(i,j) = field.p(i,j+1);
            }
        }
        else if (cell->is_border(border_position::BOTTOM))
        {
            if(cell->is_border(border_position::RIGHT))
            {
                field.p(i,j) = 0.5*(field.p(i+1,j) + field.p(i,j-1));
            }
            else if(cell->is_border(border_position::LEFT))
            {
                field.p(i,j) = 0.5*(field.p(i,j-1) + field.p(i-1,j));
             }

            else{
            field.p(i,j) = field.p(i,j-1);
            }
        }
        else if(cell->is_border(border_position::RIGHT))
        {
            if(cell->is_border(border_position::LEFT))
            {
                field.p(i,j) = 0.5*(field.p(i+1,j) + field.p(i-1,j));
            }
            else{
            field.p(i,j) = field.p(i+1,j);
            }
        }
        else if(cell->is_border(border_position::LEFT))
        {
            field.p(i,j) = field.p(i-1,j);
        }
    }
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       double wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) 
{
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

void MovingWallBoundary::apply_pressures(Fields &field)
{
    unsigned int i, j;
    for (const auto& cell: _cells)
    {
        if(cell->is_border(border_position::BOTTOM)){
            i = cell->i();
            j = cell->j();
            field.p(i,j) = field.p(i,j-1);
        }
    }
}

InFlow::InFlow(std::vector<Cell *> cells, std::map<int, double> inlet_velocity) : _cells(cells), 
                                                                            _inlet_velocity(inlet_velocity){};

InFlow::InFlow(std::vector<Cell *> cells, std::map<int, double> inlet_velocity, 
            double wall_temperature): _cells(cells), _inlet_velocity(inlet_velocity), 
                                                    _wall_temperature(wall_temperature){};

void InFlow::apply(Fields &field)
{
    unsigned int i, j;
    for(const auto& cell: _cells)
    {
        i = cell->i();
        j = cell->j();
        if(cell->is_border(border_position::RIGHT))
        {
            //assuming inlet velocity is only in u
            field.u(i, j) = _inlet_velocity[PlaneShearFlow::inflow_wall_id];
            field.v(i,j) = -field.v(i+1, j);
            field.p(i,j) = field.p(i+1,j);

        }
        if (cell->is_border(border_position::LEFT))
        {
            field.u(i-1, j) = _inlet_velocity[PlaneShearFlow::inflow_wall_id];
            field.v(i,j) = -field.v(i-1,j);
            field.p(i,j) = field.p(i-1,j);
        }
    }
}

void InFlow::apply_pressures(Fields &field)
{
    unsigned int i, j;
    for(const auto& cell: _cells)
    {
        i = cell->i();
        j = cell->j();
        if(cell->is_border(border_position::RIGHT))
        {
            //assuming inlet velocity is only in u
            field.p(i,j) = field.p(i+1,j);

        }
        if (cell->is_border(border_position::LEFT))
        {
            field.p(i,j) = field.p(i-1,j);
        }
    }
}

OutFlow::OutFlow(std::vector<Cell *>cells, double outlet_pressure):_cells(cells), _outlet_pressure(outlet_pressure) {};
OutFlow::OutFlow(std::vector<Cell *>cells, double wall_temperature, double outlet_pressure):_cells(cells), 
                                                                        _wall_temperature(wall_temperature),
                                                                        _outlet_pressure(outlet_pressure){};

void OutFlow::apply(Fields &field)
{
    //this deals with boundary sharing on only one border of the cell. 
    unsigned int i, j;
    for (const auto& cell: _cells)
    {
        i = cell->i();
        j = cell->j();
        if(cell -> is_border(border_position::LEFT))
        {
            field.u(i, j) = field.u(i-1,j);
            field.v(i, j) = field.v(i-1,j);
            field.p(i,j) = _outlet_pressure;
        }
        if(cell -> is_border(border_position::RIGHT))
        {
            field.u(i,j) = field.u(i+1,j);
            field.v(i, j) = field.v(i+1, j);
            field.p(i,j) = _outlet_pressure;
        }
        if(cell -> is_border(border_position::TOP))
        {
            field.u(i,j) = field.u(i, j+1);
            field.v(i,j) = field.v(i,j+1);
            field.p(i,j) = _outlet_pressure;
        }
        if(cell ->is_border(border_position::BOTTOM))
        {
            field.u(i,j-1) = field.u(i,j);
            field.v(i,j-1) = field.v(i,j);
            field.p(i,j) = _outlet_pressure;
        }
    }
}

void OutFlow::apply_pressures(Fields &field)
{
    //this deals with boundary sharing on only one border of the cell. 
    unsigned int i, j;
    for (const auto& cell: _cells)
    {
        i = cell->i();
        j = cell->j();
        if(cell -> is_border(border_position::LEFT))
        {
            field.p(i,j) = _outlet_pressure;
        }
        if(cell -> is_border(border_position::RIGHT))
        {
            field.p(i,j) = _outlet_pressure;
        }
        if(cell -> is_border(border_position::TOP))
        {
            field.p(i,j) = _outlet_pressure;
        }
        if(cell ->is_border(border_position::BOTTOM))
        {
            field.p(i,j) = _outlet_pressure;
        }
    }
}