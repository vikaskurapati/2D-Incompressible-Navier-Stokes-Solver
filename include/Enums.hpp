#pragma once

// If no geometry file is provided in the input file, lid driven cavity case
// will run by default. In the Grid.cpp, geometry will be created following
// PGM convention, which is:
// 0: fluid, 3: fixed wall, 4: moving wall
namespace LidDrivenCavity {
const int moving_wall_id = 8;
const int fixed_wall_id = 4;
const double wall_velocity = 1.0;
} // namespace LidDrivenCavity

namespace PlaneShearFlow {
    const int inflow_wall_id = 1;
    const int outflow_wall_id = 2;
    const int fixed_wall_id = 4;
}


enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
};

namespace border {
const int TOP = 0;
const int BOTTOM = 1;
const int LEFT = 2;
const int RIGHT = 3;
} // namespace border

enum class cell_type {

    FLUID,
    INFLOW,
    OUTFLOW,
    FIXED_WALL,
    MOVING_WALL,
    DEFAULT,
    HOT_FIXED_WALL,
    COLD_FIXED_WALL,
    ADIABATIC_FIXED_WALL
};
