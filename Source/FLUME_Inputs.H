#ifndef FLUME_INPUTS_H_
#define FLUME_INPUTS_H_

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <optional>
#include <string>
#include <vector>

namespace flume {

//! Largest flume.n_cell in x or y: 2^30, which leaves room below the largest
//! int for ghost-cell and node indices.
inline constexpr int n_cell_max = 1073741824;

//! The type of a boundary segment (flume.bc.<name>.type).
enum class BoundaryType : int { wall, outflow };

/**
 * A segment of a side that is not periodic (flume.bc.<name>.*): the boundary
 * cells cell_lo to cell_hi along the side, all with one type. read_inputs
 * resolves flume.bc.<name>.range to these cells and checks them.
 */
struct BoundarySegment
{
    std::string name;
    int direction = 0;  //!< direction of the normal to the side: 0 x, 1 y
    int high = 0;       //!< 0 the low side, 1 the high side
    BoundaryType type = BoundaryType::wall;
    int cell_lo = 0;    //!< first boundary cell along the side
    int cell_hi = -1;   //!< last boundary cell along the side
};

/**
 * The run's inputs, all in the "flume" namespace of the inputs file.
 *
 * Every key has a range check that aborts with a message naming the key. A key
 * in the "flume" namespace that no reader asked for also aborts, so a
 * misspelled key cannot be ignored silently. A value of the wrong type (text
 * for a number, 1.5 for an integer) is stopped by AMReX's parser instead.
 */
struct Inputs
{
    //! Number of cells in x and y (flume.n_cell; required; each 1 to n_cell_max)
    amrex::Array<int, 2> n_cell{{0, 0}};

    //! Lower corner of the domain in x and y, m (flume.prob_lo; default 0 0; finite)
    amrex::Array<amrex::Real, 2> prob_lo{{amrex::Real(0.0), amrex::Real(0.0)}};

    //! Upper corner of the domain in x and y, m (flume.prob_hi; required; finite;
    //! > prob_lo; (prob_hi - prob_lo) / n_cell resolvable in amrex::Real)
    amrex::Array<amrex::Real, 2> prob_hi{{amrex::Real(0.0), amrex::Real(0.0)}};

    //! Whether the domain is periodic in x and y (flume.periodic; default 0 0;
    //! each 0 or 1). A side that is not periodic is made of the segments
    //! below, and is a wall where no segment covers it.
    amrex::Array<int, 2> periodic{{0, 0}};

    //! Segments of the sides that are not periodic (flume.bc.segments and
    //! flume.bc.<name>.side, .type, .range; default none: every such side is
    //! a wall), in the order listed.
    std::vector<BoundarySegment> bc_segments;

    //! Largest box edge in cells (flume.max_grid_size; default 64; >= 1)
    int max_grid_size = 64;

    //! Plotfile name prefix (flume.plot_file; default "plt"; not empty)
    std::string plot_file = "plt";

    //! Steps between plotfiles (flume.plot_int; default 0; >= 0). The first
    //! and the last state are always written; 0 writes only those.
    int plot_int = 0;

    //! Time at which the run stops, s (flume.stop_time; no default; finite, >= 0).
    //! In double precision in every build, as the time is.
    std::optional<double> stop_time;

    //! Largest number of steps (flume.max_step; no default; >= 0). With neither
    //! stop_time nor max_step the run takes no step.
    std::optional<int> max_step;

    //! CFL number of the time step, with the rates of the two directions
    //! added (flume.cfl; default 0.5; > 0 and <= 1)
    amrex::Real cfl = amrex::Real(0.5);

    //! Order of the scheme in space and time (flume.order; default 1; 1 or 2).
    //! 1: cell values at the faces, forward Euler. 2: piecewise-linear
    //! reconstruction with the minmod limiter, two-stage step.
    int order = 1;

    //! Acceleration of gravity, m/s2 (flume.gravity; default 9.81; finite, > 0)
    amrex::Real gravity = amrex::Real(9.81);

    //! Depth at or below which a cell, or the depth a cell presents to a
    //! face, counts as dry, m (flume.dry_depth; default 1e-6; finite, >= 0)
    amrex::Real dry_depth = amrex::Real(1.0e-6);

    //! Manning coefficient of the bed, s/m^(1/3) (flume.friction.manning;
    //! default 0: no friction; finite, >= 0). The same in every cell.
    amrex::Real friction_manning = amrex::Real(0.0);

    //! Bed elevation at the lower corner of the domain, m, on the absolute
    //! vertical datum (flume.bed.elevation; default 0; finite)
    amrex::Real bed_elevation = amrex::Real(0.0);

    //! Bed slope dz_b/dx and dz_b/dy (flume.bed.slope; default 0 0; finite)
    amrex::Array<amrex::Real, 2> bed_slope{{amrex::Real(0.0), amrex::Real(0.0)}};

    //! Whether the plane of the bed continues across a periodic side
    //! (flume.bed.continue_slope; default 0; 0 or 1; 1 needs flume.periodic = 1
    //! in a direction where flume.bed.slope is not 0). 0: the bed repeats, so
    //! a sloping bed has a step of slope x domain length at the side. 1: the
    //! ghost cells beyond the side take the bed of the cell they cover plus
    //! slope x the distance between the two, an endless inclined channel.
    int bed_continue_slope = 0;

    //! Depth of a paraboloid bowl added to the bed, m (flume.bed.bowl_depth;
    //! default 0: no bowl; finite; needs flume.bed.bowl_radius): the bed moves
    //! by -bowl_depth (1 - r^2 / bowl_radius^2), r the distance from the centre
    //! of the domain. A negative value makes a mound.
    amrex::Real bed_bowl_depth = amrex::Real(0.0);

    //! Radius at which the bowl crosses the plane bed, m (flume.bed.bowl_radius;
    //! finite, > 0; required with flume.bed.bowl_depth and only read with it)
    amrex::Real bed_bowl_radius = amrex::Real(0.0);

    //! Amplitude of a rough bed, m (flume.bed.noise; default 0; finite, >= 0):
    //! each cell's bed is moved by a fixed pseudo-random amount in
    //! [-noise, noise), the same on any number of ranks.
    amrex::Real bed_noise = amrex::Real(0.0);

    //! Elevation of the initial, level water surface, m, on the same datum
    //! (flume.init.surface_elevation; no default: the domain starts dry; finite)
    std::optional<amrex::Real> init_surface_elevation;

    //! Slope of the initial water surface in x and y (flume.init.surface_slope;
    //! default 0 0; finite; needs flume.init.surface_elevation; not with a
    //! dam). The surface is surface_elevation + slope . (x - centre)
    //! + surface_curvature r^2, measured from the centre of the domain.
    amrex::Array<amrex::Real, 2> init_surface_slope{{amrex::Real(0.0), amrex::Real(0.0)}};

    //! Curvature of the initial water surface, 1/m (flume.init.surface_curvature;
    //! default 0; finite; needs flume.init.surface_elevation; not with a dam)
    amrex::Real init_surface_curvature = amrex::Real(0.0);

    //! Amplitude of a sine wave on the initial water surface, m
    //! (flume.init.wave_amplitude; default 0: no wave; finite; needs
    //! flume.init.surface_elevation and a periodic domain along the wave; not
    //! with a dam, a slope or a curvature). One wavelength spans the domain
    //! along flume.init.wave_direction, and the velocity along it is that of
    //! a simple wave: flume.init.velocity + 2 (sqrt(g h) - sqrt(g h_level)),
    //! with h_level the depth under the level surface.
    amrex::Real init_wave_amplitude = amrex::Real(0.0);

    //! Direction the wave runs along: 0 for x, 1 for y (flume.init.wave_direction;
    //! "x" or "y"; default x; needs flume.init.wave_amplitude)
    int init_wave_direction = 0;

    //! Initial velocity in x and y where there is water, m/s (flume.init.velocity;
    //! default 0 0; finite; needs flume.init.surface_elevation)
    amrex::Array<amrex::Real, 2> init_velocity{{amrex::Real(0.0), amrex::Real(0.0)}};

    //! Coordinate of a dam across the domain, m (flume.init.dam_position; no
    //! default: no dam; finite; inside the domain in the dam direction; needs
    //! flume.init.surface_elevation and flume.init.dam_surface_elevation).
    //! Cells whose centre is below it start under flume.init.surface_elevation,
    //! the others under flume.init.dam_surface_elevation.
    std::optional<amrex::Real> init_dam_position;

    //! Elevation of the initial water surface on the high-coordinate side of
    //! the dam, m (flume.init.dam_surface_elevation; finite; required with
    //! flume.init.dam_position and only read with it)
    amrex::Real init_dam_surface_elevation = amrex::Real(0.0);

    //! Initial velocity on the high-coordinate side of the dam, m/s
    //! (flume.init.dam_velocity; default: flume.init.velocity; finite; needs
    //! flume.init.dam_position)
    amrex::Array<amrex::Real, 2> init_dam_velocity{{amrex::Real(0.0), amrex::Real(0.0)}};

    //! Direction the dam position is measured along: 0 for x, 1 for y
    //! (flume.init.dam_direction; "x" or "y"; default x; needs
    //! flume.init.dam_position)
    int init_dam_direction = 0;
};

//! Read and check the "flume" namespace. Aborts on a bad value or an unknown key.
Inputs read_inputs ();

} // namespace flume

#endif
