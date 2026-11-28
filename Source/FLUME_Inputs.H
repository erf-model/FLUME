#ifndef FLUME_INPUTS_H_
#define FLUME_INPUTS_H_

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>

#include <optional>
#include <string>
#include <vector>

namespace flume {

//! Largest flume.n_cell in x or y: 2^30, which leaves room below the largest
//! int for ghost-cell and node indices.
inline constexpr int n_cell_max = 1073741824;

//! The type of a boundary segment (flume.bc.<name>.type).
enum class BoundaryType : int { wall, outflow, zero_gradient, discharge };

//! Where the rain comes from (flume.rain.type).
enum class RainType : int { none, constant, hyetograph, gridded, moving_cell, microphysics };

//! How a hyetograph's rate, or a gridded rain field, runs between its entries
//! (flume.rain.interpolation); also the discharge of an inflow hydrograph
//! (flume.bc.<name>.interpolation).
enum class RainInterpolation : int { step, linear };

//! The infiltration loss model (flume.infiltration.type).
enum class InfiltrationType : int { none, green_ampt, curve_number };

/**
 * A soil parameter of the infiltration model: one value for every cell, or
 * an ESRI ASCII grid of values (flume.infiltration.<model>.<name> or
 * flume.infiltration.<model>.<name>_file; exactly one of the two).
 */
struct SoilParameter
{
    amrex::Real value = amrex::Real(0.0);  //!< the value, when `file` is empty
    std::string file;                      //!< the grid, or empty
};

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
    //! With type = discharge: the hydrograph file, lines of a time in s and
    //! a discharge in m3/s (flume.bc.<name>.file; required; not empty)
    std::string file;
    //! With type = discharge: how the discharge runs between the entries of
    //! the file (flume.bc.<name>.interpolation; default linear; step or linear)
    RainInterpolation interpolation = RainInterpolation::linear;
};

//! A gauge point (flume.gauges.<name>.location): the time series of the
//! cell that holds it is written to gauge_<name>.txt.
struct GaugePoint
{
    std::string name;
    amrex::Array<amrex::Real, 2> location{{amrex::Real(0.0), amrex::Real(0.0)}};  //!< x y, m
};

/**
 * A storm-drain inlet (flume.drains.<name>.*): a grate in the cell that
 * holds `location`, with the weir and orifice capacity of FLUME_Drains.H.
 */
struct DrainInlet
{
    std::string name;
    amrex::Array<amrex::Real, 2> location{{amrex::Real(0.0), amrex::Real(0.0)}};  //!< x y, m
    amrex::Real perimeter = amrex::Real(0.0);  //!< weir length P, m
    amrex::Real area = amrex::Real(0.0);       //!< open area A, m2
    amrex::Real weir_coefficient = amrex::Real(1.66);    //!< C_w, m^(1/2)/s
    amrex::Real orifice_coefficient = amrex::Real(0.67); //!< C_o, no unit
};

/**
 * An inflow inside the domain (flume.inflows.<name>.*): the discharge of a
 * hydrograph file added to the cell that holds a point, or spread over the
 * cells a line crosses by the length of the line in each
 * (FLUME_InflowSource.H).
 */
struct InflowSpec
{
    std::string name;
    bool line = false;  //!< a line (flume.inflows.<name>.line), else a point (.location)
    //! The point x y, or the ends x0 y0 x1 y1 of the line, m (the point uses 0 and 1)
    amrex::Array<amrex::Real, 4> coords{{amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0),
                                         amrex::Real(0.0)}};
    //! Lines of a time in s and a discharge in m3/s (flume.inflows.<name>.file)
    std::string file;
    //! How the discharge runs between the entries of the file
    //! (flume.inflows.<name>.interpolation; default linear; step or linear)
    RainInterpolation interpolation = RainInterpolation::linear;
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
    //! flume.bc.<name>.side, .type, .range, and with type = discharge .file
    //! and .interpolation; default none: every such side is a wall), in the
    //! order listed.
    std::vector<BoundarySegment> bc_segments;

    //! Largest box edge in cells (flume.max_grid_size; default 64; >= 1)
    int max_grid_size = 64;

    //! Plotfile name prefix (flume.plot_file; default "plt"; not empty)
    std::string plot_file = "plt";

    //! Steps between plotfiles (flume.plot_int; default 0; >= 0). The first
    //! and the last state are always written; 0 writes only those.
    int plot_int = 0;

    //! Steps between lines of the water-balance log and of the gauge files
    //! (flume.log_int; default 1; >= 0). The first and the last step are
    //! always logged; 0 writes neither and keeps no balance.
    int log_int = 1;

    //! Checkpoint name prefix (flume.check_file; default "chk"; not empty)
    std::string check_file = "chk";

    //! Steps between checkpoints (flume.check_int; default 0; >= 0). With a
    //! value above 0 the last step is also written; 0 writes none.
    int check_int = 0;

    //! Checkpoint to continue from (flume.restart; default none: the run
    //! starts from the initial state; not empty)
    std::string restart;

    //! Whether the plotfiles carry the flood maps: maximum depth, maximum
    //! speed, maximum depth times speed and first arrival time
    //! (flume.flood_maps; default 0; 0 or 1)
    int flood_maps = 0;

    //! Depth above which water has arrived in a cell, m
    //! (flume.flood_maps.arrival_depth; default flume.dry_depth; finite,
    //! >= 0; needs flume.flood_maps = 1)
    amrex::Real arrival_depth = amrex::Real(0.0);

    //! Gauge points (flume.gauges.names and flume.gauges.<name>.location;
    //! default none), in the order listed. Each writes gauge_<name>.txt.
    std::vector<GaugePoint> gauges;

    //! Storm-drain inlets (flume.drains.names and flume.drains.<name>.*;
    //! default none), in the order listed:
    //!   - location: x y, m (required; inside the domain, on its edge too);
    //!   - perimeter: weir length P, m (required; finite, > 0);
    //!   - area: open area A, m2 (required; finite, > 0 and <= the area of a cell);
    //!   - weir_coefficient: C_w, m^(1/2)/s (default 1.66; finite, > 0);
    //!   - orifice_coefficient: C_o (default 0.67; finite, > 0 and <= 1).
    //! Each inlet has a cell of its own, not a building cell (checked when
    //! the inlets are placed, FLUME_Drains.H).
    std::vector<DrainInlet> drains;

    //! Inflows inside the domain (flume.inflows.names and
    //! flume.inflows.<name>.*; default none), in the order listed:
    //!   - location: x y of a point, m, or
    //!   - line: x0 y0 x1 y1 of the ends of a line, m (exactly one of the
    //!     two; inside the domain, on its edge too; a line longer than 0);
    //!   - file: lines of a time in s and a discharge in m3/s (required; not empty);
    //!   - interpolation: step or linear (default linear).
    //! No cell of an inflow is a building cell (checked when the inflows
    //! are placed, FLUME_InflowSource.H).
    std::vector<InflowSpec> inflows;

    //! Time at which the run stops, s (flume.stop_time; no default; finite, >= 0).
    //! In double precision in every build, as the time is.
    std::optional<double> stop_time;

    //! Largest number of steps (flume.max_step; no default; >= 0). With neither
    //! stop_time nor max_step the run takes no step.
    std::optional<int> max_step;

    //! Largest time step, s (flume.max_dt; no default: the step is set by the
    //! CFL number and the rain alone; finite, > 0). In double precision in
    //! every build, as the time is. A step of max_dt ends at the time plus
    //! max_dt, as a step of a host model that uses the same fixed step does.
    std::optional<double> max_dt;

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

    //! ESRI ASCII grid of the Manning coefficient, s/m^(1/3)
    //! (flume.friction.manning_file; default none; not empty; not with
    //! flume.friction.manning). Each cell takes the raster cell that holds
    //! its centre; every value a cell takes must be >= 0.
    std::string friction_manning_file;

    //! Where the rain comes from (flume.rain.type; default none; none,
    //! constant, hyetograph, gridded or moving_cell; microphysics, the surface
    //! precipitation of ERF's microphysics, only inside ERF)
    RainType rain_type = RainType::none;

    //! Rain rate, mm/h (flume.rain.rate; finite, >= 0; required with
    //! flume.rain.type = constant and read only with it)
    double rain_rate = 0.0;

    //! With flume.rain.type = hyetograph, a hyetograph file: lines of a time
    //! in s and a rate in mm/h. With flume.rain.type = gridded, an index
    //! file: lines of a time in s and an ESRI ASCII grid of the rate in mm/h
    //! (flume.rain.file; not empty; required with these two types and read
    //! only with them)
    std::string rain_file;

    //! With flume.rain.type = microphysics, the index file to write the
    //! host's rain to, as gridded rain that a stand-alone run reads back
    //! (FLUME_RainReplay.H); empty: none (flume.rain.replay_file; default
    //! none; not empty; read only with flume.rain.type = microphysics)
    std::string rain_replay_file;

    //! How the rate runs between the entries of the hyetograph or the index
    //! (flume.rain.interpolation; default step; step or linear; read only
    //! with flume.rain.type = hyetograph or gridded)
    RainInterpolation rain_interpolation = RainInterpolation::step;

    //! Peak rain rate of the moving cell, mm/h, at its centre
    //! (flume.rain.cell.peak_rate; finite, >= 0; required with
    //! flume.rain.type = moving_cell and read only with it)
    double rain_cell_peak_rate = 0.0;

    //! Standard deviation of the moving cell's Gaussian rate, m
    //! (flume.rain.cell.sigma; finite, > 0; required with flume.rain.type =
    //! moving_cell and read only with it)
    double rain_cell_sigma = 0.0;

    //! Centre of the moving cell at t = 0, m (flume.rain.cell.center; x y;
    //! finite; required with flume.rain.type = moving_cell and read only
    //! with it)
    amrex::Array<double, 2> rain_cell_center{{0.0, 0.0}};

    //! Velocity of the moving cell, m/s (flume.rain.cell.velocity; x y;
    //! default 0 0; finite; read only with flume.rain.type = moving_cell)
    amrex::Array<double, 2> rain_cell_velocity{{0.0, 0.0}};

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

    //! ESRI ASCII grid of the bed elevation, m, on the run's vertical datum
    //! (flume.bed.file; default none: the plane bed above; not empty; not
    //! with flume.bed.elevation, .slope, .continue_slope, .bowl_depth,
    //! .bowl_radius or .noise). Bilinear between the raster's cell centres at
    //! each cell centre; the raster must hold the domain.
    std::string bed_file;

    //! ESRI ASCII grid of building heights above the bed, m
    //! (flume.buildings.file; default none: no buildings; not empty). Each
    //! cell takes the raster cell that holds its centre; a value > 0 raises
    //! the bed of the cell by that height (or by flume.buildings.height), 0
    //! leaves it; a value below 0 aborts.
    std::string buildings_file;

    //! Height given to every building cell, m (flume.buildings.height;
    //! default none: the heights of the file; finite, > 0; needs
    //! flume.buildings.file). With it, the file is a footprint mask: any
    //! value > 0 is a building.
    std::optional<amrex::Real> buildings_height;

    //! The infiltration loss model (flume.infiltration.type; default none;
    //! none, green_ampt or curve_number). Building cells
    //! (flume.buildings.file) infiltrate nothing with either model.
    InfiltrationType infiltration_type = InfiltrationType::none;

    //! Green-Ampt saturated hydraulic conductivity, mm/h
    //! (flume.infiltration.green_ampt.conductivity or .conductivity_file;
    //! finite, >= 0; one of the two required with
    //! flume.infiltration.type = green_ampt and read only with it). 0 is an
    //! impervious cell.
    SoilParameter ga_conductivity;

    //! Green-Ampt suction head at the wetting front, mm
    //! (flume.infiltration.green_ampt.suction or .suction_file; finite,
    //! >= 0; one of the two required with green_ampt and read only with it)
    SoilParameter ga_suction;

    //! Green-Ampt initial moisture deficit, saturated minus initial moisture
    //! content (flume.infiltration.green_ampt.deficit or .deficit_file;
    //! finite, >= 0 and <= 1; one of the two required with green_ampt and
    //! read only with it)
    SoilParameter ga_deficit;

    //! Runoff curve number (flume.infiltration.cn.number or .number_file;
    //! finite, > 0 and <= 100; one of the two required with
    //! flume.infiltration.type = curve_number and read only with it). 100 is
    //! an impervious cell. The number must match flume.infiltration.cn.ia_ratio:
    //! handbook tables are for 0.2.
    SoilParameter cn_number;

    //! Initial abstraction over potential retention, Ia / S
    //! (flume.infiltration.cn.ia_ratio; default 0.2; finite, >= 0 and < 1;
    //! read only with flume.infiltration.type = curve_number)
    amrex::Real cn_ia_ratio = amrex::Real(0.2);

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

/**
 * Read and check the "flume" namespace. Aborts on a bad value or an unknown key.
 *
 * With `host_grid`, FLUME runs inside a host model that owns the grid, the
 * output files and the end of the run (ERF, FLUME_ERFFloodEngine.cpp): n_cell,
 * prob_lo, prob_hi and periodic are taken from `host_grid`, and the keys the
 * host sets (n_cell, prob_lo, prob_hi, periodic, max_grid_size, plot_file,
 * plot_int, check_file, check_int, restart, stop_time, max_step) abort when
 * given.
 */
Inputs read_inputs (const amrex::Geometry* host_grid = nullptr);

} // namespace flume

#endif
