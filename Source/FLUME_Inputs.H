#ifndef FLUME_INPUTS_H_
#define FLUME_INPUTS_H_

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <optional>
#include <string>

namespace flume {

//! Largest flume.n_cell in x or y: 2^30, which leaves room below the largest
//! int for ghost-cell and node indices.
inline constexpr int n_cell_max = 1073741824;

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

    //! Acceleration of gravity, m/s2 (flume.gravity; default 9.81; finite, > 0)
    amrex::Real gravity = amrex::Real(9.81);

    //! Bed elevation at the lower corner of the domain, m, on the absolute
    //! vertical datum (flume.bed.elevation; default 0; finite)
    amrex::Real bed_elevation = amrex::Real(0.0);

    //! Bed slope dz_b/dx and dz_b/dy (flume.bed.slope; default 0 0; finite)
    amrex::Array<amrex::Real, 2> bed_slope{{amrex::Real(0.0), amrex::Real(0.0)}};

    //! Elevation of the initial, level water surface, m, on the same datum
    //! (flume.init.surface_elevation; no default: the domain starts dry; finite)
    std::optional<amrex::Real> init_surface_elevation;

    //! Initial velocity in x and y where there is water, m/s (flume.init.velocity;
    //! default 0 0; finite; needs flume.init.surface_elevation)
    amrex::Array<amrex::Real, 2> init_velocity{{amrex::Real(0.0), amrex::Real(0.0)}};
};

//! Read and check the "flume" namespace. Aborts on a bad value or an unknown key.
Inputs read_inputs ();

} // namespace flume

#endif
