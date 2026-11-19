#ifndef FLUME_INPUTS_H_
#define FLUME_INPUTS_H_

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

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
};

//! Read and check the "flume" namespace. Aborts on a bad value or an unknown key.
Inputs read_inputs ();

} // namespace flume

#endif
