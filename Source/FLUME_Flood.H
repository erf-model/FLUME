#ifndef FLUME_FLOOD_H_
#define FLUME_FLOOD_H_

#include "FLUME_Inputs.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <string>

namespace flume {

/**
 * The flood model on one level: grid, state, bed, time loop and plotfiles.
 *
 * The grid is two-dimensional: nx by ny cells and one cell in z. The z extent
 * of the domain is 0 to 1 m and carries no meaning.
 *
 * The state (FLUME_State.H) and the bed elevation have no ghost cells yet:
 * nothing reads a neighbour until the flux kernel exists. Until then a step
 * advances the time and leaves the state as it is.
 */
class Flood
{
public:
    //! Build the grid and fill the bed and the initial state. Aborts if any
    //! initial value is not finite.
    explicit Flood (const Inputs& in);

    //! Write the initial plotfile, run the time loop, write the last plotfile.
    void evolve ();

    /**
     * The next time step, s: cfl / max_wave_rate(). A step that reaches
     * flume.stop_time, or falls short of it by round-off only (4 machine
     * epsilons of the stop time), is set to end on it. With no wave anywhere
     * (a dry domain) the step goes straight to flume.stop_time; without a
     * stop time it aborts. Aborts on a negative or not-finite depth or wave
     * speed, and on a step that is zero, not finite or too small to advance
     * the time. `lands` is set when the step ends at flume.stop_time.
     */
    double compute_dt (bool& lands) const;

    //! Write <plot_file><step>: depth, discharge_x, discharge_y,
    //! bed_elevation, surface_elevation (bed + depth).
    void write_plotfile () const;

    [[nodiscard]] double time () const { return m_time; }
    [[nodiscard]] int step () const { return m_step; }
    [[nodiscard]] const amrex::MultiFab& state () const { return m_state; }
    [[nodiscard]] const amrex::MultiFab& bed () const { return m_bed; }
    [[nodiscard]] const amrex::Geometry& geom () const { return m_geom; }

private:
    Inputs m_in;
    amrex::Geometry m_geom;
    amrex::BoxArray m_ba;
    amrex::DistributionMapping m_dm;
    amrex::MultiFab m_state; //!< State::ncomp components, no ghost cells
    amrex::MultiFab m_bed;   //!< bed elevation, m, no ghost cells
    //! Time, s. In double precision in every build: a single-precision sum of
    //! steps stops growing after about 2^24 steps.
    double m_time = 0.0;
    int m_step = 0;
};

} // namespace flume

#endif
