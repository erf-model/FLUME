#ifndef FLUME_RAIN_H_
#define FLUME_RAIN_H_

#include "FLUME_Forcing.H"
#include "FLUME_Source.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

namespace flume {

/**
 * Rain as a source: each step adds to the depth of every cell the rain that
 * the provider (FLUME_Forcing.H) gives for that cell from the start to the
 * end of the step. Rain falls straight down, so it carries no horizontal
 * momentum: the discharge is not changed, and the velocity falls as the
 * depth grows. Rain falls on every cell, dry or wet, building or street.
 *
 * Over a run, the depth added to a cell is the integral of the rain rate over
 * the run's time, to round-off: the depth of each step is the difference of
 * the integral at the two ends of the step, and the ends of consecutive steps
 * are the same number.
 */
class Rain : public Source
{
public:
    //! `forcing` must outlive this source.
    Rain (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
          const amrex::Geometry& geom, const Forcing& forcing);

    [[nodiscard]] std::string name () const override { return "rain"; }

    //! Add the rain from ctx.time to ctx.time_end to the depth of the valid
    //! cells of `state`.
    void apply (amrex::MultiFab& state, const SourceContext& ctx, amrex::Real dt) override;

    //! The rain depth of the last step, m, per cell (one component, no ghost
    //! cells): set by apply, read by the curve-number infiltration.
    [[nodiscard]] const amrex::MultiFab& step_depth () const { return m_depth; }

    [[nodiscard]] std::string balance_name () const override { return "rain"; }

    //! The rain added to the depth of this rank's cells since the start, m3:
    //! the cell area times the sum, over the steps, of the rain depth of each
    //! step, less the part the compensated sum still holds back (`m_carry`,
    //! at most the rounding of a depth), which is not yet in the state.
    [[nodiscard]] double balance_volume () const override;

private:
    const Forcing& m_forcing;
    amrex::MultiFab m_depth;  //!< rain depth of the step, m: one component, no ghost cells
    amrex::MultiFab m_carry;  //!< rain rounded away so far, m, added with the next step's
    double m_cell_area = 0.0; //!< m2
    double m_volume = 0.0;    //!< rain on this rank's cells since the start, m3
    double m_volume_carry = 0.0;  //!< rounding carried in the sum m_volume
};

//! Add component 0 of `depth`, plus the rain that earlier additions rounded
//! away (`carry`, updated), to the depth of the valid cells of `state` (same
//! boxes). The discharges are not read or written.
void add_rain_depth (amrex::MultiFab& state, const amrex::MultiFab& depth, amrex::MultiFab& carry);

} // namespace flume

#endif
