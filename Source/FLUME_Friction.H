#ifndef FLUME_FRICTION_H_
#define FLUME_FRICTION_H_

#include "FLUME_Source.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

namespace flume {

/**
 * Manning friction on the discharge, implicit in the discharge, one cell at
 * a time.
 *
 * With h the depth, q = (hu, hv) the discharge, n the Manning coefficient
 * and g gravity, friction alone changes the discharge by
 *
 *   dq/dt = -g n^2 |q| q / h^(7/3) .
 *
 * The step of `dt` is backward Euler in q with the depth after the step held
 * fixed: q_new + dt g n^2 |q_new| q_new / h^(7/3) = q, solved in closed form.
 * q_new has the direction of q and the size r |q|, with
 *
 *   r = 2 / (1 + sqrt(1 + 4 dt g n^2 |q| / h^(7/3)))   in (0, 1] .
 *
 * So friction never reverses the flow and needs no limit on the step. The
 * update acts once over the whole step, after the fluxes and the bed slope
 * of all its stages (Docs/sphinx_doc/Model.rst).
 *
 * The coefficient is a per-cell field (s/m^(1/3)), filled from
 * flume.friction.manning. make_sources() builds this source only when that
 * coefficient is above 0.
 */
class Friction : public Source
{
public:
    //! Allocate the coefficient field on the boxes of the state and fill it
    //! with `manning` (> 0).
    Friction (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
              amrex::Real manning, amrex::Real gravity, amrex::Real dry_depth);

    [[nodiscard]] std::string name () const override { return "manning_friction"; }

    //! One implicit friction step of `dt` on the valid cells of `state`.
    void apply (amrex::MultiFab& state, const SourceContext& ctx, amrex::Real dt) override;

    //! The coefficient field, s/m^(1/3): one component, no ghost cells.
    [[nodiscard]] const amrex::MultiFab& manning () const { return m_manning; }

private:
    amrex::MultiFab m_manning;
    amrex::Real m_gravity = amrex::Real(0.0);
    amrex::Real m_dry_depth = amrex::Real(0.0);
};

/**
 * The step of Friction::apply on the valid cells of `state`, with the
 * coefficient of each cell from `manning` (same boxes, one component).
 *
 * A cell whose depth is at or below `dry_depth` (or not finite) is not
 * touched: update_state has already set its discharge to zero. A depth whose
 * power 7/3 underflows to zero gives a zero discharge. A cell without
 * discharge keeps it at zero. The ghost cells of `state` are not read.
 */
void apply_manning_friction (amrex::MultiFab& state, const amrex::MultiFab& manning,
                             amrex::Real gravity, amrex::Real dry_depth, amrex::Real dt);

} // namespace flume

#endif
