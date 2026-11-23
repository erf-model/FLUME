#ifndef FLUME_DRAINS_H_
#define FLUME_DRAINS_H_

#include "FLUME_Inputs.H"
#include "FLUME_Source.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Extension.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <string>
#include <vector>

namespace flume {

/**
 * Storm-drain inlets: water that leaves the surface into a drain
 * (flume.drains.*, Docs/sphinx_doc/theory/SourceTerms.rst).
 *
 * Each inlet sits in the cell that holds its location and takes, at the
 * depth d of that cell, the smaller of the weir and the orifice capacity of
 * a grate in a sag (FHWA HEC-22, 3rd edition, 2009):
 *
 *   Q_w = C_w P d^(3/2),   Q_o = C_o A sqrt(2 g d),   Q = min(Q_w, Q_o),
 *
 * with P the perimeter, A the open area, C_w the weir coefficient (1.66 in
 * SI units) and C_o the orifice coefficient (0.67). The weir rate is the
 * smaller below the depth d* = C_o A sqrt(2 g) / (C_w P) and the
 * orifice rate above it. Over a step the cell depth follows
 * dh/dt = -Q(h) / (dx dy), integrated in closed form (drain_depth_after):
 * the loss does not depend on the step, and the depth never falls below 0.
 * The water leaves the model: no pipe network, no surcharge back to the
 * surface. The depth each cell has drained since the start of the run is
 * kept, written to the plotfiles as `drained_depth`, and counted in the
 * water balance as the column `drained`.
 *
 * The source acts after the rain and the infiltration, before the friction
 * (make_sources). Each inlet has its own cell; a building cell holds none.
 */
class Drains : public Source
{
public:
    //! Place the inlets of in.drains. Aborts, naming the inlet's location
    //! key, when two inlets share a cell or an inlet sits in a building cell
    //! (`buildings` above 0; may be null).
    Drains (const Inputs& in, const amrex::Geometry& geom, const amrex::BoxArray& ba,
            const amrex::DistributionMapping& dm, const amrex::MultiFab* buildings);

    [[nodiscard]] std::string name () const override { return "storm_drains"; }
    void apply (amrex::MultiFab& state, const SourceContext& ctx, amrex::Real dt) override;

    [[nodiscard]] std::vector<PlotField> plot_fields () const override
    {
        return {{"drained_depth", &m_drained}};
    }

    [[nodiscard]] std::string balance_name () const override { return "drained"; }
    [[nodiscard]] bool balance_removes () const override { return true; }

    //! Minus the drained volume of this rank's cells, m3: the cell area
    //! times the sum of the drained depth and the rounding carried in it.
    [[nodiscard]] double balance_volume () const override;

    //! The drained depth and the rounding carried in it.
    void write_checkpoint (CheckpointWriter& chk) const override;
    void read_checkpoint (const CheckpointReader& chk) override;

    //! Depth each cell has drained since the start, m: one component, no
    //! ghost cells.
    [[nodiscard]] const amrex::MultiFab& drained () const { return m_drained; }

private:
    amrex::MultiFab m_weir;     //!< a = C_w P / (dx dy), 1/(m^(1/2) s); 0 without an inlet
    amrex::MultiFab m_orifice;  //!< b = C_o A sqrt(2 g) / (dx dy), m^(1/2)/s
    amrex::MultiFab m_drained;
    amrex::MultiFab m_drained_carry;  //!< rounding carried in the sum m_drained
    amrex::Real m_dry_depth = amrex::Real(0.0);
    double m_cell_area = 0.0;  //!< m2
    //! name:i:j:P:A:C_w:C_o of each inlet, in the checkpoint; a restart
    //! with other inlets stops.
    std::vector<std::string> m_words;
};

/**
 * The depth, m, of a cell that starts a time `tau` (s) at the depth `h0`
 * (m) and loses water only to an inlet, dh/dt = -min(a h^(3/2), b h^(1/2)),
 * with a > 0 (1/(m^(1/2) s)) and b > 0 (m^(1/2)/s). The rates cross at
 * d* = b / a. Above d* (orifice): sqrt(h) = sqrt(h0) - b t / 2, which
 * reaches d* after t1 = 2 (sqrt(h0) - sqrt(d*)) / b. At or below d* (weir):
 * h = h_s / (1 + a t sqrt(h_s) / 2)^2 from the depth h_s it starts at. So
 * the result is in [0, h0]: the weir branch never empties a cell. In double
 * precision.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
drain_depth_after (double h0, double a, double b, double tau) noexcept
{
    if (!(h0 > 0.0) || !(tau > 0.0)) { return h0; }
    const double ds = b / a;
    double h = h0;
    double t = tau;
    if (h > ds) {
        const double s0 = std::sqrt(h);
        const double sd = std::sqrt(ds);
        const double t1 = 2.0 * (s0 - sd) / b;
        if (t <= t1) {
            const double s = s0 - 0.5 * b * t;
            const double r = (s > sd) ? s * s : ds;
            return (r < h0) ? r : h0;
        }
        h = ds;
        t -= t1;
    }
    const double x = 1.0 + 0.5 * a * t * std::sqrt(h);
    const double r = h / (x * x);
    return (r < h0) ? r : h0;
}

//! The step of Drains::apply on the valid cells of `state`: each cell where
//! `weir` is above 0 drains for `dt` (drain_depth_after), and the depth
//! removed is added to `drained` (rounding carried in `drained_carry`). The
//! ghost cells of `state` are not read.
void apply_drains (amrex::MultiFab& state, const amrex::MultiFab& weir,
                   const amrex::MultiFab& orifice, amrex::MultiFab& drained,
                   amrex::MultiFab& drained_carry, amrex::Real dry_depth, amrex::Real dt);

} // namespace flume

#endif
