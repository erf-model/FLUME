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
 * Each inlet sits in the cell that holds its location and takes water at a
 * rate Q set by its capacity law (flume.drains.<name>.type).
 *
 * weir_orifice (the default): the smaller of the weir and the orifice
 * capacity of a grate in a sag at the depth d of the cell (FHWA HEC-22, 3rd
 * edition, 2009):
 *
 *   Q_w = C_w P d^(3/2),   Q_o = C_o A sqrt(2 g d),   Q = min(Q_w, Q_o),
 *
 * with P the perimeter, A the open area, C_w the weir coefficient (1.66 in
 * SI units) and C_o the orifice coefficient (0.67). The weir rate is the
 * smaller below the depth d* = C_o A sqrt(2 g) / (C_w P) and the
 * orifice rate above it.
 *
 * approach_flow: the capacity of an inlet in water that flows over it, from
 * the speed u and the depth d of the cell (Xia et al. 2022, in the form of
 * Annis, Badas and Mascaro 2026, their equation 3):
 *
 *   Q = k u A Fr^alpha,   Fr = u / sqrt(g d),
 *
 * with A the inlet area, k a coefficient (0.302) and alpha an exponent
 * (-0.816). With alpha above -1 an inlet takes nothing from still water.
 *
 * Over a step the cell depth follows dh/dt = -Q / (dx dy), integrated in
 * closed form (drain_depth_after, approach_depth_after; the second holds the
 * speed of the cell fixed over the step, as the removal of water does): the
 * loss does not depend on the step, and the depth never falls below 0.
 * The water leaves the model: no pipe network, no surcharge back to the
 * surface. The depth each cell has drained since the start of the run is
 * kept, written to the plotfiles as `drained_depth`, and counted in the
 * water balance as the column `drained`.
 *
 * An inlet is given by its location, or a list file gives one inlet per
 * line, all with the parameters of their name (flume.drains.<name>.file).
 * An inlet of a location has its own cell. The inlets of one file may share
 * cells: n of them in a cell take n times the rate of one.
 *
 * The source acts after the rain and the infiltration, before the friction
 * (make_sources). A building cell holds no inlet.
 */
class Drains : public Source
{
public:
    //! Place the inlets of in.drains; the list files are read here
    //! (read_inlet_file). Aborts, naming the inlet's location key, or its
    //! file and line, when two inlets that are not of one file share a cell
    //! or an inlet sits in a building cell (`buildings` above 0; may be
    //! null).
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
    //! approach_flow: r = k A g^(-alpha/2) / (dx dy), in m^p s^(-2p) with
    //! p = -alpha/2, so that dh/dt = -r u^(1 + alpha) h^p; 0 without such
    //! an inlet
    amrex::MultiFab m_rate;
    amrex::MultiFab m_exponent;  //!< alpha of the approach_flow inlet of the cell
    bool m_any_weir = false;      //!< a weir_orifice inlet is placed
    bool m_any_approach = false;  //!< an approach_flow inlet is placed
    amrex::MultiFab m_drained;
    amrex::MultiFab m_drained_carry;  //!< rounding carried in the sum m_drained
    amrex::Real m_dry_depth = amrex::Real(0.0);
    double m_cell_area = 0.0;  //!< m2
    //! One word per name, in the checkpoint; a restart with other inlets
    //! stops. An inlet of a location: name:i:j:P:A:C_w:C_o (weir_orifice)
    //! or name:i:j:approach_flow:A:k:alpha. The inlets of a file:
    //! name:file:<inlets>:<cells>:<hash of the cells and their counts>:
    //! then P:A:C_w:C_o or approach_flow:A:k:alpha.
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

/**
 * The depth, m, of a cell that starts a time `tau` (s) at the depth `h0`
 * (m) and loses water only to an approach_flow inlet at a fixed speed,
 * dh/dt = -c h^p, with c >= 0 (m^(1-p)/s) and p in [0, 1/2]. With q = 1 - p,
 * h^q = h0^q - q c t, which reaches 0 after t = h0^q / (q c): the cell is
 * then empty. With z = q c tau / h0^q, the part of that time that has
 * passed, h = h0 (1 - z)^(1/q), taken as h0 exp(log1p(-z) / q): a loss far
 * below h0 is not the difference of two nearly equal powers. So the result
 * is in [0, h0]. In double precision. (The step, apply_approach_drains,
 * also takes a film that this leaves at or below the dry depth.)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
approach_depth_after (double h0, double c, double p, double tau) noexcept
{
    if (!(h0 > 0.0) || !(tau > 0.0) || !(c > 0.0)) { return h0; }
    const double q = 1.0 - p;
    const double z = q * c * tau / std::pow(h0, q);
    if (!(z < 1.0)) { return 0.0; }
    const double r = h0 * std::exp(std::log1p(-z) / q);
    return (r < h0) ? r : h0;
}

//! The step of Drains::apply on the valid cells of `state`: each cell where
//! `weir` is above 0 drains for `dt` (drain_depth_after), and the depth
//! removed is added to `drained` (rounding carried in `drained_carry`). The
//! ghost cells of `state` are not read.
void apply_drains (amrex::MultiFab& state, const amrex::MultiFab& weir,
                   const amrex::MultiFab& orifice, amrex::MultiFab& drained,
                   amrex::MultiFab& drained_carry, amrex::Real dry_depth, amrex::Real dt);

//! The same for the approach_flow inlets: each cell where `rate` is above 0
//! and the depth above `dry_depth` drains for `dt` at the speed it has
//! (approach_depth_after, with c = rate u^(1 + exponent) and
//! p = -exponent / 2). A cell at or below `dry_depth` holds no discharge
//! and is left.
void apply_approach_drains (amrex::MultiFab& state, const amrex::MultiFab& rate,
                            const amrex::MultiFab& exponent, amrex::MultiFab& drained,
                            amrex::MultiFab& drained_carry, amrex::Real dry_depth,
                            amrex::Real dt);

/**
 * Read a list of inlet locations (flume.drains.<name>.file, `key` in the
 * messages): lines of an x and a y in m; `#` starts a comment that runs to
 * the end of the line, and blank lines are skipped. The I/O rank reads the
 * file and sends the points to every rank. Every rank aborts with the same
 * message, naming the key, the file and the line, when the file cannot be
 * opened or is not a regular file, a line does not hold exactly two items,
 * an item is not a finite number, a point is outside the domain of `geom`
 * (its edge is inside), or the file holds no point. `line` takes the line
 * of each point.
 */
void read_inlet_file (const std::string& key, const std::string& path,
                      const amrex::Geometry& geom, std::vector<double>& x,
                      std::vector<double>& y, std::vector<int>& line);

} // namespace flume

#endif
