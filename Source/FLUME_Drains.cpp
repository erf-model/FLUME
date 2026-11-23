#include "FLUME_Drains.H"
#include "FLUME_Gauges.H"
#include "FLUME_RemoveWater.H"
#include "FLUME_State.H"
#include "FLUME_WaterBalance.H"

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_Print.H>

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

namespace flume {

Drains::Drains (const Inputs& in, const amrex::Geometry& geom, const amrex::BoxArray& ba,
                const amrex::DistributionMapping& dm, const amrex::MultiFab* buildings)
    : m_weir(ba, dm, 1, 0), m_orifice(ba, dm, 1, 0), m_drained(ba, dm, 1, 0),
      m_drained_carry(ba, dm, 1, 0), m_dry_depth(in.dry_depth), m_cell_area(cell_area(geom))
{
    m_weir.setVal(amrex::Real(0.0));
    m_orifice.setVal(amrex::Real(0.0));
    m_drained.setVal(amrex::Real(0.0));
    m_drained_carry.setVal(amrex::Real(0.0));
    const double root_2g = std::sqrt(2.0 * static_cast<double>(in.gravity));
    std::vector<amrex::IntVect> cells;
    for (std::size_t n = 0; n < in.drains.size(); ++n) {
        const DrainInlet& d = in.drains[n];
        const std::string key = "flume.drains." + d.name + ".location";
        const amrex::IntVect iv = gauge_cell(d.location, geom);
        for (std::size_t m = 0; m < cells.size(); ++m) {
            if (cells[m] == iv) {
                amrex::Abort(key + " puts inlet " + d.name + " in cell ("
                             + std::to_string(iv[0]) + ", " + std::to_string(iv[1])
                             + "), which holds inlet " + in.drains[m].name
                             + " already: a cell holds one inlet (merge them into one with the"
                               " sum of their perimeters and areas)");
            }
        }
        cells.push_back(iv);
        const amrex::Box cell(iv, iv);
        const double a = static_cast<double>(d.weir_coefficient) * static_cast<double>(d.perimeter)
                         / m_cell_area;
        const double b = static_cast<double>(d.orifice_coefficient) * static_cast<double>(d.area)
                         * root_2g / m_cell_area;
        m_weir.setVal(static_cast<amrex::Real>(a), cell, 0, 1, 0);
        m_orifice.setVal(static_cast<amrex::Real>(b), cell, 0, 1, 0);
        std::ostringstream w;
        w << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10) << d.name << ':'
          << iv[0] << ':' << iv[1] << ':' << d.perimeter << ':' << d.area << ':'
          << d.weir_coefficient << ':' << d.orifice_coefficient;
        m_words.push_back(w.str());
        // The depth at which the weir and the orifice rates are equal (in a
        // single-precision build the kernel's, from a and b rounded to
        // amrex::Real, differs from this one by about 1e-7).
        const double ds = b / a;
        amrex::Print() << "FLUME: drain inlet " << d.name << " at (" << d.location[0] << ", "
                       << d.location[1] << ") m in cell (" << iv[0] << ", " << iv[1]
                       << "): perimeter " << d.perimeter << " m, open area " << d.area
                       << " m2, weir coefficient " << d.weir_coefficient
                       << ", orifice coefficient " << d.orifice_coefficient
                       << "; orifice above a depth of " << ds << " m\n";
    }
    // No inlet in a building cell. One pass over the grid for all the
    // inlets (the weir rate is > 0 in the inlet cells only); the cells are
    // searched one by one only to name the inlet before stopping.
    if (buildings != nullptr && !cells.empty()) {
        amrex::MultiFab both(m_weir.boxArray(), m_weir.DistributionMap(), 1, 0);
        amrex::MultiFab::Copy(both, *buildings, 0, 0, 1, 0);
        amrex::MultiFab::Multiply(both, m_weir, 0, 0, 1, 0);
        if (both.max(0, 0) > amrex::Real(0.0)) {
            for (std::size_t n = 0; n < cells.size(); ++n) {
                if (buildings->max(amrex::Box(cells[n], cells[n]), 0, 0) > amrex::Real(0.0)) {
                    const std::string& g = in.drains[n].name;
                    amrex::Abort("flume.drains." + g + ".location puts inlet " + g
                                 + " in a building cell (flume.buildings.file)");
                }
            }
        }
    }
}

void
Drains::apply (amrex::MultiFab& state, const SourceContext& /*ctx*/, amrex::Real dt)
{
    apply_drains(state, m_weir, m_orifice, m_drained, m_drained_carry, m_dry_depth, dt);
}

double
Drains::balance_volume () const
{
    // As for infiltration: with the carry the volume follows the state to
    // round-off of a double, also in a single-precision build.
    return -(volume(m_cell_area, local_sum(m_drained))
             + volume(m_cell_area, local_sum(m_drained_carry)));
}

void
Drains::write_checkpoint (CheckpointWriter& chk) const
{
    chk.put_words(name() + "_inlets", m_words);
    chk.put_multifab(name() + "_depth", m_drained);
    chk.put_multifab(name() + "_depth_carry", m_drained_carry);
}

void
Drains::read_checkpoint (const CheckpointReader& chk)
{
    // The drained depth of a cell belongs to its inlet: the inlets must be
    // those the checkpoint was run with.
    if (chk.get_words(name() + "_inlets") != m_words) {
        std::string was;
        std::string now;
        for (const auto& w : chk.get_words(name() + "_inlets")) { was += (was.empty() ? "" : " ") + w; }
        for (const auto& w : m_words) { now += (now.empty() ? "" : " ") + w; }
        chk.fail("its drain inlets (" + was + ") differ from those of the inputs (" + now
                 + "): name:i:j:perimeter:area:weir_coefficient:orifice_coefficient of each,"
                   " from flume.drains.*");
    }
    chk.get_multifab(name() + "_depth", m_drained);
    chk.get_multifab(name() + "_depth_carry", m_drained_carry);
}

void
apply_drains (amrex::MultiFab& state, const amrex::MultiFab& weir, const amrex::MultiFab& orifice,
              amrex::MultiFab& drained, amrex::MultiFab& drained_carry, amrex::Real dry_depth,
              amrex::Real dt)
{
    const auto& s = state.arrays();
    const auto& aa = weir.const_arrays();
    const auto& bb = orifice.const_arrays();
    const auto& ff = drained.arrays();
    const auto& fc = drained_carry.arrays();
    // Each cell reads and writes only itself. In double precision.
    amrex::ParallelFor(state,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            const double a = static_cast<double>(aa[box_no](i, j, k));
            if (!(a > 0.0)) { return; }  // no inlet
            const auto& c = s[box_no];
            const amrex::Real h = c(i, j, k, State::h);
            // No water, or a NaN depth (left for the check of the state; an
            // infinite one gives a NaN loss below and is left too).
            if (!(h > amrex::Real(0.0))) { return; }
            const double b = static_cast<double>(bb[box_no](i, j, k));
            const double w = static_cast<double>(h);
            const double h1 = drain_depth_after(w, a, b, static_cast<double>(dt));
            // drain_depth_after returns a depth in [0, w], so w - h1 is in
            // [0, w], and so is its rounding to amrex::Real, since w = h is
            // an amrex::Real: the loss never exceeds the water of the cell.
            const amrex::Real loss = static_cast<amrex::Real>(w - h1);
            if (!(loss > amrex::Real(0.0))) { return; }
            remove_water(c, ff[box_no](i, j, k), fc[box_no](i, j, k), i, j, k, h, loss, dry_depth);
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace flume
