#ifndef FLUME_SOURCE_H_
#define FLUME_SOURCE_H_

#include "FLUME_Checkpoint.H"
#include "FLUME_Forcing.H"
#include "FLUME_Inputs.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace flume {

//! A field a source adds to the plotfiles: its name and the field (one
//! component, the boxes of the state, no ghost cells needed).
using PlotField = std::pair<std::string, const amrex::MultiFab*>;

//! What a source may read besides the state. The references live as long as
//! the call to Source::apply.
struct SourceContext
{
    const amrex::MultiFab& bed;   //!< bed elevation, m, with filled ghost cells
    const amrex::Geometry& geom;  //!< the grid
    double time;                  //!< time at the start of the step, s
    double time_end;              //!< time at the end of the step, s
};

/**
 * A source term: a change of the state that is not a flux between cells.
 *
 * Each source acts once per step, after the fluxes and the bed slope of all
 * the stages of the step, on the valid cells of the state at the end of the
 * step (Docs/sphinx_doc/theory/SourceTerms.rst). Sources act one after
 * another in the order of the list that make_sources() returns. A source
 * keeps any field it needs (a coefficient map, a running total) and switches
 * itself on from its own inputs.
 *
 * The bed-slope term is not a source in this sense: it is part of the
 * explicit stage, where it balances the flux of still water exactly
 * (FLUME_Advance.H).
 *
 * The device kernels of a source live in free functions, not in member
 * functions, so no kernel captures `this`.
 */
class Source
{
public:
    Source () = default;
    virtual ~Source () = default;
    Source (const Source&) = delete;
    Source& operator= (const Source&) = delete;
    Source (Source&&) = delete;
    Source& operator= (Source&&) = delete;

    //! The name printed at start-up and in messages, for example
    //! "manning_friction".
    [[nodiscard]] virtual std::string name () const = 0;

    //! One step of `dt`, s, on the valid cells of `state`. The ghost cells of
    //! `state` are not filled and must not be read.
    virtual void apply (amrex::MultiFab& state, const SourceContext& ctx, amrex::Real dt) = 0;

    //! The column of this source in the water-balance log
    //! (FLUME_WaterBalance.H), or empty when the source neither adds nor
    //! removes water (friction).
    [[nodiscard]] virtual std::string balance_name () const { return {}; }

    //! Whether the column counts water removed (a sink, such as
    //! infiltration) instead of water added (rain).
    [[nodiscard]] virtual bool balance_removes () const { return false; }

    //! The volume, m3, this source has added to the surface water of this
    //! rank's cells since the start of the run, in double precision; below 0
    //! for water removed. Called on every rank; the caller sums the ranks.
    [[nodiscard]] virtual double balance_volume () const { return 0.0; }

    //! The fields this source adds to the plotfiles, after the fields of the
    //! state, in this order (none by default).
    [[nodiscard]] virtual std::vector<PlotField> plot_fields () const { return {}; }

    //! Write what this source carries from step to step (none by default),
    //! under keys and field names that start with name(). Collective.
    virtual void write_checkpoint (CheckpointWriter& /*chk*/) const {}

    //! Read back what write_checkpoint wrote. Collective.
    virtual void read_checkpoint (const CheckpointReader& /*chk*/) {}
};

/**
 * The sources the inputs switch on, in the order they act: rain when
 * `forcing` is not null (flume.rain.type), then infiltration
 * (flume.infiltration.type), which so sees the rain of the step, then
 * Manning friction (flume.friction.manning > 0 or
 * flume.friction.manning_file), which sees the depth after both.
 * `forcing` must outlive the sources. `buildings` (may be null) is the
 * building height of each cell (flume.buildings.file): building cells
 * infiltrate nothing. Prints the names, or "none", at start-up.
 */
std::vector<std::unique_ptr<Source>> make_sources (const Inputs& in, const amrex::Geometry& geom,
                                                   const amrex::BoxArray& ba,
                                                   const amrex::DistributionMapping& dm,
                                                   const Forcing* forcing,
                                                   const amrex::MultiFab* buildings);

} // namespace flume

#endif
