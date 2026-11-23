#ifndef FLUME_SOURCE_H_
#define FLUME_SOURCE_H_

#include "FLUME_Forcing.H"
#include "FLUME_Inputs.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <memory>
#include <string>
#include <vector>

namespace flume {

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
};

/**
 * The sources the inputs switch on, in the order they act: rain when
 * `forcing` is not null (flume.rain.type), then Manning friction
 * (flume.friction.manning > 0 or flume.friction.manning_file), which so
 * sees the depth after the rain. `forcing` must outlive the sources. Prints
 * the names, or "none", at start-up.
 */
std::vector<std::unique_ptr<Source>> make_sources (const Inputs& in, const amrex::Geometry& geom,
                                                   const amrex::BoxArray& ba,
                                                   const amrex::DistributionMapping& dm,
                                                   const Forcing* forcing);

} // namespace flume

#endif
