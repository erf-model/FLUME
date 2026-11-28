#ifndef FLUME_FLOODMAPS_H_
#define FLUME_FLOODMAPS_H_

#include "FLUME_Checkpoint.H"
#include "FLUME_Inputs.H"
#include "FLUME_Source.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <vector>

namespace flume {

/**
 * The flood maps (flume.flood_maps = 1): per cell, the largest depth, the
 * largest speed and the largest depth times speed the run has reached, and
 * the time water first arrived (Docs/sphinx_doc/Running.rst).
 *
 * The maps follow the state at the start and at the end of every step; the
 * stages inside a step are not seen. A cell is wet when its depth is above
 * flume.dry_depth; the speed of a dry cell is 0. In a wet cell the depth
 * times the speed is the magnitude of the discharge, sqrt(hu^2 + hv^2).
 * Water arrives in a cell at the first of those times at which its depth is
 * above flume.flood_maps.arrival_depth; a cell water has not reached holds
 * -1.
 */
class FloodMaps
{
public:
    FloodMaps (const Inputs& in, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    //! m_alias holds views of m_maps, and plot_fields() pointers into m_alias.
    FloodMaps (const FloodMaps&) = delete;
    FloodMaps& operator= (const FloodMaps&) = delete;
    FloodMaps (FloodMaps&&) = delete;
    FloodMaps& operator= (FloodMaps&&) = delete;

    //! Start the maps from the state `state` at time `time`, s: the maxima
    //! are its own, and only a cell wet above arrival_depth has an arrival
    //! time (`time`).
    void start (const amrex::MultiFab& state, double time);
    //! Take the state `state` at time `time`, s, into the maps.
    void update (const amrex::MultiFab& state, double time);

    //! max_depth (m), max_speed (m/s), max_depth_speed (m2/s),
    //! arrival_time (s).
    [[nodiscard]] std::vector<PlotField> plot_fields () const;

    //! The four maps, as the field "flood_maps". Collective.
    void write_checkpoint (CheckpointWriter& chk) const;
    //! Read the maps back; aborts if the checkpoint holds none (it was
    //! written with flume.flood_maps = 0). Collective.
    void read_checkpoint (const CheckpointReader& chk);

private:
    //! max_depth, max_speed, max_depth_speed, arrival_time; no ghost cells
    amrex::MultiFab m_maps;
    amrex::Real m_dry_depth = amrex::Real(0.0);
    amrex::Real m_arrival_depth = amrex::Real(0.0);
    //! One-component aliases of the components of m_maps, for the plotfiles
    std::vector<amrex::MultiFab> m_alias;
};

//! The step of FloodMaps::update on the valid cells of `maps` (four
//! components, as above).
void update_flood_maps (amrex::MultiFab& maps, const amrex::MultiFab& state,
                        amrex::Real dry_depth, amrex::Real arrival_depth, amrex::Real time);

} // namespace flume

#endif
