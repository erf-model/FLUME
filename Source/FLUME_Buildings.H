#ifndef FLUME_BUILDINGS_H_
#define FLUME_BUILDINGS_H_

#include <AMReX_INT.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

namespace flume {

/**
 * Buildings as raised bed (flume.buildings.file): a cell whose building
 * height is above 0 has its bed raised, and water enters it only over the
 * top. No water is stored inside a building and no load on its walls is
 * computed.
 *
 * `height` has the boxes of `bed`, one component, values >= 0 in the valid
 * cells. Only the valid cells of `bed` change; its ghost cells are filled
 * afterwards. Each returns the number of cells raised, over all ranks.
 */

//! Raise the bed of each cell by its height (0 leaves it).
amrex::Long add_building_heights (amrex::MultiFab& bed, const amrex::MultiFab& height);

//! Raise the bed of each cell whose height is above 0 by `fixed` (> 0): the
//! height field is a footprint mask (flume.buildings.height).
amrex::Long add_building_footprints (amrex::MultiFab& bed, const amrex::MultiFab& height,
                                     amrex::Real fixed);

//! Number of valid cells, over all ranks, whose first component is below 0.
amrex::Long count_below_zero (const amrex::MultiFab& mf);

} // namespace flume

#endif
