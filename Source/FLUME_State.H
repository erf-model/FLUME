#ifndef FLUME_STATE_H_
#define FLUME_STATE_H_

namespace flume {

/**
 * Components of the conserved state, all cell-centred.
 *
 * The bed elevation is a separate one-component MultiFab on the same boxes.
 * Bed and water surface are elevations on one absolute vertical datum, in
 * metres; the depth is their difference.
 */
namespace State {
enum : int {
    h = 0, //!< water depth, m
    hu,    //!< discharge per unit width in x, m2/s
    hv,    //!< discharge per unit width in y, m2/s
    ncomp
};
} // namespace State

} // namespace flume

#endif
