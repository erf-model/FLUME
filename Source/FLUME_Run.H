#ifndef FLUME_RUN_H_
#define FLUME_RUN_H_

#include "FLUME_Inputs.H"

namespace flume {

/**
 * The skeleton run: build the grid the inputs describe and write one plotfile
 * that holds the water depth, zero everywhere. There is no solver yet.
 *
 * The grid is two-dimensional: nx by ny cells and one cell in z. The z extent
 * of the domain is 0 to 1 m and carries no meaning.
 */
void run (const Inputs& in);

} // namespace flume

#endif
