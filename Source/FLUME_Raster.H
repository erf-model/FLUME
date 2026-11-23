#ifndef FLUME_RASTER_H_
#define FLUME_RASTER_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace flume {

/**
 * The part of a raster of one value per cell, read from an ESRI ASCII grid
 * (.asc), that the cells of a domain read.
 *
 * The file is a header of keyword-value pairs, then nrows lines of ncols
 * values, the northernmost row first:
 *
 *   ncols         <int >= 1>
 *   nrows         <int >= 1>
 *   xllcorner     <x of the lower-left corner>   (or xllcenter: of its centre)
 *   yllcorner     <y of the lower-left corner>   (or yllcenter)
 *   cellsize      <cell edge, > 0>               (or dx and dy, as GDAL writes them)
 *   NODATA_value  <value>                        (optional)
 *
 * Keywords are not case-sensitive. The values are stored with the southern
 * row first: value(i, j) is column i from the west and row j from the south
 * of the part kept. A file has no datum and no unit: they must be those of
 * the inputs (flume.prob_lo and flume.prob_hi in metres on the run's
 * horizontal coordinates; elevations in metres on the run's vertical datum).
 */
struct Raster
{
    std::string key;    //!< the input that names the file, for messages
    std::string path;   //!< the file
    int nx = 0;         //!< columns kept
    int ny = 0;         //!< rows kept
    double x0 = 0.0;    //!< x of the west edge of the part kept, m
    double y0 = 0.0;    //!< y of the south edge of the part kept, m
    double dx = 0.0;    //!< column width, m
    double dy = 0.0;    //!< row height, m
    std::vector<double> value;  //!< nx * ny values, index i + nx * j; 0 where NODATA
    std::vector<char> nodata;   //!< 1 where the file holds NODATA_value
};

/**
 * Read an ESRI ASCII grid and keep the columns and rows the cells of `geom`
 * can read, with one more on each side where there is one. The I/O rank
 * reads and parses the file and sends the part kept to every rank, so only
 * the I/O rank holds the file. Every rank aborts with the same message,
 * naming `key` and the file, when the file cannot be opened or is not a
 * regular file, a header keyword is unknown, repeated or missing, a header
 * value is out of range, a value is not a finite number, the file holds more
 * or fewer than ncols x nrows values, the raster does not hold the domain
 * (its edges may fall on the domain's, within a millionth of a raster cell),
 * the part kept has more than 2^31 - 1 cells, or a value kept is beyond the
 * range of amrex::Real.
 */
Raster read_raster (const std::string& key, const std::string& path, const amrex::Geometry& geom);

/**
 * Fill component `comp` of the valid cells of `dst` from `r` by bilinear
 * interpolation between the raster's cell centres, at the centre of each
 * cell of `geom`. Within half a raster cell of the raster's edge the values
 * of the edge cells hold (no extrapolation). Aborts, naming the key, when a
 * cell reads a NODATA value with a weight above zero.
 */
void sample_bilinear (amrex::MultiFab& dst, int comp, const Raster& r,
                      const amrex::Geometry& geom);

/**
 * Fill component `comp` of the valid cells of `dst` with the value of the
 * raster cell that holds the centre of each cell of `geom`. A centre on the
 * edge between two raster cells, to within 16 epsilons of amrex::Real of
 * the position in raster cells, takes the one above it in x or y. Aborts,
 * naming the key, when a cell reads a NODATA value.
 */
void sample_nearest (amrex::MultiFab& dst, int comp, const Raster& r,
                     const amrex::Geometry& geom);

/**
 * Fill component `comp` of the valid cells of `dst` with the mean over each
 * cell of `geom` of the raster, taken as constant on each raster cell: the
 * sum of the raster values weighted by the area each raster cell shares with
 * the cell, divided by the cell's area. The integral of the result over the
 * domain is then the integral of the raster over the domain, to round-off,
 * whatever the two grids are. The raster's edge cells extend beyond its edge
 * (within the tolerance read_raster allows). The edges of the cells, in
 * raster cells, are computed by one formula from their index, so neighbours
 * share each edge exactly and the weights of a raster cell sum to the area
 * it shares with the domain. Aborts, naming the key, when a cell reads a
 * NODATA value with an area above zero.
 */
void sample_conservative (amrex::MultiFab& dst, int comp, const Raster& r,
                          const amrex::Geometry& geom);

} // namespace flume

#endif
