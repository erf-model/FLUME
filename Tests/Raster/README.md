# Raster input tests

## Purpose

Check that FLUME reads ESRI ASCII grids (`.asc`) for the bed
(`flume.bed.file`), building heights or footprints (`flume.buildings.file`,
`flume.buildings.height`) and the Manning coefficient
(`flume.friction.manning_file`), and that every malformed file stops the run
with a message naming the key and the file. The format and the sampling are
described in `Docs/sphinx_doc/InputFiles.rst` (Raster files) and
`Docs/sphinx_doc/theory/Terrain.rst`.

## Files

| File | What it is |
| --- | --- |
| `make_rasters.py` | Writes every `.asc` file here and in `bad/`, byte for byte as committed. |
| `check_raster.py` | Reads the plotfiles of a run and checks the bed (plane, or a file interpolated bilinearly), building cells, dry or overtopped buildings, the lake at rest, conservation and that water reached the streets. |
| `inputs_bed` | 48 x 32 cells of 1 m, no step: the bed alone. |
| `inputs_streets` | Dam break on 40 x 40 cells of 1 m with walls, run to 8 s. |
| `inputs_lake` | The same domain with the surface at 4 m and no dam, run to 4 s. |
| `inputs_manning` | The friction channel of `Tests/Friction/inputs` without `flume.friction.manning`. |
| `bed_plane.asc` | z = 2 + 0.03 x - 0.02 y + 0.001 x y at the centres of 5 m cells (`xllcenter`), larger than the domain, one NODATA cell outside the cells read. |
| `bed_coarse.asc` | 6 x 8 cells of 8 x 4 m (`dx`, `dy`) holding the domain exactly. |
| `streets_bed.asc` | A bed falling to the east: 0.075 m at the western raster centres, 0.025 m at the eastern ones. |
| `streets_heights.asc` | Blocks of 3, 5 and 8 m on a 2 m raster with its corner at (-0.75, -0.75), so a grid centre is a quarter or three quarters of a raster cell from a raster edge; NODATA -1 (not used). |
| `streets_mask.asc` | The same footprints as 0/1. |
| `manning_rows.asc` | One Manning coefficient per row: 0.02, 0.03, 0.045, 0.06 from the south. |
| `manning_uniform.asc` | 0.03 everywhere. |
| `manning_ties.asc` | Rows of 0.1 m from y = -2 m: every grid centre on a raster edge; the row above holds the grid row's coefficient, the row below 0.09. |
| `bad/*.asc` | One malformed file per check of the reader. |

The Manning tests use `Tests/check_friction.py --manning-rows`.

## Data and acknowledgements

Every grid is synthetic, made by `make_rasters.py`. No external data are used
or redistributed, and nothing needs to be downloaded.

## How to run

From the build directory:

```bash
ctest -R raster --output-on-failure
```

By hand, from an empty directory (file keys take absolute paths, or paths
relative to the directory the run starts in):

```bash
<build>/flume <repo>/Tests/Raster/inputs_streets \
    flume.bed.file=<repo>/Tests/Raster/streets_bed.asc \
    flume.buildings.file=<repo>/Tests/Raster/streets_heights.asc
python3 <repo>/Tests/Raster/check_raster.py str \
    --bed <repo>/Tests/Raster/streets_bed.asc \
    --buildings <repo>/Tests/Raster/streets_heights.asc \
    --dry-buildings --wet-beyond 14
```

To write the grids again: `python3 make_rasters.py` (then `git diff` shows
no change).

## Expected Results

| Test | Result |
| --- | --- |
| `raster_bed_plane` | largest bed error 8.9e-16 m |
| `raster_bed_coarse` | largest bed error 0 |
| `raster_streets` | 256 building cells, all dry; sum of the depth unchanged; 0.381 m in the streets beyond x = 14 m at 8 s |
| `raster_lake_at_rest_order_1`, `_order_2` | 160 dry cells; depth error 8.9e-16 m and discharge 1.6e-14 m^2/s at most |
| `raster_streets_overtop` | water over the 0.3 m buildings, 0.189 m deep at most; sum of the depth unchanged |
| `raster_manning_ties` | each row follows its own coefficient: every centre took the row above |
| `raster_manning_rows` | largest relative discharge error 5.63e-4 (order 1), 1.91e-4 (order 2) |
| `raster_manning_uniform_file` | identical to `flume.friction.manning = 0.03` |

The numbers are what the runs give, not the bounds the checks use.

## CTest

`raster_bed_plane`, `raster_bed_coarse`, `raster_bed_relative_path`,
`raster_streets`, `raster_streets_second_order`, `raster_streets_mask`,
`raster_streets_small_boxes`, `raster_lake_at_rest_order_1`,
`raster_lake_at_rest_order_2`, `raster_streets_overtop`,
`raster_manning_rows`, `raster_manning_rows_second_order`,
`raster_manning_ties`, `raster_manning_uniform_value`,
`raster_manning_uniform_file`, and the `abort_raster_*` tests. All run on one
rank. See `Docs/sphinx_doc/Verification.rst`.
