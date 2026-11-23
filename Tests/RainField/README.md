# Gridded rain and moving rain cell tests

## Purpose

Check the rain that FLUME takes from `flume.rain.type = gridded` (a time
series of ESRI ASCII grids of the rate in mm/h, listed in an index file) and
`flume.rain.type = moving_cell` (a Gaussian rain cell that moves at a
constant velocity): that each solver cell receives the integral of the input
over its area and over the step, that a closed domain receives the integral
of the input over the domain and over the run to round-off, and that every
malformed index, grid or cell input stops the run with a message naming the
key, the file and the line or entry. The model is in
`Docs/sphinx_doc/Model.rst` (Rain that varies in space) and the inputs in
`Docs/sphinx_doc/Inputs.rst` (Gridded rain index files).

## Files

| File | What it is |
| --- | --- |
| `make_rain_grids.py` | Writes every grid and index file below; running it again writes them byte for byte. |
| `check_rain_field.py` | Computes the rain from the grid files or the cell's numbers, independently of the solver (areas shared by grid and solver cells; time split at the index times; error functions and adaptive Simpson quadrature for the cell), and checks every cell after one step (`--cells`) or the water in the domain in every plotfile (`--volume`). |
| `inputs_flat` | 24 x 16 cells of 2 m in 6 boxes, walls, flat dry bed, gridded rain to 1200 s. |
| `grid_a.asc`, `grid_b.asc` | 18 x 15 cells of 3.7 x 2.9 m from (-5.3, -4.1): a smooth pattern and a peak near (30, 20), plus noise. `grid_a` has one NODATA cell outside every domain. |
| `grid_c.asc` | 50 x 32 cells of 1.3 m given by the centre of the lower-left cell; a linear ramp plus noise. |
| `grid_aligned.asc` | 17 x 17 cells of 1 m from (-1, -1), NODATA outside 0..14 m. |
| `grid_zero.asc`, `grid_burst.asc` | 11 x 8 cells of 10 m from (-20, 0): no rain, and about 3590 mm/h. |
| `index_steps.txt` | `grid_a`, `grid_b`, `grid_c` at 0, 300 and 900 s; has comments. |
| `index_offset.txt` | `grid_a` at -300 s and `grid_b` at 600 s. |
| `index_fine.txt` | `grid_c` alone. |
| `index_aligned.txt` | `grid_aligned` alone. |
| `index_burst.txt` | No rain, the burst from 100 s to 101 s, then none. |
| `bad/` | One malformed index or grid per check of the reader. |
| `RainFieldTests.cmake` | The CTest definitions, included from `Tests/CMakeLists.txt`. |

The test without a stop time uses `Tests/Skeleton/inputs`.

## Data and acknowledgements

Every grid and index file is synthetic, written by `make_rain_grids.py` for
these tests. No external data are used or redistributed, and nothing needs
to be downloaded.

## How to run

From the build directory:

```bash
ctest -R "rain_gridded|rain_cell|abort_rain_gridded|abort_rain_cell" --output-on-failure
```

By hand, from an empty directory:

```bash
mpiexec -n 2 <build>/flume <repo>/Tests/RainField/inputs_flat \
    flume.rain.file=<repo>/Tests/RainField/index_steps.txt
python3 <repo>/Tests/RainField/check_rain_field.py flat \
    --index <repo>/Tests/RainField/index_steps.txt --volume --tol-eps 4
```

```bash
<build>/flume <repo>/Tests/RainField/inputs_flat flume.rain.type=moving_cell \
    flume.rain.cell.peak_rate=120 flume.rain.cell.sigma=6 \
    "flume.rain.cell.center=40 26" flume.max_step=1 flume.plot_int=1
python3 <repo>/Tests/RainField/check_rain_field.py flat --cell 120 6 40 26 0 0 --cells 1
```

## Expected Results

| Test | Result |
| --- | --- |
| `rain_gridded_cells_linear_np1`, `_np2` | every cell within 1e-15 of the largest depth (relative) |
| `rain_gridded_cells_fine` | every cell within 5e-16 (relative) |
| `rain_gridded_cells_aligned_nodata_ring` | runs (no NODATA read through round-off); every cell within 3.7e-15 (relative) |
| `rain_gridded_volume_step`, `_linear_shifted` | water equal to the rain within 7e-16 of it at step 679 |
| `rain_gridded_dry_burst_no_stop_time` | finite steps before the burst; water equal to the rain within 3e-16 |
| `rain_cell_at_rest_np1`, `_np2` | every cell within 6e-16 (relative) |
| `rain_cell_fast_step` | every cell within 1.1e-13 (relative); 1.3e-2 with one Gauss-Legendre interval |
| `rain_cell_volume` | water equal to the rain within 1.4e-16 of it |
| `rain_cell_zero_peak_long_step` | one step of 1e7 s, no water |

The numbers are what the runs give in double precision, not the bounds the
checks use (`RainFieldTests.cmake`, `Docs/sphinx_doc/Testing.rst`).

## CTest

`rain_gridded_cells_linear_np1`, `rain_gridded_cells_linear_np2`,
`rain_gridded_cells_fine`, `rain_gridded_cells_aligned_nodata_ring`,
`rain_gridded_volume_step`,
`rain_gridded_volume_linear_shifted`, `rain_gridded_dry_burst_no_stop_time`,
`rain_cell_at_rest_np1`, `rain_cell_at_rest_np2`, `rain_cell_fast_step`,
`rain_cell_volume`, `rain_cell_zero_peak_long_step`, and the `abort_rain_gridded_*` and `abort_rain_cell_*`
tests. They run on 1 or 2 ranks. See `Docs/sphinx_doc/Testing.rst`.
