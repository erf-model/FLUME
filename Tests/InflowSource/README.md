# Inflows inside the domain: points and lines

## Purpose

Check the inflows of `flume.inflows.*`: a hydrograph added as depth to the
cell that holds a point, or to the cells a line crosses, each taking the part
of the line's length inside it, with no momentum. One step from a dry domain
checks the parts of a line against an exact clipping to each cell, and the
time step against the closed-form root of the rain-limited cubic, also for
a line along a cell face, a line through corners of the grid, and a point
with rain as well. A closed
basin checks that the volume added is the integral of each hydrograph, with
both interpolations, a late start and rain as well. A tiny inflow into deep
water checks the compensated sum. A steep channel fed by a line across it
reaches the normal depth. Further tests pin the independence of the
decomposition, restarts on the same and another rank count, and that every
malformed input stops the run with a message naming the key. The model is
in `Docs/sphinx_doc/theory/SourceTerms.rst` (Inflows inside the domain) and
the inputs in `Docs/sphinx_doc/Inputs.rst`.

## Files

| File | What it is |
| --- | --- |
| `check_inflow_source.py` | Reads the last plotfile of a one-step run from a dry, flat, closed domain fed by one steady inflow (`--q Q`, `--line X0 Y0 X1 Y1` or `--point X Y`) and checks the step against the root of the rain-limited cubic (found by bisection), each cell's depth against its part of the inflow times Q dt over the cell area (the parts from a Liang-Barsky clipping of the line to each cell), and that no cell holds a discharge. The cells are half open: a line along a face is in the cell above or to the right, and a cell the line only touches at a corner is not one of its cells. `--cells N` checks the count of cells, `--rain MMH` adds steady rain to the rate of the cubic and to each depth, and the tolerance is 1e-12 in double precision, 1e-5 in single. `--deep H0` checks instead that each cell holds H0 plus its part of Q t / A, within 4 epsilons of H0, at any step. |
| `inputs_dry` | 8 x 10 cells of 2 x 1 m in boxes of 4, closed, flat and dry, one step. The tests add the inflow. |
| `inputs_basin` | 20 x 20 cells of 2 m in 9 boxes, closed and flat, still water 0.5 m deep, to 300 s. Point `a` at (11, 21) m, in cell (5, 10); line `b` from (3.3, 5.1) to (25.7, 31.9) m, across 25 cells. |
| `inputs_deep` | One closed cell of 1 m, still water 1 m deep, a point inflow of 1e-15 m3/s, 20000 steps. |
| `inputs_channel` | 200 x 2 cells of 1 m in 4 boxes, periodic across, bed slope 0.02, Manning 0.03, dry at the start, to 600 s. Side xlo is a wall; line `up` crosses the channel at x = 0.5 m with 4 m3/s (1 m2/s). Side xhi is an outflow segment. |
| `hydrograph_triangle.txt` | 0 to 2 m3/s in 60 s, back to 0 at 200 s: 200 m3 linear in time, 280 m3 step in time. |
| `hydrograph_steps.txt` | 1 m3/s for 50 s, 3 m3/s for 70 s, then none: 260 m3 step in time, 205 m3 linear in time. |
| `hydrograph_late.txt` | 2 m3/s from 50 s on, the first time after the start of the run. |
| `hydrograph_half.txt`, `hydrograph_steady.txt`, `hydrograph_tiny.txt`, `hydrograph_zero.txt` | Steady 0.5, 4, 1e-15 and 0 m3/s. |
| `InflowSourceTests.cmake` | The CTest definitions, included from `Tests/CMakeLists.txt`. |

The channel is compared with the normal depth by
`Tests/Inflow/check_inflow.py`, the water balance by `Tests/check_balance.py`
and the restarts by `Tests/Restart/check_restart.py`. The abort tests of an
inflow in a building cell take the buildings of
`Tests/Raster/streets_heights.asc`, and that of a negative discharge the bad
file `Tests/Rain/bad/negative.txt`.

## Data and acknowledgements

No data are downloaded, and nothing needs to be. The decks and the
hydrographs are synthetic, written for these tests.

## How to run

From the build directory:

```bash
ctest -R inflow_source --output-on-failure
```

By hand, from an empty directory:

```bash
mpiexec -n 4 <build>/flume <repo>/Tests/InflowSource/inputs_dry flume.inflows.names=l "flume.inflows.l.line=1.3 0.7 14.2 8.9" flume.inflows.l.file=<repo>/Tests/InflowSource/hydrograph_half.txt
python3 <repo>/Tests/InflowSource/check_inflow_source.py dry --q 0.5 --line 1.3 0.7 14.2 8.9
```

```bash
mpiexec -n 4 <build>/flume <repo>/Tests/InflowSource/inputs_basin flume.inflows.a.file=<repo>/Tests/InflowSource/hydrograph_triangle.txt flume.inflows.b.file=<repo>/Tests/InflowSource/hydrograph_steps.txt
python3 <repo>/Tests/check_balance.py --expect-final inflow_a=200 inflow_b=260 stored=1260
```

## Expected Results

| Test | Result |
| --- | --- |
| `inflow_source_line_first_step` | 16 cells, largest part 0.121951; step 0.71564988791 s, the closed form to all printed digits; depths within 1.6e-16 of the largest |
| `inflow_source_line_reversed` | the same line from its other end: the same step, depths within 0 of the largest |
| `inflow_source_line_edge_to_edge` | 8 cells of 1/8 each; step 0.70980978329 s; depths within 3.1e-16 |
| `inflow_source_point_first_step` | cell (4, 6); step 0.355696995613 s; depth within 4.7e-16 |
| `inflow_source_line_on_face` | 10 x 10 cells on 3 m, the line along x = 0.6 m: cells (2, 1) and (2, 2), to the right of the face; step 0.0590025868352 s; depths within 1.7e-16 |
| `inflow_source_line_through_corners` | 30 x 45 cells of 0.1 m, the diagonal from (3, 0) to (0, 4.5): 60 cells (70 with the slivers of rounding at the corners); step 0.0384373602233 s |
| `inflow_source_point_with_rain` | the point and 360000 mm/h (0.1 m/s) of rain: step 0.318034064141 s, the root with the sum of the rates; depths within 1.3e-16 |
| `inflow_source_basin_volume_np4` | inflow_a 200 m3 and inflow_b 260 m3 exactly; stored 1260 m3 to 2.3e-13 |
| `inflow_source_basin_volume_swapped` | 280 m3 and 205 m3 with the other interpolations, stored 1285 m3 |
| `inflow_source_basin_late_start` | 600 m3 from a file that starts at 50 s, and the note; stored 1400 m3 |
| `inflow_source_basin_with_rain` | rain 8 m3 to 7e-15, inflows 200 and 260 m3, stored 1268 m3 |
| `inflow_source_deep_water` | 1 m plus 1.596e-12 m after 20000 steps, exact; without the compensated sum the depth stays 1 m |
| `inflow_source_basin_boxes_4_np3` | 20 plotfiles of 8000 cells, every value that of 2 ranks in 9 boxes to the last bit |
| `inflow_source_restart_np2`, `_np2_np3` | from step 100: the plotfiles to the last bit; the log byte for byte on 2 ranks and within 3.6e-16 of each column's scale on 3 |
| `inflow_source_normal_depth_order1` | the normal depth, 0.394424 m, from 40 m to 190 m within 0.77 % (bound 1 %) |
| `inflow_source_normal_depth_order2` | within 0.17 % in the depth and 0.29 % in the discharge (bound 0.5 %) |
| `abort_inflow_source_restart_moved`, `_line` | a restart from `inflow_source_basin_reference` with point `a` in cell (6, 10), or with the end of line `b` moved, stops |

The normal-depth errors are those of the discharge segment in
`inflow_normal_depth_steep_order*` (`Tests/Inflow`): the water that enters
with no momentum is set moving by the bed within the first 40 m. The
numbers are what the runs give in double precision, not the bounds the
checks use (`InflowSourceTests.cmake`).

## CTest

`inflow_source_line_first_step`, `inflow_source_line_reversed`,
`inflow_source_line_edge_to_edge`, `inflow_source_point_first_step`,
`inflow_source_line_on_face`, `inflow_source_line_through_corners`,
`inflow_source_point_with_rain`,
`inflow_source_basin_volume_np4`, `inflow_source_basin_volume_swapped`,
`inflow_source_basin_late_start`, `inflow_source_basin_with_rain`,
`inflow_source_deep_water`, `inflow_source_basin_reference`,
`inflow_source_basin_boxes_4_np3`, `inflow_source_restart_np2`,
`inflow_source_restart_np2_np3`, `inflow_source_normal_depth_order1`,
`inflow_source_normal_depth_order2`, `abort_inflow_source_restart_moved`,
`abort_inflow_source_restart_line`, and the `abort_inflow_source_*` tests.
They run on 1 to 4 ranks.
