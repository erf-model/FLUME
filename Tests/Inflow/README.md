# Inflow hydrograph tests

## Purpose

Check the inflow that FLUME takes through a boundary segment of type
`discharge` (`flume.bc.<name>.type = discharge`): a discharge in time from a
file (`flume.bc.<name>.file`), held or linear between entries, spread evenly
over the width of the segment. The tests check that the depth at the face
follows the characteristics for subcritical inflow and the critical depth
for supercritical inflow; that steady inflow down a steep slope reaches the
normal depth; that the water that enters equals the integral of the
hydrograph; that a restart continues the run to the last bit; and that every
malformed input or file stops the run with a message naming the key. The
model is in `Docs/sphinx_doc/theory/Boundaries.rst` (Inflow hydrograph) and
the inputs in `Docs/sphinx_doc/Inputs.rst` (Boundaries).

## Files

| File | What it is |
| --- | --- |
| `check_inflow.py` | Reads the last plotfile of a channel run and compares the depth and the discharge with a closed form: still uniform flow (`--mode uniform`), the plateau of a rarefaction (`plateau`), the fan from the critical depth onto a dry bed (`fan`) or the normal depth (`normal`). |
| `inputs_channel` | 400 x 2 cells over 200 x 4 m, periodic across, flat dry bed, no friction, segment `in` (discharge, xlo) and segment `out` (outflow, xhi), second order. No stop time: each test sets it. |
| `inputs_basin` | 20 x 20 cells over 40 x 40 m in boxes of 8, walls, flat bed, 0.5 m of still water, to 300 s. Segment `a` on xlo from 10 m to 30 m, segment `b` on yhi from 0 m to 10 m (held between entries). |
| `hydrograph_steady.txt` | 4 m3/s (1 m2/s over the 4 m channel). |
| `hydrograph_half.txt` | 2 m3/s (0.5 m2/s). |
| `hydrograph_triangle.txt` | 0, 2 and 0 m3/s at 0, 60 and 200 s: 200 m3 linear, 280 m3 held. |
| `hydrograph_steps.txt` | 1 m3/s from 0 s, 3 from 50 s, 0 from 120 s: 260 m3 held, 205 m3 linear. |
| `hydrograph_late.txt` | 2 m3/s from 50 s, with no line at 0 s. |
| `hydrograph_zero.txt` | 0 m3/s. |
| `InflowTests.cmake` | The CTest definitions, included from `Tests/CMakeLists.txt`. |

The abort tests reuse the malformed files of `Tests/Rain/bad`, which the same
reader takes, and the 8 x 8 m deck of the other abort tests.

## Data and acknowledgements

Every hydrograph is synthetic, written for these tests. No external data are
used or redistributed, and nothing needs to be downloaded.

## How to run

From the build directory:

```bash
ctest -R inflow --output-on-failure
```

By hand, from an empty directory (`flume.bc.<name>.file` takes an absolute
path, or a path relative to the directory the run starts in):

```bash
mpiexec -n 2 <build>/flume <repo>/Tests/Inflow/inputs_channel \
    flume.stop_time=30 flume.bc.in.file=<repo>/Tests/Inflow/hydrograph_steady.txt
python3 <repo>/Tests/Inflow/check_inflow.py chan --side xlo --q 1 --mode fan \
    --time 30 --range 0 60 --rtol 0.01
```

```bash
mpiexec -n 4 <build>/flume <repo>/Tests/Inflow/inputs_basin \
    flume.bc.a.file=<repo>/Tests/Inflow/hydrograph_triangle.txt \
    flume.bc.b.file=<repo>/Tests/Inflow/hydrograph_steps.txt
python3 <repo>/Tests/check_balance.py --expect-final in_a=200 in_b=260 stored=1260
```

## Expected Results

| Test | Result |
| --- | --- |
| `inflow_uniform_subcritical_order1`, `_order2` | uniform after 3419 steps: 250 eps in the depth, 1030 eps in the discharge |
| `inflow_rarefaction_*_order1` | plateau depth 0.869019 m: 0.022 % in the depth, 0.068 % in the discharge, on each side |
| `inflow_rarefaction_*_order2` | 0.006 % and 0.028 % |
| `inflow_dry_fan_order1`, `_order2` | 1.07 % and 0.75 % in the depth, next to the side |
| `inflow_normal_depth_steep_order1` | normal depth 0.394424 m: 0.77 % high (first-order bed source) |
| `inflow_normal_depth_steep_order2` | 0.21 % in the depth, 0.35 % in the discharge |
| `inflow_basin_volume_np4` | in_a 200, in_b 260, stored 1260 m3 (difference 2.3e-13) |
| `inflow_basin_volume_swapped` | in_a 280, in_b 205, stored 1285 m3 |
| `inflow_basin_late_start` | in_a 600 m3; the run prints a note on the first time |
| `inflow_basin_lake_at_rest` | 0.5 m and no discharge in every cell |
| `inflow_dry_no_stop_time` | three steps on a dry domain without a stop time |
| `inflow_restart_np2` | plotfiles to the last bit, log byte for byte |
| `inflow_restart_np2_np1` | plotfiles to the last bit, log within 1e-14 of each column |

The numbers are what the runs give in double precision, not the bounds the
checks use (`InflowTests.cmake`, `Docs/sphinx_doc/Verification.rst`).

## CTest

`inflow_uniform_subcritical_order1`, `_order2`,
`inflow_rarefaction_{xlo,xhi,ylo,yhi}_order{1,2}`,
`inflow_dry_fan_order1`, `_order2`, `inflow_normal_depth_steep_order1`,
`_order2`, `inflow_basin_volume_np4`, `inflow_basin_volume_swapped`,
`inflow_basin_late_start`, `inflow_basin_lake_at_rest`,
`inflow_dry_no_stop_time`, `inflow_restart_reference`,
`inflow_restart_np2`, `inflow_restart_np2_np1`, and the `abort_inflow_*` tests. They run on 1, 2 or 4
ranks. See `Docs/sphinx_doc/Verification.rst`.
