# FLUME as ERF's flood engine

## Purpose

Check FLUME run by ERF as its flood engine (`erf.flood.engine = flume`):
that ERF builds, starts, advances, checkpoints, restarts and plots it, and
that the flood it gives is the stand-alone run's, to the last bit, when each
ERF step is one FLUME step. Also check that the FLUME keys ERF sets stop the
run with a message naming ERF's key. The coupling is described in
`Docs/sphinx_doc/Running.rst` (Inside ERF).

## Files

| File | What it is |
| --- | --- |
| `flume_plane` | FLUME's half of the case, read by both runs (`FILE=<path>/flume_plane`): rain of 50 mm/h on a dry plane sloping 1 cm per metre towards an outflow side at x = 64 m, Manning n = 0.01, a wall at x = 0, `flume.max_dt = 0.25`, a log line every 64 steps (600 is not a multiple: the last line comes from the end of the run), gauges at mid-length and at the outlet. |
| `inputs_standalone` | The keys of the stand-alone run (`max_dt_plane`, registered in `Tests/CMakeLists.txt` in every build) that ERF sets in the coupled run: 64 x 8 cells of 1 m in two boxes of 32 x 8, periodic in y, a plotfile every 200 steps, 150 s. |
| `inputs_erf` | ERF with FLUME: a still, uniform anelastic atmosphere over 64 x 8 x 64 m in 32 x 4 x 8 cells (boxes of 16 x 4 x 8), periodic in y, slip walls in x and z, `erf.fixed_dt = 0.25`, 600 steps, a checkpoint every 300 steps, `erf.flood.grid_ratio = 2`, a flood plotfile every 200 steps. |
| `check_coupled.py` | Checks that each flood plotfile `floodNNNNN` holds FLUME step N at time 0.25 N exactly; with `--reference`, that the plotfiles and logs are those of the other run; with `--column-mean`, that the flood fields of ERF's plotfile are the means of FLUME's cells over each column; with `--restart`, that ERF continued from its checkpoint gives them again; with `--substeps`, that FLUME landed on the end of every ERF step after more FLUME steps than ERF steps. |
| `ERFTests.cmake` | The CTest definitions, included from `Tests/CMakeLists.txt` with `FLUME_WITH_ERF=ON`. |

## Why the two runs agree to the last bit

ERF refines its 32 x 4 columns by 2 into the flood grid: 64 x 8 cells of 1 m
on [0, 64] x [0, 8] m, in two boxes of 32 x 8 cells, which are the boxes
`flume.max_grid_size = 32` gives the stand-alone run. FLUME's step is the
shortest of `flume.max_dt`, the CFL step and the step with rain. On this
plane the CFL and rain steps stay above 0.99 s to the end (0.998 s at the
shortest, in a run with `flume.max_dt=1000`), so every step of
the stand-alone run is 0.25 s, and so is every ERF step: each ERF step is one
FLUME step. The time sums exactly in binary (0.25 s), so FLUME lands on the
same times in both runs, and every operation is the same.

## Data and acknowledgements

Every input is synthetic, written for these tests. No external data are used
or redistributed, and nothing needs to be downloaded.

## How to run

From the build directory of a `-DFLUME_WITH_ERF=ON` build:

```bash
ctest -R erf_flume --output-on-failure
```

By hand, each run from its own empty directory:

```bash
mpiexec -n 2 <build>/flume <repo>/Tests/ERF/inputs_standalone FILE=<repo>/Tests/ERF/flume_plane
```

```bash
mpiexec -n 2 <build>/ERF/Exec/erf_exec <repo>/Tests/ERF/inputs_erf FILE=<repo>/Tests/ERF/flume_plane
python3 <repo>/Tests/ERF/check_coupled.py flood --dt 0.25 --reference <stand-alone directory> \
    --log water_balance.txt gauge_mid.txt gauge_outlet.txt
```

## Expected Results

| Test | Result |
| --- | --- |
| `max_dt_plane` (stand-alone) | plotfiles at steps 0, 200, 400, 600, at 0, 50, 100, 150 s |
| `erf_flume_plane` | the same four plotfiles, every cell, and the three logs, byte for byte; the flood fields of `plt00600` within 1.1e-16 of the column means; the same files again after the restart from `chk00300` |
| `erf_flume_plane_substeps` | more FLUME steps than ERF steps (128 in 75 ERF steps of 2 s); water balance residual 8.3e-17 m3 |
| `erf_flume_plane_stop_time` | `FLUME: 174 step(s), time 52.2 s`: 174 ERF steps of 0.3 s, 1.4e-13 s short of `stop_time = 52.2`, end the run; one FLUME step each with `flume.max_dt = 0.29999999` |
| `abort_erf_flume_host_*` | each of the 12 keys stops the run, with ERF's key in the message |

## CTest

`max_dt_plane` (every build), `erf_flume_plane`, `erf_flume_plane_substeps`,
`erf_flume_plane_stop_time` (2 ranks each) and the `abort_erf_flume_host_*`
tests (1 rank). See
`Docs/sphinx_doc/TestCatalog.rst`.
