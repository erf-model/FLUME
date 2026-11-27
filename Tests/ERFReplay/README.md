# Replay of the rain from ERF's microphysics

## Purpose

Check `flume.rain.replay_file`. Inside ERF, with
`flume.rain.type = microphysics`, FLUME writes the rain it takes from ERF
as gridded rain: an index file and one ESRI ASCII grid of the rate per ERF
step with rain. A stand-alone run with `flume.rain.type = gridded` reads it
back without the atmosphere.

The checks:

- the index is closed at the end of the coupled rain by a grid of no rain,
  and its grids hold the rain of the coupled water balance (the writer on
  its own);
- the stand-alone replay gives the coupled flood to round-off (the coupling
  as a whole);
- the replay runs on after the coupled rain ends, with no more rain, while
  the water still moves (the recession);
- a restart of the coupled run from ERF's checkpoint writes the same index
  and grids, to the byte;
- a restart without the index, with an index that is not the run's, or
  from a checkpoint without what the writer holds back, stops.

The replay description is in `Docs/sphinx_doc/Running.rst` (Replay of
ERF's rain). The checks are bounded for double precision; a
single-precision build fails them.

## Files

| File | What it is |
| --- | --- |
| `inputs_replay` | The stand-alone replay: the grid ERF sets for `Tests/ERFRain/inputs_kessler` (32 x 32 cells of 50 m, periodic in x and y), in 16 boxes of 8 x 8 cells, `flume.max_dt = 0.25` (ERF's step), a plotfile every 60 steps, 30 s. It is run with `FILE=<repo>/Tests/ERFRain/flume_basin` and then `flume.rain.type=gridded` and `flume.rain.file=<index>` on the command line, so the flood is the coupled run's. |
| `check_replay.py` | The checks above, in the directory of the coupled run. It runs the replay in `replay/`, the run to 60 s in `recession/`, the three restarts that must stop in `abort_index/`, `abort_stale/` and `abort_held/`, and the restart. |
| `ERFReplayTests.cmake` | The CTest definitions, included from `Tests/CMakeLists.txt` with `FLUME_WITH_ERF=ON`. |

The coupled run uses `Tests/ERFRain/inputs_kessler`, `sounding_moist` and
`flume_basin` (see `Tests/ERFRain/README.md`), with
`amr.max_grid_size = 8 8 32` (4 boxes for 4 ranks), `erf.fixed_dt = 0.25`
(below ERF's CFL limit of about 0.34 s, so every step is 0.25 s),
`flume.max_dt = 0.25` and `flume.rain.replay_file = rain/replay.txt`.
`check_replay.py` reads the plotfiles with `Tests/check_plotfile.py` and the
water balance with `Tests/check_balance.py`.

## Data and acknowledgements

No external data are used or redistributed, and nothing needs to be
downloaded. The atmosphere is that of ERF's
`Tests/FloodEngineExample` (from the ERF repository, `Submodules/ERF`,
BSD-3-Clause-LBNL licence), through `Tests/ERFRain`, with thanks to the ERF
developers.

## How to run

From the build directory of a `-DFLUME_WITH_ERF=ON -DFLUME_MPI=ON` build:

```bash
ctest -R "rain_replay" --output-on-failure
```

By hand, from an empty directory:

```bash
mpiexec -n 4 <build>/ERF/Exec/erf_exec <repo>/Tests/ERFRain/inputs_kessler \
    FILE=<repo>/Tests/ERFRain/flume_basin \
    erf.input_sounding_file=<repo>/Tests/ERFRain/sounding_moist \
    amr.max_grid_size="8 8 32" erf.fixed_dt=0.25 flume.max_dt=0.25 \
    flume.rain.replay_file=rain/replay.txt
mkdir replay && cd replay
mpiexec -n 4 <build>/flume <repo>/Tests/ERFReplay/inputs_replay \
    FILE=<repo>/Tests/ERFRain/flume_basin flume.rain.type=gridded \
    flume.rain.file=../rain/replay.txt
```

Add `flume.stop_time=60` to the second run to carry the flood on for 30 s
after the rain.

## Expected Results

| Check | Result |
| --- | --- |
| Index | 121 lines: the grid of no rain at 0 s (step 1 has no rain), 119 grids of rain from 0.25 s to 29.75 s, the grid of no rain at 30 s. Their rain, 0.342445125773134 m3, is that of the water balance, 0.3424451257736819 m3 (difference 5.5e-13 m3: the rain still held in the cells' rounding carry, bound 5.7e-11 m3). |
| Replay | Every depth of flood00060 and flood00120 within 0.008 of 64 eps (h + R) of the coupled depth (R = 4.2e-7 m, the largest rain of a cell); both discharges within 64 eps of the largest discharge; rain 0.34244512577368186 m3 against 0.3424451257736819 m3 (1.6e-16 of it). |
| Recession | The run reaches 60 s; the rain at 60 s is that at 30 s, to the bit; a depth changes by 0.011 m from 30 s to 60 s. |
| Restart | From `chk00060`: the index and the 120 grids (119 of rain and the grid of no rain) are those of the first run, to the byte. |
| Aborts | Without the index: `flume.rain.replay_file: rain/replay.txt: cannot be opened`. With the index cut back to 7.5 s and closed there (as a run that went on from an older checkpoint and was stopped leaves it): `flume.rain.replay_file: rain/replay.txt: does not hold what the run wrote up to the checkpoint`. From a checkpoint without `rain_replay_held`: `flume.restart: ...: it was written without flume.rain.replay_file`. |

The tolerance of the replay: each step's rate is written with 17 digits,
read back and multiplied by the step, so each step's rain is that of the
coupled run within 2 eps of it, and the rain of a cell over the run within a
few eps of R. Where the rain changes a depth by less than one unit in its
last place, the two runs can round to neighbouring doubles: eps h. The flow
carries these differences on, most in thin layers of water, so 64 eps
(h + R) is a margin measured for this case (the largest difference is 0.008
of it), not a bound for every case: with ERF steps of 0.1 s over 60 s the
same comparison reaches 2.7 times it in cells 3e-6 m deep. The discharges
are compared with the largest discharge of either component, since the y
discharge of this flow is round-off.

Each check fails against these one-line changes to the writer
(`Source/FLUME_RainReplay.cpp`):

| Change | Check that fails |
| --- | --- |
| The rate in mm/(1000 h) (3.6e3 for 3.6e6) | Index (the rain of the grids) |
| No closing line | Index (the last line) |
| A step with rain does not take the place of the closing line before it | Index (two lines at one time) |
| Rates with 9 digits | Replay (the depths, the discharges and the rain) |
| A step of no rain writes a grid | Index (a grid of no rain first) |

Not reached by this test: a negative depth from ERF (Kessler gives none in
this run), a step of no rain after the first, and a step that starts an
instant before the end of the last (a single-precision ERF). The writer is
checked on these by `rain_replay_unit` (`Tests/RainReplay`), which also has
the mutants of the restart.

## CTest

`erf_flume_rain_replay` (4 ranks), `abort_erf_flume_rain_replay_directory`
and `abort_erf_flume_rain_replay_hash` (1 rank each). The stand-alone checks
of the key are `abort_rain_replay_file_with_constant`,
`abort_rain_replay_file_without_type` and `abort_rain_replay_file_empty_name`
(`Tests/Rain`). The writer on its own: `rain_replay_unit`
(`Tests/RainReplay`). See `Docs/sphinx_doc/TestCatalog.rst`.
