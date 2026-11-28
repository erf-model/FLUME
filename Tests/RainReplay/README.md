# The writer of the rain replay on its own

## Purpose

Check the writer of `flume.rain.replay_file` (`Source/FLUME_RainReplay.cpp`)
on rain a host model can hand over but the coupled test
(`Tests/ERFReplay`) does not reach. Kessler rain inside ERF is never
negative and does not stop once it starts.

- **Negative depths:** what a cell holds back is taken from its later rain,
  and only the part above 0 is written.
- **A step of no rain:** it adds no line.
- **A step that starts 1e-8 s before the end of the last:** as a
  single-precision host gives it, it takes the place of the closing line.
- **Restarts:**
  - from a checkpoint after a step of no rain;
  - from a checkpoint after a step whose closing line the next step writes
    over;
  - a restart that takes no step after it.
- **Aborts:**
  - a restart with an index that is not the run's;
  - a step that starts before the end of the last.

## Files

| File | What it is |
| --- | --- |
| `replay_unit.cpp` | The driver, linked to `flume_lib`. It hands four host steps of 0.5 s to a writer on 4 x 2 cells of 1 m (boxes of 2 x 2 cells), writing `rain/replay.txt`. Modes: `run` (with checkpoints `chk2` and `chk3` after steps 2 and 3), `restart2`, `restart3`, `resume3` (restart from `chk3`, no step) and `early` (a step from 0.25 s after one that ended at 0.5 s). |
| `check_replay_unit.py` | The checks below. It runs the driver's other modes in `restart2/`, `restart3/`, `resume3/`, `stale/` and `early/`. |
| `RainReplayTests.cmake` | The driver target and the CTest definition, included from `Tests/CMakeLists.txt` in every build. |

The depths of each step (m, cell `i + 4 j`):

| Step | Time, s | Depths | What it tests |
| --- | --- | --- | --- |
| 1 | 0 to 0.5 | cell 0: 1e-3; cell 1: -2e-3; cell 5: 5e-4 | Joins the line of no rain at 0 s. Cell 1 holds back 2e-3. |
| 2 | 0.5 to 1 | cell 1: 1e-3; cell 6: -3e-4 | No cell has rain, so no line. |
| 3 | 1 to 1.5 | cell 1: 3e-3; cell 3: 2e-3; cell 6: 1e-4 | Its own line. Cell 1 writes 2e-3; cell 6 still holds back 2e-4. |
| 4 | 1.5 - 1e-8 to 2 | cell 0: 1e-3; cell 6: 4e-4 | Joins the closing line of step 3. Cell 6 writes 2e-4. |

## Data and acknowledgements

No external data are used or redistributed, and nothing needs to be
downloaded.

## How to run

From the build directory:

```bash
ctest -R rain_replay_unit --output-on-failure
```

By hand, from an empty directory, on 2 ranks:

```bash
mpiexec -n 2 <build>/replay_unit run
python3 <repo>/Tests/RainReplay/check_replay_unit.py mpiexec -n 2 <build>/replay_unit
```

## Expected Results

| Check | Result |
| --- | --- |
| Index | `0 replay_000001.asc`, `0.5 replay_zero.asc`, `1 replay_000002.asc`, `1.5 replay_000003.asc`, `2 replay_zero.asc`, to the byte. |
| Grids | Every rate of the three grids equals, to the bit, the check's own result: depth plus what is held back, the part above 0, times 3.6e6 / (t1 - t0), in double. The grid of no rain is all 0. |
| Restarts | From `chk2` and from `chk3`, the index has a line cut short appended (as a later run that was stopped would leave it) and the later grids removed. The restart writes the index and grids of the run again, to the byte. From `chk3` without a step, the index is that of the run closed at 1.5 s. |
| Aborts | `stale/`: `flume.rain.replay_file: rain/replay.txt: does not hold what the run wrote up to the checkpoint`. `early/`: `... the host step from 0.25 s starts before the end of the last step with rain`. |

Each check fails against these one-line changes to
`Source/FLUME_RainReplay.cpp`:

| Change | Check that fails |
| --- | --- |
| A negative depth is written | Grids |
| Nothing is held back | Index (step 2 gets a line) |
| The restart does not read what is held back | Restart from `chk2` (grids 2 and 3) |
| The join tolerance is 0 (no single-precision margin) | The run stops at step 4 |
| The restart does not accept the line of the next grid where the checkpoint had its closing line | Restart from `chk3` stops |
| The restart cuts the index to its size without writing the closing line again | Restart from `chk3` without a step |

In a single-precision build (`FLUME_PRECISION=SINGLE`) the depths and what
a cell holds back are floats. CTest then passes `--single` to the check,
which rounds each depth and each sum to the nearest float, as the writer
does, so the grids still match to the bit.

## CTest

`rain_replay_unit` (2 ranks with MPI, 1 without). See
`Docs/sphinx_doc/TestCatalog.rst`.
