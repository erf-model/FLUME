#!/usr/bin/env python3
"""Check the writer of flume.rain.replay_file on its own (replay_unit.cpp).

Usage (in the directory where `replay_unit run` ran, as the check of
Tests/run_case.py):

    check_replay_unit.py [--single] <command that runs replay_unit, without its mode...>

--single: replay_unit was built in single precision (FLUME_PRECISION=SINGLE).
The depths and what a cell holds back are then floats: the check rounds each
depth and each sum to the nearest float, as the writer does, and multiplies
by 3.6e6 / (t1 - t0) in double, so the grids still match to the bit.

run: the index rain/replay.txt must be, to the byte,

    0 replay_000001.asc
    0.5 replay_zero.asc
    1 replay_000002.asc
    1.5 replay_000003.asc
    2 replay_zero.asc

The first step takes the place of the line of no rain at the start; the
second step has no rain (cell 1 makes up part of what it holds back, cell 6
holds back more) and adds no line; the third does not start where the first
ended, so it gets its own line; the fourth starts 1e-8 s before the end of the
third (within the tolerance of a single-precision host, 1e-6 of the step) and
takes the place of its closing line. Every grid must hold, to the bit, the
rate in mm/h of the depths of replay_unit.cpp with the rule of
FLUME_RainReplay.H (what a cell holds back is added to its depth; the part
above 0 is written, divided by the step and times 3.6e6; the rest is held),
computed here in the same double arithmetic; the grid of no rain is all 0.

restart: in restart2/, from the checkpoint chk2 of the run (written after the
second step, with cells 1 and 6 holding back) and the run's index with a line
cut short appended (as a later run that was stopped leaves it) and grids 2
and 3 removed, `replay_unit restart2` must write the index and every grid of
the run again, to the byte. The same in restart3/ from chk3 (after the third
step, cell 6 holding back; the index holds the line of grid 3 where chk3 had
its closing line) with grid 3 removed. In resume3/, from chk3 with the same index,
`replay_unit resume3` (no step after the restart) must leave the index of the
run as chk3 had it, closed at 1.5 s: the restart writes the closing line again
over the line of grid 3.

aborts: in stale/, chk2 with an index of the one line of no rain at 0 s
must stop with "does not hold what the run wrote up to the
checkpoint"; in early/, `replay_unit early` (a step from 0.25 s after one
that ended at 0.5 s) must stop with "starts before the end of the last step
with rain".
"""
import filecmp
import re
import shutil
import struct
import subprocess
import sys
from pathlib import Path

DEPTHS = [
    [1.0e-3, -2.0e-3, 0.0, 0.0, 0.0, 5.0e-4, 0.0, 0.0],
    [0.0, 1.0e-3, 0.0, 0.0, 0.0, 0.0, -3.0e-4, 0.0],
    [0.0, 3.0e-3, 0.0, 2.0e-3, 0.0, 0.0, 1.0e-4, 0.0],
    [1.0e-3, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0e-4, 0.0],
]
TIMES = [(0.0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5 - 1.0e-8, 2.0)]
INDEX = ("0 replay_000001.asc\n0.5 replay_zero.asc\n1 replay_000002.asc\n"
         "1.5 replay_000003.asc\n2 replay_zero.asc\n")


ROUND = lambda x: x  # noqa: E731  (float32 rounding with --single)


def to_float32(x):
    """x rounded to the nearest float, as a Python float."""
    return struct.unpack("f", struct.pack("f", x))[0]


def fail(msg):
    sys.exit(f"FAIL: {msg}")


def expected_grids():
    """The rates of each grid of rain, by cell i + 4 j, as the writer computes them."""
    held = [0.0] * 8
    grids = []
    for d, (t0, t1) in zip(DEPTHS, TIMES):
        w = [0.0] * 8
        for c in range(8):
            # The sum of two floats is exact in double, so rounding it once
            # gives the float sum.
            x = ROUND(ROUND(d[c]) + held[c])
            w[c] = max(x, 0.0)
            held[c] = min(x, 0.0)
        if max(w) > 0.0:
            scale = 3.6e6 / (t1 - t0)
            grids.append([v * scale for v in w])
    return grids


def read_grid(path):
    lines = path.read_text().splitlines()
    head = dict(line.split() for line in lines[:5])
    if head != {"ncols": "4", "nrows": "2", "xllcorner": "0", "yllcorner": "0", "cellsize": "1"}:
        fail(f"{path}: header {head}")
    rows = [[float(x) for x in line.split()] for line in lines[5:]]
    if len(rows) != 2 or any(len(r) != 4 for r in rows):
        fail(f"{path}: not 2 rows of 4 values")
    # Rows from the top: j = 1 first.
    return rows[1] + rows[0]


def run(command, where, expect_abort=None):
    print(f"in {where.name}/: {' '.join(command)}", flush=True)
    r = subprocess.run(command, cwd=where, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True)
    (where / "run.log").write_text(r.stdout)
    if expect_abort is None:
        if r.returncode != 0:
            print(r.stdout[-3000:])
            fail(f"{' '.join(command)} exited with {r.returncode} in {where}")
        return
    if r.returncode == 0:
        fail(f"{' '.join(command)} in {where} did not stop")
    if not re.search(expect_abort, r.stdout):
        print(r.stdout[-3000:])
        fail(f"{' '.join(command)} in {where} did not print {expect_abort!r}")
    print(f"PASS: {where.name}/ stops with {expect_abort!r}")


def fresh(where, chk, rain):
    if where.exists():
        shutil.rmtree(where)
    where.mkdir()
    shutil.copytree(chk, where / chk.name)
    shutil.copytree(rain, where / "rain")


def main():
    global ROUND
    driver = sys.argv[1:]
    if driver and driver[0] == "--single":
        ROUND = to_float32
        driver = driver[1:]
    if not driver:
        fail("usage: check_replay_unit.py <command that runs replay_unit...>")
    here = Path(".").resolve()
    rain = here / "rain"
    text = (rain / "replay.txt").read_text()
    if text != INDEX:
        fail(f"the index is\n{text}not\n{INDEX}")
    print("PASS: the index has the lines of the four steps, joined where they meet")
    for n, want in enumerate(expected_grids(), 1):
        got = read_grid(rain / f"replay_{n:06d}.asc")
        if got != want:
            fail(f"replay_{n:06d}.asc holds {got}, not {want}")
    if read_grid(rain / "replay_zero.asc") != [0.0] * 8:
        fail("replay_zero.asc is not all 0")
    print("PASS: the grids hold the rates of the depths with what the cells hold back, to the bit")

    names = sorted(p.name for p in rain.iterdir())
    for k, removed in ((2, (2, 3)), (3, (3,))):
        r = here / f"restart{k}"
        fresh(r, here / f"chk{k}", rain)
        for n in removed:
            (r / "rain" / f"replay_{n:06d}.asc").unlink()
        with open(r / "rain" / "replay.txt", "a") as f:
            f.write("2.5 replay_0000")
        run(driver + [f"restart{k}"], r)
        if sorted(p.name for p in (r / "rain").iterdir()) != names:
            fail(f"after the restart from chk{k} rain/ holds other files than after the run")
        _, mismatch, errors = filecmp.cmpfiles(rain, r / "rain", names, shallow=False)
        if mismatch or errors:
            fail(f"after the restart from chk{k} these files differ from those of the run: "
                 f"{' '.join(mismatch + errors)}")
        print(f"PASS: the restart from chk{k} writes the index and the grids of the run, "
              f"to the byte")

    r = here / "resume3"
    fresh(r, here / "chk3", rain)
    run(driver + ["resume3"], r)
    want = "".join(INDEX.splitlines(True)[:3]) + "1.5 replay_zero.asc\n"
    got = (r / "rain" / "replay.txt").read_text()
    if got != want:
        fail(f"after the restart from chk3 without a step the index is\n{got}not\n{want}")
    print("PASS: the restart from chk3 without a step leaves the index closed at 1.5 s")

    s = here / "stale"
    fresh(s, here / "chk2", rain)
    (s / "rain" / "replay.txt").write_text("0 replay_zero.asc\n")
    run(driver + ["restart2"], s,
        expect_abort=r"flume\.rain\.replay_file: rain/replay\.txt: does not hold what the run "
                     r"wrote up to the checkpoint")

    e = here / "early"
    if e.exists():
        shutil.rmtree(e)
    e.mkdir()
    run(driver + ["early"], e,
        expect_abort=r"flume\.rain\.replay_file: rain/replay\.txt: the host step from 0\.25 s "
                     r"starts before the end of the last step with rain")


if __name__ == "__main__":
    main()
