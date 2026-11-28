#!/usr/bin/env python3
"""Check FLUME's start inside ERF from a checkpoint written without the flood.

Usage (in the directory of the run, as the check of Tests/run_case.py):

    check_start_without_flood.py <prefix> --start-time <s> --end-time <s>
                                 --surface <m> --wet <n>
                                 --restart <checkpoint> <command...>

The run here is ERF with erf.flood.enable = false: it has written the ERF
checkpoint <checkpoint>, with no FLUME state in it, and no flood output. The
command continues from that checkpoint with the flood on and is run here. It
must exit with 0 and print ERF's "The checkpoint holds no flood state", and
FLUME then starts from its initial state (a flat water surface at <surface>
m) at ERF's time, <start-time> s. Then:

- FLUME writes its summary at ERF's last step (max_step, counted from the
  start of the first run): its output ends with "FLUME: N step(s), time
  <end-time> s".
- water_balance.txt and each gauge file start at <start-time> s, at FLUME
  step 0, and the water balance closes to round-off (Tests/check_balance.py).
- In each flood plotfile <prefix>NNNNN (at least two), the arrival time of
  each of the <n> cells wet at the start (<surface> above the bed of the
  plotfile by more than the arrival depth, flume.dry_depth = 1e-6 m by
  default, which the deck keeps) is <start-time> s exactly, and that of every other cell is -1 or
  later than <start-time> s.
"""
import re
import subprocess
import sys
from pathlib import Path

TESTS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TESTS))
import check_balance  # noqa: E402
from check_plotfile import read_cell_data, read_cell_h, read_header  # noqa: E402
from compare_plotfiles import plotfiles  # noqa: E402


# flume.flood_maps.arrival_depth: its default, flume.dry_depth's default.
ARRIVAL_DEPTH = 1.0e-6


def fail(msg):
    sys.exit(f"FAIL: {msg}")


def first_data_line(path):
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            return line.split()
    fail(f"{path}: no data line")


def main():
    args = sys.argv[1:]
    if "--restart" not in args:
        fail("--restart <checkpoint> <command...> is required")
    i = args.index("--restart")
    restart, args = args[i + 1:], args[:i]
    opts = {}
    if len(args) != 9:
        fail("usage: check_start_without_flood.py <prefix> --start-time <s> --end-time <s> "
             "--surface <m> --wet <n> --restart <checkpoint> <command...>")
    prefix = args[0]
    for k in range(1, 9, 2):
        opts[args[k]] = args[k + 1]
    try:
        t0, t_end = float(opts["--start-time"]), float(opts["--end-time"])
        surface, wet = float(opts["--surface"]), int(opts["--wet"])
    except (KeyError, ValueError):
        fail(f"cannot read --start-time, --end-time, --surface and --wet from {args[1:]}")
    if len(restart) < 2:
        fail("--restart needs a checkpoint and a command")
    chk, command = restart[0], restart[1:]
    if not (Path(chk) / "Header").is_file():
        fail(f"no ERF checkpoint {chk} here")
    if (Path(chk) / "Level_0" / "FLUME").exists() or Path("water_balance.txt").exists():
        fail(f"the run here kept the flood ({chk}/Level_0/FLUME or water_balance.txt): "
             f"it must run with erf.flood.enable = false")

    print(f"restart: {' '.join(command)}", flush=True)
    run = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         encoding="utf-8", errors="backslashreplace")
    print(run.stdout, flush=True)
    if run.returncode != 0:
        fail(f"the restart exited with {run.returncode}")
    if "The checkpoint holds no flood state" not in run.stdout:
        fail("ERF did not report a checkpoint without flood state")
    m = re.search(r"^FLUME: (\d+) step\(s\), time (\S+) s$", run.stdout, re.M)
    if not m or float(m.group(2)) != t_end:
        fail(f"no FLUME summary at {t_end} s: FLUME did not see ERF's last step")
    print(f"PASS: FLUME wrote its summary after {m.group(1)} step(s), at {t_end} s")

    logs = ["water_balance.txt"] + sorted(str(p) for p in Path(".").glob("gauge_*.txt"))
    if len(logs) < 2:
        fail(f"logs {logs}: water_balance.txt and at least one gauge file needed")
    for log in logs:
        row = first_data_line(log)
        if int(row[0]) != 0 or float(row[1]) != t0:
            fail(f"{log} starts at step {row[0]}, time {row[1]} s, not step 0 at {t0} s")
    ok, msg, _ = check_balance.check(Path("water_balance.txt"))
    if not ok:
        fail(f"water balance: {msg}")
    print(f"PASS: {', '.join(logs)} start at step 0, {t0} s; balance: {msg}")

    names = plotfiles(".", prefix)
    if len(names) < 2:
        fail(f"{len(names)} plotfiles {prefix}NNNNN here, at least 2 needed")
    for name in names:
        h = read_header(Path(name) / "Header")
        if "arrival_time" not in h["vars"] or "bed_elevation" not in h["vars"]:
            fail(f"{name} holds no arrival_time or bed_elevation (flume.flood_maps = 1 needed)")
        _, fab = read_cell_h(Path(name) / "Level_0" / "Cell_H")
        cells, _ = read_cell_data(Path(name) / "Level_0", fab)
        ia, ib = h["vars"].index("arrival_time"), h["vars"].index("bed_elevation")
        n_wet = 0
        for ij, v in cells.items():
            if surface - v[ib] > ARRIVAL_DEPTH:
                n_wet += 1
                if v[ia] != t0:
                    fail(f"{name} cell {ij}, wet at the start: arrival time {v[ia]!r} s, "
                         f"not {t0} s")
            elif v[ia] != -1.0 and not v[ia] > t0:
                fail(f"{name} cell {ij}, dry at the start: arrival time {v[ia]!r} s, "
                     f"neither -1 nor after {t0} s")
        if n_wet != wet:
            fail(f"{name}: {n_wet} cells wet at the start, {wet} expected")
    print(f"PASS: {len(names)} plotfiles; in each the {wet} cells wet at the start arrived at "
          f"{t0} s and no other cell before")


if __name__ == "__main__":
    main()
