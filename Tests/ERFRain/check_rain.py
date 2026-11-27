#!/usr/bin/env python3
"""Check the rain FLUME took from ERF's microphysics (flume.rain.type = microphysics).

Usage (in the directory of the run, as the check of Tests/run_case.py):

    check_rain.py cells <ratio> <step>... [--min-spread S] [--min-frozen F]
    check_rain.py volume <ratio> <step> [--restart <checkpoint> <command...>]

ERF writes the accumulated surface precipitation of each column, kg/m2 (=
mm of liquid water), in its 2D plotfile plt2d<step> (precip_total_accum,
liquid and frozen together, and precip_frozen_accum), and FLUME the flood in
flood<step>. Each column of ERF's level 0 covers ratio x ratio flood cells.

cells: no water moves (flume.dry_depth above any depth reached). At each
step given, every flood cell must hold, in m, the precipitation of its column
divided by 1000, to 1e-12 of the largest precipitation. The largest
precipitation over the columns must be at least S times the smallest (above
0), so that a cell given another column's rain fails; and with
--min-frozen, the frozen part must be at least F of the total over the grid,
so that a run that leaves out the frozen part fails.

volume: the rain column of the last line of water_balance.txt (time of
plt2d<step>) must be the precipitation summed over the columns, times the
column area, divided by 1000, to 1e-12 of it plus the rain each cell may
still hold in its carry (FLUME_Rain.cpp: the part of a step's rain that the
addition to the depth rounds away waits for the next step): one machine
epsilon of the largest depth per cell, times the cell area. A negative
precipitation of a step waits in the carry whole (FLUME_Rain.cpp) and is
not in the bound: Kessler gives none in these runs. The water must have moved: in some
cell the change of depth from flood00000 must differ from its column's
precipitation / 1000 by at least 1e-3 m.

With --restart, the run here has written the ERF checkpoint <checkpoint>
(ending in its step n). The flood plotfiles after step n are removed, the
water balance is kept (the restart must cut it back), the command, which
continues from that checkpoint, is run here; it must exit with 0 and print
FLUME's "restarted from <checkpoint>/Level_0/FLUME at step n". Every flood
plotfile and water_balance.txt must then be those of the first run, to the
last bit, and the volume is checked again.
"""
import filecmp
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

TESTS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TESTS))
sys.path.insert(0, str(TESTS / "ERF"))
from check_balance import read_log  # noqa: E402
from check_coupled import read_layers  # noqa: E402
from check_plotfile import read_cell_data, read_cell_h, read_header  # noqa: E402
from compare_plotfiles import plotfiles  # noqa: E402


def fail(msg):
    sys.exit(f"FAIL: {msg}")


def precip(step):
    """{(i, j): (total, frozen)} of plt2d<step>, kg/m2, and its header."""
    name = f"plt2d{step:05d}"
    h = read_header(Path(name) / "Header")
    for v in ("precip_total_accum", "precip_frozen_accum"):
        if v not in h["vars"]:
            fail(f"{name} has no {v}")
    ct, cf = h["vars"].index("precip_total_accum"), h["vars"].index("precip_frozen_accum")
    cells, _ = read_layers(Path(name) / "Level_0")
    out = {}
    for (i, j, k), v in cells.items():
        if (i, j) in out:
            fail(f"{name}: more than one layer at column {(i, j)}")
        out[(i, j)] = (v[ct], v[cf])
    return out, h


def flood(step):
    name = f"flood{step:05d}"
    h = read_header(Path(name) / "Header")
    _, fab = read_cell_h(Path(name) / "Level_0" / "Cell_H")
    cells, _ = read_cell_data(Path(name) / "Level_0", fab)
    return cells, h


def check_cells(ratio, steps, min_spread, min_frozen):
    for step in steps:
        p, hp = precip(step)
        cells, hf = flood(step)
        if hf["time"] != hp["time"]:
            fail(f"flood{step:05d} at {hf['time']!r} s, plt2d{step:05d} at {hp['time']!r} s")
        c = hf["vars"].index("depth")
        if len(cells) != ratio * ratio * len(p):
            fail(f"flood{step:05d} has {len(cells)} cells, {ratio}^2 x {len(p)} columns expected")
        top = max(v[0] for v in p.values())
        low = min(v[0] for v in p.values())
        if top <= 0.0:
            fail(f"plt2d{step:05d}: no precipitation")
        worst = 0.0
        for (i, j), (total, _) in p.items():
            for a in range(ratio):
                for b in range(ratio):
                    d = cells[(ratio * i + a, ratio * j + b)][c]
                    worst = max(worst, abs(1000.0 * d - total))
        if worst > 1e-12 * top:
            fail(f"flood{step:05d}: 1000 x depth differs from precip_total_accum of its column "
                 f"by {worst!r} kg/m2, bound {1e-12 * top!r}")
        spread = top / low if low > 0.0 else float("inf")
        if spread < min_spread:
            fail(f"plt2d{step:05d}: precipitation from {low!r} to {top!r} kg/m2, a ratio of "
                 f"{spread:.3g}, below {min_spread}: a cell given another column's rain would pass")
        frozen = sum(v[1] for v in p.values()) / sum(v[0] for v in p.values())
        if frozen < min_frozen:
            fail(f"plt2d{step:05d}: frozen part {frozen:.3g} of the total, below {min_frozen}")
        print(f"PASS: flood{step:05d}: every cell holds its column's precipitation / 1000 "
              f"(largest difference {worst:.3g} of {top:.3g} kg/m2; spread {spread:.3g}, "
              f"frozen part {frozen:.3g})")


def check_volume(ratio, step):
    p, hp = precip(step)
    dx = hp["dx"]
    expected = sum(v[0] for v in p.values()) * dx[0] * dx[1] / 1000.0
    eps, names, rows = read_log("water_balance.txt")
    last = dict(zip(names, rows[-1]))
    if last["time"] != hp["time"]:
        fail(f"the last line of water_balance.txt is at {last['time']!r} s, plt2d{step:05d} "
             f"at {hp['time']!r} s")
    cells, hf = flood(step)
    c = hf["vars"].index("depth")
    area = dx[0] * dx[1] / ratio ** 2
    bound = 1e-12 * expected + len(cells) * eps * max(v[c] for v in cells.values()) * area
    if expected <= 0.0 or abs(last["rain"] - expected) > bound:
        fail(f"rain {last['rain']!r} m3 in water_balance.txt, precipitation over the grid "
             f"{expected!r} m3, bound {bound!r}")
    start, _ = flood(0)
    moved = max(abs(cells[(ratio * i + a, ratio * j + b)][c]
                    - start[(ratio * i + a, ratio * j + b)][c] - v[0] / 1000.0)
                for (i, j), v in p.items() for a in range(ratio) for b in range(ratio))
    if moved < 1e-3:
        fail(f"flood{step:05d}: the depth of every cell changed by its column's precipitation "
             f"/ 1000 to within {moved!r} m: no water moved")
    print(f"PASS: rain {last['rain']!r} m3, precipitation over the grid {expected!r} m3 "
          f"(difference {abs(last['rain'] - expected):.3g}, bound {bound:.3g}); the water moved (by up to "
          f"{moved:.3g} m in a cell)")


def check_restart(ratio, step, chk, command):
    m = re.search(r"(\d+)$", chk)
    if not m or not (Path(chk) / "Header").is_file() or not command:
        fail(f"--restart needs a checkpoint written here (no {chk}/Header) and a command")
    n = int(m.group(1))
    ref = Path(tempfile.mkdtemp(prefix="first_run_", dir="."))
    after = [p for p in plotfiles(".", "flood") if int(p[len("flood"):]) > n]
    if not after:
        fail(f"no flood plotfile after the checkpoint's step {n}: the restart would show nothing")
    for p in plotfiles(".", "flood"):
        shutil.copytree(p, ref / p)
    shutil.copy("water_balance.txt", ref / "water_balance.txt")
    for p in after:
        shutil.rmtree(p)
    print(f"removed {' '.join(after)}; restart: {' '.join(command)}", flush=True)
    run = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         encoding="utf-8", errors="backslashreplace")
    print(run.stdout, flush=True)
    if run.returncode != 0:
        fail(f"the restart exited with {run.returncode}")
    if not re.search(rf"restarted from {re.escape(chk)}/Level_0/FLUME at step {n}\b", run.stdout):
        fail(f"FLUME did not report continuing from {chk}/Level_0/FLUME at step {n}")
    r = subprocess.run([sys.executable, str(TESTS / "compare_plotfiles.py"), "flood", str(ref),
                        "--cells"])
    if r.returncode != 0:
        fail(f"the flood plotfiles after the restart differ from those of the first run ({ref})")
    if not filecmp.cmp("water_balance.txt", ref / "water_balance.txt", shallow=False):
        fail(f"water_balance.txt after the restart differs from that of the first run ({ref})")
    print("PASS: after the restart, the flood plotfiles and water_balance.txt are those of the "
          "first run")
    check_volume(ratio, step)


def main():
    args = sys.argv[1:]
    restart = []
    if "--restart" in args:
        i = args.index("--restart")
        restart = args[i + 1:]
        args = args[:i]
    opts = {"--min-spread": 1.0, "--min-frozen": 0.0}
    for key in opts:
        if key in args:
            i = args.index(key)
            opts[key] = float(args[i + 1])
            del args[i:i + 2]
    if len(args) < 3 or args[0] not in ("cells", "volume") \
            or (args[0] == "volume" and len(args) != 3) or (restart and args[0] != "volume"):
        fail("usage: check_rain.py cells <ratio> <step>... [--min-spread S] [--min-frozen F] | "
             "volume <ratio> <step> [--restart <checkpoint> <command...>]")
    ratio, steps = int(args[1]), [int(s) for s in args[2:]]
    if args[0] == "cells":
        check_cells(ratio, steps, opts["--min-spread"], opts["--min-frozen"])
        return
    check_volume(ratio, steps[0])
    if restart:
        check_restart(ratio, steps[0], restart[0], restart[1:])


if __name__ == "__main__":
    main()
