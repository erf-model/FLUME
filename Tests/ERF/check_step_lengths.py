#!/usr/bin/env python3
"""Check that the flood inside ERF does not depend on ERF's step beyond FLUME's time error.

Usage (in any directory, as the check of a test):

    check_step_lengths.py <prefix> --coarse <dir> --fine <dir> --finer <dir> --own <dir>
                          --runs <dir>...

Each directory holds a run of the same flood (Tests/ERF/flume_plane) with a
last plotfile <prefix>NNNNN at the same time, on the same grid:

  --coarse stand-alone, every step 2 dt (flume.max_dt = 0.5 s)
  --fine   stand-alone, every step flume.max_dt = dt (0.25 s)
  --finer  stand-alone, every step dt / 2 (flume.max_dt = 0.125 s)
  --own    stand-alone, FLUME's own steps (CFL and rain limited, no cap)
  --runs   inside ERF, each with another ERF step and no cap on FLUME's

FLUME's first-order time stepping has an error that halves with the step, so
F* = 2 F(finer) - F(fine) is the flood with the time error of the fine run
removed to leading order (Richardson extrapolation), for F the depth and the
discharge in x. The space error is the same in every run (one grid), so the
distance of a run from F* is its time error.

Checked:
  1. The time error is first order, so that F* is a sound reference: the
     observed ratio |F(coarse) - F(fine)| / |F(fine) - F(finer)| (L1 norm
     over the cells) is within 2 +- 0.1. This uses three runs; a ratio of
     distances from F* would be 2 by construction.
  2. No run inside ERF is further from F* than the stand-alone run with
     FLUME's own steps, within a margin of 2 % (L1 and largest value, depth
     and discharge). Inside ERF, FLUME takes its own CFL step and shortens
     only the step that lands on the end of an ERF step, so no step is longer
     than in that run. A step that went past FLUME's own limit, rain counted
     twice or left out in a step, or a step of ERF missed or repeated, would
     leave a run further away. This is observed, not a theorem: shorter steps
     need not give a smaller error everywhere. The stand-alone run with every
     step capped at 2 s is 1.0103 times further than FLUME's own steps in the
     largest depth difference. The margin covers that.
  3. Each run inside ERF differs from the run with FLUME's own steps by at
     most half of that run's time error (L1 and largest value), so the flood
     inside ERF is the stand-alone flood within FLUME's time error.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from check_dam_break import read  # noqa: E402
from compare_plotfiles import plotfiles  # noqa: E402

COMPONENTS = ((0, "depth"), (1, "discharge_x"))
MARGIN = 1.02
ORDER_BAND = 0.1
CLOSE = 0.5


USAGE = ("usage: check_step_lengths.py <prefix> --coarse <dir> --fine <dir> --finer <dir> "
         "--own <dir> --runs <dir>...")


def fail(msg):
    sys.exit(f"FAIL: {msg}")


def last_plotfile(where, prefix):
    names = plotfiles(where, prefix)
    if not names:
        fail(f"no plotfile {prefix}NNNNN in {where}")
    h, cells, _ = read(Path(where) / names[-1])
    return h, cells, names[-1]


def errors(cells, ref, comp, scale):
    """L1 (mean) and largest |cells - ref| over the cells, divided by scale."""
    d = [abs(cells[k][comp] - ref[k]) for k in ref]
    return sum(d) / len(d) / scale, max(d) / scale


def column(cells, comp):
    return {k: v[comp] for k, v in cells.items()}


def main():
    args = sys.argv[1:]
    opts = {"--coarse": None, "--fine": None, "--finer": None, "--own": None}
    runs = []
    if not args:
        fail(USAGE)
    prefix, rest = args[0], args[1:]
    i = 0
    while i < len(rest):
        if rest[i] in opts and i + 1 < len(rest):
            opts[rest[i]] = rest[i + 1]
            i += 2
        elif rest[i] == "--runs":
            runs = rest[i + 1:]
            break
        else:
            fail(f"unexpected argument {rest[i]!r}")
    if None in opts.values() or not runs:
        fail(USAGE)

    hc, coarse, _ = last_plotfile(opts["--coarse"], prefix)
    hf, fine, _ = last_plotfile(opts["--fine"], prefix)
    hff, finer, _ = last_plotfile(opts["--finer"], prefix)
    ho, own, _ = last_plotfile(opts["--own"], prefix)
    t = hf["time"]
    if any(c.keys() != fine.keys() for c in (coarse, finer, own)):
        fail("the stand-alone runs are not on the same grid")
    for name, h in (("--coarse", hc), ("--finer", hff), ("--own", ho)):
        if h["time"] != t:
            fail(f"{name}: last plotfile at {h['time']!r} s, the fine run's at {t!r} s")
    others = []
    for d in runs:
        h, cells, name = last_plotfile(d, prefix)
        # ERF's steps add exactly in binary to the end time here.
        if h["time"] != t:
            fail(f"{d}/{name}: at {h['time']!r} s, the fine run's at {t!r} s")
        if cells.keys() != fine.keys():
            fail(f"{d}: not on the grid of the stand-alone runs")
        others.append((d, h["step"], cells))

    ok = True
    for comp, label in COMPONENTS:
        ref = {k: 2.0 * finer[k][comp] - fine[k][comp] for k in fine}
        scale = max(abs(v) for v in ref.values())
        if not scale > 0.0:
            fail(f"{label}: zero everywhere at {t} s, nothing to compare")
        d_coarse = errors(coarse, column(fine, comp), comp, scale)[0]
        d_fine = errors(fine, column(finer, comp), comp, scale)[0]
        if not d_fine > 0.0:
            fail(f"{label}: the 0.25 s and 0.125 s runs are the same, no time error to measure")
        ratio = d_coarse / d_fine
        print(f"{label} at {t} s (scale {scale:.4e}): "
              f"|F(coarse) - F(fine)| / |F(fine) - F(finer)| = {ratio:.4f}")
        if abs(ratio - 2.0) > ORDER_BAND:
            fail(f"{label}: the time error is not first order "
                 f"(ratio {ratio:.4f}, 2 +- {ORDER_BAND} expected)")
        e_own = errors(own, ref, comp, scale)
        print(f"  {'own steps':24s} {ho['step']:5d} steps  L1 {e_own[0]:.4e}  max {e_own[1]:.4e}")
        for d, steps, cells in others:
            e = errors(cells, ref, comp, scale)
            near = errors(cells, column(own, comp), comp, scale)
            good = e[0] <= MARGIN * e_own[0] and e[1] <= MARGIN * e_own[1]
            close = near[0] <= CLOSE * e_own[0] and near[1] <= CLOSE * e_own[1]
            ok = ok and good and close
            print(f"  {Path(d).name:24s} {steps:5d} steps  L1 {e[0]:.4e}  max {e[1]:.4e}"
                  f"  {'ok' if good else 'FURTHER than own steps'};"
                  f"  from own steps / its error: L1 {near[0] / e_own[0]:.3f}"
                  f"  max {near[1] / e_own[1]:.3f}  {'ok' if close else 'NOT within half'}")
    if not ok:
        fail(f"a run inside ERF is further from the time-converged flood than FLUME's own steps "
             f"(margin {MARGIN}), or differs from them by more than {CLOSE} of their time error")
    print(f"PASS: {len(others)} runs inside ERF are within FLUME's own time error at {t} s")


if __name__ == "__main__":
    main()
