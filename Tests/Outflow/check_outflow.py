#!/usr/bin/env python3
"""Check what crossed the open sides of a run from its water-balance log.

    check_outflow.py (--no-inflow | --takes-in FACTOR) [--expect-stored V]
                     [--min-out V] [--rtol R] [--factor F] [--min-lines N]
                     [<log>]

The log (default water_balance.txt) holds, on every line, the stored volume,
the volume each source has added or removed and, for every outflow or
zero_gradient segment, the volume that has left through it since step 0
(column out_<name>, below 0 for water that entered).

--no-inflow      An outflow segment lets no water in. With

                     tol = F * eps * max over the lines of
                           (stored + water added + |out_*| summed),

                 eps the machine epsilon of the state (from the header), on
                 every line each out_<name> is at least -tol, it has not
                 fallen by more than tol since the line before, and the
                 stored volume is at most the volume at step 0 plus the water
                 added, plus tol. The mass flux of a face that acts as a wall
                 is zero or the rounding of one product, which tol allows.
--takes-in       The comparison run with zero_gradient segments: on the last
                 line at least one out_<name> is below 0 and the stored
                 volume is at least FACTOR times the volume at step 0 plus
                 the water added. Fails if the segments let no water in.
--expect-stored  The stored volume of the last line is V within R (default
                 1e-12) times V, plus tol: a closed form, for a run that
                 loses no water either.
--min-out        The out_<name> columns of the last line add up to at least V:
                 the side is open, not a wall.
--min-lines      The log holds at least N lines (default 3).
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_balance as cb  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--no-inflow", action="store_true")
    mode.add_argument("--takes-in", type=float, metavar="FACTOR")
    p.add_argument("--expect-stored", type=float)
    p.add_argument("--min-out", type=float)
    p.add_argument("--rtol", type=float, default=1.0e-12)
    p.add_argument("--factor", type=float, default=64.0)
    p.add_argument("--min-lines", type=int, default=3)
    p.add_argument("log", nargs="?", default="water_balance.txt")
    a = p.parse_args()

    eps, names, rows = cb.read_log(a.log)
    added, _ = cb.read_terms(a.log)
    out = [n for n in names if n.startswith("out_")]
    if not out:
        sys.exit(f"FAIL: {a.log}: no out_<name> column: the run has no open segment")
    if len(rows) < a.min_lines:
        sys.exit(f"FAIL: {a.log}: {len(rows)} lines, at least {a.min_lines} expected")
    col = {n: names.index(n) for n in names}
    stored = [r[col["stored"]] for r in rows]
    gain = [sum(r[col[n]] for n in added) for r in rows]
    scale = max(s + g + sum(abs(r[col[n]]) for n in out) for s, g, r in zip(stored, gain, rows))
    tol = a.factor * eps * scale
    last = rows[-1]

    if a.no_inflow:
        for n in out:
            v = [r[col[n]] for r in rows]
            low = min(v)
            fall = max([x - y for x, y in zip(v, v[1:])], default=0.0)
            if low < -tol or fall > tol:
                sys.exit(f"FAIL: {n}: water entered through the segment: the column is as low"
                         f" as {low:.6g} m3 and falls by up to {fall:.6g} m3 between two lines"
                         f" (allowed {tol:.3g} m3)")
        excess = max(s - stored[0] - g for s, g in zip(stored, gain))
        if excess > tol:
            sys.exit(f"FAIL: the domain stores {excess:.6g} m3 more than it held at step 0 plus"
                     f" the water added (allowed {tol:.3g} m3)")
        print(f"no water entered through {', '.join(out)}: lowest value "
              f"{min(min(r[col[n]] for r in rows) for n in out):.3g} m3 (allowed {-tol:.3g}),"
              f" stored {stored[-1]:.9g} m3 of {stored[0] + gain[-1]:.9g} m3 held and added")
    else:
        lowest = min(last[col[n]] for n in out)
        need = a.takes_in * (stored[0] + gain[-1])
        if not lowest < 0.0 or stored[-1] < need:
            sys.exit(f"FAIL: the zero_gradient run was to take water in: its lowest out column"
                     f" ends at {lowest:.6g} m3 and it stores {stored[-1]:.6g} m3 (at least"
                     f" {need:.6g} m3 expected)")
        print(f"water entered: lowest out column {lowest:.6g} m3, stored {stored[-1]:.6g} m3"
              f" from {stored[0] + gain[-1]:.6g} m3 held and added")

    if a.expect_stored is not None:
        err = abs(stored[-1] - a.expect_stored)
        if err > a.rtol * abs(a.expect_stored) + tol:
            sys.exit(f"FAIL: stored volume {stored[-1]:.17g} m3, expected {a.expect_stored:.17g}"
                     f" m3 (difference {err:.3g})")
        print(f"stored volume {stored[-1]:.12g} m3 as expected")
    if a.min_out is not None:
        left = sum(last[col[n]] for n in out)
        if left < a.min_out:
            sys.exit(f"FAIL: {left:.6g} m3 left through {', '.join(out)}, at least"
                     f" {a.min_out:.6g} m3 expected")
        print(f"{left:.6g} m3 left through {', '.join(out)} (at least {a.min_out:.6g})")
    print("PASS")


if __name__ == "__main__":
    main()
