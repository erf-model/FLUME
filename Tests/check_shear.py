#!/usr/bin/env python3
"""Check a FLUME run of a shear layer carried by a uniform flow.

    check_shear.py <prefix> --time T --direction x|y --position X --depth H
        --speed U --along-lo A --along-hi B [--cells C] [--width-max W]

The run is periodic in both directions, on a flat bed at 0. The water is H
deep and moves at U (not 0; either sign) along the given direction. The velocity along the
other direction is A where the coordinate is below X and B beyond, so there is
a shear layer at X and another at the periodic side. Nothing but that velocity
changes, and the flow carries its two jumps a distance U T.

The run directory must hold exactly two plotfiles of the prefix: <prefix>00000
and the one of the last step. Checked in the last one (Tests/check_plotfile.py
holds the reader):

  * its time is T;
  * the depth is H and the discharge along the flow is H U in every cell,
    within 4 machine epsilons (of the precision of the data) of the value;
  * nothing varies across the flow direction;
  * the velocity across the flow stays between A and B in every cell, within
    4 epsilons: the transport makes no new extremum;
  * the sum of the discharge across the flow has not changed by more than 16
    epsilons of H (|A| + |B|) times the number of cells;
  * the velocity across the flow passes (A + B) / 2 twice, within C cells
    (default 1) of X + U T and of the periodic side moved by U T;
  * with --width-max W: along the flow, at most W cells hold a velocity
    across the flow that is more than a tenth of |B - A| from both A and B
    (the two layers together). The count is printed in every run.
"""
import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_dam_break as cd  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("prefix")
    p.add_argument("--time", type=float, required=True)
    p.add_argument("--direction", choices=["x", "y"], required=True)
    p.add_argument("--position", type=float, required=True)
    p.add_argument("--depth", type=float, required=True)
    p.add_argument("--speed", type=float, required=True)
    p.add_argument("--along-lo", type=float, required=True)
    p.add_argument("--along-hi", type=float, required=True)
    p.add_argument("--cells", type=float, default=1.0)
    p.add_argument("--width-max", type=int)
    p.add_argument("--time-rtol", type=float, default=1e-12)
    a = p.parse_args()

    found = sorted(d.name for d in Path(".").iterdir()
                   if d.is_dir() and re.fullmatch(re.escape(a.prefix) + r"\d{5,}", d.name))
    if len(found) != 2 or found[0] != f"{a.prefix}00000":
        sys.exit(f"FAIL: expected {a.prefix}00000 and one later plotfile, found {found}")
    _, first, _ = cd.read(Path(found[0]))
    head, last, real_size = cd.read(Path(found[1]))
    eps = 2.0 ** -52 if real_size == 8 else 2.0 ** -23
    errors = []
    if not abs(head["time"] - a.time) <= a.time_rtol * a.time:
        errors.append(f"time: got {head['time']!r}, expected {a.time!r}")

    d = 0 if a.direction == "x" else 1
    n = head["domain_hi"][d] + 1
    n_across = head["domain_hi"][1 - d] + 1
    dx = head["dx"][d]
    lo = head["prob_lo"][d]
    length = n * dx
    i_flow = 1 + d    # discharge along the flow
    i_cross = 2 - d   # discharge across it

    def cell(cells, m, t):
        return cells[(m, t)] if d == 0 else cells[(t, m)]

    h, q = a.depth, a.depth * a.speed
    v_min, v_max = min(a.along_lo, a.along_hi), max(a.along_lo, a.along_hi)
    v_scale = max(abs(v_min), abs(v_max))
    worst = {"depth": 0.0, "flow": 0.0, "across": 0.0, "bound": 0.0}
    for m in range(n):
        ref = cell(last, m, 0)
        for t in range(n_across):
            v = cell(last, m, t)
            if not all(math.isfinite(x) for x in v):
                sys.exit(f"FAIL: cell ({m}, {t}) holds a value that is not finite")
            worst["depth"] = max(worst["depth"], abs(v[0] - h) / h)
            worst["flow"] = max(worst["flow"], abs(v[i_flow] - q) / abs(q))
            worst["across"] = max(worst["across"], abs(v[0] - ref[0]), abs(v[i_flow] - ref[i_flow]),
                                  abs(v[i_cross] - ref[i_cross]))
            vel = v[i_cross] / v[0]
            worst["bound"] = max(worst["bound"], (vel - v_max) / v_scale, (v_min - vel) / v_scale)
    if not worst["depth"] <= 4 * eps:
        errors.append(f"the depth differs from {h} by {worst['depth']:.3g} of it")
    if not worst["flow"] <= 4 * eps:
        errors.append(f"the discharge along the flow differs from {q} by {worst['flow']:.3g} of it")
    if not worst["across"] <= 4 * eps * max(h, abs(q), h * v_scale):
        errors.append(f"the state varies across the flow direction by {worst['across']:.3g}")
    if not worst["bound"] <= 4 * eps:
        errors.append(f"the velocity across the flow leaves [{v_min}, {v_max}] by {worst['bound']:.3g}"
                      f" of {v_scale}")

    s0 = math.fsum(v[i_cross] for v in first.values())
    s1 = math.fsum(v[i_cross] for v in last.values())
    s_tol = 16 * eps * h * (abs(a.along_lo) + abs(a.along_hi)) * n * n_across
    if not abs(s1 - s0) <= s_tol:
        errors.append(f"the discharge across the flow sums to {s1!r}, at the start to {s0!r}")

    # Where the velocity across the flow passes the middle value, around the
    # periodic domain.
    mid = 0.5 * (a.along_lo + a.along_hi)
    vel = [cell(last, m, 0)[i_cross] / cell(last, m, 0)[0] - mid for m in range(n)]
    crossings = []
    for m in range(n):
        v0, v1 = vel[m], vel[(m + 1) % n]
        if v0 * v1 <= 0.0 and v0 != v1:
            crossings.append((lo + (m + 0.5 + v0 / (v0 - v1)) * dx - lo) % length + lo)
    margin = 0.1 * (v_max - v_min)
    width = sum(1 for v in vel if v_min + margin < v + mid < v_max - margin)
    if a.width_max is not None and not width <= a.width_max:
        errors.append(f"{width} cells are inside the two shear layers, at most {a.width_max}")
    want = sorted((x + a.speed * a.time - lo) % length + lo for x in (a.position, lo))
    crossings.sort()
    if len(crossings) != 2:
        errors.append(f"the velocity across the flow passes {mid} at {crossings}; expected two places, {want}")
    else:
        for got, exact in zip(crossings, want):
            off = abs(got - exact)
            off = min(off, length - off)
            if not off <= a.cells * dx:
                errors.append(f"shear layer at {got:.6g} m, expected {exact:.6g} m within {a.cells} cell(s)")

    report = (f"depth and discharge along the flow within {max(worst['depth'], worst['flow']) / eps:.3g}"
              f" epsilons; shear layers at {[round(c, 4) for c in crossings]} m (exact {want}),"
              f" {width} cells inside them")
    if errors:
        print(f"FAIL: {a.prefix}")
        for e in errors:
            print("  " + e)
        print("  (" + report + ")")
        sys.exit(1)
    print(f"PASS: {found[1]}: {report}")


if __name__ == "__main__":
    main()
