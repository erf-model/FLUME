#!/usr/bin/env python3
"""Check a FLUME run of a thin film of uneven depth down an inclined periodic
channel with stiff Manning friction against the local Manning relation.

    check_film_local_manning.py <prefix> --time T --bed-slope SX --manning N
        --rtol R [--mass-tol-eps K] [--time-rtol R]

The run has a plane bed falling along x with slope SX < 0, periodic sides
with flume.bed.continue_slope = 1, and a film whose depth varies along x
only. When the friction time scale is much shorter than the step, friction
balances the slope of the water surface in every cell, and the discharge
along x is close to the local Manning discharge

    q = h^(5/3) S^(1/2) / N ,  S = -(eta(i+1) - eta(i-1)) / (2 dx) ,

with eta the surface elevation (the bed across the seam continued by the
slope). The inertia of the film and the discretisation make up the rest.

The run directory must hold exactly two plotfiles of the prefix: <prefix>00000
and the one of the last step. Checked:

  * the time of the last plotfile is T, within R times T;
  * every depth is above 0 and every value finite, in both plotfiles;
  * the mass is the same in both, within K epsilons;
  * in every cell of the last plotfile, |hu - q| / q is at most R, and hv
    is at most R q.
"""
import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_dam_break as cd  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("prefix")
    p.add_argument("--time", type=float, required=True)
    p.add_argument("--bed-slope", type=float, required=True)
    p.add_argument("--manning", type=float, required=True)
    p.add_argument("--rtol", type=float, required=True)
    p.add_argument("--mass-tol-eps", type=float, default=64.0)
    p.add_argument("--time-rtol", type=float, default=1e-12)
    a = p.parse_args()
    if not (a.bed_slope < 0.0 and a.manning > 0.0):
        sys.exit("FAIL: the check needs a bed falling along x and a Manning coefficient above 0")

    found = sorted(d.name for d in Path(".").iterdir()
                   if d.is_dir() and re.fullmatch(re.escape(a.prefix) + r"\d{5,}", d.name))
    if len(found) != 2 or found[0] != f"{a.prefix}00000":
        sys.exit(f"FAIL: expected {a.prefix}00000 and one later plotfile, found {found}")
    _, first, real_size = cd.read(Path(found[0]))
    head, last, _ = cd.read(Path(found[1]))
    eps = 2.0 ** -52 if real_size == 8 else 2.0 ** -23
    errors = []

    if not abs(head["time"] - a.time) <= a.time_rtol * a.time:
        errors.append(f"time {head['time']!r}, expected {a.time!r}")

    for name, cells in (("first", first), ("last", last)):
        for (i, j), v in sorted(cells.items()):
            if not (all(math.isfinite(x) for x in v) and v[0] > 0.0):
                errors.append(f"{name} plotfile, cell ({i}, {j}): values {v!r}")
    if errors:
        for e in errors[:20]:
            print("FAIL: " + e)
        sys.exit(1)

    m0 = math.fsum(v[0] for v in first.values())
    m1 = math.fsum(v[0] for v in last.values())
    print(f"mass {m0!r} then {m1!r}")
    if not abs(m1 - m0) <= a.mass_tol_eps * eps * m0:
        errors.append(f"mass changed by {(m1 - m0) / m0:.3e} of itself")

    xs = sorted({i for i, _ in last})
    nx = len(xs)
    lo = xs[0]
    dx = head["dx"][0]

    def eta(i, j):
        k = (i - lo) % nx
        shift = (i - lo - k) // nx
        v = last[(lo + k, j)]
        return v[0] + v[3] + a.bed_slope * dx * nx * shift

    worst = 0.0
    for (i, j), v in sorted(last.items()):
        s = -(eta(i + 1, j) - eta(i - 1, j)) / (2.0 * dx)
        if not s > 0.0:
            errors.append(f"cell ({i}, {j}): surface slope {s!r} does not fall along x")
            continue
        q = v[0] ** (5.0 / 3.0) * math.sqrt(s) / a.manning
        rel = max(abs(v[1] - q), abs(v[2])) / q
        worst = max(worst, rel)
        if not rel <= a.rtol:
            errors.append(f"cell ({i}, {j}): discharge ({v[1]!r}, {v[2]!r}), local Manning"
                          f" {q!r}: relative difference {rel:.3e} above {a.rtol:.3e}")
    print(f"largest relative difference from the local Manning discharge: {worst:.6e}")

    if errors:
        for e in errors[:20]:
            print("FAIL: " + e)
        if len(errors) > 20:
            print(f"FAIL: ... and {len(errors) - 20} more")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
