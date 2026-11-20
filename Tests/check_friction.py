#!/usr/bin/env python3
"""Check a FLUME run of uniform flow down an inclined periodic channel with
Manning friction against its closed-form solution.

    check_friction.py <prefix> --time T --depth H --slope SX SY --manning N
        --q-rtol R [--q-rtol-eps E] [--q-rtol-min R0] [--speed0 U0] [--gravity G]
        [--depth-tol-eps K] [--time-rtol R]

The run has a plane bed with slope S = (SX, SY) (dz/dx, dz/dy), periodic
sides with flume.bed.continue_slope = 1 along every direction in which S is
not 0, and water of uniform depth H (the surface parallel to the bed). The
water starts at the speed U0 down the slope. Fluxes are then the same on
every face, the depth stays H, and the discharge q (per unit width, along
-S / |S|) obeys

    dq/dt = g H |S| - g N^2 q^2 / H^(7/3) .

With the normal discharge q_n = H^(5/3) |S|^(1/2) / N (Manning's formula)
and T_r = q_n / (g H |S|), from q(0) = q_0 = H U0:

    q_0 < q_n:  q = q_n tanh(t / T_r + atanh(q_0 / q_n)) ,
    q_0 > q_n:  q = q_n coth(t / T_r + acoth(q_0 / q_n)) .

The run directory must hold exactly two plotfiles of the prefix: <prefix>00000
and the one of the last step. Checked, in every cell:

  * the time of the last plotfile is T, within R times T;
  * the depth is H in both plotfiles, within K epsilons of H + |bed| (the
    initial depth is a surface minus a bed);
  * the first plotfile holds q_0 along the slope (within K epsilons of
    q_0 + H);
  * |q - q_exact| / q_exact is at most R + E epsilons, q_exact the vector of
    size q(T) down the slope (and the largest such error over the cells is
    at least R0 - E epsilons, when --q-rtol-min is given). E covers rounding
    in single and double precision.
"""
import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_dam_break as cd  # noqa: E402


def exact(t, q0, qn, tr):
    if q0 == qn:
        return qn
    if q0 < qn:
        return qn * math.tanh(t / tr + math.atanh(q0 / qn))
    return qn / math.tanh(t / tr + 0.5 * math.log((q0 + qn) / (q0 - qn)))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("prefix")
    p.add_argument("--time", type=float, required=True)
    p.add_argument("--depth", type=float, required=True)
    p.add_argument("--slope", type=float, nargs=2, required=True)
    p.add_argument("--manning", type=float, required=True)
    p.add_argument("--q-rtol", type=float, required=True)
    p.add_argument("--q-rtol-eps", type=float, default=0.0)
    p.add_argument("--q-rtol-min", type=float)
    p.add_argument("--speed0", type=float, default=0.0)
    p.add_argument("--gravity", type=float, default=9.81)
    p.add_argument("--depth-tol-eps", type=float, default=64.0)
    p.add_argument("--time-rtol", type=float, default=1e-12)
    a = p.parse_args()

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

    g = a.gravity
    h = a.depth
    s = math.hypot(*a.slope)
    if not (s > 0.0 and a.manning > 0.0 and h > 0.0):
        sys.exit("FAIL: the check needs a slope, a Manning coefficient and a depth above 0")
    if not a.speed0 >= 0.0:
        sys.exit("FAIL: the closed form needs --speed0 >= 0 (water starting down the slope)")
    down = [-a.slope[0] / s, -a.slope[1] / s]
    qn = h ** (5.0 / 3.0) * math.sqrt(s) / a.manning
    tr = qn / (g * h * s)
    q0 = h * a.speed0
    qt = exact(a.time, q0, qn, tr)
    print(f"normal discharge {qn:.9g} m2/s, time scale {tr:.9g} s, "
          f"q(0) {q0:.9g}, q(T) {qt:.9g} m2/s")

    q0_tol = a.depth_tol_eps * eps * (q0 + h)
    q_rtol = a.q_rtol + a.q_rtol_eps * eps
    worst = 0.0
    for name, cells, q_want in (("first", first, q0), ("last", last, qt)):
        for (i, j), v in sorted(cells.items()):
            h_tol = a.depth_tol_eps * eps * (h + abs(v[3]))
            if not abs(v[0] - h) <= h_tol:
                errors.append(f"{name} plotfile, cell ({i}, {j}): depth {v[0]!r}, expected {h!r}")
                continue
            want = (q_want * down[0], q_want * down[1])
            err = math.hypot(v[1] - want[0], v[2] - want[1])
            if name == "first":
                if not err <= q0_tol:
                    errors.append(f"first plotfile, cell ({i}, {j}): discharge ({v[1]!r}, {v[2]!r}),"
                                  f" expected ({want[0]!r}, {want[1]!r})")
                continue
            rel = err / q_want
            worst = max(worst, rel) if not math.isnan(rel) else math.inf
            if not rel <= q_rtol:
                errors.append(f"last plotfile, cell ({i}, {j}): discharge ({v[1]!r}, {v[2]!r}),"
                              f" expected ({want[0]!r}, {want[1]!r}): relative error {rel:.3e}"
                              f" above {q_rtol:.3e}")
    print(f"largest relative error of the discharge: {worst:.6e}")
    if a.q_rtol_min is not None and not worst >= a.q_rtol_min - a.q_rtol_eps * eps:
        errors.append(f"largest relative error {worst:.6e} is below the lower bound"
                      f" {a.q_rtol_min - a.q_rtol_eps * eps:.6e}")

    if errors:
        for e in errors[:20]:
            print("FAIL: " + e)
        if len(errors) > 20:
            print(f"FAIL: ... and {len(errors) - 20} more")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
