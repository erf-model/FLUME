#!/usr/bin/env python3
"""Check a FLUME run of water moving in a paraboloid bowl (Thacker's solutions).

    check_bowl.py <prefix> --shape planar|radial --bowl-depth H0 --bowl-radius A
        ( --eta E | --r0 R ) ( --periods N --l1-max EH EQ | --initial-only [--initial-phase P] )
        [--l1-min EH] [--gravity G] [--mass-tol-eps K]
        [--time-rtol R] [--max-steps S] [--dry-depth T]

The bed is z = -H0 (1 - r^2 / A^2), with r the distance from the centre of the
domain (flume.bed.bowl_depth, flume.bed.bowl_radius on a bed at 0). The water
has a moving shoreline, and the solutions are periodic in time (Delestre et
al., 2013, section 4.2.2). With X, Y measured from the centre:

  planar (--eta E): a plane surface that turns about the axis of the bowl,
      w = sqrt(2 G H0) / A
      h = max(E H0 / A^2 (2 X cos wt + 2 Y sin wt - E) - z, 0)
      u = -E w sin wt,  v = E w cos wt
  radial (--r0 R): a paraboloid surface that rises and falls,
      w = sqrt(8 G H0) / A,  B = (A^2 - R^2) / (A^2 + R^2),  D = 1 - B cos wt
      h = max(H0 (sqrt(1 - B^2) / D - 1 - r^2 / A^2 ((1 - B^2) / D^2 - 1)) - z, 0)
      u = X w B sin wt / (2 D),  v = Y w B sin wt / (2 D)

With --initial-only the run takes no step: there is one plotfile, and only the
bed and the initial state are checked (for any --bowl-depth, also negative),
against the solution at P periods (default 0).

The run directory must hold exactly two plotfiles of the prefix: <prefix>00000
and the one of the last step, whose number is not pinned. Checked:

  * the bed in both plotfiles is the bowl, within 8 epsilons of H0 times
    1 + r^2 / A^2;
  * the first plotfile is the solution at t = 0, within 8 epsilons likewise;
  * the time of the last plotfile is N periods, 2 pi N / w, within R;
  * no value is not finite and no depth is negative;
  * mass: the sum of the depth over all cells is the same in the two
    plotfiles within K epsilons times the sum;
  * L1 error: the mean of |h - h_exact| over the cells where the run or the
    solution holds water at the end, or the run at the start, divided by the mean of the initial depth there, is at
    most the first number of --l1-max (and at least --l1-min, when given: the
    run did not just keep the initial state, which is also the answer after
    whole periods); the mean of |q - q_exact|, both components, divided by the
    mean initial depth times sqrt(G H0), is at most the second;
  * with --dry-depth T: a cell of the last plotfile no deeper than T holds no
    discharge (flume.dry_depth);
  * with --max-steps S: the last step number is at most S (a guard on the
    step count).
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
    p.add_argument("--shape", choices=["planar", "radial"], required=True)
    p.add_argument("--bowl-depth", type=float, required=True)
    p.add_argument("--bowl-radius", type=float, required=True)
    p.add_argument("--eta", type=float)
    p.add_argument("--r0", type=float)
    p.add_argument("--periods", type=float)
    p.add_argument("--l1-max", type=float, nargs=2)
    p.add_argument("--initial-only", action="store_true")
    p.add_argument("--initial-phase", type=float, default=0.0)
    p.add_argument("--l1-min", type=float)
    p.add_argument("--gravity", type=float, default=9.81)
    p.add_argument("--mass-tol-eps", type=float, default=64.0)
    p.add_argument("--time-rtol", type=float, default=1e-9)
    p.add_argument("--max-steps", type=int)
    p.add_argument("--dry-depth", type=float)
    a = p.parse_args()
    if (a.shape == "planar") != (a.eta is not None) or (a.shape == "radial") != (a.r0 is not None):
        sys.exit("check_bowl.py: --shape planar needs --eta, --shape radial needs --r0")

    if a.initial_only == (a.periods is not None) or a.initial_only == (a.l1_max is not None):
        sys.exit("check_bowl.py: give --periods and --l1-max, or --initial-only")

    h0, rad, g = a.bowl_depth, a.bowl_radius, a.gravity
    if a.shape == "planar":
        w = math.sqrt(2.0 * g * abs(h0)) / rad
    else:
        w = math.sqrt(8.0 * g * abs(h0)) / rad
        big_b = (rad * rad - a.r0 * a.r0) / (rad * rad + a.r0 * a.r0)
    t_end = 0.0 if a.initial_only else a.periods * 2.0 * math.pi / w

    def bowl(x, y):
        return -h0 * (1.0 - (x * x + y * y) / (rad * rad))

    def exact(x, y, t):
        """Surface elevation (which may be below the bed) and velocity."""
        if a.shape == "planar":
            s = a.eta * h0 / (rad * rad) * (2.0 * x * math.cos(w * t) + 2.0 * y * math.sin(w * t) - a.eta)
            return s, -a.eta * w * math.sin(w * t), a.eta * w * math.cos(w * t)
        d = 1.0 - big_b * math.cos(w * t)
        s = h0 * (math.sqrt(1.0 - big_b * big_b) / d - 1.0
                  - (x * x + y * y) / (rad * rad) * ((1.0 - big_b * big_b) / (d * d) - 1.0))
        f = 0.5 * w * big_b * math.sin(w * t) / d
        return s, x * f, y * f

    found = sorted(d.name for d in Path(".").iterdir()
                   if d.is_dir() and re.fullmatch(re.escape(a.prefix) + r"\d{5,}", d.name))
    n_expected = 1 if a.initial_only else 2
    if len(found) != n_expected or found[0] != f"{a.prefix}00000":
        sys.exit(f"FAIL: expected {a.prefix}00000 and {n_expected - 1} later plotfile(s), found {found}")
    step = int(found[-1][len(a.prefix):])
    head0, first, real_size = cd.read(Path(found[0]))
    head, last, _ = cd.read(Path(found[-1]))
    eps = 2.0 ** -52 if real_size == 8 else 2.0 ** -23
    errors = []

    nx = head["domain_hi"][0] + 1
    ny = head["domain_hi"][1] + 1
    dx, dy = head["dx"][0], head["dx"][1]

    def centre(c):
        return (c[0] + 0.5 - 0.5 * nx) * dx, (c[1] + 0.5 - 0.5 * ny) * dy

    if head["step"] != step:
        errors.append(f"step: got {head['step']}, expected {step}")
    if not abs(head["time"] - t_end) <= a.time_rtol * t_end:
        errors.append(f"time: got {head['time']!r}, expected {t_end!r}")
    if a.max_steps is not None and step > a.max_steps:
        errors.append(f"the run took {step} steps, at most {a.max_steps} expected")

    worst_bed = 0.0
    worst_init = 0.0
    speed = w * max(abs(head["prob_hi"][0] - head["prob_lo"][0]), abs(head["prob_hi"][1] - head["prob_lo"][1]))
    for c, v in first.items():
        x, y = centre(c)
        z = bowl(x, y)
        scale = abs(h0) * (1.0 + (x * x + y * y) / (rad * rad))
        worst_bed = max(worst_bed, abs(v[3] - z) / scale, abs(last[c][3] - z) / scale)
        s, u, vv = exact(x, y, a.initial_phase * 2.0 * math.pi / w)
        h = max(s - v[3], 0.0)
        worst_init = max(worst_init, abs(v[0] - h) / scale, abs(v[1] - h * u) / (scale * (1.0 + speed)),
                         abs(v[2] - h * vv) / (scale * (1.0 + speed)))
    if not worst_bed <= 8.0 * eps:
        errors.append(f"the bed is off the bowl by {worst_bed / eps:.3g} epsilons, at most 8")
    if not worst_init <= 8.0 * eps:
        errors.append(f"the first plotfile is off the solution at t = 0 by {worst_init / eps:.3g}"
                      " epsilons, at most 8")
    if a.initial_only:
        if errors:
            print(f"FAIL: {a.prefix}")
            for e in errors:
                print("  " + e)
            sys.exit(1)
        print(f"PASS: {found[0]}: bed within {worst_bed / eps:.3g} and state within"
              f" {worst_init / eps:.3g} epsilons of the closed form,"
              f" {sum(1 for v in first.values() if v[0] > 0.0)} wet cells of {len(first)}")
        return

    lowest = math.inf
    for c, v in last.items():
        if not all(math.isfinite(q) for q in v):
            errors.append(f"cell {c} holds a value that is not finite")
            break
        lowest = min(lowest, v[0])
    if not lowest >= 0.0:
        errors.append(f"the lowest depth is {lowest!r}: negative")

    if a.dry_depth is not None:
        moving = [c for c, v in last.items() if v[0] <= a.dry_depth and (v[1] != 0.0 or v[2] != 0.0)]
        if moving:
            errors.append(f"{len(moving)} cell(s) no deeper than {a.dry_depth} m hold a discharge,"
                          f" e.g. cell {moving[0]}: {last[moving[0]][:3]}")

    m0 = math.fsum(v[0] for v in first.values())
    m1 = math.fsum(v[0] for v in last.values())
    if not abs(m1 - m0) <= a.mass_tol_eps * eps * m0:
        errors.append(f"mass: the depths sum to {m1!r}, at the start to {m0!r}: difference"
                      f" {m1 - m0:.3g}, tolerance {a.mass_tol_eps * eps * m0:.3g}")

    sum_h = 0.0
    sum_q = 0.0
    sum_h0 = 0.0
    n_wet = 0
    for c, v in last.items():
        x, y = centre(c)
        s, u, vv = exact(x, y, t_end)
        h = max(s - v[3], 0.0)
        if h > 0.0 or v[0] > 0.0 or first[c][0] > 0.0:
            n_wet += 1
            sum_h += abs(v[0] - h)
            sum_q += abs(v[1] - h * u) + abs(v[2] - h * vv)
            sum_h0 += first[c][0]
    l1_h = sum_h / sum_h0
    l1_q = sum_q / (sum_h0 * math.sqrt(g * h0))
    if not l1_h <= a.l1_max[0]:
        errors.append(f"L1 error of the depth {l1_h:.4g}, at most {a.l1_max[0]}")
    if not l1_q <= a.l1_max[1]:
        errors.append(f"L1 error of the discharge {l1_q:.4g}, at most {a.l1_max[1]}")
    if a.l1_min is not None and not l1_h >= a.l1_min:
        errors.append(f"L1 error of the depth {l1_h:.4g}, expected at least {a.l1_min}:"
                      " did the water move?")

    if errors:
        print(f"FAIL: {a.prefix}")
        for e in errors:
            print("  " + e)
        sys.exit(1)
    print(f"PASS: {found[1]}: {step} steps, {n_wet} wet cells, L1 error {l1_h:.4g} (depth)"
          f" {l1_q:.4g} (discharge), lowest depth {lowest:.3g}, mass change"
          f" {(m1 - m0) / (eps * m0):.3g} epsilons, bed within {worst_bed / eps:.3g} and start"
          f" within {worst_init / eps:.3g} epsilons")


if __name__ == "__main__":
    main()
