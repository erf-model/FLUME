#!/usr/bin/env python3
"""Check a FLUME run of a smooth simple wave against its exact solution.

    check_wave.py <prefix> --direction x|y --surface S --amplitude A --time T
        --l1-max EH EQ [--l1-min EH] [--velocity U V] [--gravity G]
        [--mass-tol-eps K] [--time-rtol R]

The run starts from flume.init.wave_amplitude over a flat bed, in a domain
that is periodic along the wave: with s the distance from the low side of the
domain over its length L, b the bed, H = S - b and c = sqrt(G h),

    h = H + A sin(2 pi s),   u = U + 2 (c - c_H)   along the wave,   v = V across it.

This is a simple wave: u - 2 c is the same everywhere, and each depth travels
unchanged at the speed u + c = U + 3 c - 2 c_H (the shallow-water equations in
characteristic form; Docs/sphinx_doc/Testing.rst). At the time T the depth at
x is the initial depth at x0, where x0 + (U + 3 c(x0) - 2 c_H) T = x. The
solution is smooth until the wave breaks, at

    T_b = L sqrt(H - A') / (3 pi A sqrt(G)),   A' between -A and A

(the steepest point of the front); the check stops if T is not below the
smallest such time, L sqrt(H - |A|) / (3 pi |A| sqrt(G)). x0 is found by bisection.

The run directory must hold exactly two plotfiles of the prefix: <prefix>00000
and the one of the last step. Checked:

  * the bed is flat, and the first plotfile holds the wave within 64 epsilons
    of the depth (times 1 + the speed, for the discharges);
  * the time of the last plotfile is T, within R times T;
  * the state does not vary across the wave: within 64 epsilons of H;
  * mass: the sum of the depth is the same in the two plotfiles within K
    epsilons times the sum;
  * L1 error: the mean over the cells of |h - h_exact|, divided by |A|, is at
    most the first number of --l1-max (and at least --l1-min, when given);
    the mean of |q - q_exact| along the wave, divided by |A| c_H, is at most
    the second. The discharge across the wave is h V, so its error, divided
    the same way, is at most the first number times |V| / c_H (plus 64
    epsilons). The exact values are those at the cell centres, as is the
    initial state; the mean over a cell differs by 1e-4 of A on 100 cells.
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
    p.add_argument("--direction", choices=["x", "y"], required=True)
    p.add_argument("--surface", type=float, required=True)
    p.add_argument("--amplitude", type=float, required=True)
    p.add_argument("--time", type=float, required=True)
    p.add_argument("--l1-max", type=float, nargs=2, required=True)
    p.add_argument("--l1-min", type=float)
    p.add_argument("--velocity", type=float, nargs=2, default=[0.0, 0.0])
    p.add_argument("--gravity", type=float, default=9.81)
    p.add_argument("--mass-tol-eps", type=float, default=16.0)
    p.add_argument("--time-rtol", type=float, default=1e-12)
    a = p.parse_args()

    found = sorted(d.name for d in Path(".").iterdir()
                   if d.is_dir() and re.fullmatch(re.escape(a.prefix) + r"\d{5,}", d.name))
    if len(found) != 2 or found[0] != f"{a.prefix}00000":
        sys.exit(f"FAIL: expected {a.prefix}00000 and one later plotfile, found {found}")
    step = int(found[1][len(a.prefix):])
    head, first, real_size = cd.read(Path(found[0]))
    head1, last, _ = cd.read(Path(found[1]))
    eps = 2.0 ** -52 if real_size == 8 else 2.0 ** -23
    errors = []

    d = 0 if a.direction == "x" else 1
    n = head["domain_hi"][d] + 1
    length = head["prob_hi"][d] - head["prob_lo"][d]
    g = a.gravity
    bed = first[(0, 0)][3]
    big_h = a.surface - bed
    amp = a.amplitude
    if not big_h - abs(amp) > 0.0:
        sys.exit("FAIL: the wave reaches the bed: the check needs water everywhere")
    c_ref = math.sqrt(g * big_h)
    t_break = length * math.sqrt(big_h - abs(amp)) / (3.0 * math.pi * abs(amp) * math.sqrt(g))
    if not a.time < t_break:
        sys.exit(f"FAIL: --time {a.time} is not below the breaking time {t_break:.6g} s")
    u_along = a.velocity[d]
    u_across = a.velocity[1 - d]

    def depth0(s):
        return big_h + amp * math.sin(2.0 * math.pi * s)

    def speed(s):
        return u_along + 3.0 * math.sqrt(g * depth0(s)) - 2.0 * c_ref

    def exact(s, t):
        """Depth and velocity along the wave at the fraction s of the domain."""
        if t == 0.0:
            h = depth0(s)
        else:
            # s0 + speed(s0) t / L = s, increasing in s0 before the wave breaks.
            lo = s - (u_along + 3.0 * math.sqrt(g * (big_h + abs(amp))) - 2.0 * c_ref) * t / length
            hi = s - (u_along + 3.0 * math.sqrt(g * (big_h - abs(amp))) - 2.0 * c_ref) * t / length
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if mid + speed(mid) * t / length < s:
                    lo = mid
                else:
                    hi = mid
            h = depth0(0.5 * (lo + hi))
        return h, u_along + 2.0 * (math.sqrt(g * h) - c_ref)

    if head1["step"] != step:
        errors.append(f"step: got {head1['step']}, expected {step}")
    if not abs(head1["time"] - a.time) <= a.time_rtol * a.time:
        errors.append(f"time: got {head1['time']!r}, expected {a.time!r}")

    scale_q = 1.0 + abs(u_along) + abs(u_across) + c_ref
    worst_init = 0.0
    worst_bed = 0.0
    for c, v in first.items():
        worst_bed = max(worst_bed, abs(v[3] - bed), abs(last[c][3] - bed))
        h, u = exact((c[d] + 0.5) / n, 0.0)
        q = [h * u, h * u_across] if d == 0 else [h * u_across, h * u]
        worst_init = max(worst_init, abs(v[0] - h) / big_h, abs(v[1] - q[0]) / (big_h * scale_q),
                         abs(v[2] - q[1]) / (big_h * scale_q))
    if worst_bed != 0.0:
        errors.append(f"the bed is not flat: it varies by {worst_bed!r} m")
    if not worst_init <= 64.0 * eps:
        errors.append(f"the first plotfile is off the wave by {worst_init / eps:.3g} epsilons,"
                      " at most 64")

    lowest = math.inf
    worst_across = 0.0
    sum_h = 0.0
    sum_q = 0.0
    sum_qt = 0.0
    for c, v in last.items():
        if not all(math.isfinite(x) for x in v):
            errors.append(f"cell {c} holds a value that is not finite")
            break
        lowest = min(lowest, v[0])
        row = last[(c[0], 0)] if d == 0 else last[(0, c[1])]
        worst_across = max(worst_across, abs(v[0] - row[0]) / big_h,
                           abs(v[1] - row[1]) / (big_h * scale_q),
                           abs(v[2] - row[2]) / (big_h * scale_q))
        h, u = exact((c[d] + 0.5) / n, a.time)
        sum_h += abs(v[0] - h)
        sum_q += abs(v[1 + d] - h * u)
        sum_qt += abs(v[2 - d] - h * u_across)
    if not worst_across <= 64.0 * eps:
        errors.append(f"the state varies across the wave by {worst_across / eps:.3g} epsilons,"
                      " at most 64")
    l1_h = sum_h / (len(last) * abs(amp))
    l1_q = sum_q / (len(last) * abs(amp) * c_ref)
    l1_qt = sum_qt / (len(last) * abs(amp) * c_ref)
    if not l1_h <= a.l1_max[0]:
        errors.append(f"L1 error of the depth {l1_h:.4g}, at most {a.l1_max[0]}")
    if not l1_q <= a.l1_max[1]:
        errors.append(f"L1 error of the discharge {l1_q:.4g}, at most {a.l1_max[1]}")
    max_qt = a.l1_max[0] * abs(u_across) / c_ref + 64.0 * eps
    if not l1_qt <= max_qt:
        errors.append(f"L1 error of the discharge across the wave {l1_qt:.4g},"
                      f" at most {max_qt:.4g}")
    if a.l1_min is not None and not l1_h >= a.l1_min:
        errors.append(f"L1 error of the depth {l1_h:.4g}, expected at least {a.l1_min}")

    m0 = math.fsum(v[0] for v in first.values())
    m1 = math.fsum(v[0] for v in last.values())
    if not abs(m1 - m0) <= a.mass_tol_eps * eps * m0:
        errors.append(f"mass: the depths sum to {m1!r}, at the start to {m0!r}: difference"
                      f" {m1 - m0:.3g}, tolerance {a.mass_tol_eps * eps * m0:.3g}")

    if errors:
        print(f"FAIL: {a.prefix}")
        for e in errors:
            print("  " + e)
        sys.exit(1)
    print(f"PASS: {found[1]}: {step} steps, L1 error {l1_h:.4g} (depth) {l1_q:.4g} (discharge) {l1_qt:.4g} (across)"
          f" of the amplitude, lowest depth {lowest:.6g}, mass change"
          f" {(m1 - m0) / (eps * m0):.3g} epsilons, start within {worst_init / eps:.3g} epsilons")


if __name__ == "__main__":
    main()
