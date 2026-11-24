#!/usr/bin/env python3
"""Check a FLUME plotfile of a channel fed by an inflow hydrograph (a
discharge segment) against a closed form: the last plotfile of the run,
<prefix> followed by the largest step.

    check_inflow.py <prefix> --side xlo|xhi|ylo|yhi --q Q --mode MODE
        [--range A B] [--rtol R] [--tol-eps K] [--gravity G]
        [--depth H] [--h0 H0 --u0 U0] [--time T] [--manning N --slope S]

The water enters across the side --side with the discharge Q per unit width
(m2/s). x is the distance of a cell centre from that side (m) and the normal
discharge is the discharge along the normal into the domain. Each mode gives
the depth h(x) and the normal discharge q(x):

  uniform  Still the state the run started from: depth H, discharge Q,
           everywhere. Every cell within K machine epsilons (of the precision
           of the data; default 64) of max(1, value).
  plateau  Subcritical inflow into water of depth H0 moving at U0 (into the
           domain) when the discharge drops to Q at t = 0. A rarefaction runs
           downstream; across it u - 2 c is constant (c = sqrt(g h)), so next
           to the side

               Q / h_b - 2 sqrt(g h_b) = U0 - 2 sqrt(g H0) ,  q = Q ,

           on 0 <= x < (u_b + c_b) t. Q must be below H0 U0 and h_b above the
           critical depth.
  fan      Supercritical inflow onto a dry flat bed without friction. The
           face of the side is at the critical depth h_c = (Q / sqrt(g))^(2/3)
           and moves at c_c = sqrt(g h_c); a rarefaction runs from it to the
           dry front at x = 3 c_c t, across which u + 2 c is constant:

               c = c_c - x / (3 T) ,  h = c^2 / g ,  q = h (c_c + 2 x / (3 T)) .
  normal   Steady uniform flow down a plane bed of slope S (> 0, falling away
           from the side) with Manning's coefficient N:

               h = (Q N / sqrt(S))^(3/5) ,  q = Q .

Modes other than uniform compare the cells whose x is in [A, B] (default:
all) within R (default 1e-3) of the depth and of Q. The discharge along the
side must stay within R Q (uniform: K epsilons of max(1, Q)) everywhere.
"""
import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_dam_break as cd  # noqa: E402


def plateau_depth(q, h0, u0, g):
    """The root h_b of q / h - 2 sqrt(g h) = u0 - 2 sqrt(g h0), by bisection
    (the left side falls as h grows)."""
    r = u0 - 2.0 * math.sqrt(g * h0)
    lo, hi = 1e-12, 1.0
    while q / hi - 2.0 * math.sqrt(g * hi) > r:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if q / mid - 2.0 * math.sqrt(g * mid) > r:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def last_plotfile(prefix):
    """The plotfile <prefix><step> with the largest step, in the current directory."""
    pat = re.compile(re.escape(prefix) + r"(\d{5,})$")
    found = [(int(m.group(1)), f) for f in Path(".").iterdir()
             if f.is_dir() and (m := pat.match(f.name))]
    if not found:
        sys.exit(f"FAIL: no plotfile {prefix}<step> in {Path('.').resolve()}")
    return max(found)[1]


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("prefix")
    p.add_argument("--side", choices=["xlo", "xhi", "ylo", "yhi"], required=True)
    p.add_argument("--q", type=float, required=True)
    p.add_argument("--mode", choices=["uniform", "plateau", "fan", "normal"], required=True)
    p.add_argument("--range", type=float, nargs=2, default=[-math.inf, math.inf])
    p.add_argument("--rtol", type=float, default=1e-3)
    p.add_argument("--tol-eps", type=float, default=64.0)
    p.add_argument("--gravity", type=float, default=9.81)
    p.add_argument("--depth", type=float)
    p.add_argument("--h0", type=float)
    p.add_argument("--u0", type=float)
    p.add_argument("--time", type=float)
    p.add_argument("--manning", type=float)
    p.add_argument("--slope", type=float)
    a = p.parse_args()
    g = a.gravity

    plotfile = last_plotfile(a.prefix)
    print(f"inflow: {plotfile}, time {cd.cp.read_header(plotfile / 'Header')['time']!r} s")
    h, cells, real_size = cd.read(plotfile)
    eps = 2.0 ** -52 if real_size == 8 else 2.0 ** -23
    d = 0 if a.side[0] == "x" else 1
    high = a.side.endswith("hi")
    into = -1.0 if high else 1.0
    lo, hi = h["prob_lo"][d], h["prob_hi"][d]
    dx = h["dx"][d]

    if a.mode == "uniform":
        if a.depth is None:
            sys.exit("FAIL: --mode uniform needs --depth")
        expected = lambda x: (a.depth, a.q)  # noqa: E731
    elif a.mode == "plateau":
        if a.h0 is None or a.u0 is None:
            sys.exit("FAIL: --mode plateau needs --h0 and --u0")
        hb = plateau_depth(a.q, a.h0, a.u0, g)
        hc = (a.q / math.sqrt(g)) ** (2.0 / 3.0)
        if not (a.q < a.h0 * a.u0 and hb > hc):
            sys.exit(f"FAIL: --mode plateau: Q = {a.q} must be below H0 U0 = {a.h0 * a.u0} and"
                     f" h_b = {hb} above h_c = {hc}")
        print(f"inflow: plateau depth h_b = {hb!r} m, u_b + c_b = {a.q / hb + math.sqrt(g * hb)!r} m/s")
        expected = lambda x: (hb, a.q)  # noqa: E731
    elif a.mode == "fan":
        if a.time is None or not a.time > 0.0:
            sys.exit("FAIL: --mode fan needs --time > 0")
        hc = (a.q / math.sqrt(g)) ** (2.0 / 3.0)
        cc = math.sqrt(g * hc)
        print(f"inflow: critical depth h_c = {hc!r} m, dry front at {3.0 * cc * a.time!r} m")

        if a.range[1] >= 3.0 * cc * a.time:
            sys.exit(f"FAIL: --mode fan: the range must end before the dry front,"
                     f" {3.0 * cc * a.time!r} m")

        def expected(x):
            c = max(cc - x / (3.0 * a.time), 0.0)
            return c * c / g, c * c / g * (cc + 2.0 * x / (3.0 * a.time))
    else:
        if a.manning is None or a.slope is None or not a.slope > 0.0:
            sys.exit("FAIL: --mode normal needs --manning and --slope > 0")
        hn = (a.q * a.manning / math.sqrt(a.slope)) ** 0.6
        print(f"inflow: normal depth h_n = {hn!r} m")
        expected = lambda x: (hn, a.q)  # noqa: E731

    errors = []
    worst_h = worst_q = worst_t = 0.0
    n = 0
    for (i, j), c in cells.items():
        k = (i, j)[d]
        centre = lo + (k + 0.5) * dx
        x = (hi - centre) if high else (centre - lo)
        qn = into * c[1 + d]
        qt = c[2 - d]
        if a.mode == "uniform":
            tol_t = a.tol_eps * eps * max(1.0, abs(a.q))
        else:
            tol_t = a.rtol * a.q
        worst_t = max(worst_t, abs(qt))
        if not abs(qt) <= tol_t:
            errors.append(f"cell {(i, j)}: discharge along the side {qt!r}, above {tol_t:.3e}")
        if not (a.range[0] <= x <= a.range[1]):
            continue
        n += 1
        he, qe = expected(x)
        if a.mode == "uniform":
            off_h = abs(c[0] - he) / max(1.0, abs(he))
            off_q = abs(qn - qe) / max(1.0, abs(qe))
            ok = off_h <= a.tol_eps * eps and off_q <= a.tol_eps * eps
        else:
            off_h = abs(c[0] - he) / he
            off_q = abs(qn - qe) / a.q
            ok = off_h <= a.rtol and off_q <= a.rtol
        worst_h = max(worst_h, off_h)
        worst_q = max(worst_q, off_q)
        if not ok:
            errors.append(f"cell {(i, j)} at x = {x} m: depth {c[0]!r}, discharge {qn!r};"
                          f" expected {he!r}, {qe!r}")
    if n == 0:
        sys.exit(f"FAIL: no cell has x in {a.range}")
    unit = "eps" if a.mode == "uniform" else ""
    scale = eps if a.mode == "uniform" else 1.0
    print(f"inflow: {n} cells compared, largest relative error {worst_h / scale:.3g} {unit} in the"
          f" depth and {worst_q / scale:.3g} {unit} in the discharge; along the side"
          f" {worst_t:.3e} m2/s")
    if errors:
        for e in errors[:10]:
            print(e)
        sys.exit(f"FAIL: {len(errors)} cell(s) off the {a.mode} solution")
    print("PASS")


if __name__ == "__main__":
    main()
