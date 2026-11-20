#!/usr/bin/env python3
"""Check a dry-bed dam break whose high-x side is an outflow on part of it only.

    check_partial_outflow.py <prefix> --time T --position X --depth H
        --split Y --margin M --l1-max E --bore-rtol R [--gravity G]

Still water H deep at x < X over a dry, flat bed; the side at the high end of
x is an outflow for y >= Y and a wall for y < Y (flume.bc.*). The run directory
holds <prefix>00000 and the plotfile of the last step, at time T, after the
front has reached the side (at X + 2 c1 t = the end of the domain, c1 =
sqrt(G H)). Away from y = Y, where the two halves meet, each half is a 1-D
problem of its own:

  * open rows (cell centre y >= Y + M): the exact solution of an unbounded
    domain (Ritter's; Delestre et al., 2013, section 4.1.2), as for the
    outflow over the whole side. Checked: the mean of |h - h_exact| over the
    cells of these rows divided by H, and of |q_x - q_exact| divided by H c1,
    are each at most E.
  * wall rows (cell centre y <= Y - M): the flow arriving at the wall stops in
    a bore that runs back into it. With h, u the arriving state (Ritter's at
    the side) and h2 the depth behind the bore, still water against the wall,
    the jump conditions give

        u = (h2 - h) sqrt(G (h + h2) / (2 h h2)),

    solved for h2 by bisection. The bore runs back at h u / (h2 - h), under
    1 m/s here, so the arriving state changes little while the water behind it
    piles up; the depth of the cell at the wall must be within R h2 of h2.

Also checked: the time, finite values, no negative depth, and the depths sum to
less at T than at the start (water left). Every check fails if the two halves
take the same type, if the segments are resolved to the wrong cells, or if the
pieces of the side are filled in the wrong order.
"""
import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_dam_break import read, ritter  # noqa: E402


def bore_depth(h, u, g):
    """Depth behind a bore that brings water arriving at (h, u) to rest at a wall."""
    lo, hi = h, h + 100.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (mid - h) * math.sqrt(g * (h + mid) / (2.0 * h * mid)) < u:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("prefix")
    for name in ("--time", "--position", "--depth", "--split", "--margin", "--l1-max", "--bore-rtol"):
        p.add_argument(name, type=float, required=True)
    p.add_argument("--gravity", type=float, default=9.81)
    a = p.parse_args()

    found = sorted(d.name for d in Path(".").iterdir()
                   if d.is_dir() and re.fullmatch(re.escape(a.prefix) + r"\d{5,}", d.name))
    if len(found) != 2 or found[0] != f"{a.prefix}00000":
        sys.exit(f"FAIL: expected {a.prefix}00000 and one later plotfile, found {found}")
    _, first, real_size = read(Path(found[0]))
    head, last, _ = read(Path(found[1]))
    eps = 2.0 ** -52 if real_size == 8 else 2.0 ** -23
    g = a.gravity
    h1 = a.depth
    c1 = math.sqrt(g * h1)
    t = a.time
    nx = head["domain_hi"][0] + 1
    ny = head["domain_hi"][1] + 1
    x_end = head["prob_lo"][0] + nx * head["dx"][0]
    errors = []

    if not abs(head["time"] - t) <= 1e-12 * t:
        errors.append(f"time: got {head['time']!r}, expected {t!r}")
    if not a.position + 2.0 * c1 * t > x_end:
        sys.exit("check_partial_outflow.py: the front has not reached the side at --time")
    if any(not all(math.isfinite(x) for x in v) for v in last.values()):
        errors.append("a cell holds a value that is not finite")
    if any(v[0] < 0.0 for v in last.values()):
        errors.append("a depth is negative")
    m0 = math.fsum(v[0] for v in first.values())
    m1 = math.fsum(v[0] for v in last.values())
    if not m0 - m1 > 16.0 * eps * m0:
        errors.append(f"mass: the depths sum to {m1!r}, at the start to {m0!r}; water must leave")

    def y_of(j):
        return head["prob_lo"][1] + (j + 0.5) * head["dx"][1]

    def x_of(i):
        return head["prob_lo"][0] + (i + 0.5) * head["dx"][0]

    open_rows = [j for j in range(ny) if y_of(j) >= a.split + a.margin]
    wall_rows = [j for j in range(ny) if y_of(j) <= a.split - a.margin]
    if not open_rows or not wall_rows:
        sys.exit("check_partial_outflow.py: --split and --margin leave no open or no wall row")

    # Open rows: Ritter's solution.
    err_h = 0.0
    err_q = 0.0
    for j in open_rows:
        for i in range(nx):
            h_ex, u_ex = ritter((x_of(i) - a.position) / t, h1, g)
            v = last[(i, j)]
            err_h += abs(v[0] - h_ex)
            err_q += abs(v[1] - h_ex * u_ex)
    n_open = nx * len(open_rows)
    l1_h = err_h / (n_open * h1)
    l1_q = err_q / (n_open * h1 * c1)
    if not (l1_h <= a.l1_max and l1_q <= a.l1_max):
        errors.append(f"open rows: L1 error {l1_h:.4g} (depth) and {l1_q:.4g} (discharge), at most"
                      f" {a.l1_max}")

    # Wall rows: the bore at the wall.
    h_in, u_in = ritter((x_end - a.position) / t, h1, g)
    h2 = bore_depth(h_in, u_in, g)
    worst = max(abs(last[(nx - 1, j)][0] - h2) / h2 for j in wall_rows)
    if not worst <= a.bore_rtol:
        got = [last[(nx - 1, j)][0] for j in wall_rows]
        errors.append(f"wall rows: depth at the wall {min(got):.4g} to {max(got):.4g} m, bore depth"
                      f" {h2:.4g} m (arriving h {h_in:.4g} m, u {u_in:.4g} m/s): off by {worst:.3g},"
                      f" at most {a.bore_rtol}")

    if errors:
        print(f"FAIL: {found[1]}")
        for e in errors:
            print("  " + e)
        sys.exit(1)
    print(f"PASS: {found[1]}: mass change {(m1 - m0) / m0:.4g}; open rows ({len(open_rows)}) L1"
          f" {l1_h:.4g} (depth) and {l1_q:.4g} (discharge) of Ritter's solution; wall rows"
          f" ({len(wall_rows)}) depth at the wall within {worst:.3g} of the bore depth {h2:.4g} m")


if __name__ == "__main__":
    main()
