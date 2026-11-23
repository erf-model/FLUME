#!/usr/bin/env python3
"""Check FLUME runs with gridded rain or a moving rain cell.

    check_rain_field.py <prefix> (--index FILE [--linear] | --cell P SIGMA X0 Y0 U V)
        (--cells STEP | --volume) [--tol-eps K] [--rtol R] [--n-plotfiles N]

The rain is the grids of the index FILE (lines of a time in s and the path of
an ESRI ASCII grid of the rate in mm/h, relative to FILE; `#` starts a
comment), each constant on its own cells, held from one entry to the next or,
with --linear, linear in time between entries, and held at the first and the
last grid before the first and after the last time. Or it is a cell of peak
rate P (mm/h) whose rate falls off as a Gaussian of standard deviation SIGMA
(m) from a centre at (X0 + U t, Y0 + V t) (m, m/s). Everything is computed
here from the files or the numbers, independently of the solver: the grids
by the areas each grid cell shares with a solver cell (or with the domain),
the time by splitting [0, t] at the times of the index (exact for a rate
held or linear between them), the cell by error functions in space and
adaptive Simpson quadrature in time.

--cells STEP  The plotfile of step STEP of a run on a dry domain whose
              fluxes are not computed until then (so STEP is 1): every cell
              holds the rain that fell on it from 0 to the plotfile's time,
              within K epsilons (of the precision the data are written in)
              times the largest depth, or R times it.
--volume      A closed domain with a flat bed and no water at t = 0: in every
              plotfile, the water in the domain (the sum of the depths times
              the cell area) is the integral of the rain over the domain from
              0 to the plotfile's time, within K epsilons times that integral
              times the number of the step, or R times the integral.
"""
import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_dam_break as cd  # noqa: E402

MM_PER_H = 1.0 / 3.6e6
# The overlap, as a fraction of a grid cell, below which a NODATA cell counts
# as untouched: round-off in the solver's edges. 16 epsilons of the precision
# the plotfile is written in; main() sets it from the data.
SLIVER = [16.0 * 2.0 ** -52]


def read_asc(path):
    """(x0, y0, dx, dy, values[j][i] in m/s with the southern row first, nodata[j][i])."""
    tok = Path(path).read_text().split()
    head = {}
    k = 0
    while tok[k][0].isalpha():
        head[tok[k].lower()] = float(tok[k + 1])
        k += 2
    nx, ny = int(head["ncols"]), int(head["nrows"])
    dx = head.get("cellsize", head.get("dx"))
    dy = head.get("cellsize", head.get("dy"))
    x0 = head["xllcorner"] if "xllcorner" in head else head["xllcenter"] - 0.5 * dx
    y0 = head["yllcorner"] if "yllcorner" in head else head["yllcenter"] - 0.5 * dy
    vals = [float(t) for t in tok[k:]]
    if len(vals) != nx * ny:
        sys.exit(f"FAIL: {path}: {len(vals)} values, expected {nx * ny}")
    nod = head.get("nodata_value")
    v = [[0.0] * nx for _ in range(ny)]
    bad = [[False] * nx for _ in range(ny)]
    for m, val in enumerate(vals):
        i, j = m % nx, ny - 1 - m // nx
        if nod is not None and val == nod:
            bad[j][i] = True
        else:
            v[j][i] = val * MM_PER_H
    return x0, y0, dx, dy, v, bad


def overlap(a, b, lo, hi):
    return max(0.0, min(b, hi) - max(a, lo))


def grid_integral(g, xa, xb, ya, yb):
    """The integral of the grid over the rectangle [xa, xb] x [ya, yb], m^3/s;
    fails on a NODATA cell that shares more than round-off with it."""
    x0, y0, dx, dy, v, bad = g
    ny, nx = len(v), len(v[0])
    total = 0.0
    for j in range(ny):
        ly = overlap(ya, yb, y0 + j * dy, y0 + (j + 1) * dy)
        if ly == 0.0:
            continue
        for i in range(nx):
            lx = overlap(xa, xb, x0 + i * dx, x0 + (i + 1) * dx)
            if lx == 0.0:
                continue
            if bad[j][i]:
                # A sliver of round-off in the rectangle's edges is no overlap.
                if lx <= SLIVER[0] * dx or ly <= SLIVER[0] * dy:
                    continue
                sys.exit(f"FAIL: NODATA cell ({i}, {j}) inside [{xa}, {xb}] x [{ya}, {yb}]")
            total += lx * ly * v[j][i]
    return total


def read_index(path):
    time, grids = [], []
    for line in Path(path).read_text().splitlines():
        t = line.split("#", 1)[0].split()
        if not t:
            continue
        time.append(float(t[0]))
        grids.append(read_asc(Path(path).parent / t[1]))
    return time, grids


def time_pieces(time, linear, t_end):
    """[(a, b, {k: weight at a}, {k: weight at b})]: the pieces of [0, t_end]
    between the index times, with the weight of each grid in the rate at the
    two ends of the piece. The integral over a piece is (b - a) times the
    mean of the two ends, exact for a rate held or linear on the piece."""
    def weights(t, side):
        n = len(time)
        if t < time[0] or (t == time[0] and side < 0):
            return {0: 1.0}
        if t > time[-1] or (t == time[-1] and side > 0):
            return {n - 1: 1.0}
        for k in range(n - 1):
            if time[k] <= t <= time[k + 1]:
                if not linear:
                    # held: the grid that starts at or before t, from the piece's side
                    if t == time[k + 1] and side < 0:
                        return {k: 1.0}
                    return {k + 1: 1.0} if t == time[k + 1] else {k: 1.0}
                u = (t - time[k]) / (time[k + 1] - time[k])
                return {k: 1.0 - u, k + 1: u}
        return {n - 1: 1.0}
    cuts = sorted({0.0, t_end} | {t for t in time if 0.0 < t < t_end})
    return [(a, b, weights(a, +1), weights(b, -1)) for a, b in zip(cuts, cuts[1:])]


def gridded_depths(time, grids, linear, t_end, rects):
    """For each rectangle, the rain volume from 0 to t_end, m^3."""
    pieces = time_pieces(time, linear, t_end)
    used = {k for _, _, wa, wb in pieces for k in list(wa) + list(wb)}
    out = []
    for r in rects:
        q = {k: grid_integral(grids[k], *r) for k in used}
        total = 0.0
        for a, b, wa, wb in pieces:
            ra = sum(w * q[k] for k, w in wa.items())
            rb = sum(w * q[k] for k, w in wb.items())
            total += (b - a) * 0.5 * (ra + rb)
        out.append(total)
    return out


def simpson(f, a, b, tol):
    """The integral of f over [a, b] by adaptive Simpson, to about tol."""
    def s(a, fa, m, fm, b, fb):
        return (b - a) / 6.0 * (fa + 4.0 * fm + fb)

    def rec(a, fa, b, fb, m, fm, whole, tol, depth):
        lm, rm = 0.5 * (a + m), 0.5 * (m + b)
        flm, frm = f(lm), f(rm)
        left, right = s(a, fa, lm, flm, m, fm), s(m, fm, rm, frm, b, fb)
        if depth <= 0 or abs(left + right - whole) <= 15.0 * tol:
            return left + right + (left + right - whole) / 15.0
        return (rec(a, fa, m, fm, lm, flm, left, 0.5 * tol, depth - 1)
                + rec(m, fm, b, fb, rm, frm, right, 0.5 * tol, depth - 1))

    fa, fb, m = f(a), f(b), 0.5 * (a + b)
    fm = f(m)
    return rec(a, fa, b, fb, m, fm, s(a, fa, m, fm, b, fb), tol, 40)


def cell_depths(cell, t_end, rects):
    """For each rectangle, the volume the moving cell drops from 0 to t_end, m^3."""
    p, sigma, x0, y0, u, v = cell
    s2 = math.sqrt(2.0) * sigma
    c = sigma * math.sqrt(0.5 * math.pi)

    def f(a, b, centre):
        return c * (math.erf((b - centre) / s2) - math.erf((a - centre) / s2))

    out = []
    for xa, xb, ya, yb in rects:
        g = lambda t: p * f(xa, xb, x0 + u * t) * f(ya, yb, y0 + v * t)  # noqa: E731
        # to 1e-14 of the volume the cell's peak would drop on the rectangle
        out.append(simpson(g, 0.0, t_end, 1e-14 * p * t_end * (xb - xa) * (yb - ya)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--index")
    src.add_argument("--cell", type=float, nargs=6)
    ap.add_argument("--linear", action="store_true")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--cells", type=int)
    mode.add_argument("--volume", action="store_true")
    ap.add_argument("--tol-eps", type=float, default=16.0)
    ap.add_argument("--rtol", type=float, default=0.0)
    ap.add_argument("--n-plotfiles", type=int, default=2)
    a = ap.parse_args()
    if a.linear and not a.index:
        sys.exit("FAIL: --linear goes with --index")

    found = sorted(d for d in Path(".").iterdir()
                   if d.is_dir() and re.fullmatch(re.escape(a.prefix) + r"\d{5,}", d.name))
    if len(found) < a.n_plotfiles:
        sys.exit(f"FAIL: {len(found)} plotfile(s) {a.prefix}NNNNN, expected at least {a.n_plotfiles}")
    if a.index:
        time, grids = read_index(a.index)
        depths = lambda t, rects: gridded_depths(time, grids, a.linear, t, rects)  # noqa: E731
    else:
        cell = list(a.cell)
        cell[0] *= MM_PER_H
        depths = lambda t, rects: cell_depths(cell, t, rects)  # noqa: E731

    errors = []
    if a.cells is not None:
        sel = [f for f in found if cd.read(f)[0]["step"] == a.cells]
        if len(sel) != 1:
            sys.exit(f"FAIL: no plotfile of step {a.cells}")
        head, cells, size = cd.read(sel[0])
        eps = 2.0 ** -52 if size == 8 else 2.0 ** -23
        SLIVER[0] = 16.0 * eps
        lo, dx = head["prob_lo"], head["dx"]
        keys = sorted(cells)
        rects = [(lo[0] + i * dx[0], lo[0] + (i + 1) * dx[0], lo[1] + j * dx[1], lo[1] + (j + 1) * dx[1])
                 for i, j in keys]
        area = dx[0] * dx[1]
        want = [d / area for d in depths(head["time"], rects)]
        top = max(want)
        tol = max(a.tol_eps * eps, a.rtol) * top
        worst = 0.0
        for (i, j), w in zip(keys, want):
            got = cells[(i, j)][0]
            worst = max(worst, abs(got - w))
            if abs(got - w) > tol:
                errors.append(f"cell ({i}, {j}): depth {got!r}, expected {w!r} (tol {tol:.3g})")
        print(f"{sel[0].name}: t = {head['time']} s, depths {min(want):.6e} to {top:.6e} m,"
              f" largest difference {worst:.3g} m ({worst / top:.3g} of the largest), tol {tol:.3g}")
    else:
        for f in found:
            head, cells, size = cd.read(f)
            eps = 2.0 ** -52 if size == 8 else 2.0 ** -23
            SLIVER[0] = 16.0 * eps
            lo, dx = head["prob_lo"], head["dx"]
            nx = max(k[0] for k in cells) + 1
            ny = max(k[1] for k in cells) + 1
            got = sum(v[0] for v in cells.values()) * dx[0] * dx[1]
            want = depths(head["time"], [(lo[0], lo[0] + nx * dx[0], lo[1], lo[1] + ny * dx[1])])[0]
            tol = max(a.tol_eps * eps * max(head["step"], 1), a.rtol) * abs(want)
            print(f"{f.name}: t = {head['time']:.6f} s, step {head['step']}, water {got!r} m^3,"
                  f" rain {want!r} m^3, difference {abs(got - want) / max(want, 1e-300):.3g} of it")
            if abs(got - want) > tol:
                errors.append(f"{f.name}: water {got!r} m^3, rain {want!r} m^3 (tol {tol:.3g})")
    if errors:
        print("\n".join(errors[:20]))
        print(f"FAIL: {len(errors)} error(s)")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
