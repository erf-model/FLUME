#!/usr/bin/env python3
"""Check the first step of an inflow inside the domain (flume.inflows.*).

    check_inflow_source.py <prefix> --q Q (--line X0 Y0 X1 Y1 | --point X Y)
                           [--cfl C] [--gravity G] [--dry-depth D]
                           [--rtol R] [--deep H0] [--cells N] [--rain MMH]

Reads the last plotfile <prefix>NNNNN of a run of one step from a dry, flat,
closed domain fed by one inflow of a steady discharge Q (m3/s), and checks:

* The time of the plotfile, the length of the step, against the largest dt
  with dt (R0 + S sqrt(r dt)) <= cfl, S = sqrt(g) (1/dx + 1/dy),
  R0 = sqrt(g dry_depth) (1/dx + 1/dy) (no cell is wet) and r = f_max Q / A
  the fastest rise of a depth: the root y = sqrt(dt) of the cubic
  S sqrt(r) y^3 + R0 y^2 - cfl = 0, found here by bisection.
* The depth of every cell against f_c Q dt / A, with f_c the part of the
  line's length inside cell c (1 for the cell of a point, 0 elsewhere) and A
  the cell area. The parts come from clipping the line to each cell
  (Liang-Barsky), not from the solver's cutting at grid lines. The cells
  are half open: a line along a face is in the cell above it or to its
  right, and a cell the line only touches at a corner is not one of its
  cells. --cells N checks that the line has N cells.
* With --rain MMH, steady rain of MMH mm/h on every cell as well: the rate
  of the cubic is r = MMH / 3.6e6 + f_max Q / A (the sum of the fastest
  rises, Flood::compute_dt) and each cell holds the rain's depth too.
* No discharge anywhere: the water enters with no momentum, and the stage
  of the step saw a dry domain.

Both within R relative of the largest expected value: by default 1e-12 for a
plotfile in double precision, 1e-5 (84 epsilons) in single, where the state,
the cell size and the coordinates of the line are rounded to 7 digits.

With --deep H0 the run starts from still water H0 deep in a closed domain
where nothing flows (one cell): at the last plotfile, of any step, each cell
must hold H0 plus its part of Q t / A, within 4 epsilons of H0 (of the
precision of the plotfile), and no discharge. The step is not checked.
Exit status 0 when every check holds.
"""

import argparse
import math
import struct
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_plotfile as cp  # noqa: E402


def last_plotfile(prefix):
    runs = sorted(p for p in Path(".").glob(prefix + "[0-9]*") if p.is_dir())
    if not runs:
        sys.exit(f"FAIL: no plotfile {prefix}NNNNN in {Path('.').resolve()}")
    return runs[-1]


def clip_length(x0, y0, x1, y1, xlo, ylo, xhi, yhi, open_x=False, open_y=False):
    """The length of the segment from (x0, y0) to (x1, y1) inside the box.

    A segment along the high face of x (y) belongs to the next cell, not to
    this one, when open_x (open_y): the cells are half open, as a point on a
    face goes to the cell above it or to its right."""
    t0, t1 = 0.0, 1.0
    dx, dy = x1 - x0, y1 - y0
    for p, q, is_open in ((-dx, x0 - xlo, False), (dx, xhi - x0, open_x),
                          (-dy, y0 - ylo, False), (dy, yhi - y0, open_y)):
        if p == 0.0:
            if q < 0.0 or (q == 0.0 and is_open):
                return 0.0
            continue
        r = q / p
        if p < 0.0:
            t0 = max(t0, r)
        else:
            t1 = min(t1, r)
    return max(t1 - t0, 0.0) * math.hypot(dx, dy)


def to_single(v):
    """v rounded to the nearest single-precision number."""
    return struct.unpack("f", struct.pack("f", v))[0]


def first_step(cfl, g, dry_depth, dx, dy, r):
    """sqrt of the largest dt with dt (R0 + S sqrt(r dt)) <= cfl, by bisection."""
    s = math.sqrt(g) * (1.0 / dx + 1.0 / dy)
    r0 = math.sqrt(g * dry_depth) * (1.0 / dx + 1.0 / dy)
    f = lambda y: s * math.sqrt(r) * y ** 3 + r0 * y ** 2 - cfl  # noqa: E731
    lo, hi = 0.0, 1.0
    while f(hi) < 0.0:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return (0.5 * (lo + hi)) ** 2


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("prefix")
    p.add_argument("--q", type=float, required=True)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--line", type=float, nargs=4)
    g.add_argument("--point", type=float, nargs=2)
    p.add_argument("--cfl", type=float, default=0.5)
    p.add_argument("--gravity", type=float, default=9.81)
    p.add_argument("--dry-depth", type=float, default=1e-6)
    p.add_argument("--rtol", type=float,
                   help="default 1e-12 for a double-precision plotfile, 1e-5 for single")
    p.add_argument("--deep", type=float)
    p.add_argument("--cells", type=int)
    p.add_argument("--rain", type=float, default=0.0)
    a = p.parse_args()

    plt = last_plotfile(a.prefix)
    h = cp.read_header(plt / "Header")
    if h["vars"][:3] != cp.VARS[:3]:
        sys.exit(f"FAIL: {plt}: variables {h['vars']}")
    _, fab_on_disk = cp.read_cell_h(plt / "Level_0" / "Cell_H")
    cells, real_size = cp.read_cell_data(plt / "Level_0", fab_on_disk)
    dx, dy = h["dx"][0], h["dx"][1]
    xlo, ylo = h["prob_lo"][0], h["prob_lo"][1]
    area = dx * dy
    nx, ny = h["domain_hi"][0] + 1, h["domain_hi"][1] + 1
    if a.rtol is None:
        a.rtol = 1e-12 if real_size == 8 else 1e-5
    if a.deep is None and h["step"] != 1:
        sys.exit(f"FAIL: {plt} is step {h['step']}, expected 1")

    frac = {}
    if a.line:
        # The faces as the solver places them, from the corners of the domain
        # and the count of cells; the ends of the line as the build holds
        # them (rounded to single precision in a single-precision run).
        x0, y0, x1, y1 = a.line if real_size == 8 else [to_single(v) for v in a.line]
        xhi, yhi = h["prob_hi"][0], h["prob_hi"][1]
        fx = [xlo + (xhi - xlo) * i / nx for i in range(nx + 1)]
        fy = [ylo + (yhi - ylo) * j / ny for j in range(ny + 1)]
        total = math.hypot(x1 - x0, y1 - y0)
        for j in range(ny):
            for i in range(nx):
                ln = clip_length(x0, y0, x1, y1, fx[i], fy[j], fx[i + 1], fy[j + 1],
                                 i < nx - 1, j < ny - 1)
                # A cell the line only touches, at a corner, has a length of
                # rounding.
                if ln > 1e-12 * total:
                    frac[(i, j)] = ln / total
    else:
        i = min(int(math.floor((a.point[0] - xlo) / dx)), nx - 1)
        j = min(int(math.floor((a.point[1] - ylo) / dy)), ny - 1)
        frac[(i, j)] = 1.0
    fmax = max(frac.values())
    if a.cells is not None and len(frac) != a.cells:
        sys.exit(f"FAIL: the line crosses {len(frac)} cells, expected {a.cells}")

    if a.deep is not None:
        eps = 2.0 ** -52 if real_size == 8 else 2.0 ** -23
        tol = 4.0 * eps * a.deep
        t = h["time"]
        errors = []
        worst = 0.0
        for (i, j), v in sorted(cells.items()):
            want = a.deep + frac.get((i, j), 0.0) * a.q * t / area
            worst = max(worst, abs(v[0] - want))
            if abs(v[0] - want) > tol:
                errors.append(f"cell {(i, j)}: depth {v[0]!r}, expected {want!r} (H0 + "
                              f"{want - a.deep:.6g} m), off by {abs(v[0] - want):.3e} m")
            if v[1] != 0.0 or v[2] != 0.0:
                errors.append(f"cell {(i, j)}: discharge ({v[1]!r}, {v[2]!r}), expected 0")
        if errors:
            for e in errors[:20]:
                print("FAIL:", e)
            sys.exit(1)
        print(f"PASS: {plt}: step {h['step']}, {t:.9g} s: the depths within {worst:.3e} m of"
              f" H0 plus the inflow (tolerance {tol:.3e} m)")
        return

    rain = a.rain / 3.6e6  # m/s
    dt = first_step(a.cfl, a.gravity, a.dry_depth, dx, dy, rain + fmax * a.q / area)
    errors = []
    if abs(h["time"] - dt) > a.rtol * dt:
        errors.append(f"the step is {h['time']!r} s, expected {dt!r} s "
                      f"(relative {abs(h['time'] - dt) / dt:.3e})")
    scale = (rain + fmax * a.q / area) * dt
    worst = 0.0
    for (i, j), v in sorted(cells.items()):
        want = rain * dt + frac.get((i, j), 0.0) * a.q * dt / area
        err = abs(v[0] - want)
        worst = max(worst, err / scale)
        if err > a.rtol * scale:
            errors.append(f"cell {(i, j)}: depth {v[0]!r}, expected {want!r}")
        if v[1] != 0.0 or v[2] != 0.0:
            errors.append(f"cell {(i, j)}: discharge ({v[1]!r}, {v[2]!r}), expected 0")
    if errors:
        for e in errors[:20]:
            print("FAIL:", e)
        sys.exit(1)
    print(f"PASS: {plt}: step {h['time']:.12g} s (closed form {dt:.12g} s), {len(frac)} cells"
          f" of the inflow, largest part {fmax:.6g}; depths within {worst:.3g} of the largest")


if __name__ == "__main__":
    main()
