#!/usr/bin/env python3
"""Check the water-balance log of a run (water_balance.txt).

    check_balance.py [<log>] [--factor F] [--expect-columns NAME...]
                     [--expect-final NAME=VALUE...] [--rtol R]
                     [--expect-last-step N] [--expect-lines N]

Between two lines a and b of the log the residual may change by at most

    F * (eps * n * (max(V_a, V_b) + M) + eps_d * C),

    M = sum over the other columns of |change from a to b|,
    C = |V_0| + sum over the other columns of max(|value at a|, |value at b|),

with eps the machine epsilon of the state (from the header), eps_d that of
double precision, n the number of steps between the two lines and V the
stored volume. The first term lets each step keep the water it accounts for
to F roundings of the water it holds and moves; the second is the rounding
of the residual itself, which is formed from the cumulative totals in double
precision, also in a single-precision build. A term left out of the balance fails as soon as
it moves more than that. run_case.py calls check() after every regression run
that leaves a log. It also recomputes the residual of every line from the
columns, with the signs of the '# added:', '# removed:', '# out:' and (when
there is an inflow segment) '# in:' header lines, and requires each column to
be listed once.

--expect-columns gives the names of the columns from 'stored' to the last
before 'residual'. --expect-final pins the last value of a column against a closed form, to the
relative tolerance --rtol (default 1e-12) of the larger of the value and the
largest stored volume, plus F * eps of the same. --expect-last-step and
--expect-lines give the step of the last line and the number of lines.
"""
import argparse
import math
import sys
from pathlib import Path

#: Roundings of the moved water a step may leave in the residual. The largest
#: ratio over the regression tests is printed by each check; see
#: Docs/sphinx_doc/Verification.rst for the values measured.
DEFAULT_FACTOR = 64.0

#: Machine epsilon of double precision, in which the totals are formed.
EPS_D = 2.0 ** -52


def read_log(path):
    """Return (epsilon, column names, rows of floats)."""
    eps = None
    names = None
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("# epsilon "):
            eps = float(line.split()[2])
        elif line.startswith("# step "):
            names = line[2:].split()
        elif line.strip() and not line.startswith("#"):
            rows.append([float(v) for v in line.split()])
    if eps is None or names is None:
        raise ValueError(f"{path}: no '# epsilon' or '# step' header line")
    for r in rows:
        if len(r) != len(names):
            raise ValueError(f"{path}: a line has {len(r)} values for {len(names)} columns")
    if not rows or rows[0][0] != 0:
        raise ValueError(f"{path}: the first line is not step 0")
    return eps, names, rows


def read_terms(path):
    """The columns of water added (with those of the optional '# in:' line, the
    inflow through discharge segments), and of water removed or out, from the
    header."""
    added = removed = out = None
    inflow = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("# added:"):
            added = line.split()[2:]
        elif line.startswith("# removed:"):
            removed = line.split()[2:]
        elif line.startswith("# out:"):
            out = line.split()[2:]
        elif line.startswith("# in:"):
            inflow = line.split()[2:]
    if added is None or removed is None or out is None:
        raise ValueError(f"{path}: no '# added:', '# removed:' or '# out:' header line")
    return added + inflow, removed + out


def recompute(path, eps, names, rows, factor):
    """Check the logged residual of every line against the columns: stored -
    stored at step 0 - added + removed + out - in, within factor * eps_d of the sum of
    the magnitudes. Return an error message, or None."""
    added, removed = read_terms(path)
    listed = added + removed
    columns = names[3:-1]
    if sorted(listed) != sorted(columns):
        return (f"the header lists {listed} as added, removed, out or in, but the columns are "
                f"{columns}")
    i_v = names.index("stored")
    i_r = names.index("residual")
    i_a = [names.index(c) for c in added]
    i_s = [names.index(c) for c in removed]
    v0 = rows[0][i_v]
    for r in rows:
        want = r[i_v] - v0 - sum(r[i] for i in i_a) + sum(r[i] for i in i_s)
        mag = abs(r[i_v]) + abs(v0) + sum(abs(r[i]) for i in i_a + i_s)
        if not abs(r[i_r] - want) <= factor * EPS_D * mag:
            return (f"at step {r[0]:.0f} the logged residual is {r[i_r]!r}, but the columns"
                    f" give {want!r}")
    return None


def scaled_changes(eps, names, rows):
    """For each pair of lines: (line a, line b, |change of the residual|, scale),
    with the bound factor * eps * scale (module docstring): the totals term C
    enters the scale as C * eps_d / eps."""
    i_v = names.index("stored")
    i_r = names.index("residual")
    others = [i for i in range(len(names)) if i not in (0, 1, i_v, i_r)]
    for a, b in zip(rows, rows[1:]):
        n = b[0] - a[0]
        moved = max(abs(a[i_v]), abs(b[i_v])) + sum(abs(b[i] - a[i]) for i in others)
        held = abs(rows[0][i_v]) + sum(max(abs(a[i]), abs(b[i])) for i in others)
        yield a, b, abs(b[i_r] - a[i_r]), n * moved + held * (EPS_D / eps)


def largest_ratio(path):
    """The largest change of the residual between two lines over eps * scale."""
    eps, names, rows = read_log(path)
    worst = 0.0
    for _, _, change, scale in scaled_changes(eps, names, rows):
        if change > 0.0:
            worst = max(worst, math.inf if scale == 0.0 else change / (eps * scale))
    return worst


def check(path, factor=DEFAULT_FACTOR):
    """Return (ok, message, largest ratio of a residual change to its bound)."""
    eps, names, rows = read_log(path)
    i_r = names.index("residual")
    if rows[0][i_r] != 0.0:
        return False, f"the residual of step 0 is {rows[0][i_r]!r}, not 0", math.inf
    err = recompute(path, eps, names, rows, factor)
    if err is not None:
        return False, err, math.inf
    worst = 0.0
    for a, b, change, scale in scaled_changes(eps, names, rows):
        if b[0] <= a[0]:
            return False, f"step {b[0]:.0f} follows step {a[0]:.0f}", math.inf
        bound = factor * eps * scale
        if not math.isfinite(change):
            return False, f"the residual is not finite at step {b[0]:.0f}", math.inf
        if change > bound:
            return (False,
                    f"the residual changes by {change:.3e} m3 from step {a[0]:.0f} to step "
                    f"{b[0]:.0f}, above {factor:g} (eps n (V + moved) + eps_d totals) = {bound:.3e} m3",
                    math.inf if bound == 0.0 else change / (eps * scale))
        if scale > 0.0:
            worst = max(worst, change / (eps * scale))
    last = rows[-1]
    return (True,
            f"{len(rows)} lines to step {last[0]:.0f}: residual {last[i_r]:.3e} m3, largest "
            f"change {worst:.2f} (eps n (V + moved) + eps_d totals), bound {factor:g}",
            worst)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("log", nargs="?", default="water_balance.txt")
    p.add_argument("--factor", type=float, default=DEFAULT_FACTOR)
    p.add_argument("--expect-columns", nargs="+", default=None,
                   help="the column names from 'stored' to the last before 'residual'")
    p.add_argument("--expect-final", nargs="+", default=[], metavar="NAME=VALUE")
    p.add_argument("--rtol", type=float, default=1e-12)
    p.add_argument("--expect-last-step", type=int, default=None)
    p.add_argument("--expect-lines", type=int, default=None)
    a = p.parse_args()

    ok, msg, _ = check(a.log, a.factor)
    print(f"balance: {msg}")
    if not ok:
        sys.exit(f"FAIL: {msg}")
    eps, names, rows = read_log(a.log)
    if a.expect_columns is not None:
        got = names[2:-1]
        if got != a.expect_columns:
            sys.exit(f"FAIL: the columns are {got}, expected {a.expect_columns}")
    if a.expect_last_step is not None and rows[-1][0] != a.expect_last_step:
        sys.exit(f"FAIL: the last line is step {rows[-1][0]:.0f}, expected {a.expect_last_step}")
    if a.expect_lines is not None and len(rows) != a.expect_lines:
        sys.exit(f"FAIL: the log has {len(rows)} lines, expected {a.expect_lines}")
    v_max = max(abs(r[names.index("stored")]) for r in rows)
    for item in a.expect_final:
        name, value = item.split("=", 1)
        if name not in names:
            sys.exit(f"FAIL: no column {name} in {names}")
        want = float(value)
        got = rows[-1][names.index(name)]
        scale = max(abs(want), v_max)
        tol = a.rtol * scale + a.factor * eps * scale
        print(f"balance: {name} = {got!r}, expected {want!r} (difference {got - want:.3e},"
              f" tolerance {tol:.3e})")
        if not abs(got - want) <= tol:
            sys.exit(f"FAIL: {name} = {got!r}, expected {want!r} within {tol:.3e}")
    print("PASS")


if __name__ == "__main__":
    main()
