#!/usr/bin/env python3
"""Check FLUME runs with prescribed rain against closed-form numbers.

    check_rain.py <prefix> (--rate R | --hyetograph FILE [--linear])
        (--uniform | --added-volume | --kinematic-wave N S0 L)
        [--tol-eps K] [--depth-l1-max E] [--q-out-max E] [--equilibrium-rtol E]
        [--n-plotfiles N] [--initial-depth H0] [--skip-initial]

The rain is R mm/h from t = 0 on, or the table of FILE (lines of a time in s
and a rate in mm/h; `#` starts a comment), held from one entry to the next
or, with --linear, linear between entries, and held at the end values before
the first and after the last entry. Its integral from 0 to t, I(t), is
computed here from the table, independently of the solver. Every plotfile
<prefix><step> in the directory is read (at least N with --n-plotfiles);
with --skip-initial the plotfile of step 0 is not checked (its initial
state may hold what the first step changes).

--uniform        A closed domain with a flat bed and no water at t = 0. In
                 every plotfile every cell holds the depth I(t) within K
                 epsilons (of the precision the data are written in) times
                 I(t) times the number of the step, and both discharges are
                 exactly 0. With --initial-depth the domain starts with
                 the depth H0 (m) instead: every cell holds H0 + I(t)
                 within K epsilons times H0 + I(t), not times the number of
                 the step: the solver adds the rain with a compensated sum,
                 so the error does not grow with the steps.
--added-volume   A closed domain: the sum of the depth over the cells of
                 the last plotfile minus that of the first is the number of
                 cells times I(T) - I(0), within K epsilons times the larger
                 sum times the number of the step.
--kinematic-wave Rain on a plane bed falling along x at the slope S0 with
                 Manning N, a wall at x = 0 and the outflow at x = L,
                 dry at t = 0 (Docs/sphinx_doc/Testing.rst). With
                 a = S0^(1/2) / N and i the rain rate (m/s), the depth at x
                 is i t for t < t_x = (x / (a i^(2/3)))^(3/5) and
                 (i x / a)^(3/5) after, and the discharge is a h^(5/3).
                 Checked in every plotfile after the first: the L1 error of
                 the depth over the cells relative to the L1 norm of the
                 exact depth is at most E (--depth-l1-max), and the mean
                 discharge of the last column of cells differs from the exact
                 discharge at x = L by at most E times i L (--q-out-max).
                 With --equilibrium-rtol the last plotfile must hold that
                 discharge within E times i L: at equilibrium the outlet
                 passes all the rain on the plane. Constant rain only.
"""
import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import check_dam_break as cd  # noqa: E402

MM_PER_H = 1.0 / 3.6e6


def read_table(path):
    time, rate = [], []
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].split()
        if not line:
            continue
        if len(line) != 2:
            sys.exit(f"FAIL: {path}: a line must hold a time and a rate: {line}")
        time.append(float(line[0]))
        rate.append(float(line[1]) * MM_PER_H)
    if not time or any(b <= a for a, b in zip(time, time[1:])):
        sys.exit(f"FAIL: {path}: no entry, or times not increasing")
    return time, rate


def integral(time, rate, linear, t):
    """The integral of the rain rate from time[0] to t, m."""
    if t <= time[0]:
        return rate[0] * (t - time[0])
    total = 0.0
    for k in range(len(time) - 1):
        a, b = time[k], time[k + 1]
        if t <= a:
            break
        e = min(t, b)
        if linear:
            r_e = rate[k] + (rate[k + 1] - rate[k]) * (e - a) / (b - a)
            total += 0.5 * (rate[k] + r_e) * (e - a)
        else:
            total += rate[k] * (e - a)
    if t > time[-1]:
        total += rate[-1] * (t - time[-1])
    return total


def main():
    p = argparse.ArgumentParser()
    p.add_argument("prefix")
    r = p.add_mutually_exclusive_group(required=True)
    r.add_argument("--rate", type=float)
    r.add_argument("--hyetograph")
    p.add_argument("--linear", action="store_true")
    m = p.add_mutually_exclusive_group(required=True)
    m.add_argument("--uniform", action="store_true")
    m.add_argument("--added-volume", action="store_true")
    m.add_argument("--kinematic-wave", type=float, nargs=3, metavar=("N", "S0", "L"))
    p.add_argument("--tol-eps", type=float, default=4.0)
    p.add_argument("--depth-l1-max", type=float)
    p.add_argument("--q-out-max", type=float)
    p.add_argument("--equilibrium-rtol", type=float)
    p.add_argument("--n-plotfiles", type=int, default=2)
    p.add_argument("--initial-depth", type=float)
    p.add_argument("--skip-initial", action="store_true")
    a = p.parse_args()

    if a.rate is not None:
        time, rate, linear = [0.0], [a.rate * MM_PER_H], False
    else:
        time, rate = read_table(a.hyetograph)
        linear = a.linear
    # The rain from t = 0: before the first entry the rate is the first one.
    def rain(t):
        return integral(time, rate, linear, t) - integral(time, rate, linear, 0.0)

    found = sorted(d for d in Path(".").iterdir()
                   if d.is_dir() and re.fullmatch(re.escape(a.prefix) + r"\d{5,}", d.name))
    if len(found) < a.n_plotfiles:
        sys.exit(f"FAIL: {len(found)} plotfile(s) {a.prefix}NNNNN, expected at least {a.n_plotfiles}")
    errors = []

    if a.uniform:
        worst = 0.0
        for f in found:
            head, cells, size = cd.read(f)
            eps = 2.0 ** -52 if size == 8 else 2.0 ** -23
            t, step = head["time"], head["step"]
            if a.skip_initial and step == 0:
                continue
            if a.initial_depth is None:
                want = rain(t)
                tol = a.tol_eps * eps * want * max(step, 1)
            else:
                want = a.initial_depth + rain(t)
                tol = a.tol_eps * eps * want
            for (i, j), v in cells.items():
                worst = max(worst, abs(v[0] - want))
                if abs(v[0] - want) > tol:
                    errors.append(f"{f.name} cell ({i}, {j}): depth {v[0]!r}, expected {want!r} (tol {tol:.3g})")
                if v[1] != 0.0 or v[2] != 0.0:
                    errors.append(f"{f.name} cell ({i}, {j}): discharge ({v[1]!r}, {v[2]!r}), expected 0")
            print(f"{f.name}: t = {t} s, step {step}, rain {want * 1e3:.6f} mm")
        print(f"largest depth difference {worst:.3g} m")

    elif a.added_volume:
        h0, c0, size = cd.read(found[0])
        h1, c1, _ = cd.read(found[-1])
        eps = 2.0 ** -52 if size == 8 else 2.0 ** -23
        v0 = sum(v[0] for v in c0.values())
        v1 = sum(v[0] for v in c1.values())
        want = len(c1) * (rain(h1["time"]) - rain(h0["time"]))
        tol = a.tol_eps * eps * max(v0, v1) * max(h1["step"], 1)
        print(f"sum of depth {v0!r} -> {v1!r}: added {v1 - v0!r}, expected {want!r} (tol {tol:.3g})")
        if abs(v1 - v0 - want) > tol:
            errors.append(f"added depth sum {v1 - v0!r}, expected {want!r}")

    else:
        if a.rate is None:
            sys.exit("FAIL: --kinematic-wave needs --rate")
        if a.depth_l1_max is None or a.q_out_max is None:
            sys.exit("FAIL: --kinematic-wave needs --depth-l1-max and --q-out-max")
        n, s0, length = a.kinematic_wave
        i_rate = rate[0]
        alpha = math.sqrt(s0) / n

        def depth(x, t):
            tx = (x / (alpha * i_rate ** (2.0 / 3.0))) ** 0.6
            return i_rate * t if t < tx else (i_rate * x / alpha) ** 0.6

        q_eq = i_rate * length
        l1_worst = q_worst = 0.0
        q_last = None
        for f in found:
            head, cells, _ = cd.read(f)
            t = head["time"]
            if t == 0.0:
                continue
            lo, dx = head["prob_lo"], head["dx"]
            nx = max(k[0] for k in cells) + 1
            ny = max(k[1] for k in cells) + 1
            num = den = 0.0
            for (i, j), v in cells.items():
                e = depth(lo[0] + (i + 0.5) * dx[0], t)
                num += abs(v[0] - e)
                den += e
            l1 = num / den
            q_out = sum(cells[(nx - 1, j)][1] for j in range(ny)) / ny
            q_exact = alpha * depth(length, t) ** (5.0 / 3.0)
            dq = abs(q_out - q_exact) / q_eq
            l1_worst, q_worst = max(l1_worst, l1), max(q_worst, dq)
            q_last = q_out
            print(f"{f.name}: t = {t:.6g} s, depth L1 {l1:.4f}, q_out {q_out:.6e}"
                  f" (exact {q_exact:.6e}, error {dq:.4f} of i L)")
            if l1 > a.depth_l1_max:
                errors.append(f"{f.name}: depth L1 error {l1:.4f} above {a.depth_l1_max}")
            if dq > a.q_out_max:
                errors.append(f"{f.name}: outlet discharge error {dq:.4f} of i L above {a.q_out_max}")
        print(f"largest depth L1 error {l1_worst:.4f}, largest outlet error {q_worst:.4f} of i L")
        if a.equilibrium_rtol is not None:
            rel = abs(q_last - q_eq) / q_eq
            print(f"last outlet discharge {q_last!r}, i L = {q_eq!r}, relative difference {rel:.3g}")
            if rel > a.equilibrium_rtol:
                errors.append(f"last outlet discharge {q_last!r} differs from i L = {q_eq!r} by {rel:.3g}")

    if errors:
        for e in errors[:20]:
            print("FAIL:", e)
        if len(errors) > 20:
            print(f"FAIL: ... and {len(errors) - 20} more")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
