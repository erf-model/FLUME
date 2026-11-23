#!/usr/bin/env python3
"""Check the flood maps and gauge files of a FLUME run against closed forms.

    check_flood_outputs.py <plotfile> --rain-basin RATE ARRIVAL_DEPTH
        [--bed Z] --gauge NAME I J [--gauge ...]
    check_flood_outputs.py <plotfile> --dam-break H1 H0 X0 G ARRIVAL_DEPTH
        [--gauge NAME I J ...]

<plotfile> is the prefix of the run's plotfiles (flume.plot_file); the last
of them, with the highest step, is read. It holds the maps max_depth,
max_speed, max_depth_speed and arrival_time after the state variables.

--rain-basin: rain of RATE mm/h on a closed, flat, dry domain (bed at Z m,
default 0).
Every cell holds the rain fallen by the time of the plotfile, T, in
max_depth, 0 in both speed maps, and in arrival_time the first time at
which the depth was above ARRIVAL_DEPTH m. The run must log every step:
that time is read from the gauge files, as the first line whose depth is
above ARRIVAL_DEPTH, and must lie within one step after ARRIVAL_DEPTH / rate.
Each --gauge NAME I J must have read cell (I, J) (the file's first line), and
on each of its lines the depth is the rain fallen by the line's time, the
surface elevation Z plus that, and both velocities are 0.

--dam-break: a dam break over a wet bed, still water H1 m deep at x < X0 m
and H0 m beyond, gravity G m/s2 (Tests/DamBreak/inputs), with ARRIVAL_DEPTH
between H0 and the depth behind the shock, hm. In each cell of the first row
that the shock (speed S) has passed by T, with its centre at least 8 cells
from X0 and from the shock at T, arrival_time is (x - X0) / S within one cell
divided by S plus the largest step (the run must log every step); max_depth,
max_speed and max_depth_speed are hm, um and hm um, the state behind the
shock, within 2 %. Each --gauge NAME I J (in a cell the shock has passed,
on a bed at 0) must have read cell (I, J) and end with that state: depth hm
and velocity (um, 0) within 2 %, surface elevation equal to the depth. In each
cell at least 8 cells ahead of the shock at T water has not arrived:
arrival_time is -1 and max_depth is below ARRIVAL_DEPTH (the numerical
shock runs a little ahead of the exact one, so the depth there is not
exactly H0).

Values are compared within 64 machine epsilons of the data's precision,
relative to the larger of the value and 1e-3 m (depths) unless stated.
"""
import argparse
import math
import struct
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import check_plotfile as cp  # noqa: E402
import check_dam_break as cd  # noqa: E402
import check_balance as cb  # noqa: E402

MAPS = ["max_depth", "max_speed", "max_depth_speed", "arrival_time"]
MMH = 1e-3 / 3600.0


def as_state(t, log="water_balance.txt"):
    """The time `t` as the plotfiles hold it: rounded to single precision when
    the state is single (machine epsilon of the log above 1e-10)."""
    eps, _, _ = cb.read_log(log)
    return struct.unpack("f", struct.pack("f", t))[0] if eps > 1e-10 else t


def last_plotfile(prefix):
    found = sorted(p for p in Path(".").glob(prefix + "*")
                   if p.is_dir() and p.name[len(prefix):].isdigit())
    if not found:
        sys.exit(f"FAIL: no plotfile {prefix}<step>")
    return max(found, key=lambda p: int(p.name[len(prefix):]))


def read(plotfile):
    """Header, {(i, j): {variable: value}}, machine epsilon of the data."""
    h = cp.read_header(plotfile / "Header")
    if h["vars"][:len(cp.VARS)] != cp.VARS or h["vars"][-len(MAPS):] != MAPS:
        sys.exit(f"FAIL: {plotfile}: variables {h['vars']}; expected the state first and"
                 f" {MAPS} last")
    _, fab_on_disk = cp.read_cell_h(plotfile / "Level_0" / "Cell_H")
    raw, real_size = cp.read_cell_data(plotfile / "Level_0", fab_on_disk)
    cells = {c: dict(zip(h["vars"], v)) for c, v in raw.items()}
    eps = 2.0 ** -52 if real_size == 8 else 2.0 ** -23
    return h, cells, eps


def read_gauge(path):
    """(cell (i, j), rows of step, time, depth, surface, u, v)."""
    lines = Path(path).read_text().splitlines()
    head = lines[0]
    k = head.rindex("cell (")
    i, j = (int(x) for x in head[k + 6:].rstrip(")").split(","))
    rows = [[float(x) for x in s.split()] for s in lines if s.strip() and not s.startswith("#")]
    return (i, j), rows


def rain_basin(a, h, cells, eps, errors):
    rate = a.rain_basin[0] * MMH
    arrival_depth = a.rain_basin[1]
    t_end = h["time"]
    depth = rate * t_end
    tol = 64.0 * eps
    first_wet = None
    for name, i, j in a.gauge:
        cell, rows = read_gauge(f"gauge_{name}.txt")
        if cell != (int(i), int(j)):
            errors.append(f"gauge {name} read cell {cell}, expected ({int(i)}, {int(j)})")
        if not rows or rows[0][0] != 0 or abs(rows[-1][1] - t_end) > tol * t_end:
            errors.append(f"gauge {name}: the lines do not run from step 0 to T = {t_end}")
            continue
        steps = [r[0] for r in rows]
        if steps != list(range(len(rows))):
            errors.append(f"gauge {name}: not one line per step; --rain-basin needs log_int = 1")
        for step, t, d, zs, u, v in rows:
            want = rate * t
            zs_want = a.bed + want
            if not (abs(d - want) <= tol * max(want, 1e-3)
                    and abs(zs - zs_want) <= tol * max(abs(zs_want), 1e-3)
                    and u == 0.0 and v == 0.0):
                errors.append(f"gauge {name} step {step:.0f}: depth {d!r}, surface {zs!r},"
                              f" velocity ({u!r}, {v!r}); expected {want!r}, {zs_want!r} and 0")
                break
        fw = next((r for r in rows if r[2] > arrival_depth), None)
        if fw is None:
            errors.append(f"gauge {name}: the depth never passes {arrival_depth} m")
            continue
        if first_wet is None:
            first_wet = fw
        elif fw[1] != first_wet[1]:
            errors.append(f"gauge {name}: water arrives at {fw[1]!r} s, another gauge at"
                          f" {first_wet[1]!r} s")
    if first_wet is None:
        errors.append("no gauge gives the arrival time")
        return
    k = int(first_wet[0])
    _, rows = read_gauge(f"gauge_{a.gauge[0][0]}.txt")
    t_exact = arrival_depth / rate
    if not rows[k - 1][1] <= t_exact <= first_wet[1]:
        errors.append(f"the first wet line is at {first_wet[1]!r} s, the one before at"
                      f" {rows[k - 1][1]!r} s; the rain gives {arrival_depth} m at {t_exact!r} s")
    arrival = as_state(first_wet[1])
    for c, v in cells.items():
        bad = []
        if not abs(v["max_depth"] - depth) <= tol * max(depth, 1e-3):
            bad.append(f"max_depth {v['max_depth']!r}, expected {depth!r}")
        if v["max_speed"] != 0.0 or v["max_depth_speed"] != 0.0:
            bad.append(f"max_speed {v['max_speed']!r}, max_depth_speed {v['max_depth_speed']!r}")
        if v["arrival_time"] != arrival:
            bad.append(f"arrival_time {v['arrival_time']!r}, expected {arrival!r}")
        if bad:
            errors.append(f"cell {c}: " + "; ".join(bad))
            break
    print(f"rain basin: {len(cells)} cells hold max_depth {depth!r} m and arrival"
          f" {first_wet[1]!r} s (exact {t_exact!r} s); {len(a.gauge)} gauge(s) checked")


def dam_break(a, h, cells, eps, errors):
    h1, h0, x0, g, arrival_depth = a.dam_break
    hm, um = cd.wet_bed_middle(h1, h0, g)
    if not h0 < arrival_depth < hm:
        sys.exit(f"FAIL: ARRIVAL_DEPTH {arrival_depth} is not between H0 {h0} and hm {hm}")
    s = hm * um / (hm - h0)
    t_end = h["time"]
    dx = h["dx"][0]
    xlo = h["prob_lo"][0]
    x_shock = x0 + s * t_end
    # The maps follow the state at step ends, so an arrival may be one step
    # late: the largest step, from the times of the log (one line a step).
    _, _, log_rows = cb.read_log("water_balance.txt")
    dt_max = max(b[1] - a_[1] for a_, b in zip(log_rows, log_rows[1:]))
    tol_t = dx / s + dt_max
    n_passed = n_ahead = 0
    worst_t = worst_h = worst_u = 0.0
    j = h["domain_lo"][1]
    for i in range(h["domain_lo"][0], h["domain_hi"][0] + 1):
        v = cells[(i, j)]
        x = xlo + (i + 0.5) * dx
        if x0 + 8 * dx <= x <= x_shock - 8 * dx:
            n_passed += 1
            off_t = abs(v["arrival_time"] - (x - x0) / s)
            off_h = abs(v["max_depth"] - hm) / hm
            off_u = abs(v["max_speed"] - um) / um
            off_q = abs(v["max_depth_speed"] - hm * um) / (hm * um)
            worst_t, worst_h = max(worst_t, off_t), max(worst_h, off_h)
            worst_u = max(worst_u, off_u, off_q)
            if not (off_t <= tol_t and off_h <= 0.02 and off_u <= 0.02 and off_q <= 0.02):
                errors.append(f"cell {i}: arrival {v['arrival_time']!r} s (shock {(x - x0) / s!r} s),"
                              f" max_depth {v['max_depth']!r} (hm {hm!r}), max_speed"
                              f" {v['max_speed']!r} (um {um!r}), max_depth_speed"
                              f" {v['max_depth_speed']!r} (hm um {hm * um!r})")
                break
        elif x >= x_shock + 8 * dx:
            n_ahead += 1
            if not (v["arrival_time"] == -1.0 and v["max_depth"] < arrival_depth):
                errors.append(f"cell {i} ahead of the shock: arrival {v['arrival_time']!r},"
                              f" max_depth {v['max_depth']!r}")
                break
    if n_passed < 8 or n_ahead < 8:
        errors.append(f"only {n_passed} cells behind and {n_ahead} ahead of the shock to check")
    print(f"dam break: S = {s:.6f} m/s, hm = {hm:.6f} m, um = {um:.6f} m/s; {n_passed} cells passed:"
          f" arrival within {worst_t:.4f} s ({tol_t:.4f} allowed, a cell over S plus the"
          f" largest step), max_depth within"
          f" {100 * worst_h:.3f} %, the speed maps within {100 * worst_u:.3f} %; {n_ahead} cells"
          f" ahead")
    for name, gi, gj in a.gauge:
        cell, rows = read_gauge(f"gauge_{name}.txt")
        if cell != (int(gi), int(gj)):
            errors.append(f"gauge {name} read cell {cell}, expected ({int(gi)}, {int(gj)})")
            continue
        _, t, d, zs, u, v = rows[-1]
        print(f"gauge {name} at {t!r} s: depth {d!r}, surface {zs!r}, velocity ({u!r}, {v!r})")
        if not (abs(t - t_end) <= 64 * eps * t_end and abs(d - hm) <= 0.02 * hm and zs == d
                and abs(u - um) <= 0.02 * um and v == 0.0):
            errors.append(f"gauge {name}: last line {rows[-1]}; expected T {t_end}, depth {hm},"
                          f" velocity ({um}, 0)")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("plotfile")
    p.add_argument("--rain-basin", type=float, nargs=2, metavar=("RATE", "ARRIVAL_DEPTH"))
    p.add_argument("--bed", type=float, default=0.0)
    p.add_argument("--gauge", nargs=3, action="append", default=[], metavar=("NAME", "I", "J"))
    p.add_argument("--dam-break", type=float, nargs=5, metavar=("H1", "H0", "X0", "G", "ARRIVAL_DEPTH"))
    a = p.parse_args()
    if (a.rain_basin is None) == (a.dam_break is None):
        sys.exit("FAIL: give one of --rain-basin and --dam-break")
    if a.rain_basin is not None and not a.gauge:
        sys.exit("FAIL: --rain-basin needs at least one --gauge")
    plotfile = last_plotfile(a.plotfile)
    h, cells, eps = read(plotfile)
    errors = []
    if a.rain_basin is not None:
        rain_basin(a, h, cells, eps, errors)
    else:
        dam_break(a, h, cells, eps, errors)
    if errors:
        print(f"FAIL: {plotfile}")
        for e in errors:
            print("  " + e)
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
