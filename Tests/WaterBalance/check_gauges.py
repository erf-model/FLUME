#!/usr/bin/env python3
"""Check gauge files against the cells of a plotfile written at the same step.

    check_gauges.py <plotfile> --gauge NAME I J [--gauge ...]

<plotfile> is the prefix of the run's plotfiles (flume.plot_file); the last
of them, with the highest step, is read. Each --gauge NAME I J must have read
cell (I, J) (the first line of gauge_NAME.txt), and the last line of its
file must be at the step and time of the plotfile and hold exactly that
cell's depth and surface elevation, and its discharges over the depth as the
velocity (0 at or below DRY_DEPTH, default 1e-6 m, the default of
flume.dry_depth). The division is rounded to the precision of the plotfile,
as the solver does it. On a field that varies from cell to cell, a gauge that
reads the wrong cell, or the values of another gauge, fails.
"""
import argparse
import struct
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import check_plotfile as cp  # noqa: E402
import check_flood_outputs as cf  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("plotfile")
    p.add_argument("--gauge", nargs=3, action="append", required=True,
                   metavar=("NAME", "I", "J"))
    p.add_argument("--dry-depth", type=float, default=1e-6)
    a = p.parse_args()
    plotfile = cf.last_plotfile(a.plotfile)
    h = cp.read_header(plotfile / "Header")
    _, fab_on_disk = cp.read_cell_h(plotfile / "Level_0" / "Cell_H")
    raw, real_size = cp.read_cell_data(plotfile / "Level_0", fab_on_disk)
    cells = {c: dict(zip(h["vars"], v)) for c, v in raw.items()}
    single = real_size == 4

    def rnd(x):
        return struct.unpack("f", struct.pack("f", x))[0] if single else x

    step = int(plotfile.name[len(a.plotfile):])
    errors = []
    seen = set()
    for name, gi, gj in a.gauge:
        want = (int(gi), int(gj))
        cell, rows = cf.read_gauge(f"gauge_{name}.txt")
        if cell != want:
            errors.append(f"gauge {name} read cell {cell}, expected {want}")
            continue
        c = cells[want]
        d = c["depth"]
        wet = d > a.dry_depth
        u = rnd(c["discharge_x"] / d) if wet else 0.0
        v = rnd(c["discharge_y"] / d) if wet else 0.0
        expect = [step, h["time"], d, c["surface_elevation"], u, v]
        # The file prints max_digits10 digits, which read back to the value.
        got = rows[-1][:2] + [rnd(x) for x in rows[-1][2:]]
        print(f"gauge {name} cell {cell}: {got}; plotfile {expect}")
        if got[0] != step or got[2:] != expect[2:]:
            errors.append(f"gauge {name}: last line {got}, expected {expect}")
        seen.add(tuple(expect[2:]))
    if len(seen) != len(a.gauge) and not errors:
        errors.append("two gauges read the same values: the field does not tell the cells apart")
    if errors:
        print(f"FAIL: {plotfile}")
        for e in errors:
            print("  " + e)
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
