#!/usr/bin/env python3
"""Draw the figure of the Gowanus case and print its numbers.

Usage:

    plot_gowanus.py <case directory> <run with buildings> <run without> <figure.png>

(a) the bed with the buildings and the FloodNet sensors; (b) the largest
depth of the run with the buildings; (c) the same without them; (d) the
depth at the three sensor positions, with and without the buildings.
The numbers of the README's tables are printed. Needs matplotlib and numpy.
"""
import struct
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from check_plotfile import FAB_HEAD, read_cell_h, read_header  # noqa: E402

DEEP = (0.15, 0.3, 0.5, 1.0)  # depths of the table, m


def read_plotfile(path):
    """The header and the cell values, [component][j, i] from the south-west."""
    h = read_header(path / "Header")
    boxes, fabs = read_cell_h(path / "Level_0" / "Cell_H")
    ni = h["domain_hi"][0] - h["domain_lo"][0] + 1
    nj = h["domain_hi"][1] - h["domain_lo"][1] + 1
    a = np.full((len(h["vars"]), nj, ni), np.nan)
    for name, offset in fabs:
        with open(path / "Level_0" / name, "rb") as f:
            f.seek(offset)
            m = FAB_HEAD.fullmatch(f.readline())
            if not m:
                sys.exit(f"{path}: cannot read a FAB header in {name}")
            size = int(m.group(3))
            order = [int(x) for x in m.group(4).split()]
            lo = [int(m.group(k)) for k in (5, 6, 7)]
            hi = [int(m.group(k)) for k in (8, 9, 10)]
            ncomp = int(m.group(14))
            bi, bj = hi[0] - lo[0] + 1, hi[1] - lo[1] + 1
            kind = ("<" if order[0] == size else ">") + ("f8" if size == 8 else "f4")
            v = np.frombuffer(f.read(bi * bj * ncomp * size), dtype=kind)
            a[:, lo[1]:hi[1] + 1, lo[0]:hi[0] + 1] = v.reshape(ncomp, bj, bi)
    if np.isnan(a).any():
        sys.exit(f"{path}: the boxes do not cover the domain")
    return h, a


def read_table(path):
    """The columns of a water-balance log, by name."""
    names = None
    for line in Path(path).read_text().splitlines():
        if line.startswith("# step "):
            names = line[2:].split()
    d = np.loadtxt(path, ndmin=2)
    return {n: d[:, k] for k, n in enumerate(names)}


def read_case(path):
    """case.txt: {key: values}, and the gauges as (name, x, y)."""
    out, gauges = {}, []
    for line in Path(path).read_text().splitlines():
        k, *v = line.split()
        if k == "gauge":
            gauges.append((v[0], float(v[1]), float(v[2])))
        else:
            out[k] = v
    return out, gauges


def last_plotfile(run):
    files = [p for p in Path(run).glob("flood*") if p.is_dir() and p.name[5:].isdigit()]
    return max(files, key=lambda p: int(p.name[5:]))


def main():
    if len(sys.argv) != 5:
        sys.exit("usage: plot_gowanus.py <case directory> <run with buildings> "
                 "<run without> <figure.png>")
    case, runs, png = Path(sys.argv[1]), [Path(sys.argv[2]), Path(sys.argv[3])], sys.argv[4]
    info, gauge_list = read_case(case / "case.txt")
    size, dx = float(info["size"][0]), float(info["dx"][0])
    rain_time, rain_rate = float(info["rain_time"][0]), float(info["rain_rate"][0])
    rain_m3 = float(info["rain_m3"][0])
    x0, y0 = float(info["x0"][0]), float(info["y0"][0])
    roof = np.loadtxt(case / "buildings.asc", skiprows=5)[::-1]
    solid = roof > 0.0

    labels = ["With buildings", "Without buildings"]
    heads, maps, logs, gauges = [], [], [], []
    for run in runs:
        h, a = read_plotfile(last_plotfile(run))
        heads.append(h)
        maps.append({n: a[k] for k, n in enumerate(h["vars"])})
        logs.append(read_table(run / "water_balance.txt"))
        gauges.append({g[0]: np.loadtxt(run / f"gauge_{g[0]}.txt", ndmin=2) for g in gauge_list})
    bed = maps[1]["bed_elevation"]
    # The canal: the lidar's water surface, not a street. Left out of the numbers.
    canal = bed < float(info["water_below"][0])
    street = ~solid & ~canal

    print(f"square {size:g} m, cells {dx:g} m, building cells {solid.sum()} "
          f"({100 * solid.mean():.1f} %), canal cells {canal.sum()}, rain {rain_m3:g} m3")
    print("run | steps | rain m3 | stored at end m3 | residual m3 | "
          "largest depth off the buildings and the canal m | "
          + " | ".join(f"area above {d:g} m, ha" for d in DEEP)
          + " | " + " | ".join(f"peak at {g[0]} m (time s)" for g in gauge_list))
    for label, h, m, log, gs in zip(labels, heads, maps, logs, gauges):
        deep = [1e-4 * dx * dx * np.count_nonzero(m["max_depth"][street] > d) for d in DEEP]
        peaks = [(g[:, 2].max(), g[g[:, 2].argmax(), 1]) for g in gs.values()]
        print(f"{label} | {h['step']} | {log['rain'][-1]:.1f} | {log['stored'][-1]:.1f} | "
              f"{log['residual'][-1]:.2e} | {m['max_depth'][street].max():.3f} | "
              + " | ".join(f"{d:.2f}" for d in deep) + " | "
              + " | ".join(f"{d:.3f} ({t:.0f})" for d, t in peaks))
    j, i = np.unravel_index(np.argmax(np.where(street, maps[0]["max_depth"], 0.0)), solid.shape)
    print(f"deepest cell with buildings: x0 + {(i + 0.5) * dx:g} m, "
          f"y0 + {(j + 0.5) * dx:g} m, bed {bed[j, i]:.2f} m")
    print(f"canal at the end, with buildings: surface {maps[0]['surface_elevation'][canal].mean():.2f} m "
          f"from a bed of {bed[canal].mean():.2f} m")

    fig, ax = plt.subplots(2, 2, figsize=(12.5, 10.5))
    ext = [0, size, 0, size]
    grey = np.ma.masked_where(~solid, np.ones_like(roof))
    marks = ("^", "s", "o")

    im = ax[0, 0].imshow(bed, origin="lower", extent=ext, cmap="terrain", vmin=-8, vmax=28)
    ax[0, 0].imshow(grey, origin="lower", extent=ext, cmap="Greys", vmin=0, vmax=1.6, alpha=0.75)
    fig.colorbar(im, ax=ax[0, 0], label="bed elevation, m (NAVD88)", shrink=0.85)
    ax[0, 0].set_title("(a) Bed and buildings (grey); FloodNet sensors (red)")

    for k, (a, title) in enumerate(((ax[0, 1], "(b) Largest depth, with buildings"),
                                    (ax[1, 0], "(c) Largest depth, without buildings"))):
        d = maps[k]["max_depth"]
        d = np.ma.masked_where((d < 0.02) | canal | (solid if k == 0 else False), d)
        a.set_facecolor("#f4f1ea")
        im = a.imshow(d, origin="lower", extent=ext, cmap="Blues", vmin=0, vmax=0.6)
        a.imshow(np.ma.masked_where(~canal, np.ones_like(roof)), origin="lower", extent=ext,
                 cmap="GnBu", vmin=0, vmax=1.4)
        if k == 0:
            a.imshow(grey, origin="lower", extent=ext, cmap="Greys", vmin=0, vmax=2.2)
        fig.colorbar(im, ax=a, label="largest depth, m (0.02 m and above)", shrink=0.85,
                     extend="max")
        a.set_title(title + " (canal green)")
    for a in (ax[0, 0], ax[0, 1], ax[1, 0]):
        for (name, gx, gy), mark in zip(gauge_list, marks):
            a.plot(gx - x0, gy - y0, mark, color="r", ms=6, mec="k")
        a.set_xlabel("east of x0, m")
        a.set_ylabel("north of y0, m")
        a.set_aspect("equal")

    a = ax[1, 1]
    a.axvspan(0, rain_time / 60, color="0.92", label=f"rain, {rain_rate:g} mm/h")
    for (name, _, _), color, mark in zip(gauge_list, ("C0", "C1", "C2"), marks):
        for gs, ls, label in zip(gauges, ("-", "--"), labels):
            g = gs[name]
            a.plot(g[:, 1] / 60, g[:, 2], color=color, ls=ls,
                   label=f"{name}, {label.lower()}")
        a.plot([], [], mark, color="r", mec="k", label=f"{name} on the maps")
    a.set_xlabel("time, min")
    a.set_ylabel("depth at the sensor, m")
    a.set_title("(d) Depth at the three sensor positions")
    a.legend(fontsize=7.5, ncol=2, loc="upper left")
    a.set_ylim(0, 1.3)
    a.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(png, dpi=130)
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
