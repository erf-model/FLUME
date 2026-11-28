#!/usr/bin/env python3
"""Write the rasters of the storm case: bed, curve number and Manning n.

The flood grid is ERF's domain, x from -32 to 32 km (periodic) and y from
-16 to 16 km, in cells of 250 m (ERF's 2 km columns refined by 8). The
rasters have the same cells, so each flood cell takes its raster cell.

The bed is four parallel valleys, 16 km apart. Across a valley the bed
rises 1 % from its axis to the ridges 8 km either side. Along it, a divide
at y = 4 km sends the water south of it to the ylo side and the water north
of it to the yhi side, each with a fall of 1.5 %. A channel 2 m deep and two
cells (500 m) wide runs down each axis. The catchment of the case is the
part of the valley whose axis is at x = -6 km south of the divide (x from -14 to 2 km,
y from -16 to 4 km, 320 km2), which leaves the domain through the outflow
segment ``outlet``. A town of 8 x 8 km (x from -10 to -2 km, y from -12 to
-4 km), under the storm's track, has curve number 98 and n = 0.015; the
rest of the ground has 90 and 0.04, the channels n = 0.03.

Usage: prepare_catchment.py <directory>
Writes bed.asc, cn.asc and manning.asc there.
"""
import pathlib
import sys

X_LO, X_HI, Y_LO, Y_HI = -32000.0, 32000.0, -16000.0, 16000.0
DX = 250.0
NX, NY = int((X_HI - X_LO) / DX), int((Y_HI - Y_LO) / DX)
VALLEY_AXIS, VALLEY_SPACING = -6000.0, 16000.0
CROSS_SLOPE, DOWN_SLOPE = 0.01, 0.015
Y_DIVIDE = 4000.0
CHANNEL_DEPTH = 2.0
CATCHMENT = (-14000.0, 2000.0, Y_LO, Y_DIVIDE)  # the valley at -6 km, south of the divide
TOWN = (-10000.0, -2000.0, -12000.0, -4000.0)
CN_TOWN, CN_RURAL = 98, 90
N_TOWN, N_RURAL, N_CHANNEL = 0.015, 0.04, 0.03


def axis_distance(x):
    """Distance from x to the nearest valley axis, m."""
    s = (x - VALLEY_AXIS) % VALLEY_SPACING
    return min(s, VALLEY_SPACING - s)


def in_channel(x):
    return axis_distance(x) <= 0.5 * DX


def in_town(x, y):
    return TOWN[0] <= x < TOWN[1] and TOWN[2] <= y < TOWN[3]


def bed(x, y):
    along = y - Y_LO if y < Y_DIVIDE else 2.0 * Y_DIVIDE - Y_LO - y
    z = DOWN_SLOPE * along + CROSS_SLOPE * axis_distance(x)
    return z - CHANNEL_DEPTH if in_channel(x) else z


def write(path, value, fmt):
    with open(path, "w") as f:
        f.write(f"ncols {NX}\nnrows {NY}\nxllcorner {X_LO:.1f}\nyllcorner {Y_LO:.1f}\n"
                f"cellsize {DX:.1f}\n")
        for j in reversed(range(NY)):  # the northern row first
            y = Y_LO + (j + 0.5) * DX
            f.write(" ".join(fmt.format(value(X_LO + (i + 0.5) * DX, y)) for i in range(NX)))
            f.write("\n")


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    out = pathlib.Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    write(out / "bed.asc", bed, "{:.3f}")
    write(out / "cn.asc", lambda x, y: CN_TOWN if in_town(x, y) else CN_RURAL, "{:d}")
    write(out / "manning.asc",
          lambda x, y: N_CHANNEL if in_channel(x) else (N_TOWN if in_town(x, y) else N_RURAL), "{:.3f}")
    print(f"prepare_catchment.py: {NX} x {NY} cells of {DX:g} m in {out}")


if __name__ == "__main__":
    main()
