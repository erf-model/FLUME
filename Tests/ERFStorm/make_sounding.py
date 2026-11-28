#!/usr/bin/env python3
"""Write the sounding of the storm case (Tests/ERFStorm/sounding_storm).

Thermodynamics: the analytic profile of Weisman and Klemp (1982), with the
constants of ERF's supercell case (Source/Prob/ERF_InitCustomPert_SuperCell.H):
theta_0 = 300 K, theta_tr = 343 K, z_tr = 12 km, T_tr = 213 K, relative
humidity 1 - 3/4 (z/z_tr)^(5/4), water vapour capped at 16 g/kg near the
ground (this case's choice; ERF's supercell case caps it at 14 g/kg). Wind: westerly, u = u_0 + u_s tanh(z / z_s),
v = 0, so the cell moves east.

Usage: make_sounding.py [output]   (default: sounding_storm next to it)
"""
import math
import pathlib
import sys

G, CP, RD, RV, P0 = 9.81, 1004.5, 287.0, 461.6, 1.0e5
THETA0, THETA_TR, Z_TR, T_TR = 300.0, 343.0, 12000.0, 213.0
QV_MAX = 0.016
U0, US, ZS = 6.0, 8.0, 3000.0
DZ, ZTOP = 250.0, 16500.0


def theta(z):
    if z <= Z_TR:
        return THETA0 + (THETA_TR - THETA0) * (z / Z_TR) ** 1.25
    return THETA_TR * math.exp(G * (z - Z_TR) / (CP * T_TR))


def rh(z):
    return 1.0 - 0.75 * (z / Z_TR) ** 1.25 if z <= Z_TR else 0.25


def qsat(t, p):
    # Tetens over water, Pa
    es = 610.78 * math.exp(17.27 * (t - 273.15) / (t - 35.86))
    return 0.622 * es / (p - es)


def main():
    out = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).with_name("sounding_storm")
    # Integrate the Exner function up with the virtual potential temperature.
    rows = []
    pi = 1.0  # surface pressure 1000 hPa
    z = 0.0
    qv = 0.0
    while z <= ZTOP + 1e-6:
        th = theta(z)
        p = P0 * pi ** (CP / RD)
        qv = min(QV_MAX, rh(z) * qsat(th * pi, p))
        rows.append((z, th, qv))
        # one step of the hydrostatic Exner equation with theta_v at this level
        thv = th * (1.0 + 0.61 * qv)
        z_next = z + DZ
        thv_next = theta(z_next) * (1.0 + 0.61 * qv)
        pi -= G / CP * DZ * 0.5 * (1.0 / thv + 1.0 / thv_next)
        z = z_next
    with open(out, "w") as f:
        f.write(f"1000.0  {rows[0][1]:.4f}  {1000 * rows[0][2]:.4f}\n")
        for z, th, qv in rows:
            u = U0 + US * math.tanh(z / ZS)
            f.write(f"{z:.1f}  {th:.4f}  {1000 * qv:.4f}  {u:.4f}  0.0\n")


if __name__ == "__main__":
    main()
