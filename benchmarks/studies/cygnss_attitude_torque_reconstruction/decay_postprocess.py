#!/usr/bin/env python3
"""Decay-rate calibration post-processor for run_drag_decay_calibration.jl.

Extracts the secular vis-viva SMA slope from the flight telemetry and from
each fixed-cd simulation, using the same smoothing-and-fit recipe on both so
estimator bias cancels in the comparison, and additionally subtracts the
drag-free zero reference (a conservative field cannot change SMA secularly,
so any apparent slope a drag-free run shows is measurement chain, common to
all runs). Solves the linear relation for the decay-matched cd-scale.

Usage:
  python3 decay_postprocess.py [--zero-ref errors_table.csv]

Reads data/decay_errors_cd*.csv written by the driver and the private
telemetry feather. The optional zero reference is an errors table from a
drag-free run (e.g. the +SRP rung of run_physics_ladder.jl); without it the
raw slopes are reported uncorrected, with a warning.
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd

MU = 3.98600436233e14
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
TELEMETRY = os.path.join(REPO, "data", "telemetry", "CYGNSS", "cygnss_data_48hr.feather")
EVENTS = ["state_x_time", "state_y_time", "state_z_time",
          "state_vx_time", "state_vy_time", "state_vz_time"]


def smoothed_slope(t, sma, period_s):
    ser = pd.Series(sma, index=t)
    w = int(period_s)
    sm = ser.rolling(w, center=True, min_periods=int(w * 0.9)).mean().dropna()
    a = np.vstack([sm.index.values, np.ones(len(sm))]).T
    return np.linalg.lstsq(a, sm.values, rcond=None)[0][0]


def sim_series(path, valid_by_time):
    df = pd.read_csv(path, usecols=["event", "telemetry_axis", "sim_interp_value_km"])
    df = df[df.event.isin(EVENTS)]
    piv = df.pivot_table(index="telemetry_axis", columns="event", values="sim_interp_value_km")
    t = piv.index.values
    vmask = valid_by_time.reindex(t, method="nearest").values
    r = np.sqrt(piv.state_x_time**2 + piv.state_y_time**2 + piv.state_z_time**2).values * 1e3
    s = np.sqrt(piv.state_vx_time**2 + piv.state_vy_time**2 + piv.state_vz_time**2).values * 1e3
    sma = 1.0 / (2.0 / r - s**2 / MU)
    return t[vmask], sma[vmask]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zero-ref", default=None,
                    help="errors table of a drag-free run used as the estimator-bias zero reference")
    args = ap.parse_args()

    tel = pd.read_feather(TELEMETRY)
    valid_by_time = tel.set_index("time").pvt_valid == 2
    v = (tel.pvt_valid == 2).values
    t = tel.time.values[v]
    r = np.sqrt(tel.pos_ii_1**2 + tel.pos_ii_2**2 + tel.pos_ii_3**2).values[v]
    s = np.sqrt(tel.vel_ii_1**2 + tel.vel_ii_2**2 + tel.vel_ii_3**2).values[v]
    sma_t = 1.0 / (2.0 / r - s**2 / MU)
    period = 2 * np.pi * np.sqrt(np.mean(sma_t) ** 3 / MU)
    tele_slope = smoothed_slope(t, sma_t, period)
    print(f"telemetry raw slope: {tele_slope * 86400:+.2f} m/day (orbit period {period:.0f} s)")

    bias = 0.0
    if args.zero_ref:
        tz, sz = sim_series(args.zero_ref, valid_by_time)
        bias = smoothed_slope(tz, sz, period)
        print(f"zero-reference (drag-free) apparent slope: {bias * 86400:+.2f} m/day — subtracted from all")
    else:
        print("WARNING: no --zero-ref given; slopes below are NOT bias-corrected")

    tele_corr = tele_slope - bias
    print(f"telemetry corrected decay: {tele_corr * 86400:+.2f} m/day")

    cds, slopes = [], []
    for path in sorted(glob.glob(os.path.join(HERE, "data", "decay_errors_cd*.csv"))):
        m = re.search(r"cd(\d+p?\d*)", os.path.basename(path))
        cd = float(m.group(1).replace("p", "."))
        ts, ss = sim_series(path, valid_by_time)
        sl = smoothed_slope(ts, ss, period) - bias
        cds.append(cd)
        slopes.append(sl)
        print(f"cd={cd}: corrected sim decay {sl * 86400:+.2f} m/day  ratio {sl / tele_corr:.3f}")

    if len(cds) >= 2:
        b, a = np.polyfit(cds, slopes, 1)
        cd_star = (tele_corr - a) / b
        print(f"linear model: slope = {a * 86400:+.2f} + {b * 86400:+.2f}*cd m/day")
        print(f"decay-matched cd_scale* = {cd_star:.3f}")


if __name__ == "__main__":
    main()
