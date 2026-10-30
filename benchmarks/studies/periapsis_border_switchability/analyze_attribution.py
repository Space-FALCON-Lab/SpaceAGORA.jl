"""Attribution study, Parts 1-2: variant table, overlay plots, spectral analysis. Usage: analyze_attribution.py OUT [part1|part2]"""
import sys, glob, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

import argparse
from results_root import add_results_root_argument, results_root
_ap = argparse.ArgumentParser(description="Attribution study: part1, part2, amp")
_ap.add_argument("out_dir"); _ap.add_argument("part", nargs="?", default="part1")
add_results_root_argument(_ap); _args = _ap.parse_args()
OUT = _args.out_dir; PART = _args.part; ROOT = results_root(_args.results_root)
RE = 6378137.0
VARS = {"v0": "Sun+Moon", "v1": "Moon only", "v2": "Sun only", "v3": "neither", "v4": "Sun+Moon, J3 off", "v5": "Sun+Moon, drag off",
        "v6": "Moon frozen at t0", "v7": "Sun frozen at t0", "v8": "Moon only, drag off", "v9": "Sun only, drag off"}
summ = pd.read_csv(f"{OUT}/attribution_run_summary.csv")
tags = list(dict.fromkeys(summ.tag))

def load(tag, v):
    return pd.read_csv(f"{OUT}/runs/{tag}/passes_{v}.csv")

def detrend(t, y, order):
    return y - np.polyval(np.polyfit(t, y, order), t)

def part1():
    rows = []
    for _, r in summ.iterrows():
        d = load(r.tag, r.variant); t = d.time_s / 86400
        res = detrend(t, d.radius_m, 2) if len(d) > 5 else np.zeros(len(d))
        rows.append(dict(tag=r.tag, variant=r.variant, retcode=r.retcode, passes=len(d), t_end_d=t.iloc[-1],
            ended_early=r.retcode != "Success", alt0_km=d.alt_m.iloc[0] / 1e3, alt_end_km=d.alt_m.iloc[-1] / 1e3,
            dalt_km=(d.alt_m.iloc[-1] - d.alt_m.iloc[0]) / 1e3, drad_km=(d.radius_m.iloc[-1] - d.radius_m.iloc[0]) / 1e3,
            range_rad_km=(d.radius_m.max() - d.radius_m.min()) / 1e3, pp_resid2_km=(res.max() - res.min()) / 1e3 if len(d) > 5 else np.nan,
            min_alt_km=d.alt_m.min() / 1e3))
    T = pd.DataFrame(rows); T.to_csv(f"{OUT}/attribution_table.csv", index=False)
    # attribution at 60 d (or common end): contribution of each source relative to variant v3 (neither)
    att = []
    for tag in tags:
        g = T[T.tag == tag].set_index("variant")
        def dr(v): return g.loc[v, "drad_km"] if v in g.index else np.nan
        att.append(dict(tag=tag, v0=dr("v0"), moon=dr("v1") - dr("v3"), sun=dr("v2") - dr("v3"), moon_plus_sun=dr("v1") + dr("v2") - 2 * dr("v3"),
            nonlinear_cross=dr("v0") - (dr("v1") + dr("v2") - dr("v3")), j3=dr("v0") - dr("v4"), drag=dr("v0") - dr("v5"),
            grav_only_v5=dr("v5"), moon_dragoff=dr("v8"), sun_dragoff=dr("v9"), both_dragoff=dr("v5"),
            cross_dragoff=dr("v5") - dr("v8") - dr("v9"), v0_early=g.loc["v0", "ended_early"], v5_early=g.loc["v5", "ended_early"], v5_tend=g.loc["v5", "t_end_d"]))
    A = pd.DataFrame(att); A.to_csv(f"{OUT}/attribution_contributions.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(T.pivot(index="tag", columns="variant", values="drad_km").round(1))
    print(T.pivot(index="tag", columns="variant", values="pp_resid2_km").round(1))
    print(T.pivot(index="tag", columns="variant", values="t_end_d").round(1))
    print(A.round(1))
    # overlay plots
    for sub, name in ((tags[:4], "f2"), (tags[4:8], "f1"), (tags[8:], "other")):
        fig, ax = plt.subplots(len(sub), 1, figsize=(9, 3.2 * len(sub)), sharex=True, squeeze=False)
        for a, tag in zip(ax[:, 0], sub):
            for v, c in (("v0", "k"), ("v1", "tab:blue"), ("v2", "tab:orange"), ("v5", "tab:green"), ("v3", "gray")):
                d = load(tag, v); a.plot(d.time_s / 86400, (d.radius_m - d.radius_m.iloc[0]) / 1e3, c, label=VARS[v], lw=1.3, ls="--" if v == "v3" else "-")
            a.set_title(tag, fontsize=9); a.set_ylabel("r_p - r_p0 (km)"); a.grid(alpha=.3)
        ax[0, 0].legend(fontsize=7, ncol=3); ax[-1, 0].set_xlabel("time (d)"); fig.tight_layout(); fig.savefig(f"{OUT}/fig_overlay_{name}.png", dpi=110); plt.close(fig)

CAND = {"lunar sidereal 27.32": 27.321661, "lunar synodic 29.53": 29.530589, "lunar anomalistic 27.55": 27.554550,
        "sidereal/2": 27.321661 / 2, "synodic/2": 29.530589 / 2, "anomalistic/2": 27.554550 / 2, "sun annual/2": 365.25 / 2, "sun annual/2 (182.6)": 182.625}

def fit_sines(t, y, periods):
    X = [np.ones_like(t), t, t**2]
    for P in periods: X += [np.sin(2 * np.pi * t / P), np.cos(2 * np.pi * t / P)]
    c, *_ = np.linalg.lstsq(np.array(X).T, y, rcond=None)
    return c, np.array(X).T @ c

def part2():
    from scipy.signal import lombscargle
    rows = []; fig, axs = plt.subplots(4, 3, figsize=(14, 11)); axs = axs.ravel(); k = 0
    pgs = {}
    for tag in tags:
        for v in ("v5", "v8", "v9"):
            d = load(tag, v); t = (d.time_s / 86400).values; y = d.radius_m.values / 1e3
            if t[-1] < 55: continue
            res = detrend(t, y, 2)   # quadratic trend removed (state in report)
            per = np.linspace(3, 60, 3000); w = 2 * np.pi / per
            p = lombscargle(t, res - res.mean(), w, normalize=False)
            amp = np.sqrt(4 * p / len(t))  # sinusoid amplitude estimate
            pgs[(tag, v)] = (per, amp)
            # peaks in 8-40 d
            m = (per > 8) & (per < 40); i = np.argmax(amp * m)
            # orbital-period harmonic content: spacecraft period ~ spacing of passes
            Pk = per[i]
            c, fit = fit_sines(t, y, [Pk]); a_pk = np.hypot(c[3], c[4])
            # fixed candidate fit
            cc, ffit = fit_sines(t, y, [CAND["sidereal/2"], CAND["synodic/2"], CAND["sidereal/2"] * 0 + 27.321661, 365.25 / 2]) if False else (None, None)
            rows.append(dict(tag=tag, variant=v, peak_period_d=Pk, peak_amp_km=amp[i], fit_amp_km=a_pk, resid_rms_km=np.sqrt(np.mean(res**2)),
                a_mean_km=np.mean(d.a_m) / 1e3, e_mean=np.mean(d.e), n=len(t), t_end=t[-1]))
    R = pd.DataFrame(rows); R.to_csv(f"{OUT}/spectral_peaks.csv", index=False)
    print(R.round(3).to_string())
    # periodogram figure for 4 cells x (v5,v8,v9)
    show = [tags[1], tags[5], tags[6], tags[8]]
    fig, axs = plt.subplots(2, 2, figsize=(13, 8))
    for a, tag in zip(axs.ravel(), show):
        for v, c in (("v5", "k"), ("v8", "tab:blue"), ("v9", "tab:orange")):
            if (tag, v) in pgs: a.plot(*pgs[(tag, v)], c, label=VARS[v])
        for nm, P in CAND.items():
            if P < 60: a.axvline(P, color="r", alpha=.3, lw=1)
        a.axvline(13.66, color="m", ls=":", lw=1); a.set_title(tag, fontsize=9); a.set_xlabel("period (d)"); a.set_ylabel("amplitude (km)"); a.set_xlim(3, 45); a.grid(alpha=.3)
    axs[0, 0].legend(fontsize=7); fig.suptitle("Periodogram of detrended periapsis radius (quadratic trend removed); red lines = lunar sidereal/synodic/anomalistic and halves, magenta = 13.66 d"); fig.tight_layout(); fig.savefig(f"{OUT}/fig_periodograms.png", dpi=110); plt.close(fig)

{"part1": part1, "part2": part2}.get(PART, lambda: None)()

def amp():
    """Frozen-ephemeris spectra and the derived quadrupole oscillation amplitude (see attribution_results.md for the expression)."""
    from scipy.signal import lombscargle
    mu_m, a_moon_km, sid = 4.9028005821478e12, 384400.0, 27.321661   # GM table of the study; textbook mean distance; sidereal month (d)
    J2, Re = 1.0826358e-3, 6378.1363e3; mu = 3.986004415e14
    for v in ("v6", "v7", "v1", "v2"):
        d = load("I0234_R000_W090_E1", v); t = d.time_s.values / 86400; y = d.radius_m.values / 1e3
        res = detrend(t, y, 2); per = np.linspace(8, 30, 800); p = lombscargle(t, res - res.mean(), 2*np.pi/per)
        print("frozen/live check", v, "resid p-p after quadratic %.2f km" % (res.max()-res.min()), "peak period %.2f d" % per[np.argmax(p)], "peak amp %.2f" % np.sqrt(4*p.max()/len(t)))
    rows = []
    sp = pd.read_csv(f"{OUT}/spectral_peaks.csv"); sp = sp[sp.variant == "v8"]
    for _, r in sp.iterrows():
        tag = r.tag; src = os.path.join(ROOT, "border2", "grid", tag) if not tag.startswith("I0500") else os.path.join(ROOT, "border", "grid", tag[6:])
        g = pd.read_csv(src + "/passes_3b.csv"); ip, wp = g.moon_i_p_rad.iloc[0], g.moon_omega_p_rad.iloc[0]
        a, e = r.a_mean_km * 1e3, r.e_mean; n = np.sqrt(mu / a**3); p_ = a * (1 - e**2); wdot = 0.75 * n * J2 * (Re / p_)**2 * (5*np.cos(g.i_rad.iloc[0])**2 - 1)
        c = np.cos(ip); lam = np.linspace(0, 2*np.pi, 720, endpoint=False)
        ehat = np.array([np.cos(wp), np.sin(wp)*c, np.sin(wp)*np.sin(ip)]); phat = np.array([-np.sin(wp), np.cos(wp)*c, np.cos(wp)*np.sin(ip)])
        u = np.array([np.cos(lam), np.sin(lam), 0*lam]); f = (ehat @ u) * (phat @ u); F = np.fft.rfft(f) / len(lam); c2 = 2*abs(F[2])
        lam_dot = 2*np.pi / (sid*86400); om_eff = 2 * (lam_dot - wdot); Pexp = 2*np.pi/om_eff/86400
        A = a/1e3 * 7.5 * mu_m / (n * (a_moon_km*1e3)**3) * e*np.sqrt(1-e**2) * c2 / om_eff   # amplitude of r_p oscillation, km
        rows.append(dict(tag=tag, i_p_moon_deg=np.degrees(ip), omega_p_moon_deg=np.degrees(wp), c2lambda=c2, expected_period_d=Pexp, fitted_period_d=r.peak_period_d,
            expected_amp_km=A, measured_fit_amp_km=r.fit_amp_km, ratio=r.fit_amp_km/A, a_km=r.a_mean_km, e=e))
    D = pd.DataFrame(rows); D.to_csv(f"{OUT}/quadrupole_amplitude_check.csv", index=False); print(D.round(3).to_string())

if PART == "amp": amp()
