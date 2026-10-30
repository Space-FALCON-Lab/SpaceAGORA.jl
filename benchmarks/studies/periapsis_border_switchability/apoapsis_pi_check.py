"""Pre-registered checks of the apoapsis Pi parameters (apoapsis_pi_preregistration.md, local task folder).
D: drag decay per pass against the engine truth (Mars 20x20 and Earth 20x20 per-pass tables under the results root).
B: passes and days per Odyssey phase (Smith & Bell 2005 as tabulated in Falcone & Putnam 2023, Table III) - no fitted value.
Usage: python3 -I apoapsis_pi_check.py [--results-root DIR]  -> <root>/apoapsis_pi/{dev_pairs.csv, odyssey_phases.csv, criteria.csv, diagnostics.csv}
"""
import argparse, glob, pathlib, sys, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root, REPO
import apoapsis_pi_lib as L
from apoapsis_pi_selftest import earth_h_layer, MU_MARS, H_MARS, MU_EARTH

MU_MARS50C = 4.28283703710e13      # SOURCED: Mars50c.csv header gm_m3s2 (the GM of the Mars truth elements)
OMEGA_MARS = 7.08823596e-5          # SOURCED: Mars() struct
OMEGA_EARTH = 7.2921158553e-5       # SOURCED: Earth() struct
MIN_PRED_M = 5000.0                 # pre-registered pair filter: |predicted d(r_a + r_p)| >= 5 km
EARTH_MIN_ALT_M = 130e3             # pre-registered: single fit layer (130-150 km, extended above) for the whole pass


def mars_beta():
    with open(HERE / "mars_generality_config_g20.toml", "rb") as f:
        v = tomllib.load(f)["vehicle"]
    area = v["bus_dims_m"][0] * v["bus_dims_m"][2] + 2 * v["panel_dims_m"][1] * v["panel_dims_m"][2]
    return v["cd"], area


def earth_beta(root):
    p = root / "study_config_g20.toml"
    with open(p if p.exists() else HERE / "study_config.toml", "rb") as f:
        v = tomllib.load(f)["vehicle"]
    return v["cd"], v["bus_ref_area_m2"] + v["n_wings"] * v["wing_ref_area_m2"]


def pairs(path, planet, cd, area, mu, omega, Hfun, min_alt=-np.inf):
    d = pd.read_csv(path)
    if len(d) < 3:
        return []
    out = []
    for k in range(1, len(d) - 1):           # exclude pass 1 (no earlier apoapsis) and the last row (apoapsis fallback)
        r, r0 = d.iloc[k], d.iloc[k - 1]
        H = Hfun(r.alt_m); beta = r.mass_kg / (cd * area)
        fr = L.f_rot(r.radius_m, r.a_m, r.vrel_mps, r.lat_rad, omega, mu)
        pred = 2 * L.delta_a_pass(r.a_m, r.e, r.radius_m, r.q_w_cm2 * 1e4, H, beta, mu, fr)
        obs = (r.apo_alt_m + r.alt_m) - (r0.apo_alt_m + r0.alt_m)
        out.append(dict(planet=planet, file=pathlib.Path(path).name, tag=pathlib.Path(path).parent.name, pass_=int(r["pass"]),
                        e=r.e, a_m=r.a_m, alt_m=r.alt_m, q_w_cm2=r.q_w_cm2, H_m=H, beta=beta, f_rot=fr, c=r.a_m * r.e / H,
                        pred_m=pred, obs_m=obs, obs_raw_m=r.apo_alt_m - r0.apo_alt_m,
                        keep=bool(abs(pred) >= MIN_PRED_M and r.alt_m >= min_alt)))
    return out


def frac_within(x, tol):
    return float(np.mean(np.abs(x - 1) <= tol))


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args()
    root = pathlib.Path(results_root(a.results_root)); out = root / "apoapsis_pi"; out.mkdir(parents=True, exist_ok=True)
    crit, diag = [], []
    def C(cid, stmt, value, lo, hi):
        crit.append(dict(id=cid, statement=stmt, value=value, lo=lo, hi=hi, passed=bool(lo <= value <= hi)))

    # ---- D: development, engine truth ---------------------------------------------------------
    cdm, am = mars_beta(); cde, ae = earth_beta(root)
    rows = []
    for fn in ("passes_3b.csv", "passes_deg3.csv"):
        for p in sorted(glob.glob(str(root / "mars" / "truth" / "*" / fn))):
            rows += pairs(p, "Mars", cdm, am, MU_MARS50C, OMEGA_MARS, lambda h: H_MARS)
    for p in sorted(glob.glob(str(root / "earth_g20" / "grid" / "*" / "passes_3b.csv"))):
        rows += pairs(p, "Earth", cde, ae, MU_EARTH, OMEGA_EARTH, earth_h_layer, EARTH_MIN_ALT_M)
    P = pd.DataFrame(rows); P["ratio"] = P.obs_m / P.pred_m; P["ratio_raw"] = P.obs_raw_m / P.pred_m
    P.to_csv(out / "dev_pairs.csv", index=False)
    K = P[P.keep]
    m3 = K[(K.planet == "Mars") & (K.file == "passes_3b.csv")]; md = K[(K.planet == "Mars") & (K.file == "passes_deg3.csv")]; ea = K[K.planet == "Earth"]
    C("D1", "Mars 20x20: median obs/pred", float(m3.ratio.median()), 0.97, 1.03)
    C("D2", "Mars 20x20: fraction |ratio-1|<=0.05", frac_within(m3.ratio, 0.05), 0.90, 1.0)
    C("D3", "Earth 20x20: median obs/pred", float(ea.ratio.median()), 0.97, 1.03)
    C("D4", "Earth 20x20: fraction |ratio-1|<=0.05", frac_within(ea.ratio, 0.05), 0.90, 1.0)
    for name, g in (("Mars 20x20", m3), ("Mars deg3 (secondary)", md), ("Earth 20x20", ea)):
        diag.append(dict(item=f"{name} pairs kept / cells", value=f"{len(g)} / {g.tag.nunique()}"))
        diag.append(dict(item=f"{name} ratio p05/p50/p95", value="%.4f / %.4f / %.4f" % tuple(g.ratio.quantile([.05, .5, .95]))))
        diag.append(dict(item=f"{name} raw-apoapsis ratio p05/p50/p95 (secondary)", value="%.4f / %.4f / %.4f" % tuple(g.ratio_raw.quantile([.05, .5, .95]))))
        diag.append(dict(item=f"{name} fraction within 5% (secondary for deg3)", value=f"{frac_within(g.ratio, 0.05):.4f}"))
        for lo, hi in ((0, 0.3), (0.3, 0.6), (0.6, 0.75), (0.75, 1.0)):
            s = g[(g.e >= lo) & (g.e < hi)]
            if len(s):
                diag.append(dict(item=f"{name} UNFROZEN by e in [{lo},{hi}): n, median ratio", value=f"{len(s)}, {s.ratio.median():.4f}"))

    # ---- B: blind, Odyssey flight phases (no value fitted anywhere) -------------------------------
    with open(HERE / "mars_generality_config_g20.toml", "rb") as f:
        v = tomllib.load(f)["vehicle"]
    beta_o = (v["bus_mass_kg"] + 2 * v["panel_mass_each_kg"] + v["prop_mass_kg"]) / (v["cd"] * am)   # DERIVED (Cd ASSUMED)
    rp = 3505e3          # DERIVED: Stage C calibrated Mars periapsis radius (kaula_harmonic_parameter.md section 4)
    mu = MU_MARS; H = H_MARS
    # SOURCED (Table III): apoapsis radii at phase boundaries [km], passes, days, median heat rate [W/cm2], periods [h]
    ph = pd.DataFrame(dict(phase=["Main Phase I", "Main Phase II", "Endgame", "Walkout"],
                           ra_start_km=[28533, 10038, 7408, 4904], ra_end_km=[10038, 7408, 4904, 3900],
                           passes=[108, 45, 77, 88], days=[41.8, 7.2, 8.7, 7.4], q_med=[0.24, 0.20, 0.13, 0.04],
                           T_start_h=[17.3, 4.7, 3.4, 2.3], T_end_h=[4.7, 3.4, 2.3, 1.9]))
    e_of = lambda ra: (ra - rp) / (ra + rp)
    ph["e_start"] = e_of(ph.ra_start_km * 1e3); ph["e_end"] = e_of(ph.ra_end_km * 1e3)
    ph["kQ"] = L.k_Q(rp, ph.q_med * 1e4, H, beta_o, mu)
    ph["N_pred"] = L.passes_to_go(ph.e_start, ph.e_end, ph.kQ)
    ph["days_pred"] = L.time_to_go(rp, ph.e_start, ph.e_end, ph.kQ, mu) / 86400
    ph["N_ratio"] = ph.N_pred / ph.passes; ph["days_ratio"] = ph.days_pred / ph.days
    ph["dE15"] = ph.e_start**1.5 - ph.e_end**1.5
    ph["NQ_over_dE15"] = ph.passes * ph.q_med / ph.dE15        # beta-, H-free shape quantity
    # derived descriptive quantities at the phase end (late-campaign regime)
    a_end = (ph.ra_end_km * 1e3 + rp) / 2; Vp = np.sqrt(mu * (1 + ph.e_end) / rp)
    ph["PiE_end"] = L.pi_E(a_end, ph.e_end, rp, ph.q_med * 1e4, H, beta_o, mu)
    ph["dra_per_pass_end_km"] = 2 * L.delta_a_pass(a_end, ph.e_end, rp, ph.q_med * 1e4, H, beta_o, mu) / 1e3
    ph["Lambda_end"] = L.lam_rp(ph.e_end, ph.PiE_end)
    ph["N_life_end"] = L.n_life(ph.Lambda_end, ph.q_med, 0.54)
    ph["rp_from_period_end_km"] = 2 * (mu * (ph.T_end_h * 3600 / (2 * np.pi)) ** 2) ** (1 / 3) / 1e3 - ph.ra_end_km
    ph.to_csv(out / "odyssey_phases.csv", index=False)
    for i, cid in ((0, "B1"), (1, "B2"), (2, "B3")):
        C(cid, f"{ph.phase[i]}: N_pred/N_flight", float(ph.N_ratio[i]), 0.70, 1.30)
    C("B4", "MP I + MP II + Endgame: sum N_pred / sum N_flight", float(ph.N_pred[:3].sum() / ph.passes[:3].sum()), 0.70, 1.30)
    for i, cid in ((0, "B5a"), (1, "B5b"), (2, "B5c")):
        C(cid, f"{ph.phase[i]}: days_pred/days_flight", float(ph.days_ratio[i]), 0.70, 1.30)
    s = ph.NQ_over_dE15
    C("B6a", "shape (beta-, H-free): (N Q / d e^1.5) MP I / MP II", float(s[0] / s[1]), 0.80, 1.20)
    C("B6b", "shape (beta-, H-free): (N Q / d e^1.5) MP II / Endgame", float(s[1] / s[2]), 0.80, 1.20)
    # X1 (exploratory): granularity cap at the Walkout start, delta_tol = 10 km (Falcone & Putnam reward window)
    a_w = (4904e3 + rp) / 2; e_w = e_of(4904e3)
    dra_unit = abs(2 * L.delta_a_pass(a_w, e_w, rp, 1.0, H, beta_o, mu))     # |d r_a| per pass per W/m2
    q_gran = 2 * 10e3 / dra_unit / 1e4
    C("X1", "EXPLORATORY: Q_gran(Walkout start, 10 km) / 0.04", float(q_gran / 0.04), 1 / 3, 3.0)
    # UNFROZEN diagnostics (never in a verdict)
    for rpx in (3495e3, 3515e3):
        es, ee = (ph.ra_start_km * 1e3 - rpx) / (ph.ra_start_km * 1e3 + rpx), (ph.ra_end_km * 1e3 - rpx) / (ph.ra_end_km * 1e3 + rpx)
        n = L.passes_to_go(es, ee, L.k_Q(rpx, ph.q_med * 1e4, H, beta_o, mu)) / ph.passes
        diag.append(dict(item=f"UNFROZEN Odyssey N_pred/N_flight with r_p = {rpx/1e3:.0f} km", value=" / ".join(f"{x:.3f}" for x in n)))
    implied = (ph.dE15 / (3 * ph.passes)) / L.k_Q(rp, ph.q_med * 1e4, H, beta_o, mu)
    diag.append(dict(item="UNFROZEN Odyssey implied k_Q(flight)/k_Q(model) per phase (= sqrt(H_eff/H) beta/beta_eff)", value=" / ".join(f"{x:.3f}" for x in implied)))
    diag.append(dict(item="Odyssey beta used [kg/m2]", value=f"{beta_o:.3f}"))
    diag.append(dict(item="Odyssey Q_gran at Walkout start (delta_tol 10 km) [W/cm2]", value=f"{q_gran:.4f}"))
    with open(pathlib.Path(REPO) / "benchmarks" / "studies" / "telemetry_validation" / "manifests" / "odyssey_tolson.toml", "rb") as f:
        mn = tomllib.load(f)["scenarios"][0]["maneuvers"]
    dv = np.array(mn["delta_v_mps"]); used = dv[np.array(mn["orbit_numbers"]) > 16]    # drop walk-in burns 7-16 (manifest note)
    vcp = np.sqrt(mu / rp); dvprop = vcp * (np.sqrt(1 + ph.e_start[0]) - np.sqrt(1 + ph.e_end[3]))
    diag.append(dict(item="Odyssey record: sum |dv| corridor burns (excl. walk-in, excl. final raise) / final raise [m/s]",
                     value=f"{np.abs(used[:-1]).sum():.2f} / {used[-1]:.2f}"))
    diag.append(dict(item="Odyssey propulsive dv replaced (periapsis burn, r_p fixed, 28533->3900 km) [m/s]; Pi_B = all burns / that",
                     value=f"{dvprop:.1f}; {np.abs(used).sum() / dvprop:.4f}"))
    pd.DataFrame(crit).to_csv(out / "criteria.csv", index=False); pd.DataFrame(diag).to_csv(out / "diagnostics.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
    print(pd.DataFrame(crit).to_string(index=False)); print(pd.DataFrame(diag).to_string(index=False))
    print(ph.to_string(index=False))


if __name__ == "__main__":
    main()
