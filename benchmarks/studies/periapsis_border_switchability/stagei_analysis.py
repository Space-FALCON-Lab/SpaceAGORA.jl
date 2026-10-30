"""Stage I analysis: the paper corridor (band 0.05 to 0.25 W/cm2, target 0.15; tiers low < 0.05, soft 0.25-0.30, medium 0.30-0.45, hard > 0.45; strict inequalities on the upper edge of each tier,
as Stage F's q > Q_high / q > Q_hard). Plan: tasks/2026-10-15_periapsis-region-border-targeting/stagei_plan.md (frozen before any run). Reuses stagef20_analysis.py (load_runs, build, rms, boot_ci,
summaries) with its results subdirectory redirected to the stagei tree; the band-dependent metric (tiers) is computed here.
  selftest                       exact-case checks of the tier scoring and of the criteria
  integrity   --results-root R   G1-G3 style checks: route and failures per campaign, GRAM epoch, calibration convergence to Q = 0.15, NaN scan
  params      --results-root R [--ref-root R20] [--planets mars]
                                 development refit of predictor, e_hi, e_lo from the re-flown fixed-margin dev runs (the frozen f2params procedure) -> <R>/stagei/<planet>/adapt_params.csv;
                                 prints the Stage F20 values from --ref-root next to them (information only)
  held        --results-root R [--planets mars,earth,venus]
                                 pre-registered criteria I1-I5 (Mars) and J1, J2 (Earth, Venus) and the four-tier violation counts
Usage: python3 -I stagei_analysis.py --phase ... --results-root DIR"""
import argparse, glob, json, pathlib, sys, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stagef20_analysis as sf

CI = tomllib.load(open(HERE / "stagei_config.toml", "rb"))
SF = dict(CI["stagef"]); sf.SF = SF                       # redirect sf.sfdir() to <root>/stagei
BAND, ST = CI["common"]["band"], CI["stagei"]
Q_LOW, Q_HIGH, Q_MED, Q_HARD, Q_TGT = BAND["q_low_w_cm2"], BAND["q_high_w_cm2"], BAND["q_medium_w_cm2"], BAND["q_hard_w_cm2"], BAND["q_target_w_cm2"]
AN = CI["analysis"]; X = AN["win_cost_tolerance"]
FRAC = ST["e_lo_frac"]; HCTRL = "adaptive_h" + f"{round(100 * FRAC):03d}"
ARMS = {"fixed": "fixed", "F": "adaptive", "G": "adaptive_asym", "H": HCTRL}     # arm name -> ctrl (run directory prefix)
HELD_SUMMARY = "held"
sfdir = sf.sfdir


def tiers(q):
    """Counts of the four violation tiers (plus low) of an array of periapsis heat rates (W/cm2)."""
    q = np.asarray(q, float)
    return dict(low=int((q < Q_LOW).sum()), soft=int(((q > Q_HIGH) & (q <= Q_MED)).sum()), medium=int(((q > Q_MED) & (q <= Q_HARD)).sum()), hard=int((q > Q_HARD).sum()))


def run_metrics(P, srow, dec):
    q = P.q_eff_w_cm2.values; t = P.time_s.values; apo = P.apo_alt_m.values; ok = np.isfinite(apo)
    days = (t[-1] - t[0]) / 86400.0 if len(t) > 1 else np.nan
    decay = (apo[ok][0] - apo[ok][-1]) / days if ok.sum() > 1 and days > 0 else np.nan
    tr = tiers(q); tr["over"] = tr["soft"] + tr["medium"] + tr["hard"]
    return dict(n_pass=len(P), **tr, dv=srow.get("total_dv_mps", np.nan), burns=srow.get("burns", np.nan), decay=decay, retcode=srow.get("retcode"), t_end_days=srow.get("t_end_days"),
                stopped=bool(srow.get("controller_stopped", False)), collapsed=int(srow.get("collapsed_decisions", 0)), n_dec=len(dec))


def cell_table(root, pl, arms=ARMS):
    """One row per held-out cell with run metrics of every arm; cells where any arm is missing are dropped."""
    cells = sf.cells_of(root, pl); S = sf.summaries(root, HELD_SUMMARY); rows = []
    if len(S) == 0: return pd.DataFrame()
    for t in cells[cells.split == "heldout"].tag:
        d = {}
        for k, c in arms.items():
            rd = sfdir(root) / pl / "runs" / f"{c}_pert_{t}"; sr = S[(S.planet == pl) & (S.run_id == f"{c}_pert_{t}")]
            if not (rd / "passes.csv").exists() or len(sr) == 0: d = None; break
            dec = pd.read_csv(rd / "decisions.csv") if (rd / "decisions.csv").exists() else pd.DataFrame()
            d[k] = run_metrics(pd.read_csv(rd / "passes.csv"), sr.iloc[0].to_dict(), dec)
        if d: rows.append(dict(tag=t, **{f"{k}_{m}": v for k, dd in d.items() for m, v in dd.items()}))
    return pd.DataFrame(rows)


def arm_summary(R, a):
    tot = lambda k: float(R[f"{a}_{k}"].sum()); n = tot("n_pass")
    return dict(passes=int(n), low=int(tot("low")), soft=int(tot("soft")), medium=int(tot("medium")), hard=int(tot("hard")), over=int(tot("over")), frac_low=tot("low") / n, frac_soft=tot("soft") / n,
                frac_medium=tot("medium") / n, frac_hard=tot("hard") / n, frac_over=tot("over") / n, dv_mean=float(R[f"{a}_dv"].mean()), burns_mean=float(R[f"{a}_burns"].mean()),
                decay_median_m_per_day=float(R[f"{a}_decay"].median()), t_end_days_mean=float(R[f"{a}_t_end_days"].mean()), collapsed=int(tot("collapsed")), decisions=int(tot("n_dec")),
                bad_runs=int(((R[f"{a}_retcode"] != "Success") | R[f"{a}_stopped"]).sum()))


def compare(R, S, a, rng, planet):
    """Pre-registered comparison of arm `a` with fixed. Returns the measured values and the criteria."""
    f, g = S["fixed"], S[a]
    d = (R[f"{a}_over"] / R[f"{a}_n_pass"] - R["fixed_over"] / R["fixed_n_pass"]).values; lo, hi = sf.boot_ci(d, rng, AN["bootstrap_draws"])
    out = dict(dover_mean=float(np.mean(d)), dover_ci90=[lo, hi], dv_ratio=g["dv_mean"] / f["dv_mean"] if f["dv_mean"] > 0 else float("nan"),
               decay_vs_fixed_pct=100 * (g["decay_median_m_per_day"] / f["decay_median_m_per_day"] - 1), collapse_fraction=g["collapsed"] / g["decisions"] if g["decisions"] else float("nan"))
    if planet == "mars":
        c = dict(I2_hard_le_fixed=bool(g["hard"] <= f["hard"]), I3_over_lt_fixed_ci_excl0=bool(g["frac_over"] < f["frac_over"] and hi < 0), I4_dv_le_2x_fixed=bool(g["dv_mean"] <= ST["dv_factor_mars"] * f["dv_mean"]),
                 I5_bad_runs_le_fixed=bool(g["bad_runs"] <= f["bad_runs"]))
        if a == "F": c["I1_collapse_lt_10pct"] = bool(out["collapse_fraction"] < ST["collapse_fraction_max"])
        c["all"] = all(c.values())
    else:
        loss = bool(lo > 0 or g["frac_hard"] > f["frac_hard"])
        c = dict(J1_no_loss=bool(not loss), J2_dv_le_1p1x_fixed=bool(g["dv_mean"] <= (1 + ST["dv_tolerance_earth_venus"]) * f["dv_mean"])); c["all"] = all(c.values())
    return dict(**out, criteria=c)


def phase_held(root, planets):
    rng = np.random.default_rng(AN["bootstrap_seed"]); res = {}
    for pl in planets:
        if not (sfdir(root) / pl).exists(): continue
        R = cell_table(root, pl)
        if len(R) == 0: continue
        (sfdir(root) / "analysis").mkdir(parents=True, exist_ok=True); R.to_csv(sfdir(root) / "analysis" / f"i_cells_{pl}.csv", index=False)
        S = {k: arm_summary(R, k) for k in ARMS}; cmp_ = {a: compare(R, S, a, rng, pl) for a in ("F", "G", "H")}
        if pl == "mars":
            if cmp_["F"]["criteria"]["all"]: verdict = "TWO_SIDED_CONTROLLABLE"
            elif cmp_["G"]["criteria"]["all"] or cmp_["H"]["criteria"]["all"]: verdict = "ASYMMETRIC_ONLY(" + ",".join(a for a in "GH" if cmp_[a]["criteria"]["all"]) + ")"
            else: verdict = "NO_ADAPTIVE_ARM_PASSES"
        else:
            verdict = {a: ("GENERALIZES" if cmp_[a]["criteria"]["all"] else "DEGRADES") for a in "FGH"}
        res[pl] = dict(n_cells=len(R), arms=S, vs_fixed=cmp_, verdict=verdict)
    (sfdir(root) / "analysis").mkdir(parents=True, exist_ok=True)
    (sfdir(root) / "analysis" / "i_results.json").write_text(json.dumps(res, indent=1, default=str)); print(json.dumps(res, indent=1, default=str))


def phase_params(root, planets, ref_root):
    for pl in planets:
        P = sf.load_runs(root, pl, "fixed", "pert")
        if len(P) == 0: print(pl, "no fixed pert runs"); continue
        D = sf.build(P); D = D[D.split == "dev"]
        if len(D) < 50: print(pl, "too few development scores", len(D)); continue
        sig = {p: float(sf.rms(D["e_" + p])) for p in sf.CAND}; best = min(sf.CAND, key=lambda p: sig[p]); e = D["e_" + best]
        r = dict(planet=pl, predictor=best, e_hi=float(np.quantile(e, AN["adapt_coverage"])), e_lo=float(np.quantile(e, 1 - AN["adapt_coverage"])), n_dev_scored=len(D), n_dev_cells=D.tag.nunique(), **{f"sigma_{p}": sig[p] for p in sf.CAND})
        pd.DataFrame([r]).to_csv(sfdir(root) / pl / "adapt_params.csv", index=False); print("STAGE I   ", r)
        old = pathlib.Path(ref_root) / "stagef20" / pl / "adapt_params.csv" if ref_root else None
        if old is not None and old.exists(): print("STAGE F20 ", pd.read_csv(old).iloc[0].to_dict())
        print(f"band-ln-width {np.log(Q_HIGH / Q_LOW):.3f}, spread q0.95-q0.05 {r['e_hi'] - r['e_lo']:.3f}, spread/band {(r['e_hi'] - r['e_lo']) / np.log(Q_HIGH / Q_LOW):.3f}")


def phase_integrity(root, planets):
    out = {}; cfg_epochs = {pl: sf.CFG["planets"][pl]["grid"]["epochs_utc"] for pl in planets}
    for f in sorted(glob.glob(str(sfdir(root) / "mc_*.toml"))):
        t = tomllib.load(open(f, "rb")); out[pathlib.Path(f).name] = dict(route=t["mc_route"], failed=t["mc_failed"], samples=t["mc_samples"])
    for pl in planets:
        d = sfdir(root) / pl / "cal"
        if not d.exists(): continue
        G = pd.concat([pd.read_csv(g) for g in sorted(d.glob("*/guess.csv"))], ignore_index=True)
        ep = [pd.Timestamp(x).strftime("%Y-%-m-%-d %-H:%-M:%-S") for x in cfg_epochs[pl]]
        out[pl] = dict(n_cal=len(G), all_ok=bool((G.status == "ok").all()), max_rel_q_err=float((G.q_first_w_cm2 / Q_TGT - 1).abs().max()), max_iter=int(G.n_cal_iter.max()),
                       q_target_col=sorted(set(G.q_target_w_cm2.round(12))), epochs_match=bool(all(any(str(e).startswith(x) for x in ep) for e in G.gram_core_epoch)),
                       seed_alt_m_range=[float(G.seed_alt_m.min()), float(G.seed_alt_m.max())], first_pass_vrel_range=[float(G.vrel_first_mps.min()), float(G.vrel_first_mps.max())])
        nan = 0; nruns = 0
        for p in (sfdir(root) / pl / "runs").glob("*/passes.csv"):
            T = pd.read_csv(p); nruns += 1; nan += int(T.drop(columns=["h_qhigh_m"], errors="ignore").isna().sum().sum())
        out[pl]["runs"] = nruns; out[pl]["nan_cells_excl_h_qhigh"] = nan
    (sfdir(root) / "analysis").mkdir(parents=True, exist_ok=True); (sfdir(root) / "analysis" / "i_integrity.json").write_text(json.dumps(out, indent=1, default=str)); print(json.dumps(out, indent=1, default=str))


def selftest():
    q = [0.04, 0.05, 0.2, 0.25, 0.26, 0.30, 0.31, 0.45, 0.46]; t = tiers(q)
    assert t == dict(low=1, soft=2, medium=2, hard=1), t
    assert abs(np.log(Q_HIGH / Q_LOW) - 1.6094379124341003) < 1e-12 and abs(Q_TGT - 0.15) < 1e-15
    # criteria on an exact two-cell case: arm over-fraction 0 vs fixed 0.5 in both cells, identical dv and no bad runs -> Mars I2..I5 pass, F also needs I1
    n = 10; R = pd.DataFrame({"tag": ["a", "b"]})
    for a, over in (("fixed", 5), ("F", 0)):
        for k, v in dict(n_pass=n, low=0, soft=over, medium=0, hard=0, over=over, dv=0.1, burns=1, decay=100.0, retcode="Success", t_end_days=30.0, stopped=False, collapsed=0, n_dec=20).items(): R[f"{a}_{k}"] = v
    S = {k: arm_summary(R, k) for k in ("fixed", "F")}; c = compare(R, S, "F", np.random.default_rng(1), "mars")["criteria"]
    assert c["I1_collapse_lt_10pct"] and c["I2_hard_le_fixed"] and c["I3_over_lt_fixed_ci_excl0"] and c["I4_dv_le_2x_fixed"] and c["I5_bad_runs_le_fixed"] and c["all"]   # two identical differences: the CI is a point mass at -0.5 < 0
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--phase", required=True); ap.add_argument("--results-root"); ap.add_argument("--ref-root"); ap.add_argument("--planets", default="mars,earth,venus"); a = ap.parse_args()
    pls = a.planets.split(",")
    if a.phase == "selftest": return selftest()
    {"held": lambda: phase_held(a.results_root, pls), "params": lambda: phase_params(a.results_root, pls, a.ref_root), "integrity": lambda: phase_integrity(a.results_root, pls)}[a.phase]()


if __name__ == "__main__": main()
