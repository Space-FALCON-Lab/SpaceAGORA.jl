"""Stage N-flight analysis: closed-loop flight of the chance-constrained edge rule. Plan: tasks/2026-10-15_periapsis-region-border-targeting/stagen_flight_plan.md (frozen before any flight).
Reuses stagen_analysis.py (frozen Kalman fit and filter), stagel_analysis.py (run metrics, pooling, Stage I / Stage L tables) and stagei_analysis.py unchanged. Stage I and Stage L trees are READ from <root>/stagei and <root>/stagel.
  selftest                       exact-case checks (arm code, selection rule, Clopper-Pearson, innovation join, kalman agreement with stagen_analysis)
  params   --results-root R [--n-results f.json]   filter constants and quantile multipliers from the Stage I dev fixed_pert tables -> <R>/stagen_flight/mars/kalman_params.csv (asserted equal to the offline record)
  select   --results-root R      dev half-cell selection of the band-width cap -> <R>/stagen_flight/mars/selected.csv, analysis/n_dev.json
  held     --results-root R      held-out criteria C1 to C5, comparators, closed-loop versus offline calibration -> analysis/n_held.json
  repro    --results-root R      bitwise comparison of the repro flights with the held runs (determinism) and, with --xroot, of re-flown Stage I arms with the stored Stage I tables
  integrity --results-root R     route, failures, GRAM epoch, NaN scan
Usage: python3 -I stagen_flight_analysis.py --phase ... --results-root DIR"""
import argparse, glob, json, pathlib, sys, tomllib
import numpy as np, pandas as pd
from scipy.stats import beta
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stagen_analysis as sn
import stagel_analysis as sla
import stagei_analysis as sia
import stagej_analysis as sj
sf = sia.sf
CN = tomllib.load(open(HERE / "stagen_flight_config.toml", "rb")); SP = CN["stagen"]
EPS = SP["epsilon"]; WGRID = SP["band_width_grid_m"]; DVF, DECF = SP["dv_factor"], SP["decay_fraction_min"]; OVER_FIXED = SP["over_cmp_fixed"]; CEIL = SP["dev_over_ceiling"]
PL = "mars"; QMAX, QLOW, QHARD = 0.25, 0.05, 0.45; WB, WA = sj.WB, sj.WA
LMEMBERS = ["fm04000", "fam_p080_sp000", "fam_p090_sp050", "fam_p095_sm100"]


def ndir(root): return pathlib.Path(root) / CN["stagef"]["results_subdir"]
def arm_code(eps, w): return f"cc{round(100 * eps):02d}_w" + ("none" if w <= 0 else f"{round(w):05d}")   # must equal stagen_arm_code in stagen_flight_common.jl


def n_summaries(root):
    fs = sorted(glob.glob(str(ndir(root) / "summary_*.csv")))
    if not fs: return pd.DataFrame()
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True).drop_duplicates("run_id", keep="last")


def run_n(root, tag, code, S, sub="runs"):
    d = ndir(root) / PL / sub / f"{code}_pert_{tag}"; sr = S[S.run_id == f"{code}_pert_{tag}"] if len(S) else []
    if not (d / "passes.csv").exists() or len(sr) == 0: return None
    dec = pd.read_csv(d / "decisions.csv") if (d / "decisions.csv").exists() else pd.DataFrame()
    return sia.run_metrics(pd.read_csv(d / "passes.csv"), sr.iloc[-1].to_dict(), dec)


def cp(k, n, a=0.10):
    lo = beta.ppf(a / 2, k, n - k + 1) if k > 0 else 0.0; hi = beta.ppf(1 - a / 2, k + 1, n - k) if k < n else 1.0; return [float(lo), float(hi)]


def mean_q(root, tags, ctrl_dir, ctrl):
    """Mean q_eff over ALL passes of the given cells (same definition as stagen_analysis.flown)."""
    q = [pd.read_csv(pathlib.Path(ctrl_dir) / PL / "runs" / f"{ctrl}_pert_{t}" / "passes.csv").q_eff_w_cm2.values for t in tags]; return float(np.concatenate(q).mean())


# ---- params -------------------------------------------------------------------------------------------------------------------------------
def phase_params(root, n_results=None):
    P = sj.load(root, "stagei", PL, "fixed_pert"); P = P[P.inw].reset_index(drop=True); D = P[P.split == "dev"].reset_index(drop=True)
    kp = sn.fit_kalman(D); _, S, E, _ = sn.kalman(D, kp["r"], kp["q"], kp["m0"], kp["v0"]); zd = E / np.sqrt(S)
    row = dict(planet=PL, r=kp["r"], q=kp["q"], m0=kp["m0"], v0=kp["v0"], rho=kp["rho"], n_dev_passes=len(D), n_dev_cells=int(D.tag.nunique()), **{f"t_emp_{round(100 * e):02d}": float(np.quantile(zd, 1 - e)) for e in (0.001, 0.01, 0.05)})
    if n_results:
        R = json.load(open(n_results))["I"][PL]; k0 = R["kalman_params"]; t01 = [r for r in R["rules"]["empirical"] if r["eps"] == 0.01][0]["t"]
        for k in ("r", "m0", "v0", "rho", "q"): assert abs(k0[k] - kp[k]) < 1e-12 * max(1, abs(k0[k])), (k, k0[k], kp[k])
        assert abs(t01 - row["t_emp_01"]) < 1e-9, (t01, row["t_emp_01"]); print("params equal the offline record (stagen_data/n_results.json): r, m0, v0, rho, q and t(1%)")
    (ndir(root) / PL).mkdir(parents=True, exist_ok=True); pd.DataFrame([row]).to_csv(ndir(root) / PL / "kalman_params.csv", index=False); print("KALMAN", row)


# ---- dev selection ------------------------------------------------------------------------------------------------------------------------
def select_w(T):
    """Pre-registered dev selection (plan section 4). T: one row per band-width cap with w_m, frac_over, dv_mean, decay_median, dv_ref (fixed dev dV). Feasible = frac_over <= CEIL and dv_mean <= DVF x dv_ref.
    Pick the feasible row with the largest decay (ties: larger w with 0 meaning none = infinite width). None feasible (fallback): among rows with dv <= DVF x ref the lowest frac_over, else the lowest frac_over overall."""
    T = T.copy(); T["w_key"] = np.where(T.w_m <= 0, np.inf, T.w_m); F = T[(T.frac_over <= CEIL) & (T.dv_mean <= DVF * T.dv_ref)]
    if len(F): return F.sort_values(["decay_median", "w_key"], ascending=[False, False]).iloc[0].w_m, False
    G = T[T.dv_mean <= DVF * T.dv_ref]; G = G if len(G) else T
    return G.sort_values(["frac_over", "w_key"], ascending=[True, False]).iloc[0].w_m, True


def phase_select(root):
    SN, SL = n_summaries(root), sla.l_summaries(root); tags = sla.dev_cells(root, PL); rows = []; missing = {}
    ref = [sla.run_of(root, PL, t, "fm01000", SL, "L") for t in tags]
    if any(m is None for m in ref): print("Stage L dev fixed 1000 reference missing for some half-cells"); return
    fx = sla.pool(ref)
    for w in WGRID:
        code = arm_code(EPS, w); ms = [run_n(root, t, code, SN) for t in tags]
        if any(m is None for m in ms): missing[code] = sum(m is None for m in ms); continue
        rows.append(dict(code=code, w_m=w, **sla.pool(ms), dv_ref=fx["dv_mean"]))
    T = pd.DataFrame(rows)
    if len(T) == 0: print("no complete dev arms", missing); return
    w, fb = select_w(T); sel = pd.DataFrame([dict(code=arm_code(EPS, w), w_m=float(w), fallback=fb)]); sel.to_csv(ndir(root) / PL / "selected.csv", index=False)
    (ndir(root) / "analysis").mkdir(exist_ok=True); T.to_csv(ndir(root) / "analysis" / "n_dev.csv", index=False)
    out = dict(cells=len(tags), arms_missing=missing, fixed1000_dev={k: fx[k] for k in ("frac_over", "dv_mean", "decay_median", "burns_mean", "hard", "frac_low")}, table=T.to_dict(orient="records"), selected=sel.iloc[0].to_dict())
    (ndir(root) / "analysis" / "n_dev.json").write_text(json.dumps(out, indent=1, default=float)); print(json.dumps(out, indent=1, default=float))


# ---- held-out criteria and calibration ----------------------------------------------------------------------------------------------------
def innovations(root, tag, code):
    """Per decision of one cc run: e = y(target pass) - m (the filter's predictive error), c, s, whether the pass is in the in-window rule, the realized q_eff and alt of that pass.
    meas row k (0-based) is pass k + 1; decision target_pass = periapsis_count + 1, so y of the target pass is meas row target_pass - 1."""
    d = ndir(root) / PL / "runs" / f"{code}_pert_{tag}"; D = pd.read_csv(d / "decisions.csv"); M = pd.read_csv(d / "measured.csv"); P = pd.read_csv(d / "passes.csv")
    hg = float(pd.read_csv(sla.idir(root) / PL / "cal" / tag / "guess.csv").h_guess_m[0]); out = []
    for r in D.itertuples():
        k = int(r.target_pass)
        if k - 1 >= len(M): continue
        pr = P[P["pass"] == k]
        if len(pr) == 0: continue
        pr = pr.iloc[0]; out.append(dict(tag=tag, k=k, e=float(M.y.iloc[k - 1] - r.m), s=float(r.s), c=float(r.c), alt=float(pr.alt_m), inw=bool(hg - WB <= pr.alt_m <= hg + WA), q=float(pr.q_eff_w_cm2), burned=abs(float(r.delta_v_mps)) > 0))
    return pd.DataFrame(out)


def cluster_boot(tags, num, den, rng, n):
    ut = np.array(sorted(set(tags))); idx = {t: np.where(tags == t)[0] for t in ut}; out = []
    for _ in range(n):
        i = np.concatenate([idx[t] for t in rng.choice(ut, len(ut))]); out.append(num[i].sum() / den[i].sum())
    return [float(np.percentile(out, 5)), float(np.percentile(out, 95))]


def calib(I, rng, n):
    """Closed-loop calibration on a set of decided passes I (columns e, s, c, q, tag)."""
    k = int((I.e > I.c).sum()); N = len(I); ko = int((I.q > QMAX).sum()); z = I.e / I.s; tg = I.tag.values
    return dict(n=N, innov_exceed_k=k, innov_exceed=k / N, innov_exceed_ci_binom=cp(k, N), innov_exceed_ci_cluster=cluster_boot(tg, (I.e > I.c).values.astype(float), np.ones(N), rng, n),
                realized_over_k=ko, realized_over=ko / N, realized_over_ci_binom=cp(ko, N), realized_over_ci_cluster=cluster_boot(tg, (I.q > QMAX).values.astype(float), np.ones(N), rng, n),
                z_mean=float(z.mean()), z_std=float(z.std()), z_q99=float(np.quantile(z, 0.99)), e_mean=float(I.e.mean()), e_std=float(I.e.std()), mean_c=float(I.c.mean()), mean_s=float(I.s.mean()),
                innov_exceed_contains_eps=bool(cp(k, N)[0] <= EPS <= cp(k, N)[1]), realized_over_contains_eps=bool(cp(ko, N)[0] <= EPS <= cp(ko, N)[1]))


def phase_held(root, n_results=None):
    rng = np.random.default_rng(SP["bootstrap_seed"]); NB = SP["bootstrap_draws"]; SN, SL, SI = n_summaries(root), sla.l_summaries(root), sf.summaries(root, "held")
    tags = sla.held_cells(root, PL); sel = pd.read_csv(ndir(root) / PL / "selected.csv").iloc[0]; code = str(sel.code)
    arms = {}; cells = {}
    for name, ctrl in sla.STAGE_I_ARMS.items(): cells[name] = [sla.run_of(root, PL, t, ctrl, SI, "I") for t in tags]; arms[name] = (sla.idir(root), ctrl)
    for c in LMEMBERS: cells["L_" + c] = [sla.run_of(root, PL, t, c, SL, "L") for t in tags]; arms["L_" + c] = (sla.ldir(root), c)
    cells["N_cc"] = [run_n(root, t, code, SN) for t in tags]; arms["N_cc"] = (ndir(root), code)
    cells = {k: v for k, v in cells.items() if all(m is not None for m in v)}; assert "N_cc" in cells and "I_fixed1000" in cells, "missing flown cells"
    res = dict(arm=code, selected=sel.to_dict(), n_cells=len(tags)); rows = {}
    for k, ms in cells.items():
        p = sla.pool(ms); p["mean_Q"] = mean_q(root, tags, *arms[k]); rows[k] = p
    fx = rows["I_fixed1000"]
    for k, p in rows.items(): p["dv_ratio"] = p["dv_mean"] / fx["dv_mean"]; p["decay_ratio"] = p["decay_median"] / fx["decay_median"]; p["meanQ_ratio"] = p["mean_Q"] / fx["mean_Q"]
    n = rows["N_cc"]; ko, N = n["over"], n["passes"]; ci_cp = cp(ko, N)
    over_cell = np.array([m["over"] for m in cells["N_cc"]], float); npass = np.array([m["n_pass"] for m in cells["N_cc"]], float)
    ci_cl = cluster_boot(np.array(tags), over_cell, npass, rng, NB)
    d = np.array([m["over"] / m["n_pass"] for m in cells["N_cc"]]) - np.array([m["over"] / m["n_pass"] for m in cells["I_fixed1000"]]); dlo, dhi = sf.boot_ci(d, rng, NB)
    crit = dict(
        C1a_calibration=dict(over_k=ko, n=N, over=ko / N, ci90_binom=ci_cp, ci90_cluster=ci_cl, epsilon=EPS, contains_eps_binom=bool(ci_cp[0] <= EPS <= ci_cp[1]), contains_eps_cluster=bool(ci_cl[0] <= EPS <= ci_cl[1]),
                             direction=("conservative (interval below eps)" if ci_cp[1] < EPS else "liberal (interval above eps)" if ci_cp[0] > EPS else "calibrated (interval contains eps)")),
        C1b_not_worse_than_fixed=dict(over=ko / N, fixed_over=OVER_FIXED, fixed_over_recomputed=fx["frac_over"], passes=bool(ko / N <= OVER_FIXED), paired_diff_mean=float(d.mean()), paired_diff_ci90=[dlo, dhi]),
        C2_hard_zero=dict(hard=n["hard"], passes=bool(n["hard"] == 0)),
        C3_dv=dict(dv_mean=n["dv_mean"], fixed=fx["dv_mean"], ratio=n["dv_ratio"], limit=DVF, passes=bool(n["dv_mean"] <= DVF * fx["dv_mean"])),
        C4_decay=dict(decay_median=n["decay_median"], fixed=fx["decay_median"], ratio=n["decay_ratio"], limit=DECF, passes=bool(n["decay_median"] >= DECF * fx["decay_median"])),
        C5_diagnostics=dict(meanQ_ratio=n["meanQ_ratio"], mean_Q=n["mean_Q"], fixed_mean_Q=fx["mean_Q"], frac_low=n["frac_low"], fixed_frac_low=fx["frac_low"], low_pass_pass=bool(n["frac_low"] < QLOW)))
    crit["C1_pass"] = bool(crit["C1a_calibration"]["contains_eps_binom"] and crit["C1b_not_worse_than_fixed"]["passes"])
    crit["verdict"] = dict(C1=crit["C1_pass"], C2=crit["C2_hard_zero"]["passes"], C3=crit["C3_dv"]["passes"], C4=crit["C4_decay"]["passes"]); res["criteria"] = crit
    res["table"] = {k: {kk: p[kk] for kk in ("passes", "frac_over", "over", "hard", "frac_low", "dv_mean", "dv_ratio", "burns_mean", "decay_median", "decay_ratio", "mean_Q", "meanQ_ratio", "bad_runs", "collapsed", "decisions")} for k, p in rows.items()}
    # closed-loop versus offline calibration (this is the question): innovations of the flown filter on all decided passes and on the in-window subset
    I = pd.concat([innovations(root, t, code) for t in tags], ignore_index=True); I["tag"] = I.tag.astype(str)
    cal = dict(all_decided=calib(I, rng, NB), in_window=calib(I[I.inw].reset_index(drop=True), rng, NB), burned_decisions=int(I.burned.sum()), decided=len(I))
    if n_results:
        R = json.load(open(n_results))["I"][PL]; em = [r for r in R["rules"]["empirical"] if r["eps"] == EPS][0]
        cal["offline_replay"] = dict(over=em["over"], over_k=em["over_k"], n=em["n"], ci90_binom=em["over_ci_binom"], ci90_cluster=em["over_ci_cluster"], hard=em["hard"], low=em["low"], mean_Qaim=em["mean_Qaim"], mean_Qreal=em["mean_Qreal"],
                                     standardized_innovation_lag1=R["rules"]["kalman_standardized_innov_lag1"])
        cal["offline_fixed_arm_inwindow_check"] = R["fixed_arm_replay_check"]
    # planned-versus-realized heat rate at the periapsis event: realized q against the planned edge heat rate Qmax exp(-c) on decided passes (execution + orbit-prediction + tracking error, in ln units)
    I["lnq_minus_plan"] = np.log(I.q / (QMAX * np.exp(-I.c))); res["tracking"] = dict(lnQ_real_over_Qedge_median=float(I.lnq_minus_plan.median()), lnQ_real_over_Qedge_std=float(I.lnq_minus_plan.std()),
        mean_Q_edge=float((QMAX * np.exp(-I.c)).mean()), frac_decided_passes_above_edge=float((I.q > QMAX * np.exp(-I.c)).mean()))
    (ndir(root) / "analysis").mkdir(exist_ok=True); res["calibration"] = cal; I.to_csv(ndir(root) / "analysis" / "n_innovations.csv", index=False)
    (ndir(root) / "analysis").mkdir(exist_ok=True); (ndir(root) / "analysis" / "n_held.json").write_text(json.dumps(res, indent=1, default=float)); print(json.dumps(res, indent=1, default=float))


def bit_compare(a_dir, b_dir, cols=("alt_m", "apo_alt_m", "q_eff_w_cm2", "mass_kg")):
    a, b = pd.read_csv(a_dir / "passes.csv"), pd.read_csv(b_dir / "passes.csv")
    if len(a) != len(b): return dict(identical=False, n=(len(a), len(b)))
    dmax = float(max(np.nanmax(np.abs(a[c].values - b[c].values)) for c in cols)); return dict(identical=bool(dmax == 0.0), max_abs_diff=dmax, passes=len(a))


def phase_repro(root, xroot=None):
    out = dict(determinism={}, stagei_reproduction={}); sel = pd.read_csv(ndir(root) / PL / "selected.csv").iloc[0]; code = str(sel.code)
    d = ndir(root) / PL / "repro"
    for rd in sorted(d.glob(f"{code}_pert_*")): tag = rd.name.split("_pert_")[1]; out["determinism"][tag] = bit_compare(rd, ndir(root) / PL / "runs" / rd.name)
    if True:   # <root>/stagen_flight/mars/xcheck holds fm01000, fam_p095_sp000, fam_p095_sm100 flights of one held-out cell (stagen_flight_xcheck.jl), re-flown on the current code; compare with the stored Stage I tables
        for rd in sorted((ndir(root) / PL / "xcheck").glob("*_pert_*")):
            c = rd.name.split("_pert_")[0]; tag = rd.name.split("_pert_")[1]; ci = {"fm01000": "fixed", "fam_p095_sp000": "adaptive_asym", "fam_p095_sm100": "adaptive"}[c]
            out["stagei_reproduction"][c] = dict(tag=tag, **bit_compare(rd, sla.idir(root) / PL / "runs" / f"{ci}_pert_{tag}"))
    (ndir(root) / "analysis").mkdir(exist_ok=True); (ndir(root) / "analysis" / "n_repro.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out, indent=1))


def phase_integrity(root):
    out = {}
    for f in sorted(glob.glob(str(ndir(root) / "mc_*.toml"))):
        t = tomllib.load(open(f, "rb")); out[pathlib.Path(f).name] = dict(route=t["mc_route"], failed=t["mc_failed"], samples=t["mc_samples"], host=t["host"])
    S = n_summaries(root); runs = list((ndir(root) / PL / "runs").glob("*/passes.csv")); nan = sum(int(pd.read_csv(p).drop(columns=["h_qhigh_m"], errors="ignore").isna().sum().sum()) for p in runs)
    ep = [pd.Timestamp(x).strftime("%Y-%-m-%-d %-H:%-M:%-S") for x in sf.CFG["planets"][PL]["grid"]["epochs_utc"]]
    out[PL] = dict(runs=len(runs), nan_cells_excl_h_qhigh=nan, rows=len(S), not_ok=int((S.status != "ok").sum()), epochs_match=bool(all(any(str(e).startswith(x) for x in ep) for e in S.gram_core_epoch)),
                   controller_stopped=int(S.controller_stopped.sum()), not_success=int((S.retcode != "Success").sum()))
    (ndir(root) / "analysis").mkdir(exist_ok=True); (ndir(root) / "analysis" / "n_integrity.json").write_text(json.dumps(out, indent=1, default=str)); print(json.dumps(out, indent=1, default=str))


def selftest():
    assert arm_code(0.01, 0.0) == "cc01_wnone" and arm_code(0.01, 2000.0) == "cc01_w02000"
    assert cp(5, 100)[0] < 0.05 < cp(5, 100)[1] and cp(0, 100)[0] == 0.0
    T = pd.DataFrame(dict(w_m=[0.0, 4000.0, 2000.0], frac_over=[0.005, 0.01, 0.03], dv_mean=[0.5, 0.6, 0.9], decay_median=[100.0, 110.0, 150.0], dv_ref=0.4))
    w, fb = select_w(T); assert (w, fb) == (4000.0, False), (w, fb)                      # 2000 m infeasible (over 0.03 > 0.02): best decay among feasible = 4000
    T2 = T.assign(frac_over=[0.05, 0.04, 0.03], dv_mean=[0.5, 0.9, 0.7]); w, fb = select_w(T2); assert fb and w == 2000.0, (w, fb)   # none feasible: lowest over among dV <= 2x ref (0.8): 4000 excluded (0.9), 2000 chosen
    T3 = T.assign(decay_median=[100.0, 100.0, 100.0], frac_over=0.0, dv_mean=0.5); w, fb = select_w(T3); assert w == 0.0 and not fb, w       # tie in decay: larger width (none = infinite)
    # kalman agreement: the stagen_analysis filter on a synthetic cell reproduces the hand recursion used in the Julia selftest (r 0.1, q 0.01, prior 0, 0.2; y = 1, 0)
    G = pd.DataFrame(dict(tag="a", **{"pass": [1, 2]}, y=[1.0, 0.0])); M, S, E, _ = sn.kalman(G, 0.1, 0.01, 0.0, 0.2)
    assert abs(M[0]) < 1e-15 and abs(S[0] - 0.3) < 1e-15 and abs(M[1] - 2 / 3) < 1e-15 and abs(S[1] - (1 / 15 + 0.01 + 0.1)) < 1e-15, (M, S)
    I = pd.DataFrame(dict(tag=["a"] * 50 + ["b"] * 50, e=np.r_[np.zeros(98), 1.0, 1.0], s=1.0, c=0.5, q=np.r_[np.full(99, 0.1), 0.3])); r = calib(I, np.random.default_rng(1), 200)
    assert r["innov_exceed_k"] == 2 and r["realized_over_k"] == 1 and abs(r["z_mean"] - 0.02) < 1e-12, r
    print("stagen_flight selftest ok")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--phase", required=True); ap.add_argument("--results-root"); ap.add_argument("--n-results"); ap.add_argument("--xroot"); a = ap.parse_args()
    if a.phase == "selftest": return selftest()
    root = a.results_root
    {"params": lambda: phase_params(root, a.n_results), "select": lambda: phase_select(root), "held": lambda: phase_held(root, a.n_results), "repro": lambda: phase_repro(root, a.xroot), "integrity": lambda: phase_integrity(root)}[a.phase]()


if __name__ == "__main__": main()
