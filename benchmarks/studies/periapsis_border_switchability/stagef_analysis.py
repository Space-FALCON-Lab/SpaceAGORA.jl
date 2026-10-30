"""Stage F analysis (plan: tasks/2026-10-15_periapsis-region-border-targeting/stagef_plan.md; written after the plan freeze and before any production data were analyzed, Addendum C).
Phases (python3 -I stagef_analysis.py --phase=... --results-root DIR):
  selftest   exact-case checks of the scoring code
  xhost      --root-a A --root-b B: cross-host comparison of the check runs (plan Addendum A section 3)
  f1         F1 density statistics, tests R0-R7 per planet (writes <root>/stagef/analysis/)
  f2params   development-cell predictor selection and error quantiles -> <root>/stagef/<planet>/adapt_params.csv
  f2         adaptive-versus-fixed comparison on held-out cells (writes <root>/stagef/analysis/)
All thresholds are read from stagef_config.toml [analysis]; they are identical for every planet."""
import argparse, glob, json, pathlib, sys, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import density_buffer_analysis as dba   # predict_run, lag1 (unchanged earlier code)
CFG = tomllib.load(open(HERE / "stagef_config.toml", "rb")); AN = CFG["analysis"]; SF = CFG["stagef"]
BAND = CFG["common"]["band"]; QLO, QHI, QHARD = BAND["q_low_w_cm2"], BAND["q_high_w_cm2"], BAND["q_hard_w_cm2"]; WQ = np.log(QHI / QLO)
LEVELS = AN["coverage_levels"]; KS = (1.0, 1.64, 2.33)
CAND = ["P1", "P2_3", "P2_5", "P2_10", "P3"]; ALLP = ["P0s"] + CAND
rms = dba.rms


def sfdir(root): return pathlib.Path(root) / SF["results_subdir"]


def cells_of(root, pl):
    return pd.read_csv(sfdir(root) / pl / "cells_base.csv") if (sfdir(root) / pl / "cells_base.csv").exists() else pd.read_csv(sfdir(root).parent / "cells_base.csv")


def load_runs(root, pl, ctrl, arm):
    """All in-window pass rows of one (planet, control, arm) with y, z, y_sys; one frame, column tag."""
    out = []; cells = cells_of(root, pl)
    for t in cells.tag:
        f = sfdir(root) / pl / "runs" / f"{ctrl}_{arm}_{t}" / "passes.csv"; g = sfdir(root) / pl / "cal" / t / "guess.csv"
        if not f.exists() or not g.exists(): continue
        P = pd.read_csv(f); hg = pd.read_csv(g).h_guess_m[0]
        if len(P) == 0: continue
        P["tag"] = t; P["inw"] = (P.alt_m >= hg - CFG["common"]["fit"]["window_below_m"]) & (P.alt_m <= hg + CFG["common"]["fit"]["window_above_m"])
        P["y"] = np.log(P.rho_eff / P.rho_pw); P["z"] = np.log(P.r_pert); P["ysys"] = P.y - P.z; P["split"] = cells.set_index("tag").loc[t, "split"]
        out.append(P)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def build(P):
    """Scored one-pass-ahead rows from consecutive in-window passes (predictions from earlier passes only)."""
    recs = []
    for t, T in P.groupby("tag", sort=False):
        T = T[T.inw]; brk = np.r_[True, np.diff(T["pass"].values) != 1]; rid = np.cumsum(brk)
        for r in np.unique(rid):
            R = T[rid == r]
            if len(R) < 2: continue
            pr = dba.predict_run(R.y.values, np.unwrap(R.lon_rad.values))
            for j in range(1, len(R)):
                x = R.iloc[j]; d = dict(tag=t, run=int(r), pass_=int(x["pass"]), j=j, y=x.y, z=x.z, ysys=x.ysys, H=x.h_local_m, lon=x.lon_rad, alt_m=x.alt_m, split=x.split)
                d["zprev"] = R.iloc[j - 1].z
                for p in ALLP: d["e_" + p] = x.y - pr[p][j]
                recs.append(d)
    return pd.DataFrame(recs)


def qualify(D, n=None):
    n = n or AN["min_scored_passes_per_cell"]; c = D.groupby("tag").size(); return set(c[c >= n].index)


def loco(D, col, fn):
    return {t: fn(D.loc[D.tag != t, col].values) for t in D.tag.unique()}


def score(D):
    """Per predictor: pooled sigma, LOCO sigma, empirical-quantile and Gaussian coverage, buffers in ln (= scale heights) and km."""
    rows = []; cells = list(D.tag.unique()); Hm = D.H.median()
    for p in ALLP:
        e = D["e_" + p]; sg = loco(D, "e_" + p, lambda v: np.sqrt(np.mean(v ** 2)))
        for lv in LEVELS:
            qs = loco(D, "e_" + p, lambda v, lv=lv: np.quantile(v, lv)); qv = D.tag.map(qs).values
            cov = pd.Series({t: float(np.mean(g["e_" + p].values <= qs[t])) for t, g in D.groupby("tag")})
            rows.append(dict(pred=p, kind="emp_quantile", level=lv, n=len(D), sigma=float(rms(e)), coverage=float(np.mean(e.values <= qv)), buf_ln=float(np.mean(qv)),
                             buf_km=float(np.mean(D.H.values * qv) / 1e3), percell_within=float(np.mean(np.abs(cov - lv) <= AN["percell_tolerance"])), n_cells=len(cov)))
        for k in KS:
            b = k * D.tag.map(sg).values
            rows.append(dict(pred=p, kind="gaussian_k", level=k, n=len(D), sigma=float(rms(e)), coverage=float(np.mean(e.values <= b)), buf_ln=float(np.mean(b)),
                             buf_km=float(np.mean(D.H.values * b) / 1e3), percell_within=np.nan, n_cells=len(cells)))
    return pd.DataFrame(rows)


def lag1_col(D, col):
    a, b = [], []
    for (t, r), G in D.groupby(["tag", "run"]):
        if len(G) < dba.MIN_RUN_ERR: continue
        v = G[col].values - D.loc[D.tag == t, col].mean(); a += list(v[:-1]); b += list(v[1:])
    return float(np.corrcoef(a, b)[0, 1]) if len(a) > 2 else np.nan


def f1_set(root, pl, ctrl):
    """Scored frames of the mean and perturbed arms of one control set, restricted to cells that qualify in each arm."""
    res = {}
    for arm in ("mean", "pert"):
        P = load_runs(root, pl, ctrl, arm); res[arm] = build(P) if len(P) else pd.DataFrame()
    return res


def analyze_set(root, pl, ctrl):
    S = f1_set(root, pl, ctrl); out = dict(planet=pl, ctrl=ctrl)
    if any(len(S[a]) == 0 for a in S): out["note"] = "missing arm"; return out, None, None
    qm, qp = qualify(S["mean"]), qualify(S["pert"]); out["cells_qualifying_mean"] = len(qm); out["cells_qualifying_pert"] = len(qp)
    out["R0"] = bool(len(qm) >= AN["min_cells"] and len(qp) >= AN["min_cells"])
    Dp = S["pert"][S["pert"].tag.isin(qp)].copy(); Dm = S["mean"][S["mean"].tag.isin(qm)].copy()
    if len(Dp) == 0 or len(Dm) == 0: out["note"] = "no qualifying cells"; return out, None, None
    sp, sm = score(Dp), score(Dm); out["n_scored_pert"] = len(Dp); out["n_scored_mean"] = len(Dm)
    sig = {p: float(rms(Dp["e_" + p])) for p in ALLP}; best = min(CAND, key=lambda p: sig[p]); out["best_pert"] = best; out["sigma_pert"] = sig
    out["sigma_mean"] = {p: float(rms(Dm["e_" + p])) for p in ALLP}; out["best_mean"] = min(CAND, key=lambda p: out["sigma_mean"][p])
    q = lambda D, p, lv: float(np.quantile(D["e_" + p], lv))
    eq = sp[(sp.pred == best) & (sp.kind == "emp_quantile")].set_index("level")
    out["R1"] = dict(coverage={str(l): float(eq.loc[l, "coverage"]) for l in LEVELS}, passed=bool(all(eq.loc[l, "coverage"] >= l - AN["coverage_tolerance_below"] for l in LEVELS)))
    out["R2"] = dict(fraction=float(eq.loc[0.95, "percell_within"]), n_cells=int(eq.loc[0.95, "n_cells"]), passed=bool(eq.loc[0.95, "percell_within"] >= AN["percell_fraction_required"]))
    spread = q(Dp, best, 0.95) - q(Dp, best, 0.05); out["R3"] = dict(spread_ln=spread, band_wq=float(WQ), passed=bool(spread <= WQ))
    out["R4"] = dict(ratio=sig[best] / sig["P0s"], passed=bool(sig[best] / sig["P0s"] <= AN["information_ratio_max"]))
    z = Dp.z.values; rz = lag1_col(Dp, "z"); varz = float(np.var(z)); vm = float(np.var(Dm["e_P1"]))   # mean-arm persistence error variance
    dom = bool(2 * varz * (1 - rz) >= vm); ewma = min(sig[f"P2_{n}"] for n in (3, 5, 10)); ratio = ewma / sig["P1"]
    out["R5"] = dict(random_dominated=dom, var_z=varz, lag1_z=rz, var_eP1_mean_arm=vm, ewma_over_persistence=float(ratio), ewma_beats=bool(ratio <= 1 - AN["ewma_gain_threshold"]),
                     passed=bool((ratio <= 1 - AN["ewma_gain_threshold"]) == dom))
    out["R6"] = dict(lag1_z=rz, lag1_e_P1_pert=lag1_col(Dp, "e_P1"), lag1_e_P1_mean=lag1_col(Dm, "e_P1"), lag1_y_pert=lag1_col(Dp, "y"), passed=bool(abs(rz) < AN["lag1_white_threshold"]))
    out["q95_best_ln"] = q(Dp, best, 0.95); out["q95_best_km"] = float(np.mean(Dp.H * q(Dp, best, 0.95)) / 1e3); out["median_H_km"] = float(Dp.H.median() / 1e3)
    out["decomposition"] = dict(var_y=float(np.var(Dp.y)), var_ysys=float(np.var(Dp.ysys)), var_z=varz, cov_ysys_z=float(np.cov(Dp.ysys, Dp.z)[0, 1]), std_z=float(np.std(z)),
                                sigma_P1_pert=sig["P1"], sigma_P1_mean=out["sigma_mean"]["P1"], sigma_best_mean=out["sigma_mean"][out["best_mean"]])
    pc = Dp.merge(Dm[["tag", "pass_", "e_P1"]], on=["tag", "pass_"], suffixes=("", "_mean"))   # same cell and pass index in both arms (trajectories diverge, so paired by pass index only approximately)
    out["paired_by_cell"] = [dict(tag=t, sigma_P1_pert=float(rms(g.e_P1)), sigma_P1_mean=float(rms(Dm.loc[Dm.tag == t, "e_P1"])), std_z=float(g.z.std())) for t, g in Dp.groupby("tag")]
    del pc
    return out, sp.assign(planet=pl, ctrl=ctrl, arm="pert"), sm.assign(planet=pl, ctrl=ctrl, arm="mean")


def phase_f1(root):
    A = sfdir(root) / "analysis"; A.mkdir(exist_ok=True); res = {}; scores = []
    for pl in SF["planets"]:
        if not (sfdir(root) / pl).exists(): continue
        for ctrl in ("open", "fixed"):
            o, sp, sm = analyze_set(root, pl, ctrl); res[f"{pl}_{ctrl}"] = o
            if sp is not None: scores += [sp, sm]
    # primary set per planet: open loop if R0 holds for it, else the keeper
    prim = {}
    for pl in SF["planets"]:
        o, f = res.get(f"{pl}_open", {}), res.get(f"{pl}_fixed", {})
        if not o and not f: continue
        prim[pl] = "open" if o.get("R0") else "fixed"
    gen = {}
    for arm_key in ("pert",):
        qs = {pl: res[f"{pl}_{prim[pl]}"].get("q95_best_ln") for pl in prim if res.get(f"{pl}_{prim[pl]}", {}).get("q95_best_ln") is not None}
        if len(qs) >= 2: gen = dict(q95_ln=qs, ratio_max_over_min=float(max(qs.values()) / min(qs.values())), passed=bool(max(qs.values()) / min(qs.values()) <= AN["generality_factor"]))
    res["primary_set"] = prim; res["R7"] = gen
    (A / "f1_results.json").write_text(json.dumps(res, indent=1, default=str))
    if scores: pd.concat(scores).to_csv(A / "f1_scores.csv", index=False)
    print(json.dumps({k: v for k, v in res.items()}, indent=1, default=str)[:6000])


def phase_f2params(root):
    rows = []
    for pl in SF["planets"]:
        P = load_runs(root, pl, "fixed", "pert")
        if len(P) == 0: continue
        D = build(P); D = D[D.split == "dev"]
        if len(D) < 50: print(pl, "too few development scores", len(D)); continue
        sig = {p: float(rms(D["e_" + p])) for p in CAND}; best = min(CAND, key=lambda p: sig[p])
        e = D["e_" + best]; hi, lo = float(np.quantile(e, AN["adapt_coverage"])), float(np.quantile(e, 1 - AN["adapt_coverage"]))
        r = dict(planet=pl, predictor=best, e_hi=hi, e_lo=lo, n_dev_scored=len(D), n_dev_cells=D.tag.nunique(), **{f"sigma_{p}": sig[p] for p in CAND})
        rows.append(r); pd.DataFrame([r]).to_csv(sfdir(root) / pl / "adapt_params.csv", index=False); print(r)
    return rows


def run_metrics(P, summ_row):
    if len(P) == 0: return None
    q = P.q_eff_w_cm2.values; t = P.time_s.values; apo = P.apo_alt_m.values; ok = np.isfinite(apo)
    days = (t[-1] - t[0]) / 86400.0 if len(t) > 1 else np.nan
    decay = (apo[ok][0] - apo[ok][-1]) / days if ok.sum() > 1 and days > 0 else np.nan
    return dict(n_pass=len(P), low=int((q < QLO).sum()), soft=int((q > QHI).sum()), hard=int((q > QHARD).sum()), dv=summ_row.get("total_dv_mps", np.nan), burns=summ_row.get("burns", np.nan),
                decay_m_per_day=decay, retcode=summ_row.get("retcode"), t_end_days=summ_row.get("t_end_days"), stopped=bool(summ_row.get("controller_stopped", False)),
                collapsed=summ_row.get("collapsed_decisions", 0), horizon=summ_row.get("horizon_days"))


def summaries(root, phase):
    fs = glob.glob(str(sfdir(root) / f"summary_{phase}_*.csv")); return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True) if fs else pd.DataFrame()


def boot_ci(d, rng, n, lo=5, hi=95):
    d = np.asarray(d, float); d = d[np.isfinite(d)]
    if len(d) < 2: return (np.nan, np.nan)
    m = [rng.choice(d, len(d)).mean() for _ in range(n)]; return float(np.percentile(m, lo)), float(np.percentile(m, hi))


def phase_f2(root, variant="adaptive"):
    S1, S2 = summaries(root, "f1"), summaries(root, "f2" if variant == "adaptive" else "f2track"); rng = np.random.default_rng(AN["bootstrap_seed"]); res = {}; X = AN["win_cost_tolerance"]
    for pl in SF["planets"]:
        cells = cells_of(root, pl) if (sfdir(root) / pl).exists() else None
        if cells is None: continue
        rows = []
        for t in cells[cells.split == "heldout"].tag:
            m = {}
            for ctrl, arm, S in (("fixed", "pert", S1), (variant, "pert", S2)):
                f = sfdir(root) / pl / "runs" / f"{ctrl}_{arm}_{t}" / "passes.csv"
                if not f.exists() or len(S) == 0: m = None; break
                sr = S[(S.planet == pl) & (S.run_id == f"{ctrl}_{arm}_{t}")]
                if len(sr) == 0: m = None; break
                m[ctrl] = run_metrics(pd.read_csv(f), sr.iloc[0].to_dict())
            if m and m["fixed"] and m[variant]: rows.append(dict(tag=t, **{f"{c}_{k}": v for c, d in m.items() for k, v in d.items()}))
        if not rows: continue
        R = pd.DataFrame(rows); R = R.rename(columns={c: c.replace(variant + "_", "adaptive_", 1) for c in R.columns if c.startswith(variant + "_")}) if variant != "adaptive" else R; R.to_csv(sfdir(root) / "analysis" / (f"f2_cells_{pl}.csv" if variant == "adaptive" else f"f2track_cells_{pl}.csv"), index=False)
        tot = lambda c, k: float(R[f"{c}_{k}"].sum()); frac = lambda c, k: tot(c, k) / tot(c, "n_pass")
        dsoft = (R.adaptive_soft / R.adaptive_n_pass - R.fixed_soft / R.fixed_n_pass).values; lo, hi = boot_ci(dsoft, rng, AN["bootstrap_draws"])
        bad = lambda c: int(((R[f"{c}_retcode"] != "Success") | R[f"{c}_stopped"]).sum())
        w1 = bool(frac("adaptive", "soft") < frac("fixed", "soft") and hi < 0 and frac("adaptive", "hard") <= frac("fixed", "hard"))
        dec_a, dec_f = float(R.adaptive_decay_m_per_day.median()), float(R.fixed_decay_m_per_day.median())
        w2 = bool(dec_a >= (1 - X) * dec_f and R.adaptive_dv.mean() <= (1 + X) * R.fixed_dv.mean() and bad("adaptive") <= bad("fixed"))
        loss = bool((lo > 0) or frac("adaptive", "hard") > frac("fixed", "hard"))
        cost_win = bool((not loss) and dec_a > (1 + X) * dec_f and R.adaptive_dv.mean() <= (1 + X) * R.fixed_dv.mean())
        outcome = "LOSS" if loss else "WIN" if (w1 and w2) else "SAFETY_TRADE" if w1 else "COST_WIN" if cost_win else "TIE"
        res[pl] = dict(n_cells=len(R), frac_soft=dict(fixed=frac("fixed", "soft"), adaptive=frac("adaptive", "soft")), frac_hard=dict(fixed=frac("fixed", "hard"), adaptive=frac("adaptive", "hard")),
                       frac_low=dict(fixed=frac("fixed", "low"), adaptive=frac("adaptive", "low")), dsoft_mean=float(np.mean(dsoft)), dsoft_ci90=[lo, hi], dv_mean=dict(fixed=float(R.fixed_dv.mean()), adaptive=float(R.adaptive_dv.mean())),
                       burns_mean=dict(fixed=float(R.fixed_burns.mean()), adaptive=float(R.adaptive_burns.mean())), decay_median_m_per_day=dict(fixed=dec_f, adaptive=dec_a),
                       passes_total=dict(fixed=tot("fixed", "n_pass"), adaptive=tot("adaptive", "n_pass")), bad_runs=dict(fixed=bad("fixed"), adaptive=bad("adaptive")),
                       collapsed_decisions_total=int(R.adaptive_collapsed.sum()), W1=w1, W2=w2, outcome=outcome)
    (sfdir(root) / "analysis" / ("f2_results.json" if variant == "adaptive" else "f2track_results.json")).write_text(json.dumps(res, indent=1, default=str)); print(json.dumps(res, indent=1, default=str))


def phase_xhost(a, b):
    rows = []; ok_all = True; nb = 0; nf = 0
    for fa in sorted(pathlib.Path(a).rglob("passes.csv")) + sorted(pathlib.Path(a).rglob("calibration.csv")) + sorted(pathlib.Path(a).rglob("guess.csv")) + sorted(pathlib.Path(a).rglob("nodes.csv")):
        rel = fa.relative_to(a); fb = pathlib.Path(b) / rel; nf += 1
        if not fb.exists(): rows.append((str(rel), "MISSING in B")); ok_all = False; continue
        if fa.read_bytes() == fb.read_bytes(): nb += 1; rows.append((str(rel), "bitwise identical")); continue
        A_, B_ = pd.read_csv(fa), pd.read_csv(fb); num = [c for c in A_.columns if np.issubdtype(A_[c].dtype, np.number) and c not in ("wall_s", "cal_wall_s", "run_wall_s")]
        rel_d = max(float(np.nanmax(np.abs(A_[c].values - B_[c].values) / np.maximum(np.abs(A_[c].values), 1e-300))) if len(A_) == len(B_) else np.inf for c in num)
        altd = float(np.nanmax(np.abs(A_.alt_m.values - B_.alt_m.values))) if "alt_m" in A_.columns and len(A_) == len(B_) else 0.0
        tol = 1e-6 if rel.name == "passes.csv" else 1e-9; good = rel_d <= tol and altd <= 1.0 and len(A_) == len(B_); ok_all &= good
        rows.append((str(rel), f"max rel diff {rel_d:.3e}, max alt diff {altd:.3e} m, rows {len(A_)}/{len(B_)} -> {'within tolerance' if good else 'OUT OF TOLERANCE'}"))
    for r in rows: print(*r, sep="  ")
    print(f"files compared {nf}, bitwise identical {nb}, overall {'PASS' if ok_all and nf > 0 else 'FAIL'}")
    return ok_all and nf > 0


def selftest():
    # exact case: e of persistence on y = linear ramp is constant, so sigma = step; coverage of the empirical quantile is 1
    rec = []
    for t in "ab":
        y = np.arange(30) * 0.01
        for j in range(1, 30): rec.append(dict(tag=t, run=1, pass_=j, j=j, y=y[j], z=0.0, ysys=y[j], H=5000.0, lon=0.1 * j, alt_m=1.0, split="dev", zprev=0.0, **{"e_" + p: 0.01 if p == "P1" else 0.5 for p in ALLP}))
    D = pd.DataFrame(rec); S = score(D); r = S[(S.pred == "P1") & (S.kind == "emp_quantile") & (S.level == 0.95)].iloc[0]
    assert abs(r.sigma - 0.01) < 1e-12 and r.coverage == 1.0 and abs(r.buf_ln - 0.01) < 1e-12
    # no leakage: changing the last y does not change its own prediction
    y = np.arange(10, dtype=float); lon = 0.1 * np.arange(10); p1 = dba.predict_run(y, lon); y2 = y.copy(); y2[-1] = 99; p2 = dba.predict_run(y2, lon)
    assert all(p1[k][-1] == p2[k][-1] for k in p1)
    # white-noise case: persistence variance 2 s^2, EWMA below it
    rng = np.random.default_rng(1); z = rng.normal(0, 0.3, 4000); pr = dba.predict_run(z, np.zeros(4000)); e1 = rms((z - pr["P1"])[10:]); e10 = rms((z - pr["P2_10"])[10:])
    assert abs(e1 / 0.3 - np.sqrt(2)) < 0.06 and e10 < e1
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--phase", required=True); ap.add_argument("--results-root"); ap.add_argument("--root-a"); ap.add_argument("--root-b"); a = ap.parse_args()
    if a.phase == "selftest": return selftest()
    if a.phase == "xhost": sys.exit(0 if phase_xhost(a.root_a, a.root_b) else 1)
    (sfdir(a.results_root) / "analysis").mkdir(parents=True, exist_ok=True)
    {"f1": phase_f1, "f2params": phase_f2params, "f2": phase_f2, "f2track": lambda r: phase_f2(r, "track")}[a.phase](a.results_root)


if __name__ == "__main__": main()
