"""Stage L analysis: the dV / campaign-time / violation frontier at the paper corridor. Plan: tasks/2026-10-15_periapsis-region-border-targeting/stagel_plan.md (frozen before any flight).
Reuses stagei_analysis.py (tiers, run_metrics) and stagef20_analysis.py (load_runs, build, boot_ci) unchanged; Stage I's tree is READ from <root>/stagei, Stage L's own runs live in <root>/stagel.
  selftest                       exact-case checks (arm codes, selection and dominance rules, trend, decay-vs-Q diagnostic)
  params   --results-root R      family parameters: e_hi(p) = empirical p-quantile of the Stage I development one-pass-ahead error, e_lo from the Stage I refit -> <R>/stagel/<planet>/family_params.csv
  repro    --results-root R      bit-level check of the repro flights (fm01000, fam_p095_sp000, fam_p095_sm100) against the Stage I fixed / G / F runs on the same cell
  select   --results-root R      development selection from the dev half-cell flights: frontier, M*_dV, M*_decay, held-out member list -> <R>/stagel/<planet>/selected.csv, analysis/l_dev.json
  held     --results-root R      held-out criteria L1, L2, the frontier tables and L3 trends -> analysis/l_held.json
  diag     --results-root R      per-pass apoapsis decay versus the pass's heat rate (diagnostic, not gated) -> analysis/l_diag.json
  integrity --results-root R     route, failures, GRAM epoch, NaN scan
Usage: python3 -I stagel_analysis.py --phase ... --results-root DIR"""
import argparse, glob, json, pathlib, sys, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stagei_analysis as sia
sf = sia.sf
CL = tomllib.load(open(HERE / "stagel_config.toml", "rb")); SP = CL["stagel"]
CEIL, KMAX, DVF, DECF = SP["violation_ceiling"], SP["held_max_members"], SP["dv_factor"], SP["decay_fraction_min"]
QBINS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.45, np.inf]
STAGE_I_ARMS = {"I_fixed1000": "fixed", "I_F": "adaptive", "I_G": "adaptive_asym", "I_H": sia.HCTRL}
FIXED = "fm01000"


def ldir(root): return pathlib.Path(root) / CL["stagef"]["results_subdir"]
def idir(root): return pathlib.Path(root) / "stagei"


def arm_code(kind, margin_m=None, p=None, s=None):
    """Must equal stagel_arm_code in stagel_common.jl."""
    if kind == "fixedm": return f"fm{round(margin_m):05d}"
    return f"fam_p{round(100 * p):03d}_s{'m' if s < 0 else 'p'}{round(100 * abs(s)):03d}"


def grid():
    g = [dict(kind="fixedm", margin_m=m, p_hi=np.nan, s_frac=np.nan, code=arm_code("fixedm", m)) for m in SP["fixed_margins_m"]]
    g += [dict(kind="adaptive_l", margin_m=np.nan, p_hi=p, s_frac=s, code=arm_code("adaptive_l", p=p, s=s)) for p in SP["p_hi"] for s in SP["s_frac"]]
    return g


def dev_cells(root, pl):
    """The dev half-cells: first dev cell of each (epoch, inclination) group in file order (same rule as stagel_run.jl)."""
    D = pd.read_csv(idir(root) / pl / "cells_base.csv"); D = D[D.split == "dev"]
    return D[~D.duplicated(["epoch_idx", "inc"], keep="first")].tag.tolist()


def held_cells(root, pl):
    D = pd.read_csv(idir(root) / pl / "cells_base.csv"); return D[D.split == "heldout"].tag.tolist()


def l_summaries(root):
    fs = [f for ph in ("dev", "held", "heldall") for f in glob.glob(str(ldir(root) / f"summary_{ph}_*.csv"))]
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True) if fs else pd.DataFrame()


def run_of(root, pl, tag, ctrl, S, stage):
    """Metrics of one run; stage "L" reads <root>/stagel, stage "I" reads <root>/stagei; None when the run or its summary row is missing."""
    d = (ldir(root) if stage == "L" else idir(root)) / pl / "runs" / f"{ctrl}_pert_{tag}"
    sr = S[(S.planet == pl) & (S.run_id == f"{ctrl}_pert_{tag}")] if len(S) else []
    if not (d / "passes.csv").exists() or len(sr) == 0: return None
    dec = pd.read_csv(d / "decisions.csv") if (d / "decisions.csv").exists() else pd.DataFrame()
    return sia.run_metrics(pd.read_csv(d / "passes.csv"), sr.iloc[-1].to_dict(), dec)


def pool(ms):
    """Pooled metrics of a list of per-run metric dicts (one per cell)."""
    n = sum(m["n_pass"] for m in ms); tot = lambda k: float(sum(m[k] for m in ms))
    return dict(n_cells=len(ms), passes=int(n), low=int(tot("low")), soft=int(tot("soft")), medium=int(tot("medium")), hard=int(tot("hard")), over=int(tot("over")), frac_low=tot("low") / n, frac_over=tot("over") / n,
                dv_mean=float(np.mean([m["dv"] for m in ms])), burns_mean=float(np.mean([m["burns"] for m in ms])), decay_median=float(np.median([m["decay"] for m in ms])),
                bad_runs=int(sum((m["retcode"] != "Success") or m["stopped"] for m in ms)), collapsed=int(tot("collapsed")), decisions=int(tot("n_dec")))


def nondominated(T):
    """Rows of T (columns dv_mean, decay_median) not dominated: no other row has dv <= and decay >= with one strict."""
    keep = []
    for i, a in T.iterrows():
        dom = ((T.dv_mean <= a.dv_mean) & (T.decay_median >= a.decay_median) & ((T.dv_mean < a.dv_mean) | (T.decay_median > a.decay_median))).any()
        keep.append(not dom)
    return T[keep]


def select_members(T, ceiling=CEIL, kmax=KMAX):
    """Pre-registered development selection (plan section 4). T: one row per arm with code, frac_over, dv_mean, decay_median. Returns (M_dv, M_decay, held_codes, fallback).
    Feasible = frac_over <= ceiling. M_dv = min dV feasible (ties: higher decay, then code); M_decay = max decay feasible (ties: lower dV, then code). Held set = the non-dominated feasible arms,
    thinned to kmax by evenly spaced dV rank (both extremes, hence M_dv and M_decay, always kept). No feasible arm: fallback = the arm with the lowest frac_over (ties: lower dV, then code), held set = {fallback}."""
    F = T[T.frac_over <= ceiling]
    if len(F) == 0:
        f = T.sort_values(["frac_over", "dv_mean", "code"]).iloc[0].code; return f, f, [f], True
    m_dv = F.sort_values(["dv_mean", "decay_median", "code"], ascending=[True, False, True]).iloc[0].code
    m_de = F.sort_values(["decay_median", "dv_mean", "code"], ascending=[False, True, True]).iloc[0].code
    N = nondominated(F).sort_values(["dv_mean", "code"]); codes = N.code.tolist()
    if len(codes) > kmax: codes = [codes[i] for i in sorted(set(np.round(np.linspace(0, len(codes) - 1, kmax)).astype(int)))]
    for c in (m_dv, m_de):
        if c not in codes: codes.append(c)
    return m_dv, m_de, codes, False


def trend(T, cols=("dv_mean", "decay_median", "frac_over", "burns_mean")):
    """L3: Spearman rank correlation of each metric with p_hi inside each s slice and with s inside each p slice (family rows only; slices with fewer than 3 points are skipped)."""
    Fm = T[T.kind == "adaptive_l"]; out = {"vs_p_at_s": {}, "vs_s_at_p": {}}
    for s, G in Fm.groupby("s_frac"):
        if len(G) >= 3: out["vs_p_at_s"][f"{s:+.2f}"] = {c: float(G.p_hi.corr(G[c], method="spearman")) for c in cols}
    for p, G in Fm.groupby("p_hi"):
        if len(G) >= 3: out["vs_s_at_p"][f"{p:.2f}"] = {c: float(G.s_frac.corr(G[c], method="spearman")) for c in cols}
    return out


def phase_params(root, pls):
    for pl in pls:
        P = sf.load_runs(root, pl, "fixed", "pert"); D = sf.build(P); D = D[D.split == "dev"]
        ap = pd.read_csv(idir(root) / pl / "adapt_params.csv").iloc[0]; best = str(ap.predictor); e = D["e_" + best]
        r = dict(planet=pl, predictor=best, e_lo=float(ap.e_lo), n_dev_scored=len(e), **{f"e_hi_{round(100 * p):03d}": float(np.quantile(e, p)) for p in SP["p_hi"]})
        assert len(e) == int(ap.n_dev_scored), (len(e), ap.n_dev_scored)
        assert abs(r["e_hi_095"] - float(ap.e_hi)) < 1e-9, (r["e_hi_095"], ap.e_hi)   # p = 0.95 reproduces the Stage I e_hi exactly
        (ldir(root) / pl).mkdir(parents=True, exist_ok=True); pd.DataFrame([r]).to_csv(ldir(root) / pl / "family_params.csv", index=False); print("FAMILY", r)


def phase_repro(root, pls):
    SI = sf.summaries(root, "held"); res = {}
    for pl in pls:
        d = ldir(root) / pl / "repro"
        if not d.exists(): continue
        out = {}
        for code, ci in (("fm01000", "fixed"), ("fam_p095_sp000", "adaptive_asym"), ("fam_p095_sm100", "adaptive")):
            for rd in sorted(d.glob(f"{code}_pert_*")):
                tag = rd.name.split("_pert_")[1]; ref = idir(root) / pl / "runs" / f"{ci}_pert_{tag}" / "passes.csv"
                if not ref.exists(): out[code] = "reference missing"; continue
                a, b = pd.read_csv(rd / "passes.csv"), pd.read_csv(ref)
                if len(a) != len(b): out[code] = dict(identical=False, n=(len(a), len(b))); continue
                dmax = float(max(np.nanmax(np.abs(a[c].values - b[c].values)) for c in ("alt_m", "apo_alt_m", "q_eff_w_cm2", "mass_kg")))
                out[code] = dict(identical=bool(dmax == 0.0), max_abs_diff=dmax, passes=len(a), tag=tag)
        res[pl] = out
    print(json.dumps(res, indent=1)); (ldir(root) / "analysis").mkdir(exist_ok=True); (ldir(root) / "analysis" / "l_repro.json").write_text(json.dumps(res, indent=1))


def dev_table(root, pl, S):
    tags = dev_cells(root, pl); rows = []; missing = {}
    for a in grid():
        ms = [run_of(root, pl, t, a["code"], S, "L") for t in tags]
        if any(m is None for m in ms): missing[a["code"]] = sum(m is None for m in ms); continue
        rows.append(dict(**a, **pool(ms)))
    return pd.DataFrame(rows), missing


def phase_select(root, pls):
    S = l_summaries(root); res = {}
    for pl in pls:
        T, miss = dev_table(root, pl, S)
        if len(T) == 0 or FIXED not in set(T.code): print(pl, "no complete dev arms / fixed 1000 m missing", miss); continue
        m_dv, m_de, codes, fb = select_members(T); fx = T[T.code == FIXED].iloc[0]
        sel = pd.DataFrame([dict(code=c, role=("M_dv" if c == m_dv else "") + ("+" if c == m_dv and c == m_de else "") + ("M_decay" if c == m_de else "") or "frontier", kind=T[T.code == c].iloc[0].kind,
                                 margin_m=T[T.code == c].iloc[0].margin_m, p_hi=T[T.code == c].iloc[0].p_hi, s_frac=T[T.code == c].iloc[0].s_frac) for c in codes])
        sel.to_csv(ldir(root) / pl / "selected.csv", index=False)
        (ldir(root) / "analysis").mkdir(exist_ok=True); T.to_csv(ldir(root) / "analysis" / f"l_dev_{pl}.csv", index=False)
        res[pl] = dict(n_cells=int(T.n_cells.iloc[0]), arms_complete=len(T), arms_missing=miss, fixed1000=fx[["frac_over", "dv_mean", "decay_median", "burns_mean", "hard"]].to_dict(), ceiling=CEIL, feasible=int((T.frac_over <= CEIL).sum()),
                       M_dv=m_dv, M_decay=m_de, held_codes=codes, fallback=fb, trend=trend(T))
        print(pl, json.dumps({k: v for k, v in res[pl].items() if k != "trend"}, default=str))
    (ldir(root) / "analysis").mkdir(exist_ok=True); (ldir(root) / "analysis" / "l_dev.json").write_text(json.dumps(res, indent=1, default=str))


def criteria(R, a, fx_over, fx_dv, fx_dec, rng):
    """L1 / L2 conditions for arm a against the Stage I fixed 1000 m arm: over-0.25 fraction <= fixed's, mean dV <= DVF x fixed's, median decay >= DECF x fixed's (plan section 5)."""
    c = dict(over_le_fixed=bool(a["frac_over"] <= fx_over), dv_le_2x=bool(a["dv_mean"] <= DVF * fx_dv), decay_ge_90pct=bool(a["decay_median"] >= DECF * fx_dec))
    c["all"] = all(c.values()); return c


def phase_held(root, pls):
    SL = l_summaries(root); SI = sf.summaries(root, "held"); rng = np.random.default_rng(SP["bootstrap_seed"]); res = {}
    for pl in pls:
        tags = held_cells(root, pl); sel_path = ldir(root) / pl / "selected.csv"
        sel = pd.read_csv(sel_path) if sel_path.exists() else pd.DataFrame(columns=["code", "role"])
        rows = []; per_cell = {}
        for name, ctrl in STAGE_I_ARMS.items():
            ms = [run_of(root, pl, t, ctrl, SI, "I") for t in tags]
            if all(m is not None for m in ms): rows.append(dict(code=name, kind="stageI", margin_m=np.nan, p_hi=np.nan, s_frac=np.nan, **pool(ms))); per_cell[name] = ms
        for a in grid():
            ms = [run_of(root, pl, t, a["code"], SL, "L") for t in tags]
            if all(m is not None for m in ms): rows.append(dict(**a, **pool(ms))); per_cell[a["code"]] = ms
        T = pd.DataFrame(rows)
        if "I_fixed1000" not in set(T.code): print(pl, "Stage I fixed arm missing"); continue
        fx = T[T.code == "I_fixed1000"].iloc[0]; T["dev_selected"] = T.code.isin(sel.code); T["meets_ceiling"] = T.frac_over <= CEIL
        T["dv_ratio"] = T.dv_mean / fx.dv_mean; T["decay_ratio"] = T.decay_median / fx.decay_median
        (ldir(root) / "analysis").mkdir(exist_ok=True); T.to_csv(ldir(root) / "analysis" / f"l_held_{pl}.csv", index=False)
        out = dict(n_cells=len(tags), fixed1000=fx[["frac_over", "hard", "dv_mean", "decay_median", "burns_mean"]].to_dict(), arms_complete=len(T))
        crit = {}
        for role in ("M_dv", "M_decay"):
            c = sel[sel.role.str.contains(role)].code.tolist()
            if not c or c[0] not in set(T.code): crit[role] = "not flown"; continue
            a = T[T.code == c[0]].iloc[0]; d = (np.array([m["over"] / m["n_pass"] for m in per_cell[c[0]]]) - np.array([m["over"] / m["n_pass"] for m in per_cell["I_fixed1000"]]))
            lo, hi = sf.boot_ci(d, rng, SP["bootstrap_draws"])
            crit[role] = dict(code=c[0], **a[["frac_over", "hard", "dv_mean", "dv_ratio", "decay_median", "decay_ratio", "burns_mean", "frac_low"]].to_dict(), dover_mean=float(d.mean()), dover_ci90=[lo, hi],
                              criteria=criteria(None, a, fx.frac_over, fx.dv_mean, fx.decay_median, rng), meets_ceiling_held=bool(a.frac_over <= CEIL), strict_dv_le_1p1x=bool(a.dv_ratio <= 1 + SP["dv_tolerance_strict"]))
        out["selected"] = crit
        # descriptive frontier on the held-out cells (all complete arms; NOT used for selection)
        Fz = T[T.meets_ceiling]; out["frontier_held_descriptive"] = dict(ceiling=CEIL, nondominated=nondominated(Fz).code.tolist() if len(Fz) else [],
                                                                       min_dv_meeting_ceiling=(Fz.sort_values("dv_mean").iloc[0].code if len(Fz) else None))
        out["trend_held"] = trend(T)
        out["table"] = T[["code", "kind", "p_hi", "s_frac", "margin_m", "dev_selected", "frac_over", "frac_low", "hard", "dv_mean", "dv_ratio", "burns_mean", "decay_median", "decay_ratio", "bad_runs"]].to_dict(orient="records")
        res[pl] = out
    verdict = {}
    if "mars" in res: m = res["mars"]["selected"].get("M_dv"); verdict["L1_mars"] = (m["criteria"]["all"] if isinstance(m, dict) else "not flown")
    ev = [res[p]["selected"].get("M_dv") for p in ("earth", "venus") if p in res]
    if ev: verdict["L2_earth_venus"] = {p: (res[p]["selected"]["M_dv"]["criteria"]["all"] if isinstance(res[p]["selected"].get("M_dv"), dict) else "not flown") for p in ("earth", "venus") if p in res}
    res["verdict"] = verdict; (ldir(root) / "analysis").mkdir(exist_ok=True)
    (ldir(root) / "analysis" / "l_held.json").write_text(json.dumps(res, indent=1, default=str)); print(json.dumps(res, indent=1, default=str))


def decay_pairs(P, dec):
    """(q_prev, q_next, d_apo) per consecutive pass pair with finite apoapsis and no burn decided for the second pass (a burn before that pass changes the orbit the diagnostic is meant to see);
    d_apo = apo_k - apo_(k+1) > 0 is decay."""
    P = P.sort_values("pass"); burn = set()
    if len(dec) and "target_pass" in dec: burn = set(dec[dec.delta_v_mps.abs() > 0].target_pass.astype(int))
    k = P["pass"].values; apo = P.apo_alt_m.values; q = P.q_eff_w_cm2.values; out = []
    for i in range(len(P) - 1):
        if k[i + 1] - k[i] == 1 and np.isfinite(apo[i]) and np.isfinite(apo[i + 1]) and int(k[i + 1]) not in burn: out.append((q[i], q[i + 1], apo[i] - apo[i + 1]))
    return out


def diag_from_pairs(pairs):
    A = np.array(pairs, float)
    if len(A) == 0: return {}
    bins = []
    for lo, hi in zip(QBINS[:-1], QBINS[1:]):
        m = (A[:, 1] >= lo) & (A[:, 1] < hi)
        bins.append(dict(q_lo=lo, q_hi=None if np.isinf(hi) else hi, n=int(m.sum()), median_decay_m_per_pass=float(np.median(A[m, 2])) if m.any() else None))
    slope = lambda x, y: float((x * y).sum() / (x * x).sum())
    return dict(n=len(A), bins=bins, slope_decay_per_q_next=slope(A[:, 1], A[:, 2]), slope_decay_per_q_prev=slope(A[:, 0], A[:, 2]),
                corr_next=float(np.corrcoef(A[:, 1], A[:, 2])[0, 1]), corr_prev=float(np.corrcoef(A[:, 0], A[:, 2])[0, 1]))


def phase_diag(root, pls):
    res = {}
    for pl in pls:
        pairs = []
        for base in (ldir(root), idir(root)):
            for rd in sorted((base / pl / "runs").glob("*_pert_*")):
                if not (rd / "passes.csv").exists(): continue
                dec = pd.read_csv(rd / "decisions.csv") if (rd / "decisions.csv").exists() else pd.DataFrame()
                pairs += decay_pairs(pd.read_csv(rd / "passes.csv"), dec)
        res[pl] = diag_from_pairs(pairs)
    (ldir(root) / "analysis").mkdir(exist_ok=True); (ldir(root) / "analysis" / "l_diag.json").write_text(json.dumps(res, indent=1, default=str)); print(json.dumps(res, indent=1, default=str))


def phase_integrity(root, pls):
    out = {}
    for f in sorted(glob.glob(str(ldir(root) / "mc_*.toml"))):
        t = tomllib.load(open(f, "rb")); out[pathlib.Path(f).name] = dict(route=t["mc_route"], failed=t["mc_failed"], samples=t["mc_samples"])
    S = l_summaries(root)
    for pl in pls:
        runs = list((ldir(root) / pl / "runs").glob("*/passes.csv")) if (ldir(root) / pl / "runs").exists() else []
        nan = sum(int(pd.read_csv(p).drop(columns=["h_qhigh_m"], errors="ignore").isna().sum().sum()) for p in runs)
        ep = [pd.Timestamp(x).strftime("%Y-%-m-%-d %-H:%-M:%-S") for x in sf.CFG["planets"][pl]["grid"]["epochs_utc"]] if pl in sf.CFG["planets"] else []
        sp = S[S.planet == pl] if len(S) else S
        out[pl] = dict(runs=len(runs), nan_cells_excl_h_qhigh=nan, rows=len(sp), not_ok=int((sp.status != "ok").sum()) if len(sp) else 0, epochs_match=bool(all(any(str(e).startswith(x) for x in ep) for e in sp.gram_core_epoch)) if len(sp) and ep else None)
    (ldir(root) / "analysis").mkdir(exist_ok=True); (ldir(root) / "analysis" / "l_integrity.json").write_text(json.dumps(out, indent=1, default=str)); print(json.dumps(out, indent=1, default=str))


def selftest():
    assert arm_code("fixedm", 500.0) == "fm00500" and arm_code("adaptive_l", p=0.8, s=-1.0) == "fam_p080_sm100" and arm_code("adaptive_l", p=0.99, s=0.5) == "fam_p099_sp050" and arm_code("adaptive_l", p=0.95, s=0.0) == "fam_p095_sp000"
    g = grid(); assert len(g) == len(SP["fixed_margins_m"]) + len(SP["p_hi"]) * len(SP["s_frac"]) and len({a["code"] for a in g}) == len(g) and FIXED in {a["code"] for a in g}
    # selection on an exact table: B and C feasible and non-dominated, D feasible but dominated by C, A infeasible; M_dv = B (min dV), M_decay = C (max decay)
    T = pd.DataFrame(dict(code=list("ABCD"), kind="adaptive_l", frac_over=[0.05, 0.005, 0.008, 0.009], dv_mean=[0.1, 0.5, 0.9, 1.0], decay_median=[100.0, 80.0, 90.0, 85.0], p_hi=0.9, s_frac=0.0, margin_m=np.nan, burns_mean=1.0))
    m_dv, m_de, codes, fb = select_members(T); assert (m_dv, m_de, fb) == ("B", "C", False) and sorted(codes) == ["B", "C"], (m_dv, m_de, codes)
    f = select_members(T.assign(frac_over=[0.05, 0.03, 0.02, 0.04])); assert f == ("C", "C", ["C"], True), f      # none feasible: lowest over-fraction
    big = pd.DataFrame(dict(code=[f"X{i}" for i in range(8)], kind="x", frac_over=0.0, dv_mean=np.arange(8.0), decay_median=np.arange(8.0) * 10, p_hi=0.9, s_frac=0.0, margin_m=np.nan, burns_mean=1.0))
    mdv, mde, codes, fb = select_members(big); assert len(codes) == KMAX and "X0" in codes and "X7" in codes, codes    # all 8 non-dominated: thinned to KMAX incl. both extremes
    # trend: dV increasing in p -> rho +1; decay decreasing in p -> -1
    F = pd.DataFrame([dict(kind="adaptive_l", p_hi=p, s_frac=s, dv_mean=p * 10 + s, decay_median=-p * 5, frac_over=-p, burns_mean=p) for p in (0.8, 0.9, 0.95, 0.99) for s in (-1.0, 0.0, 0.5)])
    tr = trend(F); assert abs(tr["vs_p_at_s"]["+0.00"]["dv_mean"] - 1) < 1e-12 and abs(tr["vs_p_at_s"]["+0.00"]["decay_median"] + 1) < 1e-12 and abs(tr["vs_s_at_p"]["0.90"]["dv_mean"] - 1) < 1e-12
    # decay diagnostic on an exact case: decay per pass = 1000 * q_next, passes 1..5, the burn-decided pass 4 pair (3 -> 4) is excluded
    q = np.array([0.1, 0.2, 0.1, 0.3, 0.15]); apo = [1e6]
    for i in range(1, 5): apo.append(apo[-1] - 1000.0 * q[i])
    P = pd.DataFrame(dict(**{"pass": [1, 2, 3, 4, 5]}, apo_alt_m=apo, q_eff_w_cm2=q)); dec = pd.DataFrame(dict(target_pass=[4], delta_v_mps=[0.5]))
    pr = decay_pairs(P, dec); assert len(pr) == 3, pr
    d = diag_from_pairs(pr); assert abs(d["slope_decay_per_q_next"] - 1000.0) < 1e-9 and d["n"] == 3
    assert abs(float(np.quantile(np.arange(101) / 100.0, 0.95)) - 0.95) < 1e-12
    print("stagel selftest ok")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--phase", required=True); ap.add_argument("--results-root"); ap.add_argument("--planets", default="mars,earth,venus"); a = ap.parse_args()
    if a.phase == "selftest": return selftest()
    pls = a.planets.split(","); root = a.results_root
    {"params": phase_params, "repro": phase_repro, "select": phase_select, "held": phase_held, "diag": phase_diag, "integrity": phase_integrity}[a.phase](root, pls)


if __name__ == "__main__": main()
