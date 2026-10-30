"""Stage N2 analysis: edge-aimed controller with a calibrated aim offset. Plan: tasks/2026-10-15_periapsis-region-border-targeting/stagen2_plan.md (frozen before any N2 flight).
Reuses stagen_flight_analysis.py (run metrics, innovations join, CP, bootstrap), stagen_analysis.py (frozen filter) and stagel_analysis.py unchanged.
  selftest | calibrate --results-root R | select | held | repro | integrity"""
import argparse, glob, json, pathlib, sys, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import stagen_flight_analysis as nf
import stagen_analysis as sn
import stagel_analysis as sla
import stagei_analysis as sia
import stagej_analysis as sj
sf = sia.sf
CN = tomllib.load(open(HERE / "stagen2_config.toml", "rb")); SP = CN["stagen2"]
PL = "mars"; QMAX, QLOW, QHARD = 0.25, 0.05, 0.45; EPS = SP["epsilon"]; N_ARM = "cc01_w04000"; W_CAP = SP["w_cap_m"]; LMEMBERS = nf.LMEMBERS
def n2dir(root): return pathlib.Path(root) / CN["stagef"]["results_subdir"]


def flight_table(root, tags, base, code):
    """Per decided pass of the runs <base>/mars/runs/<code>_pert_<tag>: predictive error e, c, s, the realized q_eff, the controller's predicted position pos = ln(Q_pred / Q_edge) (Q_pred = q_eff exp(-e), Q_edge = Q_max exp(-c)),
    band width and mid-band offset in ln Q units (width / h_local), burn flag and signed dV."""
    out = []
    for tag in tags:
        d = pathlib.Path(base) / PL / "runs" / f"{code}_pert_{tag}"; D = pd.read_csv(d / "decisions.csv"); M = pd.read_csv(d / "measured.csv"); P = pd.read_csv(d / "passes.csv")
        for r in D.itertuples():
            k = int(r.target_pass); pr = P[P["pass"] == k]
            if k - 1 >= len(M) or len(pr) == 0: continue
            pr = pr.iloc[0]; e = float(M.y.iloc[k - 1] - r.m); q = float(pr.q_eff_w_cm2); h = float(pr.h_local_m)
            out.append(dict(tag=tag, k=k, e=e, s=float(r.s), c=float(r.c), q=q, pos=float(np.log(q) - e - np.log(QMAX) + r.c), width_ln=float((r.upper_m - r.lower_m) / h), g0=float((r.upper_m - r.lower_m) / 2 / h),
                            burned=abs(float(r.delta_v_mps)) > 0, dv=float(r.delta_v_mps)))
    return pd.DataFrame(out)


def pass1_over(root, tags):
    q = [pd.read_csv(sla.idir(root) / PL / "runs" / f"fixed_pert_{t}" / "passes.csv").q_eff_w_cm2.iloc[0] for t in tags]; q = np.array(q); return float((q > QMAX).mean()), q


def ref_z(root):
    """Reference standardized innovations: Stage I DEV fixed-arm in-window passes through the frozen filter (the sample the frozen t was fit on)."""
    P = sj.load(root, "stagei", PL, "fixed_pert"); P = P[P.inw].reset_index(drop=True); D = P[P.split == "dev"].reset_index(drop=True); kp = sn.fit_kalman(D)
    _, S, E, _ = sn.kalman(D, kp["r"], kp["q"], kp["m0"], kp["v0"]); return np.sort(E / np.sqrt(S))


def coast_increments(T):
    """Step in the controller's predicted position rho (= pos) between consecutive decided passes when NO burn was decided for the later one (545 steps on the dev half-cells)."""
    inc = []
    for t, G in T.sort_values(["tag", "k"]).groupby("tag"):
        p, k, b = G.pos.values, G.k.values, G.burned.values
        inc += [p[j + 1] - p[j] for j in range(len(G) - 1) if k[j + 1] == k[j] + 1 and not b[j + 1]]
    return np.array(inc)


def burn_kappa(T):
    """dV per unit of ln-Q position moved by a burn: |dv| of the flown dev burns after the first decision over the distance they cover, taken as g0 (raise: edge to mid-band) or width - g0 (lower: far edge to mid-band)."""
    B = T[T.burned & (T.k > 2)]; d = np.where(B.dv > 0, B.g0, B.width_ln - B.g0)
    return float(B.dv.abs().sum() / d.sum()), int(len(B))


def simulate(T, g_a, t, z, q1, inc, kappa, n_sim, seed, n_pass=49):
    """Monte Carlo of the controlled position rho (plan section 3) with common random numbers (seed): rho_2 = -g_a (the first decision burns to the aim, as in every dev cell), then rho += an empirical coast step;
    rho > 0 triggers a raise burn back to -g_a, rho < -width a lower burn; a burn costs kappa x distance (the first burn: kappa x width / 2). Q_k = 0.25 exp(-t s + rho + s z) with s the median predictive sd, z resampled from the
    reference standardized innovations; the first pass is the undecided calibration pass (q1 resampled). Rates are over all passes."""
    rng = np.random.default_rng(seed); s = float(T.s.median()); c = t * s; W = float(T.width_ln.median()); zz = np.asarray(z)
    rho = np.full(n_sim, -g_a); q = rng.choice(q1, n_sim); over = (q > QMAX).astype(float); low = (q < QLOW).astype(float); sq = q.copy(); n = np.ones(n_sim); burns = np.ones(n_sim); dv = np.full(n_sim, kappa * W / 2)
    for _ in range(n_pass - 2):
        rho = rho + rng.choice(inc, n_sim); up = rho > 0; dn = rho < -W
        burns += up | dn; dv += kappa * np.where(up, rho + g_a, 0.0) + kappa * np.where(dn, -g_a - rho, 0.0); rho = np.where(up | dn, -g_a, rho)
        q = QMAX * np.exp(-c + rho + s * rng.choice(zz, n_sim)); over += q > QMAX; low += q < QLOW; sq += q; n += 1
    return dict(g_a=g_a, t=t, over_all=float(over.sum() / n.sum()), low_all=float(low.sum() / n.sum()), meanQ_all=float(sq.sum() / n.sum()), burns_per_run=float(burns.mean()), dv_per_run=float(dv.mean()))


def solve_t(T, g_a, target, z, q1, inc, kappa, n_sim, seed):
    """Edge multiplier t with predicted overall over-0.25 = target (bisection; the prediction is monotone decreasing in t, CRN)."""
    lo, hi = 0.3, 3.0
    for _ in range(30):
        mid = (lo + hi) / 2
        if simulate(T, g_a, mid, z, q1, inc, kappa, n_sim, seed)["over_all"] > target: lo = mid
        else: hi = mid
    return (lo + hi) / 2


def phase_calibrate(root):
    dev = sla.dev_cells(root, PL); held = sla.held_cells(root, PL); base = pathlib.Path(root) / "stagen_flight"
    Td = flight_table(root, dev, base, N_ARM); z = ref_z(root); dev_all = sorted(set(sj.load(root, "stagei", PL, "fixed_pert").query("split == 'dev'").tag))
    _, q1 = pass1_over(root, dev_all)       # undecided first pass: the 24 Stage I dev cells (identical first pass in every arm)
    fixed_dev = float(np.concatenate([pd.read_csv(sla.ldir(root) / PL / "runs" / f"fm01000_pert_{t}" / "passes.csv").q_eff_w_cm2.values for t in dev]).mean())
    inc = coast_increments(Td); kap, nb = burn_kappa(Td); N = SP["mc_draws"]; seed = SP["mc_seed"]; g0 = float(Td.g0.median()); t0 = SP["t_mult"]
    res = dict(dev_cells=len(dev), n_dev_decided=len(Td), n_coast_steps=len(inc), coast_step_mean=float(inc.mean()), coast_step_sd=float(inc.std()), kappa_mps_per_lnQ=kap, kappa_n_burns=nb, g0_median=g0, width_ln_median=float(Td.width_ln.median()),
               s_median=float(Td.s.median()), pass1_n=len(q1), pass1_over=float((q1 > QMAX).mean()), pass1_meanQ=float(q1.mean()), fixed_dev_meanQ=fixed_dev, ref_z_n=len(z), t_stageN=t0)
    flown = dict(over_decided=float((Td.q > QMAX).mean()), meanQ_decided=float(Td.q.mean()), burns_per_run=float(Td.burned.sum() / len(dev)))
    v = simulate(Td, g0, t0, z, q1, inc, kap, N, seed); res["validation_stageN_dev"] = dict(model=v, model_meanQ_ratio=v["meanQ_all"] / fixed_dev, flown=flown, flown_dv_per_run=0.5538, flown_meanQ_ratio_decided=flown["meanQ_decided"] / fixed_dev)
    Th = flight_table(root, held, base, N_ARM); res["held_out_check_only_stageN_flown"] = dict(over_decided=float((Th.q > QMAX).mean()), meanQ_decided=float(Th.q.mean()), burns_per_run=float(Th.burned.sum() / len(held)))
    res["stageN_edge_unchanged"] = {f"g{g:.2f}": simulate(Td, g, t0, z, q1, inc, kap, N, seed) for g in SP["g_grid"]}
    arms = []
    for g in SP["g_grid"]:
        t1 = solve_t(Td, g, 0.01, z, q1, inc, kap, N, seed); r = simulate(Td, g, t1, z, q1, inc, kap, N, seed); r["meanQ_ratio_vs_fixed_dev"] = r["meanQ_all"] / fixed_dev; arms.append(dict(label=f"eps1_g{round(100 * g):03d}", eligible=True, **r))
    g5 = SP["g_diag"]; t5 = solve_t(Td, g5, 0.05, z, q1, inc, kap, N, seed); r = simulate(Td, g5, t5, z, q1, inc, kap, N, seed); r["meanQ_ratio_vs_fixed_dev"] = r["meanQ_all"] / fixed_dev
    arms.append(dict(label=f"eps5_g{round(100 * g5):03d}_diagnostic", eligible=False, **r)); res["arms"] = arms
    n2dir(root).mkdir(parents=True, exist_ok=True); (n2dir(root) / PL).mkdir(parents=True, exist_ok=True); (n2dir(root) / "analysis").mkdir(exist_ok=True)
    G = pd.DataFrame([dict(label=a["label"], g_a=round(a["g_a"], 4), t_mult=round(a["t"], 3), eligible=a["eligible"]) for a in arms]); G.to_csv(n2dir(root) / PL / "aim_grid.csv", index=False)
    (n2dir(root) / "analysis" / "n2_calibration.json").write_text(json.dumps(res, indent=1, default=float)); print(json.dumps(res, indent=1, default=float))



# ---- run metrics, dev selection ---------------------------------------------------------------------------------------------------------
def n2_summaries(root):
    fs = sorted(glob.glob(str(n2dir(root) / "summary_*.csv")))
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True).drop_duplicates("run_id", keep="last") if fs else pd.DataFrame()


def run_n2(root, tag, code, S, sub="runs"):
    d = n2dir(root) / PL / sub / f"{code}_pert_{tag}"; sr = S[S.run_id == f"{code}_pert_{tag}"] if len(S) else []
    if not (d / "passes.csv").exists() or len(sr) == 0: return None
    dec = pd.read_csv(d / "decisions.csv") if (d / "decisions.csv").exists() else pd.DataFrame()
    return sia.run_metrics(pd.read_csv(d / "passes.csv"), sr.iloc[-1].to_dict(), dec)


def grid_codes(root):
    G = pd.read_csv(n2dir(root) / PL / "aim_grid.csv"); G["code"] = [arm_code(EPS, W_CAP, r.g_a, r.t_mult) for r in G.itertuples()]; return G


def arm_code(eps, w, g, t): return f"cc{round(100 * eps):02d}_a{round(100 * g):03d}_t{round(100 * t):03d}_w{round(w):05d}"   # must equal stagen2_arm_code in stagen2_common.jl


def select_arm(T):
    """Pre-registered dev selection (plan section 4). T rows: eligible, g_a, frac_over, dv_mean, dv_ref, decay_median. Feasible = eligible and frac_over <= CEIL and dv_mean <= DVF x dv_ref; pick the highest decay (ties: larger g_a).
    None feasible (FALLBACK, flagged): the eligible arm with the lowest frac_over among dv <= DVF x ref, else lowest frac_over."""
    E = T[T.eligible]; F = E[(E.frac_over <= SP["dev_over_ceiling"]) & (E.dv_mean <= SP["dv_factor"] * E.dv_ref)]
    if len(F): return F.sort_values(["decay_median", "g_a"], ascending=[False, False]).iloc[0], False
    G = E[E.dv_mean <= SP["dv_factor"] * E.dv_ref]; G = G if len(G) else E
    return G.sort_values(["frac_over", "g_a"], ascending=[True, False]).iloc[0], True


def phase_select(root):
    S = n2_summaries(root); SL = sla.l_summaries(root); tags = sla.dev_cells(root, PL); G = grid_codes(root)
    ref = [sla.run_of(root, PL, t, "fm01000", SL, "L") for t in tags]; fx = sla.pool(ref); fxq = float(np.concatenate([pd.read_csv(sla.ldir(root) / PL / "runs" / f"fm01000_pert_{t}" / "passes.csv").q_eff_w_cm2.values for t in tags]).mean())
    base = pathlib.Path(root) / "stagen_flight"; nS = nf.n_summaries(root); nm = [nf.run_n(root, t, N_ARM, nS) for t in tags]; rows = []
    for r in G.itertuples():
        ms = [run_n2(root, t, r.code, S) for t in tags]
        if any(m is None for m in ms): print("incomplete", r.code); continue
        Q = float(np.concatenate([pd.read_csv(n2dir(root) / PL / "runs" / f"{r.code}_pert_{t}" / "passes.csv").q_eff_w_cm2.values for t in tags]).mean())
        rows.append(dict(label=r.label, code=r.code, eligible=bool(r.eligible), g_a=r.g_a, t_mult=r.t_mult, **sla.pool(ms), dv_ref=fx["dv_mean"], meanQ_ratio=Q / fxq))
    T = pd.DataFrame(rows); sel, fb = select_arm(T)
    pd.DataFrame([dict(code=sel.code, g_a=sel.g_a, t_mult=sel.t_mult, label=sel.label, fallback=fb)]).to_csv(n2dir(root) / PL / "selected.csv", index=False)
    (n2dir(root) / "analysis").mkdir(exist_ok=True); T.to_csv(n2dir(root) / "analysis" / "n2_dev.csv", index=False)
    out = dict(fixed1000_dev={k: fx[k] for k in ("frac_over", "dv_mean", "decay_median", "burns_mean", "hard", "frac_low")}, fixed_dev_meanQ=fxq, stageN_arm_dev=sla.pool(nm) if all(m is not None for m in nm) else None,
               table=T.to_dict(orient="records"), selected=dict(code=sel.code, label=sel.label, fallback=fb))
    (n2dir(root) / "analysis" / "n2_dev.json").write_text(json.dumps(out, indent=1, default=float)); print(json.dumps(out, indent=1, default=float))


# ---- held-out criteria --------------------------------------------------------------------------------------------------------------------
def nondominated(rows):
    return [a["name"] for a in rows if not any((b["dv"] <= a["dv"] and b["decay"] >= a["decay"] and (b["dv"] < a["dv"] or b["decay"] > a["decay"])) for b in rows)]


def phase_held(root):
    rng = np.random.default_rng(SP["bootstrap_seed"]); NB = SP["bootstrap_draws"]; S2, SN, SL, SI = n2_summaries(root), nf.n_summaries(root), sla.l_summaries(root), sf.summaries(root, "held")
    tags = sla.held_cells(root, PL); sel = pd.read_csv(n2dir(root) / PL / "selected.csv").iloc[0]; code = str(sel.code); cal = json.load(open(n2dir(root) / "analysis" / "n2_calibration.json"))
    cells, dirs = {}, {}
    for name, ctrl in sla.STAGE_I_ARMS.items(): cells[name] = [sla.run_of(root, PL, t, ctrl, SI, "I") for t in tags]; dirs[name] = (sla.idir(root), ctrl)
    for c in LMEMBERS: cells["L_" + c] = [sla.run_of(root, PL, t, c, SL, "L") for t in tags]; dirs["L_" + c] = (sla.ldir(root), c)
    cells["N_cc"] = [nf.run_n(root, t, N_ARM, SN) for t in tags]; dirs["N_cc"] = (pathlib.Path(root) / "stagen_flight", N_ARM)
    cells["N2"] = [run_n2(root, t, code, S2) for t in tags]; dirs["N2"] = (n2dir(root), code)
    assert all(m is not None for m in cells["N2"]) and all(m is not None for m in cells["I_fixed1000"]), "missing flown cells"
    rows = {}
    for k, ms in cells.items():
        if any(m is None for m in ms): continue
        p = sla.pool(ms); p["mean_Q"] = nf.mean_q(root, tags, *dirs[k]); rows[k] = p
    fx = rows["I_fixed1000"]
    for p in rows.values(): p["dv_ratio"] = p["dv_mean"] / fx["dv_mean"]; p["decay_ratio"] = p["decay_median"] / fx["decay_median"]; p["meanQ_ratio"] = p["mean_Q"] / fx["mean_Q"]
    n = rows["N2"]; ko, N = n["over"], n["passes"]; ci_cp = nf.cp(ko, N)
    oc = np.array([m["over"] for m in cells["N2"]], float); npass = np.array([m["n_pass"] for m in cells["N2"]], float); ci_cl = nf.cluster_boot(np.array(tags), oc, npass, rng, NB)
    d = np.array([m["over"] / m["n_pass"] for m in cells["N2"]]) - np.array([m["over"] / m["n_pass"] for m in cells["I_fixed1000"]]); dlo, dhi = sf.boot_ci(d, rng, NB)
    crit = dict(C1a=dict(over_k=ko, n=N, over=ko / N, ci90_binom=ci_cp, ci90_cluster=ci_cl, epsilon=EPS, contains_eps=bool(ci_cp[0] <= EPS <= ci_cp[1]), contains_eps_cluster=bool(ci_cl[0] <= EPS <= ci_cl[1]),
                         direction=("conservative (interval below eps)" if ci_cp[1] < EPS else "liberal (interval above eps)" if ci_cp[0] > EPS else "calibrated (interval contains eps)")),
                C1b=dict(over=ko / N, fixed_over=SP["over_cmp_fixed"], fixed_over_recomputed=fx["frac_over"], passes=bool(ko / N <= SP["over_cmp_fixed"]), paired_diff_mean=float(d.mean()), paired_diff_ci90=[dlo, dhi]),
                C2=dict(hard=n["hard"], passes=bool(n["hard"] == 0)), C3=dict(dv_mean=n["dv_mean"], fixed=fx["dv_mean"], ratio=n["dv_ratio"], limit=SP["dv_factor"], passes=bool(n["dv_mean"] <= SP["dv_factor"] * fx["dv_mean"])),
                C4=dict(decay_median=n["decay_median"], fixed=fx["decay_median"], ratio=n["decay_ratio"], limit=SP["decay_fraction_min"], passes=bool(n["decay_median"] >= SP["decay_fraction_min"] * fx["decay_median"])),
                C5=dict(meanQ_ratio=n["meanQ_ratio"], mean_Q=n["mean_Q"], fixed_mean_Q=fx["mean_Q"], frac_low=n["frac_low"], fixed_frac_low=fx["frac_low"], burns_mean=n["burns_mean"], fixed_burns=fx["burns_mean"]))
    crit["C1"] = bool(crit["C1a"]["contains_eps"] and crit["C1b"]["passes"]); crit["verdict"] = dict(C1=crit["C1"], C2=crit["C2"]["passes"], C3=crit["C3"]["passes"], C4=crit["C4"]["passes"])
    res = dict(arm=code, selected=sel.to_dict(), criteria=crit, table={k: {kk: p[kk] for kk in ("passes", "over", "frac_over", "hard", "frac_low", "dv_mean", "dv_ratio", "burns_mean", "decay_median", "decay_ratio", "mean_Q", "meanQ_ratio", "bad_runs")} for k, p in rows.items()})
    # frontier (plan section 5): among arms with held over fraction <= 1% (the Stage L ceiling), non-dominated in (dV, decay); and over all arms
    arms = [dict(name=k, dv=p["dv_mean"], decay=p["decay_median"], over=p["frac_over"]) for k, p in rows.items()]; q = [a for a in arms if a["over"] <= EPS]
    nd_q = nondominated(q); nd_wo = nondominated([a for a in q if a["name"] != "N2"]); n2a = [a for a in arms if a["name"] == "N2"][0]
    dominated_by_n2 = [a["name"] for a in arms if a["name"] != "N2" and n2a["dv"] <= a["dv"] and n2a["decay"] >= a["decay"] and (n2a["dv"] < a["dv"] or n2a["decay"] > a["decay"])]
    res["frontier"] = dict(ceiling=EPS, qualifying=[a["name"] for a in q], nondominated_with_N2=nd_q, nondominated_without_N2=nd_wo, N2_on_frontier=bool("N2" in nd_q), N2_dominates=dominated_by_n2, N2_dominated_by=[a["name"] for a in arms if a["name"] != "N2" and a["dv"] <= n2a["dv"] and a["decay"] >= n2a["decay"] and (a["dv"] < n2a["dv"] or a["decay"] > n2a["decay"])])
    res["frontier"]["moved"] = bool(res["frontier"]["N2_on_frontier"] and (set(nd_wo) - set(nd_q) or len(dominated_by_n2) > 0))
    # prediction versus flight (calibration plan) and closed-loop innovations
    T = flight_table(root, tags, n2dir(root).parent / "stagen2", code); a = [x for x in cal["arms"] if abs(x["g_a"] - sel.g_a) < 1e-3 and abs(x["t"] - sel.t_mult) < 1e-3][0]
    ic = dict(n=len(T), innov_exceed_k=int((T.e > T.c).sum()), innov_exceed=float((T.e > T.c).mean()), innov_exceed_ci=nf.cp(int((T.e > T.c).sum()), len(T)), z_std=float((T.e / T.s).std()), z_mean=float((T.e / T.s).mean()),
              lnQ_over_Qedge_median=float(np.log(T.q / (QMAX * np.exp(-T.c))).median()), frac_above_edge=float((T.q > QMAX * np.exp(-T.c)).mean()), pos_median=float(T.pos.median()), pos_q95=float(T.pos.quantile(.95)))
    res["prediction_vs_flight"] = dict(predicted=a, flown=dict(over_all=ko / N, meanQ_ratio=n["meanQ_ratio"], low=n["frac_low"], burns_per_run=n["burns_mean"], dv_per_run=n["dv_mean"]), closed_loop=ic)
    (n2dir(root) / "analysis").mkdir(exist_ok=True); T.to_csv(n2dir(root) / "analysis" / "n2_innovations.csv", index=False)
    (n2dir(root) / "analysis" / "n2_held.json").write_text(json.dumps(res, indent=1, default=float)); print(json.dumps(res, indent=1, default=float))


def phase_repro(root):
    sel = pd.read_csv(n2dir(root) / PL / "selected.csv").iloc[0]; code = str(sel.code); out = {}
    for rd in sorted((n2dir(root) / PL / "repro").glob(f"{code}_pert_*")): out[rd.name.split("_pert_")[1]] = nf.bit_compare(rd, n2dir(root) / PL / "runs" / rd.name)
    (n2dir(root) / "analysis").mkdir(exist_ok=True); (n2dir(root) / "analysis" / "n2_repro.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out, indent=1))


def phase_integrity(root):
    out = {}
    for f in sorted(glob.glob(str(n2dir(root) / "mc_*.toml"))):
        t = tomllib.load(open(f, "rb")); out[pathlib.Path(f).name] = dict(route=t["mc_route"], failed=t["mc_failed"], samples=t["mc_samples"], host=t["host"])
    S = n2_summaries(root); runs = list((n2dir(root) / PL / "runs").glob("*/passes.csv")); nan = sum(int(pd.read_csv(p).drop(columns=["h_qhigh_m"], errors="ignore").isna().sum().sum()) for p in runs)
    ep = [pd.Timestamp(x).strftime("%Y-%-m-%-d %-H:%-M:%-S") for x in sf.CFG["planets"][PL]["grid"]["epochs_utc"]]
    out[PL] = dict(runs=len(runs), nan_cells_excl_h_qhigh=nan, rows=len(S), not_ok=int((S.status != "ok").sum()), epochs_match=bool(all(any(str(e).startswith(x) for x in ep) for e in S.gram_core_epoch)), controller_stopped=int(S.controller_stopped.sum()), not_success=int((S.retcode != "Success").sum()))
    (n2dir(root) / "analysis").mkdir(exist_ok=True); (n2dir(root) / "analysis" / "n2_integrity.json").write_text(json.dumps(out, indent=1, default=str)); print(json.dumps(out, indent=1, default=str))


def selftest():
    assert arm_code(0.01, 4000.0, 0.1, 1.5196) == "cc01_a010_t152_w04000"
    T = pd.DataFrame(dict(eligible=[True, True, True, False], g_a=[0.05, 0.1, 0.2, 0.1], frac_over=[0.01, 0.0, 0.03, 0.0], dv_mean=[0.5, 0.5, 0.5, 0.5], dv_ref=0.3, decay_median=[120.0, 110.0, 150.0, 200.0]))
    s, fb = select_arm(T); assert s.g_a == 0.05 and not fb      # 0.2 infeasible (over), the diagnostic is not eligible, best decay among feasible
    s, fb = select_arm(T.assign(frac_over=0.05)); assert fb and s.g_a == 0.2                  # none feasible: ties in lowest over -> larger g_a
    s, fb = select_arm(T.assign(decay_median=100.0, frac_over=0.0)); assert s.g_a == 0.2 and not fb
    assert nondominated([dict(name="a", dv=1, decay=1, over=0), dict(name="b", dv=2, decay=2, over=0), dict(name="c", dv=3, decay=1, over=0)]) == ["a", "b"]
    # simulate: zero steps => rho stays at -g_a, no burns beyond the first; q = 0.25 exp(-t s - g_a + s z); with z = 0 and q1 = 0.15 -> over only if exp(-t s - g_a) > 1: never
    Tt = pd.DataFrame(dict(s=[0.34], width_ln=[0.63], g0=[0.31])); r = simulate(Tt, 0.1, 1.5, np.zeros(10), np.array([0.15]), np.zeros(5), 0.56, 50, 1)
    assert r["over_all"] == 0.0 and abs(r["burns_per_run"] - 1.0) < 1e-12 and abs(r["meanQ_all"] - (0.15 + 47 * 0.25 * np.exp(-1.5 * 0.34 - 0.1)) / 48) < 1e-12, r
    # raise burn: constant step +0.08 from -0.05 crosses 0 on the 1st step (rho=0.03>0): burn of kappa*(0.03+0.05), back to -0.05, every pass
    r = simulate(Tt, 0.05, 1.5, np.zeros(10), np.array([0.15]), np.array([0.08]), 1.0, 20, 1); assert abs(r["burns_per_run"] - 48.0) < 1e-9 and abs(r["dv_per_run"] - (0.63 / 2 + 47 * 0.08)) < 1e-9, r
    T2 = solve_t(Tt, 0.1, 0.0, np.array([-1.0, 0.0, 1.0]), np.array([0.15]), np.array([0.0]), 0.5, 100, 3); print("stagen2 selftest ok")



def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--phase", required=True); ap.add_argument("--results-root"); a = ap.parse_args()
    if a.phase == "selftest": return selftest()
    r = a.results_root; {"calibrate": phase_calibrate, "select": phase_select, "held": phase_held, "repro": phase_repro, "integrity": phase_integrity}[a.phase](r)


if __name__ == "__main__": main()
