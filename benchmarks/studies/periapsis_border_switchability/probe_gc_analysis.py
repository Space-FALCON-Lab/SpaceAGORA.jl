"""Probe GC analysis: gravity degree and GRAM look-ahead cache vs direct 20x20 (plan: tasks/2026-10-15_periapsis-region-border-targeting/probe_gc_plan.md).

Usage: python3 -I probe_gc_analysis.py --results-root DIR [--out FILE]      reads DIR/probe_gc/{summary_main_*.csv,runs/*}
       python3 -I probe_gc_analysis.py --selftest                           synthetic data, no flights
Thresholds are the frozen ones in the plan; this script only evaluates them.
"""
import argparse, csv, glob, math, os, sys, tempfile, tomllib
import numpy as np

WINDOW_ALT_M = 250e3          # ASSUMED dense window for within-pass structure (plan section 4)
P2_ALT_RMS_M, P2_Q_REL_RMS = 150.0, 0.05
P3_RATIO = (0.8, 1.25); P3_LAG = 0.1; P1_MIN_RATIO = 1.5; D0_TOL_M = 1e-6


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def load_run(d):
    """Per-pass arrays from passes.csv and, if present, the perturbation log (kind 3 knots) and summary."""
    rows = read_csv(os.path.join(d, "passes.csv"))
    col = lambda k: np.array([float(r[k]) for r in rows])
    run = dict(alt=col("alt_m"), q=col("q_eff_w_cm2"), apo=col("apo_alt_m"), r=col("r_pert"), t=col("time_s"))
    log = os.path.join(d, "gram_density_perturbation_log.csv")
    run["knots"] = None; run["base_seeds"] = None
    if os.path.isfile(log):
        k = []
        with open(log) as f:
            next(f)
            for line in f:
                s = line.split(",")
                if s[0] == "3":
                    k.append((int(s[2]), float(s[3]), float(s[4]), float(s[5])))
        run["knots"] = np.array(k) if k else None
        if os.path.isfile(log + ".summary.toml"):
            with open(log + ".summary.toml", "rb") as f:
                run["base_seeds"] = tomllib.load(f).get("base_seeds")
    return run


def decay(apo):
    return apo[:-1] - apo[1:]


def diffstats(a, b, rel=False):
    n = min(len(a), len(b)); d = (a[:n] - b[:n]) / (b[:n] if rel else 1.0)
    return dict(n=n, mean=float(np.mean(d)), rms=float(np.sqrt(np.mean(d * d))), max=float(np.max(np.abs(d))))


def lag1(x):
    x = np.asarray(x, float); x = x - x.mean()
    den = float(np.sum(x * x))
    return float(np.sum(x[:-1] * x[1:]) / den) if len(x) > 2 and den > 0 else float("nan")


def pert_stats(run):
    """ln r at periapsis per pass, pooled within-pass variance/lag-1 on the dense window, variance of per-pass decay."""
    lr = np.log(run["r"])
    out = dict(var_peri=float(np.var(lr, ddof=1)), lag1_peri=lag1(lr), var_decay=float(np.var(decay(run["apo"]), ddof=1)))
    kn = run["knots"]
    if kn is not None:
        ss = 0.0; dof = 0; num = 0.0; den = 0.0
        for p in np.unique(kn[:, 0]):
            m = (kn[:, 0] == p) & (kn[:, 2] <= WINDOW_ALT_M)
            x = np.log(kn[m, 3])
            if len(x) < 3: continue
            x = x - x.mean(); ss += float(np.sum(x * x)); dof += len(x) - 1
            num += float(np.sum(x[:-1] * x[1:])); den += float(np.sum(x * x))
        out["var_within"] = ss / dof if dof else float("nan"); out["lag1_within"] = num / den if den else float("nan")
    return out


def trend_std(alt):
    t = np.arange(len(alt)); c = np.polyfit(t, alt, 1)
    return float(np.std(alt - np.polyval(c, t), ddof=2))


def wall_per_pass(row):
    return float(row["run_wall_s"]) / max(int(row["passes"]), 1)


def ratio_ok(x, lo, hi):
    return bool(lo <= x <= hi)


def analyze(runs, summ):
    """runs: run_id -> run dict; summ: run_id -> summary row. Returns (report dict, text)."""
    rep = {}; L = []
    key = lambda c, d, a, r, s="": f"{c}{d}_{a}_r{r}{s}"
    cfgs = [("direct", 20), ("cache", 20), ("direct", 10), ("cache", 10)]
    # D0 determinism
    d0 = {}
    for c, d in cfgs:
        for a in ("mean", "pert"):
            k1, k2 = key(c, d, a, 1), key(c, d, a, 2)
            if k1 in runs and k2 in runs:
                n = min(len(runs[k1]["alt"]), len(runs[k2]["alt"])); d0[f"{c}{d}_{a}"] = float(np.max(np.abs(runs[k1]["alt"][:n] - runs[k2]["alt"][:n])))
    rep["D0_max_alt_diff_m"] = d0; rep["D0_pass"] = bool(d0) and all(v <= D0_TOL_M for v in d0.values())
    L.append(f"D0 replicate determinism, max |alt r1 - alt r2| (m): {d0}  -> {'PASS' if rep['D0_pass'] else 'FAIL (repeat spread is the noise floor)'}")
    # S1
    seeds = {k: runs[k]["base_seeds"] for k in runs if "_pert_r1" in k and "_alt" not in k and runs[k]["base_seeds"] is not None}
    rep["S1_pass"] = len({str(v) for v in seeds.values()}) == 1 and len(seeds) > 0
    L.append(f"S1 base_seeds equal across perturbed primary flights: {rep['S1_pass']} ({seeds})")
    base = {a: key("direct", 20, a, 1) for a in ("mean", "pert")}
    # M1 / P1
    p1 = {}; allfaster = True; ratios = []
    for a in ("mean", "pert"):
        for d in (20, 10):
            for r in (1, 2):
                kd, kc = key("direct", d, a, r), key("cache", d, a, r)
                if kd in summ and kc in summ:
                    ratio = wall_per_pass(summ[kd]) / wall_per_pass(summ[kc]); p1[f"{a}{d}_r{r}"] = ratio; ratios.append(ratio); allfaster &= ratio > 1.0
    med = float(np.median(ratios)) if ratios else float("nan")
    rep["P1"] = dict(ratios=p1, median=med, cache_faster_everywhere=bool(allfaster), worthwhile=bool(allfaster and med >= P1_MIN_RATIO))
    L.append(f"P1 wall/pass direct/cache: {p1} median {med:.3f}; faster everywhere {allfaster}; worthwhile (>= {P1_MIN_RATIO}): {rep['P1']['worthwhile']}")
    wp = {k: wall_per_pass(v) for k, v in summ.items() if v.get("warm", "false") in ("false", "False")}
    L.append("M1 wall s per pass: " + ", ".join(f"{k}={v:.3f}" for k, v in sorted(wp.items())))
    L.append("M2 steps (n_sol_t, nf, naccept, density_calls): " + ", ".join(f"{k}=({v['n_sol_t']},{v['nf']},{v['naccept']},{v['density_calls']})" for k, v in sorted(summ.items())))
    # M3-M5, P2, P4
    P = {}
    for a in ("mean", "pert"):
        if base[a] not in runs: continue
        b = runs[base[a]]
        for c, d in cfgs[1:]:
            k = key(c, d, a, 1)
            if k not in runs: continue
            x = runs[k]
            P[f"{c}{d}_{a}"] = dict(alt_m=diffstats(x["alt"], b["alt"]), q_rel=diffstats(x["q"], b["q"], rel=True), decay_m=diffstats(decay(x["apo"]), decay(b["apo"])))
    rep["metrics_vs_direct20"] = P
    for k, v in P.items(): L.append(f"M3-5 {k}: alt {v['alt_m']}  q_rel {v['q_rel']}  decay {v['decay_m']}")
    k = "cache20_mean"
    if k in P:
        rep["P2"] = dict(alt_rms_m=P[k]["alt_m"]["rms"], q_rel_rms=P[k]["q_rel"]["rms"], passed=bool(P[k]["alt_m"]["rms"] <= P2_ALT_RMS_M and P[k]["q_rel"]["rms"] <= P2_Q_REL_RMS))
        L.append(f"P2 (cache20 vs direct20, mean arm): alt rms {rep['P2']['alt_rms_m']:.2f} m (<= {P2_ALT_RMS_M}), q rel rms {rep['P2']['q_rel_rms']:.4f} (<= {P2_Q_REL_RMS}) -> {'PASS' if rep['P2']['passed'] else 'FAIL'}")
    if "direct10_mean" in P and base["mean"] in runs:
        ref = trend_std(runs[base["mean"]]["alt"]); rep["P4"] = dict(alt_rms_m=P["direct10_mean"]["alt_m"]["rms"], trend_std_m=ref, ratio=P["direct10_mean"]["alt_m"]["rms"] / ref if ref > 0 else float("nan"))
        L.append(f"P4 (10x10 vs 20x20, direct, mean; no gate): alt rms {rep['P4']['alt_rms_m']:.2f} m vs trend-removed std of direct20 periapsis altitude {ref:.2f} m (ratio {rep['P4']['ratio']:.3f}); q_rel {P['direct10_mean']['q_rel']}")
    # M6 / P3
    S = {k: pert_stats(runs[k]) for k in runs if "_pert_" in k}
    rep["pert_stats"] = S
    for k, v in sorted(S.items()): L.append(f"M6 {k}: {v}")
    kd, kc, ka, kac = key("direct", 20, "pert", 1), key("cache", 20, "pert", 1), key("direct", 20, "pert", 1, "_alt"), key("cache", 20, "pert", 1, "_alt")
    if kd in S and kc in S:
        rat = lambda a, b, f: S[a][f] / S[b][f] if S[b][f] else float('nan')
        p3 = dict(var_ratio_peri=rat(kc, kd, "var_peri"), lag1_change_peri=S[kc]["lag1_peri"] - S[kd]["lag1_peri"], var_ratio_decay=rat(kc, kd, "var_decay"))
        if "var_within" in S[kd]: p3.update(var_ratio_within=rat(kc, kd, "var_within"), lag1_change_within=S[kc]["lag1_within"] - S[kd]["lag1_within"])
        ok = ratio_ok(p3["var_ratio_peri"], *P3_RATIO) and abs(p3["lag1_change_peri"]) <= P3_LAG
        if "var_ratio_within" in p3: ok = ok and ratio_ok(p3["var_ratio_within"], *P3_RATIO) and abs(p3["lag1_change_within"]) <= P3_LAG
        else: ok = False; p3["note"] = "within-pass knots missing"
        p3["passed"] = bool(ok)
        if ka in S:   # realization spread reference: direct primary vs alternative seed
            p3["ref_seed_spread_var_ratio_peri"] = S[ka]["var_peri"] / S[kd]["var_peri"]; p3["ref_seed_spread_var_ratio_decay"] = S[ka]["var_decay"] / S[kd]["var_decay"]
        if kac in S and ka in S:
            p3["alt_seed_cache_over_direct_var_ratio_decay"] = S[kac]["var_decay"] / S[ka]["var_decay"]
        rep["P3"] = p3; L.append(f"P3 (cache20 vs direct20, perturbed): {p3} -> {'PASS' if ok else 'FAIL'}")
        if base["pert"] in runs:
            cs = {}
            for k in (key("cache", 20, "pert", 1), key("direct", 10, "pert", 1), key("cache", 10, "pert", 1)):
                if k in runs and len(runs[k]["r"]) > 2:
                    n = min(len(runs[k]["r"]), len(runs[kd]["r"])); cs[k] = float(np.corrcoef(np.log(runs[k]["r"][:n]), np.log(runs[kd]["r"][:n]))[0, 1])
            rep["S1_corr_ln_r_vs_direct20"] = cs; L.append(f"S1 correlation of per-pass ln r at periapsis vs direct20: {cs}")
    rep["usable_cache"] = bool(rep["P1"]["worthwhile"] and rep.get("P2", {}).get("passed") and rep.get("P3", {}).get("passed"))
    L.append(f"Pre-registered reading: cache usable (P1 worthwhile and P2 and P3) = {rep['usable_cache']}")
    return rep, "\n".join(L)


def load_all(root):
    pg = os.path.join(root, "probe_gc"); summ = {}
    for f in sorted(glob.glob(os.path.join(pg, "summary_*_*.csv"))):
        for r in read_csv(f):
            if r["status"] == "ok": summ[r["run_id"]] = r
    runs = {k: load_run(os.path.join(pg, "runs", k)) for k in summ if os.path.isfile(os.path.join(pg, "runs", k, "passes.csv"))}
    summ = {k: v for k, v in summ.items() if k in runs and v.get("warm") not in ("true", "True") and not k.startswith(("pilot", "calls"))}
    runs = {k: v for k, v in runs.items() if k in summ}
    return runs, summ


def selftest():
    rng = np.random.default_rng(1); n = 40
    def mk(alt0, scale, lnr, extra=0.0):
        alt = alt0 + scale * np.cumsum(rng.normal(size=n)) * 0 + extra * np.arange(n); apo = 2.5e7 - 5e3 * np.arange(n) - extra * 10 * np.arange(n) + 50 * np.sin(np.arange(n) ** 2 + 3 * extra)
        r = np.exp(lnr); q = 0.15 * (1 + 0.001 * np.arange(n)) * r
        kn = [(p, float(s), 120e3 + abs(s - 50) * 10e3, math.exp(0.1 * math.sin(0.3 * s + p))) for p in range(1, n + 1) for s in range(0, 100, 1)]
        return dict(alt=alt, q=q, apo=apo, r=r, t=np.arange(n) * 6e4, knots=np.array(kn), base_seeds=[7])
    lnr = rng.normal(0, 0.3, n); runs = {}; summ = {}
    for c, d in [("direct", 20), ("cache", 20), ("direct", 10), ("cache", 10)]:
        for a in ("mean", "pert"):
            for r in (1, 2):
                e = 0.0 if c == "direct" and d == 20 else 0.5
                x = mk(110e3, 1, lnr if a == "pert" else np.zeros(n), e); runs[f"{c}{d}_{a}_r{r}"] = x
                summ[f"{c}{d}_{a}_r{r}"] = dict(run_wall_s=str((100 if c == "direct" else 40) * (1 + 0.01 * r)), passes=str(n), n_sol_t="1", nf="1", naccept="1", density_calls="-1", warm="false")
    for c in ("direct", "cache"):
        runs[f"{c}20_pert_r1_alt"] = mk(110e3, 1, rng.normal(0, 0.3, n)); summ[f"{c}20_pert_r1_alt"] = dict(run_wall_s="50", passes=str(n), n_sol_t="1", nf="1", naccept="1", density_calls="-1", warm="false")
    rep, text = analyze(runs, summ)
    assert rep["D0_pass"] and rep["S1_pass"], text
    assert abs(rep["P1"]["median"] - 100 / 40) < 0.05, rep["P1"]
    assert rep["P1"]["worthwhile"], rep["P1"]
    assert abs(rep["P2"]["alt_rms_m"] - 0.5 * math.sqrt(np.mean(np.arange(n) ** 2))) < 1e-6 and rep["P2"]["passed"]
    assert rep["P3"]["passed"] and abs(rep["P3"]["var_ratio_peri"] - 1) < 1e-12
    assert abs(lag1(np.sin(0.2 * np.arange(200))) - math.cos(0.2)) < 0.03
    # file round trip: write one run dir and re-read it
    with tempfile.TemporaryDirectory() as t:
        d = os.path.join(t, "probe_gc", "runs", "direct20_mean_r1"); os.makedirs(d)
        with open(os.path.join(d, "passes.csv"), "w") as f:
            f.write("pass,time_s,alt_m,q_eff_w_cm2,apo_alt_m,r_pert\n" + "\n".join(f"{i+1},{i*6e4},{1.1e5+i},{0.15},{2.5e7-i*1e3},1.0" for i in range(5)) + "\n")
        with open(os.path.join(d, "gram_density_perturbation_log.csv"), "w") as f:
            f.write("kind,sat,pass,t_s,alt_m,r,sigma_frac,mean_density_kgm3,aux\n3,1,1,0.0,120000.0,0.9,0,0,0\n3,1,1,1.0,119000.0,1.1,0,0,0\n3,1,1,2.0,118000.0,1.0,0,0,0\n")
        with open(os.path.join(d, "gram_density_perturbation_log.csv.summary.toml"), "w") as f:
            f.write("base_seeds = [5]\n")
        with open(os.path.join(t, "probe_gc", "summary_main_x.csv"), "w") as f:
            f.write("run_id,status,passes,run_wall_s,warm,n_sol_t,nf,naccept,density_calls\ndirect20_mean_r1,ok,5,10.0,false,1,1,1,-1\nwarm,ok,3,1.0,true,1,1,1,-1\n")
        rr, ss = load_all(t)
        assert list(rr) == ["direct20_mean_r1"] and rr["direct20_mean_r1"]["base_seeds"] == [5] and len(rr["direct20_mean_r1"]["knots"]) == 3
    print(text); print("SELFTEST PASS")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--results-root"); ap.add_argument("--out"); ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest: selftest(); sys.exit(0)
    runs, summ = load_all(os.path.abspath(a.results_root)); rep, text = analyze(runs, summ); print(text)
    if a.out: open(a.out, "w").write(text + "\n")
