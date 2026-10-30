"""Stage O analysis (WP-O3): the periapsis-corridor keeper against the flown Mars Odyssey main phase. Plan (frozen): tasks/2026-10-15_periapsis-region-border-targeting/odyssey_comparison_plan.md
(sha256 f38bfde788ab702360489f6c0306e336e4b2d8d8679181c4a5ab8ae34ae9a68e), sections 4 to 5 and Revision 1 (R1.2 to R1.4). Reuses stagen2_analysis (Monte Carlo map, solve), stagen_flight_analysis (CP, cluster
bootstrap), stagei_analysis.run_metrics and stagef20_analysis.boot_ci unchanged.
  selftest   [--n2-inputs F] [--write-n2-inputs --n2-root R]   synthetic rows + the Stage N2 t = 1.520 reproduction (cached inputs in stagen2_data)
  gate       --fast-run D --ref-run D                          G1 integrator gate (1.0 s run against the 0.2 s run) and G2 per flight
  ens, ensn  --results-root R                                  step-size ensemble (F-O, N2-T): paired 0.5 s against 1.0 s, verdict NOISE or BIAS (stageo_ensemble_plan.md)
  calibrate  --results-root R [--n-pass N]                     re-solve t for 1% and for the iso-risk target from the dev1 (F-O, N2-T) tables -> o_calibration.json
  devcheck   --results-root R                                  dev1/dev2 tables, no criteria
  held       --results-root R                                  O0 to O7 on the held seeds (24), both tier sets, T1 to T4
  report     --results-root R                                  reported-not-gated items (KS, knots, burns per run, days with more than one burn, flown inside 5-95% range)
  integrity  --results-root R                                  route, failures, NaN, epoch, hosts, stop-rule flags
Layout read (the Stage N2 layout under <R>/stageo): summary_<phase>_*.csv with run_id = <code>_pert_<seed>, columns of the N2 summary plus stop_reason and t_end_days;
<R>/stageo/mars/runs/<run_id>/{passes,decisions,measured}.csv. Phases: dev1 (F-O, N2-T), dev2 (N2-1, N2-I), held (F-O, N2-1, N2-I; host 1), heldp (F-P, N2-P; host 2; never pooled with held).
passes.csv: the N2 columns plus rho_peak_kg_m3, rp_m, ra_m (WP-O2). Arm -> rows: the summary `arm` key column (FO, N2T, N21, N2I, FP, N2P; label F-O etc. normalized), else the run-id prefix; stop_reason end_radius = reached.
Flown reference: stageo_data/flown_reference.json (WP-O1, key "metrics"), read defensively; the plan's own numbers are the fallback (labeled), the flown decay bins have none (criterion O4 is then UNAVAILABLE).
Python only, no Julia. Usage: python3 -I stageo_analysis.py --phase ..."""
import argparse, contextlib, glob, json, math, pathlib, re, sys, tempfile, tomllib
import numpy as np, pandas as pd
from scipy.stats import ks_2samp
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import stagen2_analysis as n2
import stagen_flight_analysis as nf
import stagei_analysis as sia
import stagem_lib as sm
sf = sia.sf

# ---- constants: SOURCED = plan section / R1 table, DERIVED = from them, ASSUMED = named choice -----------------------------------------------------------------------------------------------------
MU = 42828.375e9                      # SOURCED R1.2 (km^3/s^2 converted to m^3/s^2)
UP, LO, ALLOW, QUAL = 0.32, 0.14, 0.54, 0.65   # SOURCED 2.2 (Odyssey corridor upper/lower, flight allowable, qualification)
RA_TARGET = 7.408e6                   # SOURCED Table III (R1.1)
FLOWN_PLAN = dict(passes=153, days=49.0, dv_mps=2.924, mean_q=0.2238, noise_std=0.424, iso_eps=18 / 153, hi_iso=0.183, lo_iso=0.052, hi_045=0.033, hi_054=0.007)   # SOURCED/DERIVED plan 2.1, R1.3
DECAY_GAP = 0.13                      # SOURCED 2.3 item 3
PASS_LIMIT, DAY_LIMIT = 173, 55.4     # 153 x 1.13 and 49.0 x 1.13 as the plan states them (O3)
DECAY_FRAC = 1 - DECAY_GAP            # 0.87 (O4)
NOISE_BOUNDS, NOISE_STOP = (0.8, 1.25), (0.5, 2.0)   # ASSUMED 2.6 / stop rule section 9
EPS1 = 0.01                           # N2 design point
BOOT_N, BOOT_SEED = 10000, 20261101   # ASSUMED section 5
N_PASS_CAL = 173                      # ASSUMED: passes per simulated flight in the Monte Carlo map (O3 pass limit); the N2 map used 49 for its 30 d horizon
G_A = 0.10                            # frozen N2 aim offset (2.5)
BINS = [("28.5-20", 20e6, 1e9), ("20-15", 15e6, 20e6), ("15-10", 10e6, 15e6), ("10-7.4", 7.4e6, 10e6)]   # lo < ra <= hi (as stageo_reference.py)
PL = "mars"; EPOCH_PREFIX = "2001-11-6 19:0:32"
PHASE_ARMS = dict(dev1=("F-O", "N2-T"), dev2=("N2-1", "N2-I"), held=("F-O", "N2-1", "N2-I"), heldp=("F-P", "N2-P"))


def norm(s): return re.sub(r"[^A-Za-z0-9]", "", str(s)).upper()        # arm label F-O <-> the driver's arm key FO (summary column `arm`)


def arm_rows(S, arm):
    """Summary rows of one arm: by the `arm` key column (stageo_run.jl writes FO, N2T, N21, N2I, FP, N2P), else by the run-id prefix."""
    if not len(S): return S
    return S[S.arm.map(norm) == norm(arm)] if "arm" in S else S[S.run_id.str.upper().str.startswith(norm(arm) + "_")]


def code_of(S, arm):
    """Run-directory code of an arm (the run_id before _pert_)."""
    return arm_rows(S, arm).run_id.iloc[0].split("_pert_")[0]
def sd(root): return pathlib.Path(root) / "stageo"
def adir(root): d = sd(root) / "analysis"; d.mkdir(parents=True, exist_ok=True); return d
def dump(root, name, res): (adir(root) / name).write_text(json.dumps(res, indent=1, default=float)); print(json.dumps(res, indent=1, default=float))


def task_dir():
    for p in HERE.parents:
        t = p / "tasks" / "2026-10-15_periapsis-region-border-targeting"
        if t.exists(): return t
    return None


# ---- flown reference (read defensively) ---------------------------------------------------------------------------------------------------------------------------------------------------
def load_flown(path=None):
    """Plan numbers as the fallback; flown_reference.json["metrics"] (WP-O1) overrides. Tier fractions: {threshold: (nominal, lo, hi)}. Bins: {label: (T2 median m/pass, T3 median m/day)} or {}."""
    f = dict(FLOWN_PLAN); f["source"] = "plan constants"; f["bins"] = {}; f["mean_t2"] = f["mean_t3"] = None
    f["gt"] = {0.25: (0.327, 0.235, 0.399), 0.30: (0.157, 0.078, 0.255), 0.32: (0.118, 0.052, 0.183), 0.45: (0.013, 0.007, 0.033), 0.54: (0.0, 0.0, 0.007), 0.65: (0.0, 0.0, 0.0)}
    f["lt"] = {0.14: (0.124, 0.078, 0.157), 0.05: (0.0, 0.0, 0.0)}
    if path is None and task_dir() is not None: path = task_dir() / "stageo_data" / "flown_reference.json"
    try: m = json.load(open(path))["metrics"]
    except Exception as e: f["source"] = f"plan constants (flown_reference.json unavailable: {type(e).__name__})"; return f
    f["source"] = str(pathlib.Path(path).name); g = lambda k, d: m.get(k, d)
    f["passes"], f["dv_mps"], f["mean_q"] = g("passes", f["passes"]), g("main_dv_mps", f["dv_mps"]), g("mean", f["mean_q"]); f["noise_std"] = g("lnq_std_mp", f["noise_std"]); f["iso_eps"] = g("iso_risk_eps", f["iso_eps"])
    f["days_table3"] = m.get("T1_kernel", {}).get("table3_days", f["days"]); f["days"] = f["days_table3"]; f["kernel_T1"] = m.get("T1_kernel")
    for key, src in (("gt", "gt_025_to_065_table"), ("lt", "lt_table")):
        for t, v in m.get(src, {}).items():
            try: f[key][float(t)] = (v["fraction"], v["interval"][0], v["interval"][1])
            except Exception: pass
    for k, v in m.get("T2_T3_bins", {}).items():
        nums = [float(x) for x in re.findall(r"[\d.]+", k.replace("Mm", ""))]
        lab = next((b[0] for b in BINS if (len(nums) == 1 and b[2] > 1e8 and abs(nums[0] - b[1] / 1e6) < 0.1) or (len(nums) == 2 and abs(nums[0] - b[1] / 1e6) < 0.1 and abs(nums[1] - b[2] / 1e6) < 0.1)), None)
        if lab: f["bins"][lab] = (v["T2_median_km_per_pass"] * 1e3, v["T3_median_km_per_day"] * 1e3)
    c = m.get("T2_T3_campaign_mean")
    if c: f["mean_t2"], f["mean_t3"] = c["T2_km_per_pass"] * 1e3, c["T3_km_per_day"] * 1e3
    return f


def flown_q_values(path=None):
    """Flown main-phase Q_s values (P19-171) from stageo_data/q_per_pass_flown.csv (columns P, q_flown_w_cm2 per stageo_reference.py), else None."""
    d = task_dir()
    try:
        T = pd.read_csv((d / "stageo_data" / "q_per_pass_flown.csv") if path is None else path); T = T[(T.P >= 19) & (T.P <= 171)]; return T.q_flown_w_cm2.values
    except Exception: return None


# ---- scoring heat rate and per-run metrics (R1.2, R1.4) ---------------------------------------------------------------------------------------------------------------------------------
def qs_of(rho, rp, ra):
    """Q_s = 0.5 rho_peak V_p^3 / 1e4 (W/cm2), V_p the inertial vis-viva speed at the pass periapsis radius (SI inputs)."""
    return 0.5 * np.asarray(rho, float) * (MU * (2 / np.asarray(rp, float) - 2 / (np.asarray(rp, float) + np.asarray(ra, float)))) ** 1.5 / 1e4


def apsides(P):
    if "ra_m" in P and "rp_m" in P: return P.rp_m.values.astype(float), P.ra_m.values.astype(float)
    return (P.a_m * (1 - P.e)).values, (P.a_m * (1 + P.e)).values


def tiers_odyssey(q):
    q = np.asarray(q, float); return dict(low=int((q < LO).sum()), inband=int(((q >= LO) & (q <= UP)).sum()), over=int((q > UP).sum()), allow=int((q > ALLOW).sum()), qual=int((q > QUAL).sum()))


def tiers_paper(q):
    t = sia.tiers(q); t["over"] = t["soft"] + t["medium"] + t["hard"]; t["inband"] = len(q) - t["low"] - t["over"]; return t


def run_o(P, dec, srow, ra_target=RA_TARGET):
    """Metrics of one flight. The pass set runs from the start to the first pass whose apoapsis radius is <= ra_target (inclusive); a flight that never gets there is flagged not reached (T1 = inf)."""
    P = P.sort_values("pass").reset_index(drop=True); P = P[np.isfinite(P.rho_peak_kg_m3) & np.isfinite(apsides(P)[1])].reset_index(drop=True)   # a flight cut mid-pass logs the open pass with ra = NaN
    rp, ra = apsides(P); hit = np.where(ra <= ra_target)[0]
    reached = bool(len(hit)) or str(srow.get("stop_reason", "")) in ("end_radius", "ra_target"); n = hit[0] + 1 if len(hit) else len(P)
    P = P.iloc[:n].reset_index(drop=True); rp, ra = rp[:n], ra[:n]; t = P.time_s.values; qs = qs_of(P.rho_peak_kg_m3.values, rp, ra); qe = P.q_eff_w_cm2.values
    base = sia.run_metrics(P, srow, dec); td = srow.get("t_end_days", np.nan)
    days = float(td) if reached and td is not None and np.isfinite(td) else ((t[-1] - t[0]) / 86400.0 if n > 1 else np.nan)
    bd = dec[np.abs(dec.delta_v_mps.values) > 0] if len(dec) else dec; burned = set(bd.target_pass.astype(int)) if len(bd) else set()
    k = P["pass"].values; lq = [math.log(qs[j] / qs[j - 1]) for j in range(1, n) if k[j] == k[j - 1] + 1 and k[j] not in burned and qs[j] > 0 and qs[j - 1] > 0]
    dra = ra[:-1] - ra[1:]; dd = np.diff(t) / 86400.0; t3 = dra / dd; a = 0.5 * (ra[:-1] + rp[:-1]); e = (ra[:-1] - rp[:-1]) / (ra[:-1] + rp[:-1])
    kobs = dra / np.array([sm.drop_per_pass(a[j], e[j], rp[j], qs[j] * 1e4, 1.0, MU) for j in range(len(dra))]) if len(dra) else np.array([])
    bins = {lab: dict(n=int(m.sum()), T2=float(np.median(dra[m])) if m.any() else None, T3=float(np.median(t3[m])) if m.any() else None) for lab, lo, hi in BINS for m in [(ra[:-1] > lo) & (ra[:-1] <= hi)]}
    rq = qs / qe; days_b = [int(x // 86400) for x in (bd.request_time_s.values if len(bd) else [])]
    return dict(base=base, n=n, reached=reached, T1_passes=n if reached else math.inf, T1_days=days if reached else math.inf, qs=qs, qe=qe, ratio=rq, lnq=lq, bins=bins, ra0=float(ra[0]), ra_end=float(ra[-1]),
                T2_mean=float(dra.mean()) if len(dra) else np.nan, T3_camp=float((ra[0] - ra[-1]) / ((t[-1] - t[0]) / 86400.0)) if n > 1 else np.nan, T4_median=float(np.median(kobs)) if len(kobs) else np.nan,
                dv=float(srow.get("total_dv_mps", np.nan)), burns=base["burns"], q_first=float(qe[0]), tiers_o=tiers_odyssey(qs), tiers_p=tiers_paper(qs), qe_over=int((qe > UP).sum()),
                any_gt054=bool((qs > ALLOW).any()), any_gt045=bool((qs > 0.45).any()), multi_burn_days=len(days_b) - len(set(days_b)), bad=bool(srow.get("retcode") not in ("Success", "Terminated") or srow.get("controller_stopped", False) or srow.get("status", "ok") != "ok"),
                stop_reason=str(srow.get("stop_reason", "")))


# ---- loading -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def summaries(root, phase):
    fs = sorted(glob.glob(str(sd(root) / f"summary_{phase}_*.csv")))
    return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True).drop_duplicates("run_id", keep="last") if fs else pd.DataFrame()


def load_arm(root, phase, arm, S=None):
    """{seed tag: run_o dict} of every complete run of `arm` in `phase`."""
    S = summaries(root, phase) if S is None else S; out = {}
    for r in arm_rows(S, arm).itertuples():
        tag = r.run_id.split("_pert_", 1)[1]; d = sd(root) / PL / "runs" / r.run_id
        if (d / "passes.csv").exists(): out[tag] = run_o(pd.read_csv(d / "passes.csv"), pd.read_csv(d / "decisions.csv") if (d / "decisions.csv").exists() else pd.DataFrame(), S[S.run_id == r.run_id].iloc[-1].to_dict())
    return out


def common_seeds(arms):
    """Missing-data rule: a seed missing from any arm that has data is dropped from every arm."""
    have = [set(r) for r in arms.values() if r]; s = set.intersection(*have) if have else set()
    return {a: {t: r[t] for t in sorted(s)} for a, r in arms.items()}, {a: sorted(set(r) - s) for a, r in arms.items() if r}


# ---- statistics ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def frac_block(k, n, tags, rng):
    """Pooled fraction with Clopper-Pearson 90% and the cluster bootstrap (by seed) 90%."""
    k, n = np.asarray(k, float), np.asarray(n, float); K, N = int(k.sum()), int(n.sum())
    return dict(k=K, n=N, frac=K / N if N else np.nan, cp90=nf.cp(K, N), cluster90=nf.cluster_boot(np.array(tags), k, n, rng, BOOT_N) if N else [np.nan, np.nan])


def med(x): x = [v for v in x if v is not None and not (isinstance(v, float) and np.isnan(v))]; return float(np.median(x)) if x else None


def arm_table(runs, rng):
    """Everything reported for one arm: both tier sets (pooled fractions + CP90 + cluster90), T1-T4, dV, noise std, first-pass Q."""
    tags = sorted(runs); R = [runs[t] for t in tags]; n = [r["n"] for r in R]; col = lambda s, k: [r[s][k] for r in R]
    o = {k: frac_block(col("tiers_o", k), n, tags, rng) for k in ("low", "over", "allow", "qual")}; p = {k: frac_block(col("tiers_p", k), n, tags, rng) for k in ("low", "soft", "medium", "hard", "over")}
    lq = np.concatenate([r["lnq"] for r in R]) if R else np.array([])
    return dict(n_runs=len(R), passes_total=int(sum(n)), odyssey=o, paper=p, qeff_over_032=frac_block([r["qe_over"] for r in R], n, tags, rng), seeds_gt054=int(sum(r["any_gt054"] for r in R)), seeds_gt045=int(sum(r["any_gt045"] for r in R)),
                T1_passes_median=med([r["T1_passes"] for r in R]), T1_days_median=med([r["T1_days"] for r in R]), not_reached=int(sum(not r["reached"] for r in R)), dv_median=med([r["dv"] for r in R]), dv_mean=float(np.nanmean([r["dv"] for r in R])),
                burns_mean=float(np.nanmean([r["burns"] for r in R])), meanQs=float(np.mean(np.concatenate([r["qs"] for r in R]))), T2_mean_median=med([r["T2_mean"] for r in R]), T3_camp_mean=float(np.nanmean([r["T3_camp"] for r in R])),
                T4_median=med([r["T4_median"] for r in R]), bins={lab: dict(T2=med([r["bins"][lab]["T2"] for r in R]), T3=med([r["bins"][lab]["T3"] for r in R]), n_seeds=int(sum(r["bins"][lab]["n"] > 0 for r in R))) for lab, _, _ in BINS},
                noise_std=float(lq.std(ddof=1)) if len(lq) > 2 else np.nan, noise_pairs=int(len(lq)), ratio_qs_over_qeff_median=float(np.median(np.concatenate([r["ratio"] for r in R]))), q_first_range=[float(min(r["q_first"] for r in R)), float(max(r["q_first"] for r in R))],
                bad_runs=int(sum(r["bad"] for r in R)), burns_per_day_gt1=int(sum(r["multi_burn_days"] for r in R)))


def paired_vs(runs, ref, rng):
    """Per-seed paired differences to the reference arm: over-0.32 fraction (Q_s), over-0.25 fraction (Q_s), T3 campaign mean, dV; bootstrap 90% of the mean."""
    tags = sorted(set(runs) & set(ref)); out = {}
    for name, f in (("d_over_032", lambda r: r["tiers_o"]["over"] / r["n"]), ("d_over_025", lambda r: r["tiers_p"]["over"] / r["n"]), ("d_T3_camp", lambda r: r["T3_camp"]), ("d_dv", lambda r: r["dv"]), ("d_T1_passes", lambda r: r["T1_passes"])):
        d = np.array([f(runs[t]) - f(ref[t]) for t in tags], float); d = d[np.isfinite(d)]; out[name] = dict(mean=float(d.mean()) if len(d) else None, ci90=list(sf.boot_ci(d, rng, BOOT_N)) if len(d) > 1 else None)
    return out


def seed_range(vals, flown):
    v = np.asarray([x for x in vals if np.isfinite(x)], float)
    if not len(v) or flown is None: return None
    lo, hi = np.percentile(v, [5, 95]); return dict(p5=float(lo), p95=float(hi), flown=float(flown), inside=bool(lo <= flown <= hi))


def tf(x): return None if x is None else bool(x)


# ---- criteria -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def criteria(T, R, fl):
    """T: arm tables (F-O, N2-1, N2-I), R: arm runs (seed-aligned). Pass values are True/False, None when an input is unavailable."""
    c = {}; fo, n1, ni = T.get("F-O"), T.get("N2-1"), T.get("N2-I")
    if fo:
        r = fo["noise_std"] / fl["noise_std"]; c["O0"] = dict(ratio=r, pooled_std=fo["noise_std"], flown_std=fl["noise_std"], pairs=fo["noise_pairs"], passes=bool(NOISE_BOUNDS[0] <= r <= NOISE_BOUNDS[1]), stop_rule_outside_0p5_2p0=bool(not NOISE_STOP[0] <= r <= NOISE_STOP[1]),
                                  note="if False, O3 and O4 are reported as noise-mismatched, not as verdicts")
    if n1:
        hi_o, ci = fl["gt"][UP][1], n1["odyssey"]["over"]["cp90"]; c["O1"] = dict(cp90_over_032=ci, flown_interval_low=hi_o, zero_gt_054=n1["odyssey"]["allow"]["k"] == 0, passes=bool(ci[1] < hi_o and n1["odyssey"]["allow"]["k"] == 0))
        q = n1["qeff_over_032"]; c["O2"] = dict(cp90_qeff_over_032=q["cp90"], cluster90=q["cluster90"], frac=q["frac"], passes=bool(q["cp90"][0] <= EPS1 <= q["cp90"][1]))
    if ni:
        c["O3"] = dict(T1_passes_median=ni["T1_passes_median"], T1_days_median=ni["T1_days_median"], limits=[PASS_LIMIT, DAY_LIMIT], not_reached=ni["not_reached"], passes=bool(ni["T1_passes_median"] <= PASS_LIMIT and ni["T1_days_median"] <= DAY_LIMIT))
        if fl["bins"]:
            rows = {}
            for lab, _, _ in BINS:
                b, fb = ni["bins"][lab], fl["bins"].get(lab)
                rows[lab] = dict(T2=b["T2"], T3=b["T3"], flown_T2=fb[0] if fb else None, flown_T3=fb[1] if fb else None,
                                 passes=None if fb is None else bool(b["T2"] is not None and b["T3"] is not None and b["T2"] >= DECAY_FRAC * fb[0] and b["T3"] >= DECAY_FRAC * fb[1]))
            ok = [v["passes"] for v in rows.values()]; c["O4"] = dict(bins=rows, factor=DECAY_FRAC, passes=None if any(v is None for v in ok) else bool(all(ok)))
        else: c["O4"] = dict(passes=None, note="UNAVAILABLE: flown decay-bin medians missing from the flown reference")
        o, p = ni["odyssey"]["over"], ni["paper"]["hard"]; q45 = frac_block([r["qs"].__gt__(0.45).sum() for r in R["N2-I"].values()], [r["n"] for r in R["N2-I"].values()], sorted(R["N2-I"]), np.random.default_rng(BOOT_SEED))
        c["O5"] = dict(frac_over_032=o["frac"], cp90=o["cp90"], flown_interval=[fl["gt"][UP][1], fl["gt"][UP][2]], zero_gt_054=ni["odyssey"]["allow"]["k"] == 0, frac_gt_045=q45["frac"], limit_045=fl["hi_045"],
                       passes=bool(o["frac"] <= fl["gt"][UP][2] and ni["odyssey"]["allow"]["k"] == 0 and q45["frac"] <= fl["hi_045"]))
        c["O6"] = dict(dv_median=ni["dv_median"], flown=fl["dv_mps"], passes=bool(ni["dv_median"] <= fl["dv_mps"]), note="manifest reference (D7); the literature totals disagree, see plan 2.3")
    if n1 and fo:
        t = sorted(set(R["N2-1"]) & set(R["F-O"])); a = np.array([R["N2-1"][s]["dv"] for s in t]); b = np.array([R["F-O"][s]["dv"] for s in t]); x = np.array([R["N2-1"][s]["T3_camp"] for s in t]); y = np.array([R["F-O"][s]["T3_camp"] for s in t])
        c["O7"] = dict(dv_mean=float(a.mean()), dv_fixed=float(b.mean()), dv_ratio=float(a.mean() / b.mean()) if b.mean() > 0 else None, T3_mean=float(x.mean()), T3_fixed=float(y.mean()), T3_ratio=float(x.mean() / y.mean()), n_paired=len(t),
                       passes=bool(a.mean() <= 2.0 * b.mean() and x.mean() >= 0.90 * y.mean()))
    return c


# ---- phases ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def collect(root, phases, arms_wanted=None):
    runs = {}
    for ph in phases:
        S = summaries(root, ph)
        for arm in PHASE_ARMS[ph]:
            if arms_wanted and arm not in arms_wanted: continue
            r = load_arm(root, ph, arm, S)
            if r: runs[arm] = r
    return runs


def phase_held(root, flown=None):
    rng = np.random.default_rng(BOOT_SEED); fl = flown or load_flown(); runs, dropped = common_seeds(collect(root, ("held", "heldp")))
    T = {a: arm_table(r, rng) for a, r in runs.items() if r}; res = dict(flown_source=fl["source"], arms=sorted(T), seeds=sorted(next(iter(runs.values()))) if runs else [], dropped_seeds=dropped, tables=T,
                                                                            paired_vs_FO={a: paired_vs(runs[a], runs["F-O"], rng) for a in runs if a != "F-O" and "F-O" in runs}, criteria=criteria(T, runs, fl))
    res["verdict"] = {k: v.get("passes") for k, v in res["criteria"].items()}; res["flown"] = {k: fl[k] for k in ("passes", "days", "dv_mps", "mean_q", "noise_std", "iso_eps")}
    dump(root, "o_held.json", res); return res


def phase_devcheck(root, flown=None):
    rng = np.random.default_rng(BOOT_SEED); fl = flown or load_flown(); runs = collect(root, ("dev1", "dev2")); T = {a: arm_table(r, rng) for a, r in runs.items()}; out = dict(flown_source=fl["source"], tables=T)
    if "F-O" in T: out["O0_dev"] = T["F-O"]["noise_std"] / fl["noise_std"]
    if "N2-1" in T: q = T["N2-1"]["qeff_over_032"]; out["N2-1_target_1pct"] = dict(frac=q["frac"], cp90=q["cp90"], contains=bool(q["cp90"][0] <= EPS1 <= q["cp90"][1]))
    if "N2-I" in T: o = T["N2-I"]["odyssey"]["over"]; out["N2-I_target_iso"] = dict(frac_Qs=o["frac"], cp90=o["cp90"], target=fl["iso_eps"], contains=bool(o["cp90"][0] <= fl["iso_eps"] <= o["cp90"][1]))
    for a in ("N2-T", "N2-1", "N2-I"):
        if a in runs:
            with band():
                Tb = n2.flight_table(root, sorted(runs[a]), sd(root), code_of(summaries(root, "dev1" if a == "N2-T" else "dev2"), a)); out.setdefault("innovation_exceedance", {})[a] = dict(n=len(Tb), frac=float((Tb.e > Tb.c).mean()), burned_frac=float(Tb.burned.mean()))
    dump(root, "o_dev.json", out); return out


@contextlib.contextmanager
def band(qmax=UP, qlow=LO):
    """Run the Stage N2 map (flight_table, simulate) on the Odyssey band: n2 reads QMAX and QLOW as module globals."""
    old = n2.QMAX, n2.QLOW; n2.QMAX, n2.QLOW = qmax, qlow
    try: yield
    finally: n2.QMAX, n2.QLOW = old


def solve_t(Td, g_a, target, z, q1, inc, kappa, n_sim, seed, n_pass, lo=-2.0, hi=4.0):
    """The Stage N2 solve_t bisection (same simulate, same 30 steps) with wider bounds, for targets above 1%, and n_pass passes per flight."""
    for _ in range(30):
        mid = (lo + hi) / 2
        if n2.simulate(Td, g_a, mid, z, q1, inc, kappa, n_sim, seed, n_pass)["over_all"] > target: lo = mid
        else: hi = mid
    return (lo + hi) / 2


def calib_core(Td, z, q1, g_a, target, n_pass, ratio=1.0, n_sim=None, seed=None):
    """t for P(Q > edge-threshold) = target on the decided-pass table Td. ratio = Q_s / q_eff (1 = the keeper's own variable): over(Q_s > edge) <=> over(q_eff > edge / ratio) is the plain map at t' = t - ln(ratio)/s, so t = t' + ln(ratio)/s."""
    inc, kap = n2.coast_increments(Td), n2.burn_kappa(Td)[0]; tp = solve_t(Td, g_a, target, z, q1, inc, kap, n_sim or n2.SP["mc_draws"], seed or n2.SP["mc_seed"], n_pass)
    return tp + math.log(ratio) / float(Td.s.median()), tp


def phase_calibrate(root, n_pass=N_PASS_CAL, flown=None):
    fl = flown or load_flown(); runs = collect(root, ("dev1",)); assert "N2-T" in runs and "F-O" in runs, "calibrate needs the dev1 F-O and N2-T runs"
    tags = sorted(runs["N2-T"])
    with band(): Td = n2.flight_table(root, tags, sd(root), code_of(summaries(root, "dev1"), "N2-T"))
    z = (Td.e / Td.s).values; q1 = np.array([r["q_first"] for r in runs["F-O"].values()]); ratio = float(np.median(np.concatenate([r["ratio"] for a in ("F-O", "N2-T") for r in runs[a].values()])))
    t1, t1p = calib_core(Td, z, q1, G_A, EPS1, n_pass); ti, tip = calib_core(Td, z, q1, G_A, fl["iso_eps"], n_pass, ratio)
    with band(): inc, kap = n2.coast_increments(Td), n2.burn_kappa(Td)[0]; v = n2.simulate(Td, G_A, 1.5196033513639122, z, q1, inc, kap, n2.SP["mc_draws"], n2.SP["mc_seed"], n_pass)
    nT = arm_table(runs["N2-T"], np.random.default_rng(BOOT_SEED))
    res = dict(n_dev_decided=len(Td), n_pass_map=n_pass, g_a=G_A, ratio_Qs_over_qeff_median=ratio, ratio_note="Q_s / q_eff measured on the dev1 passes (F-O and N2-T pooled); assumed constant across arms",
               iso_eps_Qs=fl["iso_eps"], iso_q_eff_threshold=UP / ratio, t_1pct=t1, t_1pct_plain=t1p, t_iso=ti, t_iso_plain=tip, t_1pct_rounded=round(t1, 3), t_iso_rounded=round(ti, 3), s_median=float(Td.s.median()), z_source="standardized innovations e/s of the dev N2-T decided passes",
               map_at_N2_t=v, N2T_realized=dict(qeff_over_032=nT["qeff_over_032"], Qs_over_032=nT["odyssey"]["over"], meanQs=nT["meanQs"]), note="t_iso is the q_eff-variable edge multiplier whose Q_s violation rate is fl.iso_eps; logged per plan R1.2")
    (sd(root) / PL).mkdir(parents=True, exist_ok=True); pd.DataFrame(dict(arm=["N2-1", "N2-I"], g_a=[G_A, G_A], t_mult=[round(t1, 3), round(ti, 3)])).to_csv(sd(root) / PL / "calibrated_arms.csv", index=False)   # read by stageo_run.jl (dev2, held, repro, xhost)
    dump(root, "o_calibration.json", res); return res


def phase_gate(fast, ref):
    """G1: `fast` = the 1.0 s run directory, `ref` = the 0.2 s run directory (identical otherwise)."""
    def ld(d):
        d = pathlib.Path(d); P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); dec = pd.read_csv(d / "decisions.csv") if (d / "decisions.csv").exists() else pd.DataFrame(columns=["target_pass", "delta_v_mps"])
        rp, ra = apsides(P); q = qs_of(P.rho_peak_kg_m3.values, rp, ra) if "rho_peak_kg_m3" in P else P.q_eff_w_cm2.values
        burned = {int(k) for k, v in zip(dec.target_pass, dec.delta_v_mps) if abs(v) > 0}; return P["pass"].values, q, ra, burned
    ka, qa, ra_a, ba = ld(fast); kb, qb, ra_b, bb = ld(ref); n = min(len(ka), len(kb)); first = next((int(ka[j]) for j in range(n) if (ka[j] in ba) != (kb[j] in bb)), None)
    m = n if first is None else int(np.searchsorted(ka[:n], first)); da, db = (ra_a[:m - 1] - ra_a[1:m]), (ra_b[:m - 1] - ra_b[1:m])
    rq = float(np.sqrt(np.mean(((qa[:m] - qb[:m]) / qb[:m]) ** 2))) if m else np.nan; rd = float(np.sqrt(np.mean((da - db) ** 2)) / np.sqrt(np.mean(db ** 2))) if m > 1 else np.nan
    dp, dq = abs(len(ka) - len(kb)) / len(kb), abs(np.nanmean(qa) - np.nanmean(qb)) / np.nanmean(qb)   # nanmean: a flight cut by the horizon has a last open pass with ra = NaN
    res = dict(passes_fast=len(ka), passes_ref=len(kb), first_different_burn_decision_pass=first, passes_before_divergence=m, a_rel_Q_rms=rq, a_drop_rms_rel=rd, b_passes_rel_diff=dp, b_meanQ_rel_diff=dq,
               G1_a=bool(rq <= 0.01 and rd <= 0.01), G1_b=bool(dp <= 0.02 and dq <= 0.02), thresholds=dict(a=0.01, b=0.02), note="drop rms is rms(difference)/rms(reference drop); 'different burn decision' = burn versus no burn at the same target pass")
    res["G1"] = res["G1_a"] and res["G1_b"]
    res["G2"] = {nm: noise_check(d) for nm, d in (("fast", fast), ("ref", ref))}   # G2: pass-to-pass ln Q_s std over the flown 0.424, per flight
    print(json.dumps(res, indent=1, default=float)); return res


def run_dir_o(d, srow=None):
    """run_o of one flight directory (passes.csv, decisions.csv); srow defaults to a clean row (stop_reason unknown, resolved from the apoapsis radius)."""
    d = pathlib.Path(d); dec = pd.read_csv(d / "decisions.csv") if (d / "decisions.csv").exists() else pd.DataFrame(columns=["target_pass", "delta_v_mps", "request_time_s"])
    return run_o(pd.read_csv(d / "passes.csv"), dec, srow if srow is not None else dict(stop_reason="", retcode="Success", status="ok"))


def noise_check(d, srow=None):
    """G2 for one flight: std of ln(Q_s,k / Q_s,k-1) over consecutive pairs without a burn at the later pass, divided by the flown 0.424; in-window = O0 window [0.8, 1.25], inside-stop-rule = [0.5, 2.0]."""
    r = run_dir_o(d, srow); sdv = float(np.std(r["lnq"], ddof=1)); ratio = sdv / FLOWN_PLAN["noise_std"]
    return dict(pairs=len(r["lnq"]), std=sdv, ratio=ratio, in_window=bool(NOISE_BOUNDS[0] <= ratio <= NOISE_BOUNDS[1]), inside_stop_rule=bool(NOISE_STOP[0] <= ratio <= NOISE_STOP[1]))


# ---- step-size ensemble (stageo_ensemble_plan.md; pre-registered before any flight) ---------------------------------------------------------------------------------------------------------
ENS_TOL = dict(mean_Qs=0.02, passes=0.01, days=0.02, T2_mean=0.02)   # ASSUMED, grounded in the G1 limits (Q 2%, passes 2% tightened to 1%, days and decay 2%); see the ensemble plan
ENS_T90 = {1: 6.314, 2: 2.920, 3: 2.353, 4: 2.132, 5: 2.015, 6: 1.943}   # two-sided 90% t quantile (0.95 one-sided) by degrees of freedom


def first_diff_pass(da, db):
    """First pass at which the two flights disagree on burn versus no burn (None if they never do within the shorter flight)."""
    ba = {int(k) for k, v in zip(da.target_pass, da.delta_v_mps) if abs(v) > 0}; bb = {int(k) for k, v in zip(db.target_pass, db.delta_v_mps) if abs(v) > 0}
    n = min(len(da), len(db)); ks = sorted(set(int(k) for k in da.target_pass) | set(int(k) for k in db.target_pass))
    return next((k for k in ks if (k in ba) != (k in bb)), None)


def phase_ens(root, arm="FO", fast_dt="dt0p5", ref_dt="dt1p0", phase="ens"):
    """Paired per-seed relative differences (fast_dt run against ref_dt run) in mean Q_s, passes to target, days to target and mean decay per pass; verdict NOISE or BIAS by the pre-registered rule."""
    S = summaries(root, phase); S = S[S.arm.map(norm) == norm(arm)] if len(S) else S; rows = {}
    for r in S.itertuples():
        m = re.match(r".*_pert_(ODY_s\d+)_(dt\w+)$", r.run_id)
        if m: rows.setdefault(m.group(1), {})[m.group(2)] = (r, run_dir_o(sd(root) / PL / "ensemble" / r.run_id, S[S.run_id == r.run_id].iloc[-1].to_dict()), sd(root) / PL / "ensemble" / r.run_id)
    per = {}; excluded = {}
    for seed, d in sorted(rows.items()):
        if ref_dt not in d or fast_dt not in d: excluded[seed] = "missing a step size"; continue
        (ra_, A, da), (rb_, B, db) = d[ref_dt], d[fast_dt]
        if not (A["reached"] and B["reached"]): excluded[seed] = f"not reached: ref {A['reached']} fast {B['reached']}"; continue
        rel = lambda x, y: float((x - y) / y)
        decA = pd.read_csv(da / "decisions.csv"); decB = pd.read_csv(db / "decisions.csv")
        per[seed] = dict(mean_Qs=rel(B["qs"].mean(), A["qs"].mean()), passes=rel(B["T1_passes"], A["T1_passes"]), days=rel(B["T1_days"], A["T1_days"]), T2_mean=rel(B["T2_mean"], A["T2_mean"]),
                         first_diff_burn_pass=first_diff_pass(decA, decB), passes_ref=A["T1_passes"], passes_fast=B["T1_passes"], days_ref=A["T1_days"], days_fast=B["T1_days"], meanQs_ref=float(A["qs"].mean()), meanQs_fast=float(B["qs"].mean()),
                         T2_ref=A["T2_mean"], T2_fast=B["T2_mean"], G2_ref=noise_check(da)["ratio"], G2_fast=noise_check(db)["ratio"])
    n = len(per); res = dict(arm=arm, fast=fast_dt, ref=ref_dt, seeds=list(per), excluded=excluded, per_seed=per, tolerances=ENS_TOL, metrics={})
    for k, tol in ENS_TOL.items():
        x = np.array([per[s][k] for s in per]); 
        if n < 2: res["metrics"][k] = dict(n=n, note="fewer than 2 seeds"); continue
        mean, sdv = float(x.mean()), float(x.std(ddof=1)); hw = ENS_T90[n - 1] * sdv / math.sqrt(n)
        res["metrics"][k] = dict(mean=mean, sd=sdv, t90=[mean - hw, mean + hw], within_tol=bool(abs(mean) <= tol), interval_contains_zero=bool(mean - hw <= 0 <= mean + hw), same_sign_all=bool(np.all(x > 0) or np.all(x < 0)), per_seed=x.tolist())
        res["metrics"][k]["pass"] = bool(res["metrics"][k]["within_tol"] and res["metrics"][k]["interval_contains_zero"] and not res["metrics"][k]["same_sign_all"])
    ok = n == 4 and not excluded and all(v.get("pass") for v in res["metrics"].values())
    res["verdict"] = "NOISE" if ok else ("BIAS" if n >= 2 and not excluded else f"INVALID (seeds excluded or fewer than 4: {excluded})")
    res["G2_per_seed"] = {s: dict(ref=per[s]["G2_ref"], fast=per[s]["G2_fast"]) for s in per}
    (adir(root) / f"o_ensemble_{arm}.json").write_text(json.dumps(res, indent=1, default=float)); print(json.dumps(res, indent=1, default=float)); return res


def phase_report(root, flown=None):
    rng = np.random.default_rng(BOOT_SEED); fl = flown or load_flown(); runs, _ = common_seeds(collect(root, ("held", "heldp"))); fq = flown_q_values(); out = dict(flown_source=fl["source"], arms={})
    for a, r in runs.items():
        R = list(r.values()); q = np.concatenate([x["qs"] for x in R]); T = arm_table(r, rng); e = dict(burns_per_run=T["burns_mean"], burns_flown=13, days_with_gt1_burn_total=T["burns_per_day_gt1"], flown_burns_note="diagnostic only (plan 4.1)")
        if fq is not None: e["KS_vs_flown_Qs"] = dict(zip(("D", "p"), map(float, ks_2samp(q, fq)))); e["quantiles_arm"] = np.percentile(q, [5, 25, 50, 75, 95]).tolist(); e["quantiles_flown"] = np.percentile(fq, [5, 25, 50, 75, 95]).tolist()
        e["Qs_over_qeff"] = dict(median=T["ratio_qs_over_qeff_median"], p5_p95=np.percentile(np.concatenate([x["ratio"] for x in R]), [5, 95]).tolist())
        e["flown_inside_5_95"] = {k: seed_range(v, f) for k, v, f in (("T1_passes", [x["T1_passes"] for x in R], fl["passes"]), ("T1_days", [x["T1_days"] for x in R], fl["days"]), ("dv", [x["dv"] for x in R], fl["dv_mps"]),
                                  ("meanQs", [x["qs"].mean() for x in R], fl["mean_q"]), ("frac_over_032", [x["tiers_o"]["over"] / x["n"] for x in R], fl["gt"][UP][0]), ("T3_camp", [x["T3_camp"] for x in R], fl["mean_t3"]), ("T2_mean", [x["T2_mean"] for x in R], fl["mean_t2"]))}
        e["flown_tiers"] = {f"gt_{k}": v for k, v in fl["gt"].items()} | {f"lt_{k}": v for k, v in fl["lt"].items()}; e["occupancy_vs_knots"] = "UNAVAILABLE: knots not in the flown reference (flat bounds are the primary definition, plan section 3)"; e["T4_median_kobs"] = T["T4_median"]; out["arms"][a] = e
    dump(root, "o_report.json", out); return out


def phase_integrity(root):
    out = {}
    for f in sorted(glob.glob(str(sd(root) / "mc_*.toml"))):
        t = tomllib.load(open(f, "rb")); out[pathlib.Path(f).name] = dict(route=t["mc_route"], failed=t["mc_failed"], samples=t["mc_samples"], host=t.get("host"), route_is_process=bool(t["mc_route"] == "process"))
    S = pd.concat([summaries(root, ph).assign(phase=ph) for ph in PHASE_ARMS], ignore_index=True) if any(len(summaries(root, ph)) for ph in PHASE_ARMS) else pd.DataFrame()
    runs = list((sd(root) / PL / "runs").glob("*/passes.csv")); nan = sum(int(pd.read_csv(p).drop(columns=["h_qhigh_m"], errors="ignore").isna().sum().sum()) for p in runs)
    out[PL] = dict(runs=len(runs), nan_cells_excl_h_qhigh=nan, rows=len(S))
    if len(S):
        out[PL].update(not_ok=int((S.status != "ok").sum()), not_success=int((~S.retcode.isin(["Success", "Terminated"])).sum()), controller_stopped=int(S.controller_stopped.sum()), epochs_match=bool(all(str(e).startswith(EPOCH_PREFIX) for e in S.gram_core_epoch)) if "gram_core_epoch" in S else None,
                       hosts_by_phase={ph: sorted(set(g.host.astype(str))) for ph, g in S.groupby("phase")} if "host" in S else None, stop_reasons=S.stop_reason.value_counts().to_dict() if "stop_reason" in S else None,
                       held_and_heldp_share_a_host=bool(set(S[S.phase == "held"].host) & set(S[S.phase == "heldp"].host)) if "host" in S else None)
        q1 = [pd.read_csv(sd(root) / PL / "runs" / r / "passes.csv").q_eff_w_cm2.iloc[0] for r in S.run_id if (sd(root) / PL / "runs" / r / "passes.csv").exists()]
        out[PL]["first_pass_q_range"] = [float(min(q1)), float(max(q1))] if q1 else None; out[PL]["stop_rule_first_pass_outside_0p05_0p54"] = bool(q1 and (min(q1) < 0.05 or max(q1) > 0.54))
    dump(root, "o_integrity.json", out); return out


# ---- selftest ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def export_n2_inputs(n2_root, out):
    """Cache the Stage N2 calibration inputs (dev decided-pass table, reference z, undecided first-pass q) from a Stage N2 results root."""
    dev = n2.sla.dev_cells(n2_root, n2.PL); base = pathlib.Path(n2_root) / "stagen_flight"; Td = n2.flight_table(n2_root, dev, base, n2.N_ARM); z = n2.ref_z(n2_root)
    dev_all = sorted(set(n2.sj.load(n2_root, "stagei", n2.PL, "fixed_pert").query("split == 'dev'").tag)); _, q1 = n2.pass1_over(n2_root, dev_all)
    json.dump(dict(Td=Td[["tag", "k", "s", "width_ln", "g0", "pos", "burned", "dv"]].to_dict(orient="list"), z=list(map(float, z)), q1=list(map(float, q1))), open(out, "w")); print("wrote", out)


def n2_reproduction(path):
    d = json.load(open(path)); Td = pd.DataFrame(d["Td"]); z, q1 = np.array(d["z"]), np.array(d["q1"]); res = {}
    with band(0.25, 0.05):
        for g, tgt, want in ((0.10, 0.01, 1.520), (0.10, 0.05, 1.012)):
            t, _ = calib_core(Td, z, q1, g, tgt, 49); res[f"g{g}_eps{tgt}"] = t; assert abs(t - want) < 1e-3, (g, tgt, t, want)
    return res


def _synth(root, phase, arm, tag, qs, dv, rng_ra=(28.5e6, 7.3e6), stop="end_radius", host="h"):
    """Write one synthetic flight: Q_s given per pass, rp 3.5 Mm, ra linear, one burn request at target pass 5."""
    n = len(qs); ra = np.linspace(rng_ra[0], rng_ra[1], n); rp = np.full(n, 3.5e6); V = np.sqrt(MU * (2 / rp - 2 / (rp + ra))); rho = 2 * np.asarray(qs) * 1e4 / V ** 3; a = 0.5 * (ra + rp)
    t = np.concatenate([[3e4], 3e4 + np.cumsum(0.5 * 2 * np.pi * np.sqrt(a[:-1] ** 3 / MU))]); rid = f"{norm(arm)}_c_pert_{tag}"; d = sd(root) / PL / "runs" / rid; d.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dict(**{"pass": np.arange(1, n + 1)}, time_s=t, apo_alt_m=ra - 3.39e6, q_eff_w_cm2=np.asarray(qs) / 1.02, rho_peak_kg_m3=rho, rp_m=rp, ra_m=ra)).to_csv(d / "passes.csv", index=False)
    pd.DataFrame(dict(request_time_s=[t[3]], target_pass=[5], delta_v_mps=[dv])).to_csv(d / "decisions.csv", index=False)
    row = pd.DataFrame([dict(run_id=rid, arm=norm(arm), status="ok", retcode="Success", controller_stopped=False, total_dv_mps=dv, burns=1, t_end_days=float((t[-1] - 3e4) / 86400), stop_reason=stop, host=host, gram_core_epoch=EPOCH_PREFIX + ".0", passes=n)])
    f = sd(root) / f"summary_{phase}_mars_x.csv"; row.to_csv(f, mode="a", header=not f.exists(), index=False)


def selftest(n2_inputs=None):
    # formula and tiers
    assert abs(qs_of(2e-8, 3.5e6, 20e6) - 0.5 * 2e-8 * (MU * (2 / 3.5e6 - 2 / 23.5e6)) ** 1.5 / 1e4) < 1e-18
    q = [0.13, 0.14, 0.32, 0.33, 0.54, 0.55, 0.65, 0.66]; assert tiers_odyssey(q) == dict(low=1, inband=2, over=5, allow=3, qual=1), tiers_odyssey(q)
    assert tiers_paper([0.04, 0.05, 0.25, 0.26, 0.30, 0.31, 0.45, 0.46]) == dict(low=1, soft=2, medium=2, hard=1, over=5, inband=2)
    assert nf.cp(0, 3600)[1] < 0.0014 and abs(nf.cp(18, 153)[0] - 0.0775) < 0.01
    # iso-risk conversion: ratio 1.1 shifts t by ln(1.1)/s relative to the plain map (toy decided-pass table)
    r0 = np.random.default_rng(2); k = np.arange(2, 42); Tt = pd.DataFrame(dict(tag="a", k=k, s=0.34, width_ln=0.63, g0=0.31, pos=np.cumsum(r0.normal(0, 0.07, len(k))) * 0.2 - 0.1, burned=(k % 13 == 0), dv=np.where(k % 13 == 0, 0.05, 0.0)))
    z = r0.normal(0, 1, 400); q1 = np.array([0.2])
    with band(): ti, tp = calib_core(Tt, z, q1, 0.1, 0.118, 60, 1.1, 4000, 5); t1, _ = calib_core(Tt, z, q1, 0.1, 0.118, 60, 1.0, 4000, 5)
    assert abs(ti - (tp + math.log(1.1) / 0.34)) < 1e-12 and abs(t1 - tp) < 1e-12 and ti > t1
    # synthetic campaign: 24 seeds x 150 passes
    rng = np.random.default_rng(7); tags = [f"ODY_s{k}" for k in range(13, 37)]; n = 150
    with tempfile.TemporaryDirectory() as tmp:
        root = pathlib.Path(tmp); lnq = {}
        for i, tg in enumerate(tags):
            q = 0.2 * np.exp(rng.normal(0, 0.3, n)); _synth(root, "held", "F-O", tg, q, 1.0); lnq[tg] = q                                      # F-O: noisy, no threshold content needed
            a = np.full(n, 0.2); a[: 2 if i < 12 else 1] = 0.4; _synth(root, "held", "N2-1", tg, a, 1.5)                                              # 36 passes over 0.32 in 3600 = 1.0% (qeff 0.392 > 0.32)
            b = np.full(n, 0.2); b[:18] = 0.4; _synth(root, "held", "N2-I", tg, b, 2.0)                                                                 # 0.12 over 0.32, none over 0.45
        _synth(root, "held", "F-O", "ODY_s99", lnq[tags[0]], 1.0)                                                                                        # seed missing from the other arms -> dropped
        # expected O0 from the arrays directly (burn at pass 5 removes the pair 4->5)
        d = np.concatenate([np.delete(np.diff(np.log(lnq[tg])), 3) for tg in tags]); want_r = d.std(ddof=1) / 0.424
        fl = load_flown(path="__none__"); r1 = run_o(pd.read_csv(sd(root) / PL / "runs" / "FO_c_pert_ODY_s13" / "passes.csv"), pd.read_csv(sd(root) / PL / "runs" / "FO_c_pert_ODY_s13" / "decisions.csv"), summaries(root, "held").iloc[0].to_dict())
        assert r1["n"] == 150 and r1["reached"] and r1["ra_end"] <= RA_TARGET and len(r1["lnq"]) == 148
        # decay bins: flown bin medians = 0.8x ours (pass) and 1.2x ours (fail the 0.87 factor)
        res0 = phase_held(root, fl); assert res0["dropped_seeds"]["F-O"] == ["ODY_s99"] and len(res0["seeds"]) == 24
        assert abs(res0["criteria"]["O0"]["ratio"] - want_r) < 1e-9 and res0["criteria"]["O4"]["passes"] is None
        ni = res0["tables"]["N2-I"]; fl["bins"] = {lab: (ni["bins"][lab]["T2"] * 0.8, ni["bins"][lab]["T3"] * 0.8) for lab, _, _ in BINS}; c = phase_held(root, fl)["criteria"]
        assert c["O4"]["passes"] is True and all(c[k]["passes"] for k in ("O1", "O2", "O3", "O5", "O6", "O7")), {k: v.get("passes") for k, v in c.items()}
        assert abs(c["O5"]["frac_over_032"] - 0.12) < 1e-12 and c["O5"]["frac_gt_045"] == 0 and c["O1"]["cp90_over_032"][1] < 0.052 and c["O2"]["frac"] == 0.01 and c["O3"]["T1_passes_median"] == 150
        fl["bins"] = {lab: (ni["bins"][lab]["T2"] * 1.2, ni["bins"][lab]["T3"] * 1.2) for lab, _, _ in BINS}; assert phase_held(root, fl)["criteria"]["O4"]["passes"] is False
        h = phase_held(root, fl); assert h["tables"]["N2-I"]["paper"]["over"]["frac"] == 0.12 and h["paired_vs_FO"]["N2-I"]["d_dv"]["mean"] == 1.0 and h["tables"]["N2-I"]["odyssey"]["allow"]["k"] == 0
        # O1 fails when the flown lower interval end is below the N2-1 CP bound; O3 fails on a late flight; O7 fails on dV
        fl2 = dict(fl); fl2["gt"] = dict(fl["gt"]); fl2["gt"][UP] = (0.118, 0.005, 0.183); assert phase_held(root, fl2)["criteria"]["O1"]["passes"] is False
        _synth(root, "held", "N2-I", "ODY_s13", np.full(n + 40, 0.2), 2.0, rng_ra=(28.5e6, 9.0e6), stop="horizon"); assert not run_o(pd.read_csv(sd(root) / PL / "runs" / "N2I_c_pert_ODY_s13" / "passes.csv"), pd.DataFrame(), dict(stop_reason="horizon", t_end_days=60.0))["reached"]
        # gate: identical runs pass, a 3% mean-Q change fails (b)
        g = phase_gate(sd(root) / PL / "runs" / "FO_c_pert_ODY_s14", sd(root) / PL / "runs" / "FO_c_pert_ODY_s14"); assert g["G1"] is True and g["first_different_burn_decision_pass"] is None
        _synth(root, "held", "F-O", "slow", lnq["ODY_s14"] * 1.03, 1.0); g = phase_gate(sd(root) / PL / "runs" / "FO_c_pert_slow", sd(root) / PL / "runs" / "FO_c_pert_ODY_s14"); assert g["G1_a"] is False and g["G1_b"] is False
        phase_report(root, fl); phase_integrity(root)
    # CP, paired bootstrap: exact cases
    assert sf.boot_ci(np.full(5, -0.5), np.random.default_rng(1), 100) == (-0.5, -0.5) and nf.cp(5, 10) == nf.cp(5, 10)
    # Stage N2 reproduction on N2's own dev tables
    p = pathlib.Path(n2_inputs) if n2_inputs else (task_dir() / "stagen2_data" / "n2_calibration_inputs.json" if task_dir() else None)
    if p is not None and p.exists(): print("N2 reproduction (t, tolerance 1e-3):", n2_reproduction(p))
    else: print("N2 reproduction SKIPPED: no cached inputs (run --phase selftest --write-n2-inputs --n2-root R once)")
    print("stageo selftest ok")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--phase", required=True); ap.add_argument("--results-root"); ap.add_argument("--n-pass", type=int, default=N_PASS_CAL); ap.add_argument("--fast-run"); ap.add_argument("--ref-run")
    ap.add_argument("--n2-inputs"); ap.add_argument("--n2-root"); ap.add_argument("--write-n2-inputs", action="store_true"); ap.add_argument("--flown"); a = ap.parse_args(); r = a.results_root
    if a.phase == "selftest":
        if a.write_n2_inputs: export_n2_inputs(a.n2_root, a.n2_inputs or (task_dir() / "stagen2_data" / "n2_calibration_inputs.json"))
        return selftest(a.n2_inputs)
    if a.phase == "gate": return phase_gate(a.fast_run, a.ref_run)
    if a.phase == "ens": return phase_ens(r, "FO", phase="ens")
    if a.phase == "ensn": return phase_ens(r, "N2T", phase="ensn")
    fl = load_flown(a.flown)
    {"calibrate": lambda: phase_calibrate(r, a.n_pass, fl), "devcheck": lambda: phase_devcheck(r, fl), "held": lambda: phase_held(r, fl), "report": lambda: phase_report(r, fl), "integrity": lambda: phase_integrity(r)}[a.phase]()


if __name__ == "__main__": main()
