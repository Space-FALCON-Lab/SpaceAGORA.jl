"""Stage C (Mars generality), Python half: A1 oracle labels on the Mars truth, the frozen estimator's predictions scored against them,
pre-registered pass/fail table (mars_generality_plan.md), baselines, error tables, classification, dV error, quiet-set test, figures.

Usage (numpy, pandas, scipy, matplotlib):
  python3 -I analyze_mars_generality.py [--results-root DIR] [--phase labels|score]

Reads <root>/mars/{cells.csv, truth_run_summary.csv, truth/<tag>/{calibration,passes_3b}.csv, pred*/<tag>.csv, hand/<tag>.csv, pred_cells.csv}.
Reuses analyze_selection_skill.py (metrics, policy on a history) and oracle_lib.py (DP, policies, Keplerian lever) UNCHANGED; only the central GM of
the lever (oracle_lib.MU) and the margin are set from here, from the Mars config.
"""
import argparse, importlib.util, math, pathlib, sys, tomllib
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
spec = importlib.util.spec_from_file_location("analyze_selection_skill", HERE / "analyze_selection_skill.py")
an = importlib.util.module_from_spec(spec); sys.modules["analyze_selection_skill"] = an; spec.loader.exec_module(an)
lib = an.lib
MC = tomllib.load(open(HERE / "mars_generality_config.toml", "rb"))
DAY = 86400.0
Q_HARD = MC["band"]["q_hard_w_cm2"]
MARGIN = MC["analysis"]["margin_km"] * 1000.0
an.MARGIN = MARGIN                       # policy_costs in the shared module reads this global
IMPACT_ALT = MC["analysis"]["impact_altitude_m"]
an.IMPACT_ALT = IMPACT_ALT
KFRAC = MC["analysis"]["top_fraction"]
CRIT = MC["analysis"]["critical_inclinations_deg"]
EARTH_QUIET = MC["analysis"]["earth_quiet_set"]
BINS = MC["analysis"]["lookahead_bin_edges_days"]
HALVES = MC["analysis"]["raan_halves"]
METHODS = ["v2", "v2_selfcal", "v1", "B1_hand", "B2c_critical", "B2e_earth_set", "B3_random"]
# Predictor folders (default: the FROZEN estimator). `--variant j3` scores the UNFROZEN J3 diagnostic instead; its outputs carry the suffix _j3diag.
VAR = dict(v2_folder="pred", v2sc_folder="pred_selfcal", pc_file="pred_cells.csv", suffix="")


def mars_mu(root):
    pc = pd.read_csv(root / "mars" / "planet_constants.csv").set_index("planet")
    return float(pc.loc["Mars", "mu_m3s2"])


# ------------------------------------------------------------------------------------------------ labels (A1 pipeline on the Mars truth)

def build_labels(root, margin=MARGIN):
    m = root / "mars"
    cells = pd.read_csv(m / "cells.csv"); rs = pd.read_csv(m / "truth_run_summary.csv").set_index("tag")
    rows = []
    for c in cells.itertuples():
        if c.tag not in rs.index or rs.loc[c.tag, "status"] != "ok":
            rows.append(dict(tag=c.tag, inc_deg=c.inc, raan_deg=c.raan, omega_deg=c.omega, epoch=int(c.epoch_idx), status="missing")); continue
        df = pd.read_csv(m / "truth" / c.tag / "passes_3b.csv")
        retcode = rs.loc[c.tag, "retcode"]; status = "completed" if retcode == "Success" else "decayed"
        x = df.alt_m.values; lo = df.h_qhigh_m.values + margin; hi = df.h_qlow_m.values - margin
        a = df.a_m.values; e = df.e.values; L = lib.lever_kepler(a, e); q = df.q_w_cm2.values
        row = dict(tag=c.tag, inc_deg=c.inc, raan_deg=c.raan, omega_deg=c.omega, epoch=int(c.epoch_idx), status=status, retcode=retcode,
                   t_end_days=df.time_s.iloc[-1] / DAY, n_passes=len(df), a0_km=a[0] / 1e3, e0=e[0], band_width_km=float((df.h_qlow_m - df.h_qhigh_m).iloc[0] / 1e3),
                   x0_in_shrunk_band=bool(lo[0] <= x[0] <= hi[0]), passive_q_max_w_cm2=float(q.max()), passive_exceeds_hard=bool((q > Q_HARD).any()),
                   passive_in_shrunk_band_all=bool(np.all((x >= lo) & (x <= hi))), passive_alt_min_km=float(x.min() / 1e3), passive_alt_max_km=float(x.max() / 1e3),
                   passive_excursion_km=float((x.max() - x.min()) / 1e3), lever_kepler_first_mps_per_km=1000 * float(L[0]), lever_kepler_last_mps_per_km=1000 * float(L[-1]))
        d = lib.dp_schedule(lo - x, hi - x, L, None, tie=an.ORC["dp"]["tie_lambda_mps"])
        row["dp_dv"] = np.nan if d is None else d["dv"]; row["dp_burns"] = -1 if d is None else d["burns"]
        c0, n0, _ = lib.lazy(x, lo, hi, L, s=0.0); row["lazy_s0_dv"] = c0; row["lazy_s0_burns"] = n0
        cn, nn, _ = lib.near_passive(x, lo, hi, L); row["near_passive_dv"] = cn; row["near_passive_burns"] = nn
        cr, nr, _ = lib.reset_opposite(x, lo, hi, L); row["reset_opposite_dv"] = cr
        ch, nh, _ = lib.hold(x, lo, hi, L, "nearest"); row["hold_nearest_dv"] = ch; row["hold_nearest_burns"] = nh
        row["true_cost_label_mps"] = row["dp_dv"] if status == "completed" else np.nan
        row["cost_to_termination_mps"] = row["dp_dv"] if status == "decayed" else np.nan
        row["true_class"] = "decayed" if status == "decayed" else ("breach" if row["passive_exceeds_hard"] else "safe")
        # dimensionless cost: band widths of correction (dV / (lever x shrunk band width))
        w = float(np.mean(hi - lo)); row["shrunk_band_width_m"] = w
        row["cost_band_widths"] = row["dp_dv"] / (float(np.median(L)) * w) if np.isfinite(row["dp_dv"]) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ predicted tables

def build_cell_table(root, labels, margin=MARGIN):
    m = root / "mars"; pc = pd.read_csv(m / VAR["pc_file"]).set_index("tag"); pcf = pd.read_csv(m / "pred_cells.csv").set_index("tag")   # pcf: frozen v1 flags
    T = []
    for c in labels.itertuples():
        if c.status == "missing" or c.tag not in pc.index:
            continue
        r = dict(tag=c.tag, inc_deg=c.inc_deg, raan_deg=c.raan_deg, omega_deg=c.omega_deg, epoch=c.epoch, status=c.status, true_class=c.true_class,
                 true_cost_mps=c.true_cost_label_mps, true_lazy_s0_mps=c.lazy_s0_dv, true_near_passive_mps=c.near_passive_dv, true_dp_mps=c.dp_dv,
                 true_exceeds_hard=bool(c.passive_exceeds_hard), true_t_end_days=c.t_end_days, cost_band_widths=c.cost_band_widths)
        for key, folder, kind, dcol in (("v2", VAR["v2_folder"], "est", "decayed_pred"), ("v2_selfcal", VAR["v2sc_folder"], "est", "decayed_pred_selfcal"),
                                        ("v1", "pred_v1", "est", "decayed_pred_v1"), ("B1_hand", "hand", "hand", None)):
            f = m / folder / f"{c.tag}.csv"
            df = pd.read_csv(f) if f.exists() else pd.DataFrame(columns=["t_s", "alt_m", "q_w_cm2", "h_qhigh_m", "h_qlow_m", "a_m", "e"])
            row, dec = an.predicted_row(df, kind)
            if kind == "est":
                dec = bool(pcf.loc[c.tag, dcol] if key == "v1" else pc.loc[c.tag, dcol])
            for k_, v_ in row.items():
                r[f"{key}_{k_}"] = v_
            r[f"{key}_decayed"] = dec
            r[f"{key}_breach"] = bool((not dec) and row["q_max"] > Q_HARD)
            r[f"{key}_class"] = "decayed" if dec else ("breach" if r[f"{key}_breach"] else "safe")
        r["seed_truth_m"] = pc.loc[c.tag, "seed_truth_m"]; r["seed_selfcal_m"] = pc.loc[c.tag, "seed_selfcal_m"]
        r["t_decay_pred_days"] = pc.loc[c.tag, "t_decay_pred_days"]
        T.append(r)
    return pd.DataFrame(T)


def crit_flag(T):
    return np.array([any(abs(i - ci) < 0.05 for ci in CRIT) for i in T.inc_deg])


def earth_quiet_flag(T):
    return np.array([any(abs(i - qi) < 0.05 and abs(r - qr) < 0.5 for qi, qr in EARTH_QUIET) for i, r in zip(T.inc_deg, T.raan_deg)])


def splits_of(T):
    out = {"pooled": np.ones(len(T), bool)}
    for e in sorted(T.epoch.unique()):
        out[f"epoch_{e}_held_out"] = (T.epoch == e).values
    if set(T.raan_deg) <= set(HALVES[0]) | set(HALVES[1]):     # the 128-cell grid: the configured RAAN halves
        out[f"raan_{int(HALVES[0][0])}-{int(HALVES[0][-1])}"] = T.raan_deg.isin(HALVES[0]).values
        out[f"raan_{int(HALVES[1][0])}-{int(HALVES[1][-1])}"] = T.raan_deg.isin(HALVES[1]).values
    else:                                                      # off-grid RAAN values (blind sets): halves [0, 180) and [180, 360)
        out["raan_0-179"] = (T.raan_deg < 180).values; out["raan_180-359"] = (T.raan_deg >= 180).values
    return out


def evaluate(T, qmargin=0.0):
    rng = np.random.default_rng(MC["analysis"]["random_seed"]); draws = MC["analysis"]["random_draws"]
    true_c = T.true_cost_mps.values.astype(float)
    true_unsafe = ((T.status == "decayed") | T.true_exceeds_hard).values; true_dec = (T.status == "decayed").values
    flags = {"B2c_critical": crit_flag(T), "B2e_earth_set": earth_quiet_flag(T)}
    rows = []
    for sname, mask in splits_of(T).items():
        idx = np.where(mask)[0]
        for mth in METHODS:
            if mth in ("v2", "v2_selfcal", "v1", "B1_hand"):
                score = np.where(T[f"{mth}_decayed"].values, np.inf, T[f"{mth}_lazy_dv"].values.astype(float))
                pcm = ~T[f"{mth}_decayed"].values
                ps = pcm & (T[f"{mth}_q_max"].values <= Q_HARD / (1 + qmargin))
                res = an.split_metrics(true_c[idx], score[idx], pcm[idx], ps[idx], true_unsafe[idx], true_dec[idx], KFRAC)
            else:
                base = dict(true_c=true_c[idx], pred_complete=np.ones(len(idx), bool), pred_safe=np.ones(len(idx), bool), true_unsafe=true_unsafe[idx], true_decayed=true_dec[idx], k_frac=KFRAC)
                if mth in flags:
                    q = flags[mth]; mk = lambda r, q=q: (~q).astype(float) + 1e-3 * r.random(len(T))
                else:
                    mk = lambda r: r.random(len(T))
                res = an.stochastic_metrics(mk, idx, base, draws, rng); res["n"] = len(idx)
            res.update(split=sname, method=mth); rows.append(res)
    return pd.DataFrame(rows)


def false_safe_table(T):
    """Controlled-outcome (decay only) and literal (decay or passive Q > 0.54) false-safe counts for the predicted-safe set at several Q margins."""
    true_dec = (T.status == "decayed").values; true_unsafe = true_dec | T.true_exceeds_hard.values
    rows = []
    for sname, mask in splits_of(T).items():
        for key in ("v2", "v2_selfcal", "v1"):
            pcm = ~T[f"{key}_decayed"].values; q = T[f"{key}_q_max"].values.astype(float)
            for m in (0.0, 0.05, 0.10, 0.15):
                ps = pcm & (q <= Q_HARD / (1 + m)) & mask
                rows.append(dict(split=sname, method=key, q_margin=m, n_pred_safe=int(ps.sum()), false_safe_decay_only=int((ps & true_dec).sum()), false_safe_literal=int((ps & true_unsafe).sum()),
                                 rate_decay_only=float((ps & true_dec).sum() / max(1, ps.sum())), rate_literal=float((ps & true_unsafe).sum() / max(1, ps.sum())),
                                 truly_safe_kept=int((ps & ~true_unsafe).sum()), truly_safe=int((~true_unsafe & mask).sum())))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ error tables

def lookahead_errors(T, root, pred_folder=None):
    pred_folder = pred_folder or VAR["v2_folder"]
    m = root / "mars"; recs = []
    for c in T.itertuples():
        tr = pd.read_csv(m / "truth" / c.tag / "passes_3b.csv"); f = m / pred_folder / f"{c.tag}.csv"
        if not f.exists():
            continue
        pr = pd.read_csv(f); n = min(len(tr), len(pr))
        if n < 1:
            continue
        e_abs = pr.alt_m.values[:n] - tr.alt_m.values[:n]
        e_chg = (pr.alt_m.values[:n] - pr.alt_m.values[0]) - (tr.alt_m.values[:n] - tr.alt_m.values[0])
        e_q = pr.q_w_cm2.values[:n] / tr.q_w_cm2.values[:n] - 1
        for k in range(n):
            recs.append(dict(tag=c.tag, status=c.status, k=k + 1, t_days=tr.time_s.values[k] / DAY, dt_s=pr.t_s.values[k] - tr.time_s.values[k],
                             err_alt_m=e_abs[k], err_change_m=e_chg[k], err_q_rel=e_q[k]))
    D = pd.DataFrame(recs); D["bin"] = pd.cut(D.t_days, BINS, right=True, include_lowest=True)
    return D


def lookahead_table(D, label, band_km):
    g = D.groupby("bin", observed=True)
    t = g.agg(n=("err_alt_m", "size"), rms_alt_km=("err_alt_m", lambda v: np.sqrt(np.mean(v ** 2)) / 1e3), max_alt_km=("err_alt_m", lambda v: np.max(np.abs(v)) / 1e3),
              rms_change_km=("err_change_m", lambda v: np.sqrt(np.mean(v ** 2)) / 1e3), rms_q_rel=("err_q_rel", lambda v: np.sqrt(np.mean(v ** 2))),
              pass_time_err_h_rms=("dt_s", lambda v: np.sqrt(np.mean(v ** 2)) / 3600.0)).reset_index()
    t["rms_alt_over_band"] = t.rms_alt_km / band_km
    t.insert(0, "set", label); return t


def dv_error(T, key):
    g = T[(T.status == "completed") & ~T[f"{key}_decayed"]]
    if len(g) < 3:
        return dict(method=key, n=len(g))
    p = g[f"{key}_lazy_dv"].values.astype(float); t = g.true_cost_mps.values.astype(float); err = p - t
    pos = t > 0; ratio = p[pos] / t[pos]   # ratios only where the true label is positive (19 of 77 Mars cells have a zero label); counted in n_true_zero
    return dict(method=key, n=len(g), n_true_zero=int((~pos).sum()), spearman=stats.spearmanr(p, t)[0], median_ratio=float(np.median(ratio)), q1_ratio=float(np.quantile(ratio, 0.25)), q3_ratio=float(np.quantile(ratio, 0.75)),
                rms_err_mps=float(np.sqrt(np.mean(err ** 2))), median_abs_err_mps=float(np.median(np.abs(err))), max_abs_err_mps=float(np.max(np.abs(err))),
                median_true_mps=float(np.median(t)), rms_err_over_median_true=float(np.sqrt(np.mean(err ** 2)) / np.median(t)))


def confusion(T, key):
    cls = ["safe", "breach", "decayed"]
    return pd.crosstab(T.true_class, T[f"{key}_class"]).reindex(index=cls, columns=cls, fill_value=0)


def pass1_offsets(T, root):
    m = root / "mars"; rows = []
    for c in T.itertuples():
        tr = pd.read_csv(m / "truth" / c.tag / "passes_3b.csv")
        r = dict(tag=c.tag, status=c.status, inc_deg=c.inc_deg, truth_alt_m=tr.alt_m.iloc[0])
        for key, folder in (("v2", VAR["v2_folder"]), ("v2_selfcal", VAR["v2sc_folder"]), ("v1", "pred_v1")):
            f = m / folder / f"{c.tag}.csv"
            pr = pd.read_csv(f) if f.exists() else pd.DataFrame()
            r[f"{key}_pass1_offset_m"] = float(pr.alt_m.iloc[0] - tr.alt_m.iloc[0]) if len(pr) else np.nan
            r[f"{key}_alt1_m"] = float(pr.alt_m.iloc[0]) if len(pr) else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def decay_time(T):
    g = T[(T.status == "decayed") & T.v2_decayed]
    err = g.t_decay_pred_days - g.true_t_end_days
    return pd.DataFrame([dict(n=len(g), mean_err_days=float(err.mean()) if len(g) else np.nan, rms_err_days=float(np.sqrt((err ** 2).mean())) if len(g) else np.nan,
                              min_err_days=float(err.min()) if len(g) else np.nan, max_err_days=float(err.max()) if len(g) else np.nan,
                              median_true_end_days=float(g.true_t_end_days.median()) if len(g) else np.nan)])


# ------------------------------------------------------------------------------------------------ thresholds

def verdict_table(T, M, FS, LA, DV, CM, P1, band_km, earth_bw):
    rows = []; add = lambda crit, scope, value, thr, ok, note="": rows.append(dict(criterion=crit, scope=scope, value=value, threshold=thr, passed=ok, note=note))
    for sname in M.split.unique():
        g = M[M.split == sname].set_index("method"); v = g.loc["v2"]
        b1 = g.loc["B1_hand"].spearman; b2c = g.loc["B2c_critical"].spearman; b2e = g.loc["B2e_earth_set"].spearman
        thr1 = max(b1, np.nan_to_num(b2c, nan=-1))
        add("T1 Spearman > max(B1, B2c)", sname, v.spearman, f"> {thr1:.3f} (B1 {b1:.3f}, B2c {b2c:.3f}, B2e {b2e:.3f})", bool(v.spearman > thr1))
        add("T2 top-10% recall > k/n", sname, v.recall_top, f"> {v.random_recall:.3f}", bool(v.recall_top > v.random_recall))
        add("T3 regret <= 0.76 m/s (predicted-complete pick)", sname, v.regret_complete_mps, "<= 0.76", bool(v.regret_complete_mps <= 0.76) if np.isfinite(v.regret_complete_mps) else False,
            f"band-scaled threshold {0.76 * band_km / earth_bw:.3f} m/s: {'PASS' if (np.isfinite(v.regret_complete_mps) and v.regret_complete_mps <= 0.76 * band_km / earth_bw) else 'FAIL'}; 0.72 low bound: {'PASS' if (np.isfinite(v.regret_complete_mps) and v.regret_complete_mps <= 0.72) else 'FAIL'}")
        f = FS[(FS.split == sname) & (FS.method == "v2") & (FS.q_margin == 0.10)].iloc[0]
        add("T4 false-safe, controlled outcome (decay only), Q margin 10%", sname, f.rate_decay_only, "= 0", bool(f.false_safe_decay_only == 0), f"{int(f.false_safe_decay_only)}/{int(f.n_pred_safe)} predicted safe")
        f0 = FS[(FS.split == sname) & (FS.method == "v2") & (FS.q_margin == 0.0)].iloc[0]
        add("T4b false-safe, literal (decay or passive Q>0.54), margin 0", sname, f0.rate_literal, "= 0", bool(f0.false_safe_literal == 0), f"{int(f0.false_safe_literal)}/{int(f0.n_pred_safe)}; at 10%: {int(f.false_safe_literal)}/{int(f.n_pred_safe)}")
    p1 = P1[P1.status == "completed"].v2_pass1_offset_m.dropna()
    rms1 = float(np.sqrt((p1 ** 2).mean())) / 1e3
    add("T5a pass-1 offset RMS (completed, truth seed)", "pooled", rms1, "<= 0.25 km", bool(rms1 <= 0.25), f"{len(p1)} cells; v1 RMS {float(np.sqrt((P1[P1.status == 'completed'].v1_pass1_offset_m.dropna() ** 2).mean())) / 1e3:.3f} km")
    for _, r in LA.iterrows():
        add("T5b RMS altitude error per lookahead bin", str(r.bin), r.rms_alt_km, "<= 1.0 km", bool(r.rms_alt_km <= 1.0), f"n={int(r.n)}, max {r.max_alt_km:.2f} km, /band {r.rms_alt_over_band:.3f}")
        add("T5c RMS relative Q error per lookahead bin", str(r.bin), r.rms_q_rel, "<= 0.10", bool(r.rms_q_rel <= 0.10))
    d = DV[DV.method == "v2"].iloc[0]
    if "median_ratio" in d and np.isfinite(d.get("median_ratio", np.nan)):
        add("T6a dV median predicted/true", "completed", d.median_ratio, "in [0.90, 1.10]", bool(0.90 <= d.median_ratio <= 1.10), f"q1-q3 {d.q1_ratio:.3f}-{d.q3_ratio:.3f}")
        add("T6b dV RMS error / median true label", "completed", d.rms_err_over_median_true, "<= 0.10", bool(d.rms_err_over_median_true <= 0.10), f"RMS {d.rms_err_mps:.3f} m/s, median true {d.median_true_mps:.3f}")
        add("T6c dV Spearman", "completed", d.spearman, ">= 0.95", bool(d.spearman >= 0.95))
    cm = CM.values; tot = cm.sum(); acc = np.trace(cm) / tot
    add("T7a 3-class accuracy", "pooled", acc, ">= 0.95", bool(acc >= 0.95))
    nd = cm[2].sum(); add("T7b decay recall", "pooled", cm[2, 2] / nd if nd else np.nan, ">= 0.95", bool(nd > 0 and cm[2, 2] / nd >= 0.95), f"{cm[2, 2]}/{nd}" if nd else "no decayed cells")
    ncomp = cm[:2].sum(); fd = cm[:2, 2].sum() / ncomp if ncomp else np.nan
    add("T7c false-decay rate", "pooled", fd, "<= 0.05", bool(fd <= 0.05), f"{cm[:2, 2].sum()}/{ncomp}")
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ figures

def fig_scatter(T, out):
    sp = splits_of(T); names = list(sp)
    fig, axes = plt.subplots(1, len(names), figsize=(4.2 * len(names), 4.2), sharex=True, sharey=True)
    for ax, nm in zip(np.atleast_1d(axes), names):
        g = T[sp[nm]]; ok = np.isfinite(g.true_cost_mps) & ~g.v2_decayed; miss = np.isfinite(g.true_cost_mps) & g.v2_decayed
        ax.scatter(g.true_cost_mps[ok], g.v2_lazy_dv[ok], s=14, c="#0b6e99", alpha=0.75, label="estimator v2 (predicted complete)")
        ax.scatter(g.true_cost_mps[miss], g.v2_lazy_dv[miss], s=22, facecolors="none", edgecolors="#b3261e", label="predicted decayed (truly completed)")
        okh = np.isfinite(g.true_cost_mps) & ~g.B1_hand_decayed
        ax.scatter(g.true_cost_mps[okh], g.B1_hand_lazy_dv[okh], s=10, c="#c4612b", alpha=0.4, marker="x", label="B1 double-averaged")
        lim = [0.05, 50]; ax.plot(lim, lim, color="k", lw=0.8); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_title(f"{nm} (n={len(g)})", fontsize=9); ax.grid(alpha=0.2); ax.set_xlabel("true dV label (m/s)")
    np.atleast_1d(axes)[0].set_ylabel("predicted dV (lazy s=0, m/s)"); np.atleast_1d(axes)[0].legend(fontsize=6, loc="upper left")
    fig.tight_layout(); fig.savefig(out / "fig_pred_vs_true_dv.png", dpi=130); plt.close(fig)


def fig_lookahead(LA_c, LA_a, out, band_km):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for LA, lab, col in ((LA_c, "completed cells", "#0b6e99"), (LA_a, "all cells to termination", "#c4612b")):
        x = [(b.left + b.right) / 2 for b in LA.bin]
        ax[0].plot(x, LA.rms_alt_km, "o-", color=col, label=lab + ", RMS"); ax[0].plot(x, LA.max_alt_km, "s--", color=col, alpha=0.5, label=lab + ", max")
    ax[0].axhline(1.0, color="k", lw=0.8, ls=":"); ax[0].axhline(band_km, color="gray", lw=0.8, ls="--"); ax[0].set_yscale("log"); ax[0].set_xlabel("lookahead (days)"); ax[0].set_ylabel("periapsis altitude error (km)")
    ax[0].legend(fontsize=7); ax[0].grid(alpha=0.2); ax[0].set_title("Mars: error vs lookahead (dotted: 1 km margin; dashed: band width)", fontsize=9)
    x = [(b.left + b.right) / 2 for b in LA_c.bin]; ax[1].plot(x, LA_c.rms_q_rel, "o-", color="#0b6e99"); ax[1].axhline(0.10, color="k", ls=":", lw=0.8)
    ax[1].set_xlabel("lookahead (days)"); ax[1].set_ylabel("RMS relative Q error"); ax[1].grid(alpha=0.2)
    fig.tight_layout(); fig.savefig(out / "fig_error_vs_lookahead.png", dpi=130); plt.close(fig)


def fig_confusion(T, out):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8)); cls = ["safe", "breach", "decayed"]
    for ax, key, ttl in zip(axes, ("v2", "v1", "B1_hand"), ("estimator v2", "estimator v1", "B1 double-averaged")):
        cm = confusion(T, key).values; ax.imshow(cm, cmap="Blues")
        for i in range(3):
            for j in range(3):
                ax.text(j, i, cm[i, j], ha="center", va="center", color="white" if cm[i, j] > cm.max() / 2 else "black")
        ax.set_xticks(range(3)); ax.set_xticklabels(cls); ax.set_yticks(range(3)); ax.set_yticklabels(cls); ax.set_xlabel("predicted (passive 60 d)"); ax.set_ylabel("truth"); ax.set_title(ttl, fontsize=10)
    fig.tight_layout(); fig.savefig(out / "fig_confusion.png", dpi=130); plt.close(fig)


# ------------------------------------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); ap.add_argument("--phase", default="all")
    ap.add_argument("--variant", default="", help="'j3' scores the UNFROZEN J3 diagnostic predictions (pred_j3) instead of the frozen estimator (outputs <root>/mars/j3diag); 'v3' scores estimator v3 (pred_v3, outputs <root>/mars/v3)")
    a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); m = root / "mars"
    out = m
    if a.variant == "j3":
        VAR.update(v2_folder="pred_j3", v2sc_folder="pred_j3_selfcal", pc_file="pred_cells_j3.csv", suffix="_j3diag"); out = m / "j3diag"; out.mkdir(exist_ok=True)
    elif a.variant.startswith(("v3", "v4", "v5", "num")):   # "num..." = the numerical predictor N(field) (run_numerical_predictor.jl) and its hybrids; estimator v3 (conditional zonal terms; "v3a" = the per-degree criterion): the "v2" columns of the tables then hold that predictor
        VAR.update(v2_folder=f"pred_{a.variant}", v2sc_folder=f"pred_{a.variant}_selfcal", pc_file=f"pred_cells_{a.variant}.csv", suffix=f"_{a.variant}"); out = m / a.variant; out.mkdir(exist_ok=True)
    mu = mars_mu(root); lib.MU = mu
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    L = build_labels(root); (out == m) and L.to_csv(m / "mars_cell_labels.csv", index=False)
    print(f"labels: {len(L)} cells; status counts:\n{L.status.value_counts().to_string()}")
    print("true class counts:", L.true_class.value_counts().to_dict() if "true_class" in L else {})
    # margin-scaled secondary labels (diagnostic only)
    L2 = build_labels(root, margin=MC["analysis"]["margin_km_secondary"] * 1000.0); (out == m) and L2.to_csv(m / "mars_cell_labels_margin_secondary.csv", index=False)
    if a.phase == "labels":
        return
    T = build_cell_table(root, L); T.to_csv(out / "selection_cell_table.csv", index=False)
    band_km = float(L.band_width_km.dropna().median()); earth_bw = 12.09
    M = evaluate(T); M.to_csv(out / "selection_metrics.csv", index=False)
    FS = false_safe_table(T); FS.to_csv(out / "false_safe_table.csv", index=False)
    P1 = pass1_offsets(T, root); P1.to_csv(out / "pass1_offsets.csv", index=False)
    D_all = lookahead_errors(T, root); D_comp = D_all[D_all.status == "completed"]
    LA_c = lookahead_table(D_comp, "completed", band_km); LA_a = lookahead_table(D_all, "all", band_km)
    pd.concat([LA_c, LA_a]).to_csv(out / "lookahead_error.csv", index=False)
    DV = pd.DataFrame([dv_error(T, k) for k in ("v2", "v2_selfcal", "v1", "B1_hand")]); DV.to_csv(out / "dv_error_table.csv", index=False)
    CMs = {k: confusion(T, k) for k in ("v2", "v2_selfcal", "v1", "B1_hand")}
    for k, c in CMs.items():
        c.to_csv(out / f"confusion_{k}.csv")
    DT = decay_time(T); DT.to_csv(out / "decay_time_error.csv", index=False)
    VT = verdict_table(T, M, FS, LA_c, DV, CMs["v2"], P1, band_km, earth_bw); VT.to_csv(out / "verdict_table.csv", index=False)
    # secondary: margin 0.43 km labels (diagnostic) -> metrics with the same predictions are NOT recomputed (predictions use the 1 km margin); reported label differences only
    fig_scatter(T, out); fig_lookahead(LA_c, LA_a, out, band_km); fig_confusion(T, out)
    cols = ["split", "method", "n", "n_complete_true", "spearman", "kendall", "k", "recall_top", "random_recall", "regret_complete_mps", "regret_safe_mps", "false_safe_rate", "false_decay_rate", "n_pred_safe"]
    print(M[cols].round(3).to_string())
    print(FS[FS.method == "v2"].round(3).to_string())
    print(LA_c.round(3).to_string()); print(LA_a.round(3).to_string())
    print(DV.round(3).to_string()); print(DT.round(3).to_string())
    for k in ("v2", "v1", "B1_hand"):
        print(k); print(CMs[k])
    r1 = P1[P1.status == "completed"]
    print("pass-1 offset RMS (m), completed:", {k: float(np.sqrt((r1[f"{k}_pass1_offset_m"].dropna() ** 2).mean())) for k in ("v1", "v2", "v2_selfcal")})
    print(VT.to_string())


if __name__ == "__main__":
    main()
