"""Task A2, Python half: policy cost on the predicted histories, baselines, selection-skill metrics, thresholds, figures.

Usage (numpy, pandas, scipy, matplotlib):
  python3 -I analyze_selection_skill.py [--results-root DIR] [--out DIR] [--phase all|cells|scan]

Reads <out>/pred, pred_selfcal, hand, pred_cells.csv, attrib_series.csv, scan_* (written by run_selection_skill.jl),
<root>/oracle/oracle_cell_labels.csv (A1 labels) and the truth passes under <root>/border and <root>/border2.
Writes the tables and PNG figures under <out> (default <root>/selection).
"""
import argparse, importlib.util, math, os, pathlib, sys, tomllib
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
spec = importlib.util.spec_from_file_location("oracle_lib", HERE / "oracle_lib.py")
lib = importlib.util.module_from_spec(spec); sys.modules["oracle_lib"] = lib; spec.loader.exec_module(lib)

DAY = 86400.0
CFG = tomllib.load(open(HERE / "selection_skill_config.toml", "rb"))
ORC = tomllib.load(open(HERE / "oracle_config.toml", "rb"))
Q_HARD = ORC["band"]["q_hard_w_cm2"]
MARGIN = CFG["estimator"]["margin_km"] * 1000.0
IMPACT_ALT = CFG["estimator"]["impact_altitude_m"]
PALETTE = {"v1": "#0b6e99", "v1_selfcal": "#5aa9c9", "B1_hand": "#c4612b", "B2_quiet": "#7f7f7f", "B3_random": "#b8b8b8"}
METHODS = ["v1", "v1_selfcal", "B1_hand", "B2_quiet", "B3_random"]


ROOT_LABELS = None


def old_dir(tag):
    return tag[6:]


def truth_path(tag, inc, root):
    return root / "border" / "grid" / old_dir(tag) / "passes_3b.csv" if inc == 50.0 else root / "border2" / "grid" / tag / "passes_3b.csv"


# ------------------------------------------------------------------------------------------------ policy on a history

def policy_costs(x, hq, hl, a, e):
    """Costs and burns of the reference label policy (lazy-to-edge s = 0), near_passive and the exact DP on a band-relative history."""
    if len(x) == 0:
        return dict(lazy_dv=np.nan, lazy_burns=0, near_dv=np.nan, near_burns=0, dp_dv=np.nan, dp_burns=0)
    lo = hq + MARGIN; hi = hl - MARGIN; L = lib.lever_kepler(a, e)
    c, n, _ = lib.lazy(x, lo, hi, L, s=0.0)
    cn, nn, _ = lib.near_passive(x, lo, hi, L)
    d = lib.dp_schedule(lo - x, hi - x, L, None, tie=ORC["dp"]["tie_lambda_mps"])
    return dict(lazy_dv=c, lazy_burns=n, near_dv=cn, near_burns=nn, dp_dv=np.nan if d is None else d["dv"], dp_burns=-1 if d is None else d["burns"])


def predicted_row(df, kind):
    """Status and costs from a predicted history CSV (columns of the Julia writer)."""
    x = df.alt_m.values.astype(float)
    if kind == "hand":
        low = np.where(x <= IMPACT_ALT)[0]
        decayed = len(low) > 0
        if decayed:
            df = df.iloc[:low[0]]; x = x[:low[0]]
        qmax = float(df.q_w_cm2.max()) if len(df) else np.nan
        t_end = float(df.t_s.iloc[-1] / DAY) if len(df) else 0.0
    else:
        decayed = None
        qmax = float(df.q_w_cm2.max()) if len(df) else np.nan
        t_end = float(df.t_s.iloc[-1] / DAY) if len(df) else 0.0
    row = policy_costs(x, df.h_qhigh_m.values, df.h_qlow_m.values, df.a_m.values, df.e.values)
    row.update(q_max=qmax, n_passes=len(df), t_end_days=t_end, alt_min_km=float(x.min() / 1e3) if len(x) else np.nan)
    return row, decayed


# ------------------------------------------------------------------------------------------------ metrics

def rank_metrics(true_c, score):
    ok = np.isfinite(true_c)
    if ok.sum() < 3:
        return np.nan, np.nan
    s = pd.Series(score[ok]).rank(method="average").values
    t = pd.Series(true_c[ok]).rank(method="average").values
    if np.ptp(s) == 0:
        return 0.0, 0.0
    return stats.spearmanr(s, t)[0], stats.kendalltau(score[ok], true_c[ok])[0]


def split_metrics(true_c, score, pred_complete, pred_safe, true_unsafe, true_decayed, k_frac):
    """All metrics of one score vector on one fold (arrays restricted to the fold)."""
    n = len(true_c); ok = np.isfinite(true_c)
    sp, kt = rank_metrics(true_c, score)
    k = max(1, int(math.ceil(k_frac * n)))
    top_true = set(np.argsort(np.where(ok, true_c, np.inf), kind="stable")[:min(k, ok.sum())])
    top_pred = set(np.argsort(score, kind="stable")[:k])
    recall = len(top_true & top_pred) / max(1, len(top_true)) if top_true else np.nan
    best_true = np.nanmin(true_c) if ok.any() else np.nan

    def regret(mask):
        if not mask.any():
            return np.nan, np.nan   # no cell predicted (complete|safe)
        idx = np.where(mask)[0]; pick = idx[np.argmin(score[idx])]
        return (true_c[pick] - best_true if ok[pick] else np.inf), pick
    reg_c, pick_c = regret(pred_complete); reg_s, pick_s = regret(pred_safe)
    fs = (true_unsafe[pred_safe].mean() if pred_safe.any() else np.nan)
    fd = (true_decayed[pred_complete].mean() if pred_complete.any() else np.nan)
    return dict(n=n, n_complete_true=int(ok.sum()), spearman=sp, kendall=kt, k=k, recall_top=recall, random_recall=k / n,
                regret_complete_mps=reg_c, regret_safe_mps=reg_s, picked_decayed_complete=bool(np.isinf(reg_c)) if not np.isnan(reg_c) else False,
                false_safe_rate=fs, false_decay_rate=fd, n_pred_safe=int(pred_safe.sum()), n_pred_complete=int(pred_complete.sum()))


def stochastic_metrics(make_score, fold_idx, base, draws, rng):
    """Average of split_metrics over random score draws (B2 tie-break, B3)."""
    acc = []
    for _ in range(draws):
        sc = make_score(rng)[fold_idx]
        acc.append(split_metrics(base["true_c"], sc, base["pred_complete"], base["pred_safe"], base["true_unsafe"], base["true_decayed"], base["k_frac"]))
    out = {}
    for key in acc[0]:
        vals = np.array([a[key] for a in acc], dtype=float)
        fin = vals[np.isfinite(vals)]
        out[key] = float(fin.mean()) if len(fin) else (np.inf if np.isinf(vals).any() else np.nan)
        if key.startswith("regret"):
            out[key + "_p_pick_decayed"] = float(np.isinf(vals).mean())
    return out


# ------------------------------------------------------------------------------------------------ main flow

def build_cell_table(root, out):
    labels = pd.read_csv(root / "oracle" / "oracle_cell_labels.csv")
    pc = pd.read_csv(out / "pred_cells.csv").set_index("tag")
    rows = []
    for c in labels.itertuples():
        r = dict(tag=c.tag, inc_deg=c.inc_deg, raan_deg=c.raan_deg, omega_deg=c.omega_deg, epoch=int(c.epoch), status=c.status,
                 true_cost_mps=c.true_cost_label_mps, true_lazy_s0_mps=c.lazy_s0_dv, true_near_passive_mps=c.near_passive_dv, true_dp_mps=c.dp_dv,
                 true_exceeds_hard=bool(c.passive_exceeds_hard), true_t_end_days=c.t_end_days)
        for key, folder, kind in (("v1", "pred", "est"), ("v1_selfcal", "pred_selfcal", "est"), ("B1_hand", "hand", "hand")):
            df = pd.read_csv(out / folder / f"{c.tag}.csv")
            row, dec = predicted_row(df, kind)
            if kind == "est":
                dec = bool(pc.loc[c.tag, "decayed_pred" if key == "v1" else "decayed_pred_selfcal"])
            for k_, v_ in row.items():
                r[f"{key}_{k_}"] = v_
            r[f"{key}_decayed"] = dec
            r[f"{key}_breach"] = bool((not dec) and row["q_max"] > Q_HARD)
            r[f"{key}_class"] = "decayed" if dec else ("breach" if r[f"{key}_breach"] else "safe")
        r["seed_truth_m"] = pc.loc[c.tag, "seed_truth_m"]; r["seed_selfcal_m"] = pc.loc[c.tag, "seed_selfcal_m"]
        r["t_decay_pred_days"] = pc.loc[c.tag, "t_decay_pred_days"]
        r["true_class"] = "decayed" if c.status == "decayed" else ("breach" if c.passive_exceeds_hard else "safe")
        rows.append(r)
    return pd.DataFrame(rows)


def quiet_flag(T):
    q = CFG["baselines"]["quiet_set"]
    return np.array([any(abs(i - qi) < 0.05 and abs(r - qr) < 0.5 for qi, qr in q) for i, r in zip(T.inc_deg, T.raan_deg)])


def splits_of(T):
    out = {"pooled": np.ones(len(T), bool)}
    for e in sorted(T.epoch.unique()):
        out[f"epoch_{e}_held_out"] = (T.epoch == e).values
    h = CFG["splits"]["raan_halves"]
    out["raan_0-135"] = T.raan_deg.isin(h[0]).values
    out["raan_180-315"] = T.raan_deg.isin(h[1]).values
    return out


def evaluate(T):
    rng = np.random.default_rng(CFG["baselines"]["random_seed"]); draws = CFG["baselines"]["random_draws"]
    kfrac = CFG["thresholds"]["top_fraction"]
    quiet = quiet_flag(T)
    true_c = T.true_cost_mps.values.astype(float)
    true_unsafe = ((T.status == "decayed") | T.true_exceeds_hard).values
    true_dec = (T.status == "decayed").values
    rows = []
    for sname, mask in splits_of(T).items():
        idx = np.where(mask)[0]
        for m in METHODS:
            if m in ("v1", "v1_selfcal", "B1_hand"):
                score = np.where(T[f"{m}_decayed"].values, np.inf, T[f"{m}_lazy_dv"].values.astype(float))
                pc = ~T[f"{m}_decayed"].values; ps = pc & ~T[f"{m}_breach"].values
                res = split_metrics(true_c[idx], score[idx], pc[idx], ps[idx], true_unsafe[idx], true_dec[idx], kfrac)
            else:
                base = dict(true_c=true_c[idx], pred_complete=np.ones(len(idx), bool), pred_safe=np.ones(len(idx), bool), true_unsafe=true_unsafe[idx], true_decayed=true_dec[idx], k_frac=kfrac)
                if m == "B2_quiet":
                    mk = lambda r: (~quiet).astype(float) + 1e-3 * r.random(len(T))
                else:
                    mk = lambda r: r.random(len(T))
                res = stochastic_metrics(mk, idx, base, draws, rng)
                res["n"] = len(idx)
            res.update(split=sname, method=m)
            rows.append(res)
    return pd.DataFrame(rows)


def threshold_table(M):
    th = CFG["thresholds"]; rows = []
    for sname in M.split.unique():
        g = M[M.split == sname].set_index("method")
        v = g.loc["v1"]
        fs = v.false_safe_rate
        rows.append(dict(split=sname, metric="false_safe_rate (predicted safe = completes 60 d and Q<=0.54)", value=fs, threshold=f"<= {th['max_false_safe_rate']}", passed=bool(fs <= th["max_false_safe_rate"]) if np.isfinite(fs) else False))
        rows.append(dict(split=sname, metric="false_decay_rate (predicted completes, truly decayed)", value=v.false_decay_rate, threshold="reported", passed=np.nan))
        for tag, reg in (("regret_complete_mps", v.regret_complete_mps), ("regret_safe_mps", v.regret_safe_mps)):
            rows.append(dict(split=sname, metric=f"{tag} (<= 0.76; 0.72 low)", value=reg, threshold=f"<= {th['max_selection_regret_mps']}", passed=bool(reg <= th["max_selection_regret_mps"]) if np.isfinite(reg) else False))
        b1, b2 = g.loc["B1_hand"].spearman, g.loc["B2_quiet"].spearman
        rows.append(dict(split=sname, metric="spearman > B1 and B2", value=v.spearman, threshold=f"> max(B1={b1:.3f}, B2={b2:.3f})", passed=bool(v.spearman > max(b1, np.nan_to_num(b2, nan=-1)))))
        rows.append(dict(split=sname, metric="top-10% recall > k/N", value=v.recall_top, threshold=f"> {v.random_recall:.3f}", passed=bool(v.recall_top > v.random_recall)))
    return pd.DataFrame(rows)


def confusion(T, key):
    cls = ["safe", "breach", "decayed"]
    cm = pd.crosstab(T.true_class, T[f"{key}_class"]).reindex(index=cls, columns=cls, fill_value=0)
    return cm


# ------------------------------------------------------------------------------------------------ lookahead error

def lookahead_errors(T, root, out):
    bins = CFG["lookahead"]["bin_edges_days"]; recs = []
    for c in T.itertuples():
        tr = pd.read_csv(truth_path(c.tag, c.inc_deg, root)); pr = pd.read_csv(out / "pred" / f"{c.tag}.csv")
        n = min(len(tr), len(pr))
        if n < 2:
            continue
        e_abs = pr.alt_m.values[:n] - tr.alt_m.values[:n]
        e_chg = (pr.alt_m.values[:n] - pr.alt_m.values[0]) - (tr.alt_m.values[:n] - tr.alt_m.values[0])
        e_q = pr.q_w_cm2.values[:n] / tr.q_w_cm2.values[:n] - 1
        for k in range(n):
            recs.append(dict(tag=c.tag, status=c.status, k=k + 1, t_days=tr.time_s.values[k] / DAY, dt_s=pr.t_s.values[k] - tr.time_s.values[k],
                             err_alt_m=e_abs[k], err_change_m=e_chg[k], err_q_rel=e_q[k]))
    D = pd.DataFrame(recs)
    D["bin"] = pd.cut(D.t_days, bins, right=True, include_lowest=True)
    return D


def lookahead_table(D, label):
    g = D.groupby("bin", observed=True)
    t = g.agg(n=("err_alt_m", "size"), rms_alt_km=("err_alt_m", lambda v: np.sqrt(np.mean(v ** 2)) / 1e3), max_alt_km=("err_alt_m", lambda v: np.max(np.abs(v)) / 1e3),
              rms_change_km=("err_change_m", lambda v: np.sqrt(np.mean(v ** 2)) / 1e3), max_change_km=("err_change_m", lambda v: np.max(np.abs(v)) / 1e3),
              rms_q_rel=("err_q_rel", lambda v: np.sqrt(np.mean(v ** 2))), pass_time_err_h_rms=("dt_s", lambda v: np.sqrt(np.mean(v ** 2)) / 3600.0)).reset_index()
    t.insert(0, "set", label)
    return t


def attrib_lookahead(out):
    A = pd.read_csv(out / "attrib_series.csv"); bins = CFG["lookahead"]["bin_edges_days"]; rows = []
    for v, g in A.groupby("variant"):
        g = g.copy(); g["t_days"] = g.t_truth_s / DAY
        g["err_change"] = (g.rp_pred_m - g.groupby("tag").rp_pred_m.transform("first")) - (g.rp_truth_m - g.groupby("tag").rp_truth_m.transform("first"))
        g["bin"] = pd.cut(g.t_days, bins, include_lowest=True)
        t = g.groupby("bin", observed=True).agg(n=("err_change", "size"), rms_change_km=("err_change", lambda x: np.sqrt(np.mean(x ** 2)) / 1e3), max_change_km=("err_change", lambda x: np.abs(x).max() / 1e3)).reset_index()
        t.insert(0, "variant", v); rows.append(t)
    return pd.concat(rows)


def attrib_tables(out):
    A = pd.read_csv(out / "attrib_series.csv"); rows = []
    for (v, tag), h in A.groupby(["variant", "tag"], sort=False):
        e = (h.rp_pred_m - h.rp_pred_m.iloc[0]) - (h.rp_truth_m - h.rp_truth_m.iloc[0])
        ea = h.rp_pred_m - h.rp_truth_m
        rows.append(dict(variant=v, tag=tag, passes=len(h), n_truth=int(h.n_truth.iloc[0]), n_pred=int(h.n_pred.iloc[0]),
                         rms_change_km=float(np.sqrt((e ** 2).mean()) / 1e3), max_change_km=float(e.abs().max() / 1e3),
                         rms_rp_km=float(np.sqrt((ea ** 2).mean()) / 1e3), max_rp_km=float(ea.abs().max() / 1e3),
                         end_time_err_h=float((h.t_pred_s.iloc[-1] - h.t_truth_s.iloc[-1]) / 3600.0)))
    return pd.DataFrame(rows)


def dv_error_table(T):
    rows = []
    for key in ("v1", "v1_selfcal", "B1_hand"):
        g = T[(T.status == "completed") & ~T[f"{key}_decayed"]]
        for lab, pcol, tcol in (("lazy_s0 vs label (DP)", f"{key}_lazy_dv", "true_cost_mps"), ("near_passive vs true near_passive", f"{key}_near_dv", "true_near_passive_mps")):
            ratio = g[pcol] / g[tcol]; err = g[pcol] - g[tcol]
            rows.append(dict(method=key, comparison=lab, n=len(g), spearman=stats.spearmanr(g[pcol], g[tcol])[0], median_ratio=ratio.median(), q1_ratio=ratio.quantile(0.25), q3_ratio=ratio.quantile(0.75),
                             rms_err_mps=float(np.sqrt((err ** 2).mean())), max_abs_err_mps=float(err.abs().max()), median_abs_err_mps=float(err.abs().median())))
        g2 = g
        rows.append(dict(method=key, comparison="burn count lazy_s0 (predicted vs true)", n=len(g2), spearman=np.nan, median_ratio=np.nan, q1_ratio=np.nan, q3_ratio=np.nan,
                         rms_err_mps=float(np.sqrt(((g2[f"{key}_lazy_burns"] - pd.read_csv(ROOT_LABELS).set_index("tag").loc[g2.tag, "lazy_s0_burns"].values) ** 2).mean())), max_abs_err_mps=np.nan, median_abs_err_mps=np.nan))
    return pd.DataFrame(rows)


def decay_time_table(T):
    g = T[(T.status == "decayed") & T.v1_decayed]
    err = g.t_decay_pred_days - g.true_t_end_days
    return pd.DataFrame([dict(n=len(g), mean_err_days=err.mean(), rms_err_days=float(np.sqrt((err ** 2).mean())), min_err_days=err.min(), max_err_days=err.max(),
                              median_true_end_days=g.true_t_end_days.median(), min_true_end_days=g.true_t_end_days.min())])


def margin_cv(T):
    """Predicted-Q safety factor fit inside the training folds only: flag a cell safe only if q_pred <= Q_hard/(1+m); m = smallest value that leaves no
    truly-unsafe predicted-complete training cell in the safe set. Test on the held-out fold. (v1 itself has no fitted parameter.)"""
    rows = []
    for key in ("v1", "v1_selfcal"):
        pc = ~T[f"{key}_decayed"].values; q = T[f"{key}_q_max"].values; unsafe = ((T.status == "decayed") | T.true_exceeds_hard).values
        for sname, mask in splits_of(T).items():
            train = ~mask if sname != "pooled" else mask
            need = Q_HARD / q[train & pc & unsafe] - 1
            m = float(max(0.0, need.max())) + 1e-9 if len(need) else 0.0
            safe = pc & (q <= Q_HARD / (1 + m)) & mask
            base_safe = pc & (q <= Q_HARD) & mask
            true_safe = ~unsafe & mask
            rows.append(dict(method=key, split=sname, margin_m=m, trained_on="in-sample (pooled)" if sname == "pooled" else "other folds", n_pred_safe=int(safe.sum()), false_safe=int((safe & unsafe).sum()),
                             false_safe_rate=float((safe & unsafe).sum() / max(1, safe.sum())), n_pred_safe_without_margin=int(base_safe.sum()), false_safe_without_margin=int((base_safe & unsafe).sum()),
                             recall_of_truly_safe=float((safe & true_safe).sum() / max(1, true_safe.sum()))))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ figures

def fig_scatter(T, out):
    sp = splits_of(T); names = list(sp)
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.2), sharex=True, sharey=True)
    for ax, nm in zip(axes.ravel(), names):
        g = T[sp[nm]]
        ok = np.isfinite(g.true_cost_mps) & ~g.v1_decayed
        miss = np.isfinite(g.true_cost_mps) & g.v1_decayed
        ax.scatter(g.true_cost_mps[ok], g.v1_lazy_dv[ok], s=14, c=PALETTE["v1"], alpha=0.75, label="estimator v1 (completes)")
        ax.scatter(g.true_cost_mps[miss], g.v1_lazy_dv[miss], s=22, facecolors="none", edgecolors="#b3261e", label="predicted decayed (truly completed)")
        okh = np.isfinite(g.true_cost_mps) & ~g.B1_hand_decayed
        ax.scatter(g.true_cost_mps[okh], g.B1_hand_lazy_dv[okh], s=10, c=PALETTE["B1_hand"], alpha=0.4, marker="x", label="B1 double-averaged")
        lim = [0.3, 20]; ax.plot(lim, lim, color="k", lw=0.8)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_title(f"{nm} (n={len(g)})", fontsize=10); ax.grid(alpha=0.2)
    for ax in axes[1]: ax.set_xlabel("true dV label (m/s, exact DP on truth)")
    for ax in axes[:, 0]: ax.set_ylabel("predicted dV (lazy s=0, m/s)")
    axes[0, 0].legend(fontsize=7, loc="upper left")
    fig.suptitle("Predicted vs true maintenance dV per split (truly completed cells)")
    fig.tight_layout(); fig.savefig(out / "fig_pred_vs_true_dv.png", dpi=130); plt.close(fig)


def fig_lookahead(D_all, D_comp, A, out):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for D, lab, col in ((D_comp, "completed cells (203)", PALETTE["v1"]), (D_all, "all 336 cells (to termination)", "#c4612b")):
        t = lookahead_table(D, lab); x = [(b.left + b.right) / 2 for b in t.bin]
        ax[0].plot(x, t.rms_alt_km, "o-", color=col, label=lab + ", RMS"); ax[0].plot(x, t.max_alt_km, "s--", color=col, alpha=0.5, label=lab + ", max")
    ax[0].set_xlabel("lookahead (days; bin center)"); ax[0].set_ylabel("periapsis altitude error (km)"); ax[0].set_yscale("log"); ax[0].grid(alpha=0.2); ax[0].legend(fontsize=7)
    ax[0].set_title("Estimator error vs lookahead, drag on (absolute altitude)")
    cols = {"v0_drag_on": PALETTE["v1"], "v5_drag_off": "#2a7f3f", "v3_no_bodies_drag_on": "#8a5ca8"}
    sub = A[A.variant.isin(cols)]
    for v, g in sub.groupby("variant"):
        ax[1].plot(range(len(g)), g.rms_change_km.values, "o", color=cols[v], label=v)
    ax[1].set_xticks(range(12)); ax[1].set_xticklabels(list(sub[sub.variant == "v0_drag_on"].tag), rotation=90, fontsize=7)
    ax[1].set_ylabel("RMS error of periapsis radius change (km)")
    ax[1].grid(alpha=0.2); ax[1].legend(fontsize=7); ax[1].set_title("12 attribution cells, RMS of r_p change vs truth", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "fig_error_vs_lookahead.png", dpi=130); plt.close(fig)


def fig_confusion(T, out):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2)); cls = ["safe", "breach", "decayed"]
    for ax, key, ttl in zip(axes, ("v1", "v1_selfcal", "B1_hand"), ("estimator v1 (truth seed)", "estimator v1 (self-calibrated seed)", "B1 double-averaged hand model")):
        cm = confusion(T, key).values
        ax.imshow(cm, cmap="Blues")
        for i in range(3):
            for j in range(3):
                ax.text(j, i, cm[i, j], ha="center", va="center", color="white" if cm[i, j] > cm.max() / 2 else "black")
        ax.set_xticks(range(3)); ax.set_xticklabels(cls); ax.set_yticks(range(3)); ax.set_yticklabels(cls)
        ax.set_xlabel("predicted (passive 60 d)"); ax.set_ylabel("truth"); ax.set_title(ttl, fontsize=10)
    fig.suptitle("Passive-outcome classes: safe = completes 60 d with peak Q <= 0.54; breach = completes, Q > 0.54; decayed = impact before 60 d")
    fig.tight_layout(); fig.savefig(out / "fig_confusion.png", dpi=130); plt.close(fig)


def fig_metrics(M, out):
    sp = list(M.split.unique()); fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
    w = 0.16; cap = 25.0
    for j, (met, ttl) in enumerate((("spearman", "Spearman (completed cells)"), ("recall_top", "top-10% recall"), ("regret_complete_mps", "regret of predicted-complete pick (m/s)"))):
        for i, m in enumerate(METHODS):
            raw = [M[(M.split == s) & (M.method == m)][met].values[0] for s in sp]
            v = [cap if np.isinf(x) else (x if np.isfinite(x) else 0) for x in raw]
            bars = ax[j].bar(np.arange(len(sp)) + i * w, v, w, color=PALETTE[m], label=m)
            for b, x in zip(bars, raw):
                if np.isinf(x):
                    b.set_hatch("///"); b.set_edgecolor("k")
        ax[j].set_xticks(np.arange(len(sp)) + 2 * w); ax[j].set_xticklabels(sp, rotation=30, ha="right", fontsize=8); ax[j].set_title(ttl, fontsize=9); ax[j].grid(alpha=0.2, axis="y")
    ax[2].axhline(CFG["thresholds"]["max_selection_regret_mps"], color="k", ls="--", lw=0.8); ax[2].set_yscale("symlog", linthresh=1); ax[2].set_ylim(0, cap * 1.2)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, ncol=5, loc="upper center", fontsize=8); fig.text(0.99, 0.01, "hatched regret bar = the best predicted-complete cell truly decayed (regret undefined, drawn at the axis cap)", ha="right", fontsize=7)
    fig.tight_layout(rect=(0, 0.03, 1, 0.94)); fig.savefig(out / "fig_metrics_by_split.png", dpi=130); plt.close(fig)


def fig_seed_time(T, out):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    d = (T.seed_selfcal_m - T.seed_truth_m)
    ax[0].hist(d, bins=40, color=PALETTE["v1"]); ax[0].set_xlabel("self-calibrated seed - truth seed (m)"); ax[0].set_title("Seed calibration offset"); ax[0].grid(alpha=0.2)
    c = T[T.status == "completed"]
    ax[1].scatter(c.true_cost_mps, c.v1_selfcal_lazy_dv, s=12, c=PALETTE["v1_selfcal"], label="self-calibrated seed")
    ax[1].scatter(c.true_cost_mps, c.v1_lazy_dv, s=10, c=PALETTE["v1"], alpha=0.6, label="truth seed")
    ax[1].plot([0.3, 20], [0.3, 20], "k", lw=0.8); ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.2)
    ax[1].set_xlabel("true dV (m/s)"); ax[1].set_ylabel("predicted dV (m/s)")
    fig.tight_layout(); fig.savefig(out / "fig_seed_calibration.png", dpi=130); plt.close(fig)


# ------------------------------------------------------------------------------------------------ scan

def scan_analysis(root, out, T):
    S = pd.read_csv(out / "scan_summary.csv"); H = pd.read_csv(out / "scan_hist.csv.gz")
    key = ["inc_deg", "raan_deg", "omega_deg", "epoch"]
    rows = []
    for gk, h in H.groupby(key, sort=False):
        pol = policy_costs(h.alt_m.values, h.h_qhigh_m.values, h.h_qlow_m.values, h.a_m.values, h.e.values)
        q = h.q_w_cm2.max()
        margin = float(np.min(np.minimum(h.alt_m.values - (h.h_qhigh_m.values + MARGIN), (h.h_qlow_m.values - MARGIN) - h.alt_m.values))) if len(h) else np.nan
        rows.append(dict(zip(key, gk), lazy_dv=pol["lazy_dv"], lazy_burns=pol["lazy_burns"], near_dv=pol["near_dv"], q_max=q, min_band_margin_m=margin))
    P = pd.DataFrame(rows); S = S.merge(P, on=key, how="left")
    S["pred_breach"] = (~S.decayed_pred) & (S.q_max > Q_HARD)
    S["pred_class"] = np.where(S.decayed_pred, "decayed", np.where(S.pred_breach, "breach", "safe"))
    Tm = T[key[:3] + ["epoch", "status", "true_cost_mps", "true_class"]]
    S = S.merge(Tm, on=key, how="left")
    S.to_csv(out / "scan_candidates.csv", index=False)
    tops = []
    topn = CFG["scan"]["top_n"]
    for e, g in S.groupby("epoch"):
        comp = g[~g.decayed_pred].sort_values(["lazy_dv", "min_band_margin_m"], ascending=[True, False])
        truth_e = T[(T.epoch == e) & (T.status == "completed")]
        tb = truth_e.sort_values("true_cost_mps").iloc[0]
        for r, c in enumerate(comp.head(topn).itertuples(), 1):
            tops.append(dict(epoch=e, rank=r, inc_deg=c.inc_deg, raan_deg=c.raan_deg, omega_deg=c.omega_deg, pred_dv_mps=c.lazy_dv, pred_burns=c.lazy_burns, pred_class=c.pred_class,
                             truth_exists=bool(np.isfinite(c.true_cost_mps) or isinstance(c.status, str)), truth_status=c.status if isinstance(c.status, str) else "no truth cell",
                             true_dv_mps=c.true_cost_mps, true_rank_in_epoch=(int((truth_e.true_cost_mps < c.true_cost_mps).sum() + 1) if np.isfinite(c.true_cost_mps) else np.nan),
                             n_truth_completed_in_epoch=len(truth_e), truth_best_epoch=f"{tb.tag}", truth_best_dv_mps=tb.true_cost_mps))
    TP = pd.DataFrame(tops); TP.to_csv(out / "scan_top5.csv", index=False)
    # the same ranking restricted to candidates that have a truth cell (verifiable subset): predicted top-5 vs true top-5 per epoch
    ver = []
    for e, g in S[S.status.notna()].groupby("epoch"):
        pc_ = g[~g.decayed_pred]; tc = g[g.status == "completed"]
        pt = pc_.sort_values(["lazy_dv", "min_band_margin_m"], ascending=[True, False]).head(topn); tt = tc.sort_values("true_cost_mps").head(topn)
        k_ = [tuple(x) for x in pt[["inc_deg", "raan_deg", "omega_deg"]].values]; kt_ = [tuple(x) for x in tt[["inc_deg", "raan_deg", "omega_deg"]].values]
        ver.append(dict(epoch=e, n_truth_cells_in_scan=len(g), pred_top5=";".join(f"{a:g}/{b:g}/{c:g}" for a, b, c in k_), true_top5=";".join(f"{a:g}/{b:g}/{c:g}" for a, b, c in kt_),
                        overlap=len(set(k_) & set(kt_)), pred_best_true_dv=float(pt.true_cost_mps.iloc[0]) if len(pt) else np.nan, true_best_dv=float(tt.true_cost_mps.iloc[0]),
                        true_dv_of_pred_top5=";".join(f"{v:.2f}" for v in pt.true_cost_mps.values)))
    VER = pd.DataFrame(ver); VER.to_csv(out / "scan_verifiable_top5.csv", index=False)
    # agreement where truth exists
    both = S[np.isfinite(S.true_cost_mps) & ~S.decayed_pred]
    sp = stats.spearmanr(both.lazy_dv, both.true_cost_mps)[0] if len(both) > 3 else np.nan
    fig, axes = plt.subplots(3, 4, figsize=(17, 11), sharex=True, sharey=True)
    cmap = plt.get_cmap("viridis").copy(); cmap.set_bad("#d9d9d9")
    vmax = np.nanpercentile(S.lazy_dv[~S.decayed_pred], 95)
    for i, e in enumerate(sorted(S.epoch.unique())):
        for j, inc in enumerate(sorted(S.inc_deg.unique())):
            g = S[(S.epoch == e) & (S.inc_deg == inc)]
            Z = g.assign(v=np.where(g.decayed_pred, np.nan, g.lazy_dv)).pivot(index="omega_deg", columns="raan_deg", values="v")
            ax = axes[i, j]; im = ax.imshow(Z.values, origin="lower", aspect="auto", cmap=cmap, vmin=0, vmax=vmax,
                                          extent=[Z.columns.min() - 7.5, Z.columns.max() + 7.5, Z.index.min() - 7.5, Z.index.max() + 7.5])
            gt = g[np.isfinite(g.true_cost_mps)]
            ax.scatter(gt.raan_deg, gt.omega_deg, s=26, facecolors="none", edgecolors="white", lw=0.8)
            gd = g[g.decayed_pred]
            ax.scatter(gd.raan_deg, gd.omega_deg, s=5, c="#b3261e", marker="x", lw=0.6)
            best = g[~g.decayed_pred].sort_values("lazy_dv").head(topn)
            ax.scatter(best.raan_deg, best.omega_deg, s=60, marker="*", c="#ffcc00", edgecolors="k", lw=0.5)
            ax.set_title(f"inc {inc:g}, epoch {e}", fontsize=9)
    for ax in axes[-1]: ax.set_xlabel("RAAN (deg)")
    for ax in axes[:, 0]: ax.set_ylabel("omega (deg)")
    fig.colorbar(im, ax=axes, shrink=0.6, label="predicted dV (lazy s=0, m/s); grey = predicted decayed")
    fig.suptitle("Candidate scan, predicted dV: stars = best 5 per epoch, white rings = candidates with a truth cell, red x = predicted decayed", y=0.93)
    fig.savefig(out / "fig_scan_heatmaps.png", dpi=120, bbox_inches="tight"); plt.close(fig)
    return S, TP, sp, VER


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap)
    ap.add_argument("--out", default=None); ap.add_argument("--phase", default="all")
    a = ap.parse_args()
    global ROOT_LABELS
    root = pathlib.Path(results_root(a.results_root)); out = pathlib.Path(a.out) if a.out else root / "selection"
    ROOT_LABELS = root / "oracle" / "oracle_cell_labels.csv"
    if a.phase in ("all", "cells"):
        T = build_cell_table(root, out); T.to_csv(out / "selection_cell_table.csv", index=False)
        M = evaluate(T); M.to_csv(out / "selection_metrics.csv", index=False)
        TH = threshold_table(M); TH.to_csv(out / "threshold_table.csv", index=False)
        for key in ("v1", "v1_selfcal", "B1_hand"):
            confusion(T, key).to_csv(out / f"confusion_{key}.csv")
        D_all = lookahead_errors(T, root, out); D_comp = D_all[D_all.status == "completed"]
        LA = pd.concat([lookahead_table(D_comp, "completed"), lookahead_table(D_all, "all_cells")]); LA.to_csv(out / "lookahead_error.csv", index=False)
        A = attrib_tables(out); A.to_csv(out / "attribution_cells_error.csv", index=False)
        DV = dv_error_table(T); DV.to_csv(out / "dv_error_table.csv", index=False)
        DT = decay_time_table(T); DT.to_csv(out / "decay_time_error.csv", index=False)
        MC = margin_cv(T); MC.to_csv(out / "safety_margin_cv.csv", index=False)
        AL = attrib_lookahead(out); AL.to_csv(out / "lookahead_error_attribution_cells.csv", index=False)
        fig_scatter(T, out); fig_lookahead(D_all, D_comp, A, out); fig_confusion(T, out); fig_metrics(M, out); fig_seed_time(T, out)
        pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
        cols = ["split", "method", "n", "n_complete_true", "spearman", "kendall", "k", "recall_top", "random_recall", "regret_complete_mps", "regret_safe_mps", "false_safe_rate", "false_decay_rate", "n_pred_safe"]
        print(M[cols].round(3).to_string()); print(TH.round(3).to_string()); print(LA.round(3).to_string())
        for key in ("v1", "v1_selfcal", "B1_hand"):
            print(key); print(confusion(T, key))
        print(DV.round(3).to_string()); print(DT.round(3).to_string()); print(MC.round(3).to_string()); print(AL.round(3).to_string())
        print(A.round(2).to_string())
    if a.phase in ("all", "scan"):
        T = pd.read_csv(out / "selection_cell_table.csv")
        S, TP, sp, VER = scan_analysis(root, out, T)
        print(TP.round(3).to_string()); print(VER.to_string()); print("pred dV==0 candidates:", int((S.lazy_dv == 0).sum()), "of", len(S), " predicted decayed:", int(S.decayed_pred.sum()), "pred breach:", int(S.pred_breach.sum()))
        print("scan Spearman (scan candidates with truth, predicted-complete):", sp)


if __name__ == "__main__":
    main()
