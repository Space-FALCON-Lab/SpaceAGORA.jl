"""Estimator v2 vs v1: pass-1 offset, short-period amplitude/phase tests, error vs lookahead, selection metrics, scaling collapse, planets.

Usage (numpy, pandas, scipy, matplotlib):
  python3 -I analyze_estimator_v2.py [--results-root DIR] [--phase all|score|compare|osc|blind] [--blind-dir DIR]

Reads <root>/selection (v1 outputs of run_selection_skill.jl + analyze_selection_skill.py), <root>/estimator_v2 (outputs of run_estimator_v2.jl),
the truth under <root>/border, <root>/border2, <root>/border3 and the A1 labels. Writes tables and PNG figures under <root>/estimator_v2.
`score` runs analyze_selection_skill.py's cell analysis on the v2 predictions (same folds, thresholds and baselines as v1).
"""
import argparse, importlib.util, math, pathlib, subprocess, sys, tomllib
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); sys.modules[name] = mod; spec.loader.exec_module(mod)
    return mod


SK = load_module("analyze_selection_skill", HERE / "analyze_selection_skill.py")
CFG = SK.CFG
DAY = 86400.0
C1, C2, CGREY = "#c4612b", "#0b6e99", "#7f7f7f"      # v1, v2, neutral
Q_HARD = SK.Q_HARD


def rms(x):
    x = np.asarray(x, float)
    return float(np.sqrt(np.mean(x ** 2))) if len(x) else np.nan


# ---------------------------------------------------------------------------------------------- selection / pass-1

def score(root, out, variants=("",)):
    """Run the A2 cell analysis (metrics, folds, thresholds, lookahead) on a directory of predictions."""
    for v in variants:
        d = root / ("estimator_v2" + (("_" + v) if v else ""))
        if not (d / "pred_cells.csv").exists():
            continue
        subprocess.run([sys.executable, "-I", str(HERE / "analyze_selection_skill.py"), "--results-root", str(root), "--out", str(d), "--phase", "cells"],
                       check=True, stdout=open(d / "analyze_cells.log", "w"), stderr=subprocess.STDOUT)
        print("scored", d)


def pass1_offsets(root, dirs):
    """Pass-1 periapsis altitude error (predicted - truth, m) for every cell and each prediction directory."""
    labels = pd.read_csv(root / "oracle" / "oracle_cell_labels.csv"); rows = []
    for c in labels.itertuples():
        tr = pd.read_csv(SK.truth_path(c.tag, c.inc_deg, root))
        r = dict(tag=c.tag, inc_deg=c.inc_deg, epoch=int(c.epoch), status=c.status, truth_alt0_m=tr.alt_m.iloc[0])
        for name, d in dirs.items():
            p = pd.read_csv(d / "pred" / f"{c.tag}.csv"); r[f"{name}_off_m"] = p.alt_m.iloc[0] - tr.alt_m.iloc[0]
            r[f"{name}_t1_err_s"] = p.t_s.iloc[0] - tr.time_s.iloc[0]
            ps = pd.read_csv(d / "pred_selfcal" / f"{c.tag}.csv"); r[f"{name}_selfcal_off_m"] = ps.alt_m.iloc[0] - tr.alt_m.iloc[0]
        rows.append(r)
    return pd.DataFrame(rows)


def fig_pass1(P, out, names):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    bins = np.linspace(-1.5, 1.5, 61)
    for n, col in names:
        ax[0].hist(P[f"{n}_off_m"] / 1e3, bins=bins, alpha=0.55, color=col, label=f"{n}: RMS {rms(P[f'{n}_off_m'])/1e3:.3f} km, mean {P[f'{n}_off_m'].mean()/1e3:+.3f} km")
    ax[0].set_xlabel("pass-1 periapsis altitude error, predicted - truth (km)"); ax[0].set_ylabel("cells (336)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.2)
    ax[0].set_title("Pass-1 offset, truth seed", fontsize=10)
    corr = (P["v2_off_m"] * 0 + (P["v1_off_m"] - P["v2_off_m"])) / 1e3   # predicted correction of pass-1 altitude (v1 minus v2 = what v2 removes)
    ax[1].scatter(corr, P["v1_off_m"] / 1e3, s=9, color=C1, alpha=0.6)
    lim = 1.5; ax[1].plot([-lim, lim], [-lim, lim], "k", lw=0.8); ax[1].set_xlim(-lim, lim); ax[1].set_ylim(-lim, lim)
    sl = np.polyfit(corr, P["v1_off_m"] / 1e3, 1)
    ax[1].set_xlabel("what v2 changes at pass 1: v1 - v2 altitude (km)"); ax[1].set_ylabel("v1 offset vs truth (km)")
    ax[1].set_title(f"v1 error explained by the short-period correction (slope {sl[0]:.2f}, R2 {np.corrcoef(corr, P['v1_off_m'])[0, 1] ** 2:.2f})", fontsize=9); ax[1].grid(alpha=0.2)
    ax[2].scatter(P["v1_off_m"] / 1e3, P["v2_off_m"] / 1e3, s=9, color=C2, alpha=0.6); ax[2].plot([-lim, lim], [-lim, lim], "k", lw=0.8)
    ax[2].set_xlim(-lim, lim); ax[2].set_ylim(-lim, lim); ax[2].set_xlabel("v1 offset (km)"); ax[2].set_ylabel("v2 offset (km)"); ax[2].grid(alpha=0.2)
    ax[2].set_title("Per cell: v1 vs v2", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "fig_pass1_offset.png", dpi=130); plt.close(fig)


# ---------------------------------------------------------------------------------------------- false-safe definitions

def false_safe_table(T, label):
    """Three false-safe definitions on the same folds: literal (decayed or passive Q > 0.54), controlled-outcome (decayed only),
    and literal with a Q margin m: predicted safe iff completes and Q_pred <= 0.54/(1+m)."""
    rows = []
    unsafe_lit = ((T.status == "decayed") | T.true_exceeds_hard).values
    unsafe_dec = (T.status == "decayed").values
    for key in ("v1", "v1_selfcal"):
        pc = ~T[f"{key}_decayed"].values; q = T[f"{key}_q_max"].values
        for sname, mask in SK.splits_of(T).items():
            for m in (0.0, 0.05, 0.10, 0.15):
                ps = pc & (q <= Q_HARD / (1 + m)) & mask
                true_safe = ~unsafe_lit & mask
                for dname, unsafe in (("literal", unsafe_lit), ("controlled_outcome(decay_only)", unsafe_dec)):
                    if dname.startswith("controlled") and m != 0.0:
                        continue
                    nfs = int((ps & unsafe).sum())
                    rows.append(dict(est=label, seed="truth" if key == "v1" else "selfcal", split=sname, definition=dname, q_margin=m, n_pred_safe=int(ps.sum()), n_false_safe=nfs,
                                     false_safe_rate=nfs / max(1, ps.sum()), kept_truly_safe=int((ps & true_safe).sum()), n_truly_safe=int(true_safe.sum())))
    return pd.DataFrame(rows)


def metrics_side_by_side(M1, M2):
    keep = ["split", "method", "n", "spearman", "kendall", "k", "recall_top", "random_recall", "regret_complete_mps", "regret_safe_mps", "false_safe_rate", "false_decay_rate", "n_pred_safe"]
    a = M1[M1.method.isin(["v1", "v1_selfcal", "B1_hand", "B2_quiet"])][keep].copy(); b = M2[M2.method.isin(["v1", "v1_selfcal"])][keep].copy()
    a["method"] = a.method.replace({"v1": "estimator v1", "v1_selfcal": "estimator v1 selfcal"}); b["method"] = b.method.replace({"v1": "estimator v2", "v1_selfcal": "estimator v2 selfcal"})
    return pd.concat([a, b])


def fig_metrics(S, out):
    splits = list(S.split.unique()); fig, ax = plt.subplots(1, 3, figsize=(15, 4.4)); w = 0.2
    for j, (met, ttl) in enumerate((("spearman", "Spearman"), ("recall_top", "top-10% recall"), ("false_safe_rate", "false-safe rate (literal)"))):
        for i, (m, col) in enumerate((("estimator v1", C1), ("estimator v2", C2), ("estimator v1 selfcal", "#e3a78a"), ("estimator v2 selfcal", "#5aa9c9"))):
            v = [S[(S.split == s) & (S.method == m)][met].values[0] for s in splits]
            ax[j].bar(np.arange(len(splits)) + i * w, v, w, color=col, label=m)
        ax[j].set_xticks(np.arange(len(splits)) + 1.5 * w); ax[j].set_xticklabels(splits, rotation=30, ha="right", fontsize=8); ax[j].set_title(ttl, fontsize=10); ax[j].grid(alpha=0.2, axis="y")
    ax[2].axhline(0, color="k", lw=0.8); ax[0].set_ylim(0.98, 1.0005); ax[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out / "fig_selection_metrics_v1_v2.png", dpi=130); plt.close(fig)


def fig_lookahead(L1, L2, out):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    for L, n, col in ((L1, "v1", C1), (L2, "v2", C2)):
        g = L[L.set == "completed"]; x = [float(str(b).split(",")[0].strip("([")) + 0 for b in g.bin]; x = [(float(str(b).strip("([]").split(",")[0]) + float(str(b).strip("([]").split(",")[1])) / 2 for b in g.bin]
        ax[0].plot(x, g.rms_alt_km, "o-", color=col, label=f"{n} RMS"); ax[0].plot(x, g.max_alt_km, "s--", color=col, alpha=0.5, label=f"{n} max")
        ax[1].plot(x, g.rms_change_km, "o-", color=col, label=f"{n} RMS of change from pass 1"); ax[1].plot(x, g.max_change_km, "s--", color=col, alpha=0.5, label=f"{n} max")
    for a_, t in zip(ax, ("periapsis altitude error", "error of the change from pass 1")):
        a_.set_xlabel("lookahead (days; bin center)"); a_.set_ylabel("km"); a_.set_title(t + ", 203 completed cells, drag on", fontsize=10); a_.grid(alpha=0.2); a_.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out / "fig_error_vs_lookahead_v1_v2.png", dpi=130); plt.close(fig)


# ---------------------------------------------------------------------------------------------- oscillation / amplitude / phase

def design(t, periods, order=2):
    X = [t ** k for k in range(order + 1)]
    for P in periods:
        X += [np.sin(2 * np.pi * t / P), np.cos(2 * np.pi * t / P)]
    return np.array(X).T


def fit_osc(t_d, y, periods):
    """Quadratic + sinusoids at the given periods; returns amplitude and time of maximum (days, in [0, P)) per period, and the residual RMS."""
    X = design(t_d, periods); c, *_ = np.linalg.lstsq(X, y, rcond=None)
    out = []
    for j, P in enumerate(periods):
        al, be = c[3 + 2 * j], c[4 + 2 * j]
        out.append((float(np.hypot(al, be)), float((np.arctan2(al, be) * P / (2 * np.pi)) % P)))
    return out, rms(y - X @ c)


def wrap(x, P):
    return (x + P / 2) % P - P / 2


def osc_analysis(root, est_dir, out):
    S = pd.read_csv(est_dir / "osc_series.csv"); H = pd.read_csv(est_dir / "osc_harmonics.csv")
    sp = pd.read_csv(root / "border3" / "spectral_peaks.csv"); rows = []; sun_rows = []; both_rows = []
    for h in H[H.body == "Moon"].itertuples():
        g = S[(S.tag == h.tag) & (S.set == "moon")]
        t_tr = g.t_truth_s.values / DAY; t_pr = g.t_pred_s.values / DAY
        ytr = g.rp_truth_m.values / 1e3; ypr = (g.rbar_m.values + g.sp_m.values) / 1e3
        periods = [h.P2 / DAY, h.P1 / DAY]
        ftr, rtr = fit_osc(t_tr, ytr, periods); fpr, rpr = fit_osc(t_pr, ypr, periods)
        peak = sp[(sp.tag == h.tag) & (sp.variant == "v8")]
        P2d, P1d = h.P2 / DAY, h.P1 / DAY
        rows.append(dict(tag=h.tag, i_p_deg=np.degrees(h.i_p), omega_p_deg=np.degrees(h.omega_p), e=h.e, a_km=h.a / 1e3, Pi_3B=h.Pi, a_over_b=h.aob, nu=h.nu,
                         R2=h.R2, psidot2_deg_per_day=np.degrees(h.psidot2) * DAY,
                         P2_pred_d=P2d, P2_pred_frozen_d=h.P2f / DAY, P2_free_fit_d=float(peak.peak_period_d.iloc[0]) if len(peak) else np.nan,
                         A2_closed_km=h.A2 / 1e3, A2_closed_frozen_km=h.A2f / 1e3, A2_series_km=fpr[0][0], A2_truth_km=ftr[0][0],
                         ratio_closed=ftr[0][0] / (h.A2 / 1e3), ratio_closed_frozen=ftr[0][0] / (h.A2f / 1e3), ratio_series=ftr[0][0] / fpr[0][0],
                         tmax2_closed_d=(h.tmax2 / DAY) % P2d, tmax2_series_d=fpr[0][1], tmax2_truth_d=ftr[0][1],
                         phase_err_closed_d=wrap(ftr[0][1] - (h.tmax2 / DAY) % P2d, P2d), phase_err_closed_frozen_d=wrap(ftr[0][1] - (h.tmax2f / DAY) % (h.P2f / DAY), h.P2f / DAY),
                         phase_err_series_d=wrap(ftr[0][1] - fpr[0][1], P2d),
                         P1_pred_d=P1d, A1_closed_km=h.A1 / 1e3, A1_closed_frozen_km=h.A1f / 1e3, A1_series_km=fpr[1][0], A1_truth_km=ftr[1][0],
                         phase_err1_closed_d=wrap(ftr[1][1] - (h.tmax1 / DAY) % P1d, P1d), phase_err1_series_d=wrap(ftr[1][1] - fpr[1][1], P1d),
                         resid_rms_truth_km=rtr, resid_rms_series_km=rpr,
                         rms_series_vs_truth_km=rms((ypr - ypr[0]) - (ytr - ytr[0])), rms_truth_change_km=rms(ytr - ytr[0]), rms_v2_est_vs_truth_km=rms((g.rp_v2_m.values - g.rp_v2_m.values[0]) / 1e3 - (ytr - ytr[0])),
                         rms_series_nosp_vs_truth_km=rms((g.rbar_m.values) / 1e3 - (ytr - ytr[0]))))
    O = pd.DataFrame(rows); O.to_csv(out / "osc_amplitude_phase.csv", index=False)
    for h in H[H.body == "Sun"].itertuples():
        g = S[(S.tag == h.tag) & (S.set == "sun")]; ytr = g.rp_truth_m.values / 1e3; ypr = (g.rbar_m.values + g.sp_m.values) / 1e3
        dtr = ytr - ytr[0]; dpr = ypr - ypr[0]
        sun_rows.append(dict(tag=h.tag, n=len(g), Pi_3B=h.Pi, nu=h.nu, lamdot_deg_per_day=np.degrees(h.lamdot) * DAY, A2_closed_km=h.A2 / 1e3, P2_pred_d=h.P2 / DAY,
                             truth_change_km=dtr[-1], pred_change_km=dpr[-1], truth_range_km=np.ptp(dtr), pred_range_km=np.ptp(dpr),
                             rms_pred_vs_truth_km=rms(dpr - dtr), rms_truth_change_km=rms(dtr), slope_through_origin=float((dpr @ dtr) / (dtr @ dtr)) if (dtr @ dtr) > 0 else np.nan,
                             corr=float(np.corrcoef(dpr, dtr)[0, 1]) if np.std(dtr) > 0 and np.std(dpr) > 0 else np.nan))
    SU = pd.DataFrame(sun_rows); SU.to_csv(out / "osc_sun_only.csv", index=False)
    for tag, g in S[S.set == "both"].groupby("tag"):
        ytr = g.rp_truth_m.values / 1e3; dtr = ytr - ytr[0]; dpr = (g.rbar_m.values + g.sp_m.values) / 1e3; dpr = dpr - dpr[0]
        dv2 = (g.rp_v2_m.values - g.rp_v2_m.values[0]) / 1e3
        both_rows.append(dict(tag=tag, n=len(g), truth_change_km=dtr[-1], rms_closed_form_vs_truth_km=rms(dpr - dtr), rms_v2_estimator_vs_truth_km=rms(dv2 - dtr), rms_truth_change_km=rms(dtr)))
    BO = pd.DataFrame(both_rows); BO.to_csv(out / "osc_both_bodies.csv", index=False)
    return O, SU, BO, S, H


def fig_osc(O, S, out):
    fig, ax = plt.subplots(2, 3, figsize=(16, 8.4))
    x = np.arange(len(O)); lab = [t.replace("_E", "\nE") for t in O.tag]
    ax[0, 0].bar(x - 0.2, O.A2_truth_km, 0.2, color="k", label="truth fit (V8)"); ax[0, 0].bar(x, O.A2_closed_km, 0.2, color=C2, label="closed form (Pi form)"); ax[0, 0].bar(x + 0.2, O.A2_series_km, 0.2, color="#5aa9c9", label="fit of predicted series")
    ax[0, 0].set_xticks(x); ax[0, 0].set_xticklabels(lab, fontsize=6); ax[0, 0].set_ylabel("2-lambda amplitude (km)"); ax[0, 0].legend(fontsize=7); ax[0, 0].set_title("Lunar oscillation (half-sidereal) amplitude", fontsize=10); ax[0, 0].grid(alpha=0.2, axis="y")
    ax[0, 1].scatter(O.A2_closed_km, O.A2_truth_km, color=C2, s=30); lim = [0, max(O.A2_closed_km.max(), O.A2_truth_km.max()) * 1.1]; ax[0, 1].plot(lim, lim, "k", lw=0.8)
    for r in O.itertuples(): ax[0, 1].annotate(r.tag[:5] + r.tag[6:10] + r.tag[11:15], (r.A2_closed_km, r.A2_truth_km), fontsize=5)
    ax[0, 1].set_xlabel("predicted amplitude, closed form (km)"); ax[0, 1].set_ylabel("fitted truth amplitude (km)"); ax[0, 1].set_title("Amplitude: predicted vs fitted", fontsize=10); ax[0, 1].grid(alpha=0.2)
    ax[0, 2].bar(x - 0.15, O.phase_err_closed_d, 0.3, color=C2, label="closed-form phase"); ax[0, 2].bar(x + 0.15, O.phase_err_series_d, 0.3, color="#5aa9c9", label="predicted-series phase")
    ax[0, 2].set_xticks(x); ax[0, 2].set_xticklabels(lab, fontsize=6); ax[0, 2].set_ylabel("phase error, truth - predicted (days)"); ax[0, 2].axhline(0, color="k", lw=0.8); ax[0, 2].legend(fontsize=7); ax[0, 2].grid(alpha=0.2, axis="y")
    ax[0, 2].set_title("Phase error of the 2-lambda term (time of maximum)", fontsize=10)
    ax[1, 0].scatter(O.P2_pred_d, O.P2_free_fit_d, color=C2, s=30); ax[1, 0].plot([13.4, 14.2], [13.4, 14.2], "k", lw=0.8); ax[1, 0].set_xlabel("predicted period (d)"); ax[1, 0].set_ylabel("free-fit period (d), attribution study")
    ax[1, 0].set_title("Period (Rayleigh resolution about 3 d)", fontsize=10); ax[1, 0].grid(alpha=0.2)
    for j, tag in enumerate(["I0234_R000_W090_E1", "I0634_R270_W090_E2"]):
        g = S[(S.tag == tag) & (S.set == "moon")]; a_ = ax[1, 1 + j]
        a_.plot(g.t_truth_s / DAY, (g.rp_truth_m - g.rp_truth_m.iloc[0]) / 1e3, "k", label="truth (Moon only, drag off)")
        a_.plot(g.t_pred_s / DAY, (g.rbar_m + g.sp_m - (g.rbar_m + g.sp_m).iloc[0]) / 1e3, color=C2, label="Pi-form: mean + short period"); a_.plot(g.t_pred_s / DAY, g.rbar_m / 1e3, color="#5aa9c9", ls="--", lw=1, label="Pi-form mean only")
        a_.set_title(tag, fontsize=9); a_.set_xlabel("time (d)"); a_.set_ylabel("r_p - r_p(pass 1) (km)"); a_.grid(alpha=0.2); a_.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out / "fig_oscillation_amp_phase.png", dpi=130); plt.close(fig)


def fig_collapse(O, H, SU, P, out):
    """Amplitude normalized by a Pi_3B R2 vs nu: theory A/(a Pi R2) = 1/(2 nu); Earth data (Moon fits) and predictions for Earth Sun, Mars, Venus."""
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    nu = np.logspace(-3.2, -0.6, 100); ax[0].plot(nu, 1 / (2 * nu), "k", lw=1, label="theory 1/(2 nu)")
    Hm = H[H.body == "Moon"].set_index("tag")
    for r in O.itertuples():
        h = Hm.loc[r.tag]; R2 = np.hypot(h.c2, h.s2)
        if R2 > 1.0:
            ax[0].scatter(r.nu, r.A2_truth_km * 1e3 / (h.a * h.Pi * R2), color=C2, s=36, zorder=3)
    ax[0].scatter([], [], color=C2, label="Earth-Moon, fitted truth amplitude (cells with 2-lambda harmonic R2 > 1)")
    Hs = H[H.body == "Sun"]
    for h in Hs.itertuples():
        ax[0].scatter(abs(h.nu), 1 / (2 * abs(h.nu)), facecolors="none", edgecolors=C1, s=30, zorder=2)
    ax[0].scatter([], [], facecolors="none", edgecolors=C1, label="Earth-Sun, theory only (period > record: not fit)")
    mk = {"Earth": "s", "Mars": "^", "Venus": "v"}
    for r in P.itertuples():
        ax[0].scatter(r.nu, 1 / (2 * r.nu), marker=mk[r.planet], color="#2a7f3f", s=22, alpha=0.8, zorder=2)
    ax[0].scatter([], [], marker="^", color="#2a7f3f", label="planets table (Earth, Mars, Venus, theory)")
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("nu = lambda_dot / n  (perturber-to-spacecraft rate)"); ax[0].set_ylabel("A / (a Pi_3B R2)")
    ax[0].set_title("Amplitude collapse: A = a Pi_3B R2 / (2 nu)", fontsize=10); ax[0].legend(fontsize=6, loc="upper right"); ax[0].grid(alpha=0.2, which="both")
    # (b) all passes, Moon-only: normalized predicted vs measured change
    S = pd.read_csv(out / "osc_series.csv")
    for body, tv, col in (("moon", "Moon", C2), ("sun", "Sun", C1)):
        g = S[S.set == body]
        for tag, gg in g.groupby("tag"):
            h = H[(H.tag == tag) & (H.body == tv)].iloc[0]; norm_ = h.a * h.Pi
            d_tr = (gg.rp_truth_m - gg.rp_truth_m.iloc[0]).values / norm_; d_pr = ((gg.rbar_m + gg.sp_m) - (gg.rbar_m + gg.sp_m).iloc[0]).values / norm_
            ax[1].plot(d_pr, d_tr, ".", ms=3, color=col, alpha=0.5)
    ax[1].plot([], [], ".", color=C2, label="Moon only (V8), all passes of 12 cells"); ax[1].plot([], [], ".", color=C1, label="Sun only (V9)")
    lim = np.array(ax[1].get_xlim()); lim = [min(lim[0], ax[1].get_ylim()[0]), max(lim[1], ax[1].get_ylim()[1])]; ax[1].plot(lim, lim, "k", lw=0.8)
    ax[1].set_xlabel("predicted  Delta r_p / (a Pi_3B)"); ax[1].set_ylabel("measured  Delta r_p / (a Pi_3B)"); ax[1].set_title("Response collapse, both perturbers (Pi_3B differs 2x, nu 13x)", fontsize=10)
    ax[1].legend(fontsize=7); ax[1].grid(alpha=0.2)
    # (c) planets: amplitude scale in km
    for r in P.itertuples():
        pass
    lab = [f"{r.planet}-{r.perturber}\na={r.a_km:.0f}" for r in P.itertuples()]
    ax[2].barh(np.arange(len(P)), P.A2max_km, color="#2a7f3f"); ax[2].set_yticks(np.arange(len(P))); ax[2].set_yticklabels(lab, fontsize=5); ax[2].set_xscale("log")
    ax[2].set_xlabel("maximum 2-lambda amplitude over geometry (km), A = a Pi R2max / (2 nu)"); ax[2].grid(alpha=0.2, axis="x"); ax[2].set_title("Predicted amplitude scales (illustrative orbits)", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "fig_scaling_collapse.png", dpi=130); plt.close(fig)


def attrib_compare(est_dir, out):
    A = pd.read_csv(est_dir / "attrib_series_all.csv"); rows = []; expl = []
    for (v, tag), g in A.groupby(["variant", "tag"], sort=False):
        g1 = g[g.est == "v1"].reset_index(drop=True); g2 = g[g.est == "v2"].reset_index(drop=True)
        n = min(len(g1), len(g2))
        e1 = (g1.rp_pred_m[:n] - g1.rp_truth_m[:n]).values; e2 = (g2.rp_pred_m[:n] - g2.rp_truth_m[:n]).values
        c1 = ((g1.rp_pred_m[:n] - g1.rp_pred_m[0]) - (g1.rp_truth_m[:n] - g1.rp_truth_m[0])).values; c2 = ((g2.rp_pred_m[:n] - g2.rp_pred_m[0]) - (g2.rp_truth_m[:n] - g2.rp_truth_m[0])).values
        rows.append(dict(variant=v, tag=tag, passes=n, rms_rp_v1_km=rms(e1) / 1e3, rms_rp_v2_km=rms(e2) / 1e3, max_rp_v1_km=np.abs(e1).max() / 1e3, max_rp_v2_km=np.abs(e2).max() / 1e3,
                         rms_change_v1_km=rms(c1) / 1e3, rms_change_v2_km=rms(c2) / 1e3, max_change_v1_km=np.abs(c1).max() / 1e3, max_change_v2_km=np.abs(c2).max() / 1e3,
                         pass1_v1_km=e1[0] / 1e3, pass1_v2_km=e2[0] / 1e3))
        corr = (g2.rp_pred_m[:n] - g1.rp_pred_m[:n]).values; resid = -e1       # truth - v1
        if np.std(corr) > 0:
            sl = float((corr - corr.mean()) @ (resid - resid.mean()) / ((corr - corr.mean()) @ (corr - corr.mean())))
            expl.append(dict(variant=v, tag=tag, slope_resid_on_correction=sl, corr=float(np.corrcoef(corr, resid)[0, 1]), rms_v1_resid_km=rms(resid - resid.mean()) / 1e3, rms_after_km=rms(resid - corr - (resid - corr).mean()) / 1e3))
    R = pd.DataFrame(rows); X = pd.DataFrame(expl)
    R.to_csv(out / "attribution_cells_v1_v2.csv", index=False); X.to_csv(out / "sp_explained_variance.csv", index=False)
    return R, X


def blind_score(root, out):
    """Score v2 (self-calibrated seed, no tuning) on the 45 blind truth cells and compare with v1's picks (columns pred_* of blind_cell_table.csv)."""
    BT = pd.read_csv(root / "blind" / "blind_cell_table.csv"); pc = pd.read_csv(out / "blind_pred_cells.csv").set_index("tag"); rows = []
    for c in BT.itertuples():
        df = pd.read_csv(out / "blind_pred" / f"{c.tag}.csv"); row, _ = SK.predicted_row(df, "est")
        dec = bool(pc.loc[c.tag, "decayed_pred"]); cls = "decayed" if dec else ("breach" if row["q_max"] > Q_HARD else "safe")
        rows.append(dict(tag=c.tag, group=c.group, true_class=c.true_class, v1_class=c.pred_class, v2_class=cls, v1_q=c.pred_q_max, v2_q=row["q_max"], true_q=c.true_q_max,
                         v1_dv=c.pred_dv_mps, v2_dv=row["lazy_dv"], true_dv=c.true_lazy_s0_mps if c.status == "completed" else np.nan, true_cost=c.true_cost_mps, status=c.status,
                         v1_tdec=c.pred_t_decay_days, v2_tdec=pc.loc[c.tag, "t_decay_pred_days"], true_tend=c.t_end_days, true_unsafe=bool(c.true_unsafe)))
    D = pd.DataFrame(rows); D.to_csv(out / "blind_v1_v2.csv", index=False)
    res = []
    for e in ("v1", "v2"):
        pred_safe = D[f"{e}_class"] == "safe"; comp = (D.status == "completed") & (D[f"{e}_class"] != "decayed")
        err = (D[f"{e}_dv"] - D.true_dv)[comp]
        qr = (D[f"{e}_q"] / D.true_q)[(D.status == "completed")]
        sp = stats.spearmanr(D[f"{e}_dv"][comp], D.true_cost[comp])[0]; kt = stats.kendalltau(D[f"{e}_dv"][comp], D.true_cost[comp])[0]
        fs_lit = int((pred_safe & D.true_unsafe).sum()); fs_dec = int((pred_safe & (D.status == "decayed")).sum())
        marg = {f"false_safe_margin_{int(m*100)}": int(((D[f"{e}_class"] != "decayed") & (D[f"{e}_q"] <= Q_HARD / (1 + m)) & D.true_unsafe).sum()) for m in (0, .05, .10, .15)}
        res.append(dict(est=e, n=len(D), classes_correct=int((D[f"{e}_class"] == D.true_class).sum()), n_pred_safe=int(pred_safe.sum()), false_safe_literal=fs_lit, false_safe_decay_only=fs_dec,
                        n_completed_pred_complete=int(comp.sum()), dv_rms_err_mps=rms(err), dv_median_abs_err_mps=float(np.median(np.abs(err))), dv_max_abs_err_mps=float(np.abs(err).max()),
                        spearman=sp, kendall=kt, q_ratio_median=float(qr.median()), q_ratio_min=float(qr.min()), q_ratio_max=float(qr.max()),
                        decay_time_rms_d=rms((D[f"{e}_tdec"] - D.true_tend)[D.status == "decayed"].dropna()), **marg))
    R = pd.DataFrame(res); R.to_csv(out / "blind_scores_v1_v2.csv", index=False); print(R.round(4).T.to_string())
    print(pd.crosstab(D.true_class, D.v2_class)); print(D[D.v1_class != D.v2_class][["tag", "true_class", "v1_class", "v2_class"]])
    return D, R


# ---------------------------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap)
    ap.add_argument("--phase", default="all"); ap.add_argument("--blind-dir", default=None)
    a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root)); v2 = root / "estimator_v2"; v1 = root / "selection"
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_rows", 400)
    SK.ROOT_LABELS = root / "oracle" / "oracle_cell_labels.csv"
    if a.phase in ("all", "score"):
        score(root, v2, variants=("", "quadrupole", "j2quad"))
    if a.phase in ("all", "compare"):
        dirs = {"v1": v1, "v2": v2}
        P = pass1_offsets(root, dirs); P.to_csv(v2 / "pass1_offsets.csv", index=False)
        summ = []
        for n in ("v1", "v2"):
            for k in ("off_m", "selfcal_off_m"):
                x = P[f"{n}_{k}"].values
                summ.append(dict(est=n, seed="truth seed" if k == "off_m" else "self-calibrated seed", mean_km=x.mean() / 1e3, rms_km=rms(x) / 1e3, std_km=x.std() / 1e3, min_km=x.min() / 1e3, max_km=x.max() / 1e3, abs_gt_0p5_km=int((np.abs(x) > 500).sum()), abs_gt_1_km=int((np.abs(x) > 1000).sum())))
        PS = pd.DataFrame(summ); PS.to_csv(v2 / "pass1_offset_summary.csv", index=False); print(PS.round(4).to_string())
        # by inclination
        print(P.groupby("inc_deg")[["v1_off_m", "v2_off_m"]].agg(lambda x: rms(x) / 1e3).round(3))
        fig_pass1(P, v2, (("v1", C1), ("v2", C2)))
        T1 = pd.read_csv(v1 / "selection_cell_table.csv"); T2 = pd.read_csv(v2 / "selection_cell_table.csv")
        M1 = pd.read_csv(v1 / "selection_metrics.csv"); M2 = pd.read_csv(v2 / "selection_metrics.csv")
        SBS = metrics_side_by_side(M1, M2); SBS.to_csv(v2 / "selection_metrics_v1_v2.csv", index=False); print(SBS.round(3).to_string()); fig_metrics(SBS, v2)
        FS = pd.concat([false_safe_table(T1, "v1"), false_safe_table(T2, "v2")]); FS.to_csv(v2 / "false_safe_definitions_v1_v2.csv", index=False)
        print(FS[(FS.seed == "truth")].pivot_table(index=["split", "definition", "q_margin"], columns="est", values=["n_false_safe", "n_pred_safe", "false_safe_rate"], aggfunc="first").round(3).to_string())
        TH1 = pd.read_csv(v1 / "threshold_table.csv"); TH2 = pd.read_csv(v2 / "threshold_table.csv")
        TH = TH1.merge(TH2, on=["split", "metric"], suffixes=("_v1", "_v2")); TH.to_csv(v2 / "threshold_table_v1_v2.csv", index=False); print(TH[["split", "metric", "value_v1", "passed_v1", "value_v2", "passed_v2"]].round(4).to_string())
        L1 = pd.read_csv(v1 / "lookahead_error.csv"); L2 = pd.read_csv(v2 / "lookahead_error.csv")
        LL = L1.merge(L2, on=["set", "bin"], suffixes=("_v1", "_v2")); LL.to_csv(v2 / "lookahead_error_v1_v2.csv", index=False); print(LL.round(3).to_string()); fig_lookahead(L1, L2, v2)
        DV1 = pd.read_csv(v1 / "dv_error_table.csv"); DV2 = pd.read_csv(v2 / "dv_error_table.csv"); print(DV1.round(3).to_string()); print(DV2.round(3).to_string())
        DT1 = pd.read_csv(v1 / "decay_time_error.csv"); DT2 = pd.read_csv(v2 / "decay_time_error.csv"); print(DT1.round(3).to_string()); print(DT2.round(3).to_string())
        for k in ("v1", "v1_selfcal"):
            print(k, "v1 confusion\n", pd.read_csv(v1 / f"confusion_{k}.csv", index_col=0), "\nv2 confusion\n", pd.read_csv(v2 / f"confusion_{k}.csv", index_col=0))
        R, X = attrib_compare(v2, v2); print(R.round(3).to_string()); print(X.round(3).to_string())
        # variants
        for var in ("quadrupole", "j2quad"):
            d = root / f"estimator_v2_{var}"
            if (d / "selection_metrics.csv").exists():
                M = pd.read_csv(d / "selection_metrics.csv"); Pv = pass1_offsets(root, {"v2": d, "v1": v1}); print(var, "pass-1 RMS km", rms(Pv["v2_off_m"]) / 1e3, "mean", Pv["v2_off_m"].mean() / 1e3)
                print(M[(M.method == "v1") & (M.split == "pooled")][["spearman", "kendall", "recall_top", "regret_complete_mps", "false_safe_rate"]].round(4).to_string())
                print(false_safe_table(pd.read_csv(d / "selection_cell_table.csv"), var).query("seed=='truth' and split=='pooled' and q_margin==0").to_string())
    if a.phase in ("all", "blind") and (v2 / "blind_pred_cells.csv").exists():
        blind_score(root, v2)
    if a.phase in ("all", "osc"):
        O, SU, BO, S, H = osc_analysis(root, v2, v2)
        print(O.round(3).to_string()); print(SU.round(3).to_string()); print(BO.round(3).to_string())
        fig_osc(O, S, v2)
        Pl = pd.read_csv(v2 / "planets_table.csv") if (v2 / "planets_table.csv").exists() else pd.DataFrame(columns=["planet", "nu", "a_km", "perturber", "A2max_km"])
        fig_collapse(O, H, SU, Pl, v2)


if __name__ == "__main__":
    main()
