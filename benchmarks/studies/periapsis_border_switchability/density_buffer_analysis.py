"""Density-buffer analysis (pre-registration: tasks/2026-10-15_periapsis-region-border-targeting/density_buffer_preregistration.md).
Can a k-sigma buffer from previously measured densities replace a fixed periapsis-altitude margin? Venus-GRAM MEAN density truth passes; Python only.
Usage: python3 -I density_buffer_analysis.py --results-root GRAM_ROOT [--selftest] -> <GRAM_ROOT>/venus/density_buffer/*.csv and a printed summary.
"""
import argparse, pathlib, sys, tomllib
from statistics import NormalDist
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
FIT = tomllib.load(open(HERE / "venus_gram_config.toml", "rb"))["fit"]
BAND = tomllib.load(open(HERE / "study_config.toml", "rb"))["band"]
Q_LO, Q_HI = BAND["q_low_w_cm2"], BAND["q_high_w_cm2"]; W_Q = np.log(Q_HI / Q_LO)
KS = (1.0, 1.64, 2.33); NOM = {k: NormalDist().cdf(k) for k in KS}
MARGIN_M = 1000.0; EWMA_N = (3, 5, 10); MIN_RUN_ERR = 5
PREDS = ["P0s", "P1"] + [f"P2_{n}" for n in EWMA_N] + ["P3"]


def predict_run(y, lon):
    """One-pass-ahead predictions for targets j = 1..len-1 of one run from passes < j only. Returns dict name -> array (index j; NaN at j = 0)."""
    n = len(y); out = {p: np.full(n, np.nan) for p in PREDS}
    dl = np.angle(np.exp(1j * np.diff(lon)))   # wrapped longitude increments, dl[k-1] = lon[k]-lon[k-1]
    ew = {N: y[0] for N in EWMA_N}
    for j in range(1, n):
        out["P0s"][j] = 0.0; out["P1"][j] = y[j - 1]
        for N in EWMA_N:
            a = 2.0 / (N + 1); ew[N] = a * y[j - 1] + (1 - a) * ew[N] if j > 1 else y[0]; out[f"P2_{N}"][j] = ew[N]
        # P3: pairs k = max(1, j-10)..j-1 use passes <= j-1 only
        ks = np.arange(max(1, j - 10), j)
        if len(ks) >= 4:
            dy = y[ks] - y[ks - 1]; dx = dl[ks - 1]; m = float(dy @ dx / (dx @ dx)) if dx @ dx > 0 else 0.0
            out["P3"][j] = y[j - 1] + m * dl[j - 2]
        else: out["P3"][j] = y[j - 1]
    return out


def load(root, window):
    v = pathlib.Path(root) / "venus"; cells = pd.read_csv(v / "cells.csv").set_index("tag"); rs = pd.read_csv(v / "truth_run_summary.csv").set_index("tag"); rows = []
    for t in cells.index:
        T = pd.read_csv(v / "truth" / t / "passes_3b.csv"); hg = cells.loc[t, "h_guess_m"]
        T["tag"] = t; T["done"] = rs.loc[t, "retcode"] == "Success"
        T["inw"] = (T.alt_m >= hg - FIT["window_below_m"]) & (T.alt_m <= hg + FIT["window_above_m"])
        T["y"] = np.log(T.rho / T.rho_pw); T["qmod"] = 0.5 * T.rho_pw * T.vrel_mps ** 3 / 1e4; rows.append(T)
    return pd.concat(rows, ignore_index=True)


def build(P, cellset, use_window):
    recs = []
    for t, T in P.groupby("tag", sort=False):
        if t not in cellset: continue
        T = T[T.inw] if use_window else T
        brk = np.r_[True, np.diff(T["pass"].values) != 1]; rid = np.cumsum(brk)
        for r in np.unique(rid):
            R = T[rid == r]
            if len(R) < 2: continue
            pr = predict_run(R.y.values, np.unwrap(R.lon_rad.values))
            for j in range(1, len(R)):
                x = R.iloc[j]; d = dict(tag=t, run=int(r), pass_=int(x["pass"]), j=j, y=x.y, H=x.h_local_m, qmod=x.qmod, q=x.q_w_cm2, lon=x.lon_rad, alt_m=x.alt_m)
                for p in PREDS: d["yh_" + p] = pr[p][j]; d["e_" + p] = x.y - pr[p][j]
                recs.append(d)
    return pd.DataFrame(recs)


def rms(x): x = np.asarray(x, float); return float(np.sqrt(np.mean(x ** 2))) if len(x) else np.nan


def score(D):
    """Per predictor and per sigma estimator: sigma, coverage, alarm rate, buffer, decay factor. Returns long DataFrame (cell = 'POOLED' or tag)."""
    out = []; cells = list(D.tag.unique()); D = D.sort_values(["tag", "pass_"]).copy()
    sig = {}
    for p in PREDS:
        e = "e_" + p
        D["sl_" + p] = [rms(D.loc[D.tag != t, e]) for t in D.tag]    # leave-one-cell-out sigma
        run = np.full(len(D), np.nan)
        for t in cells:
            idx = np.where(D.tag.values == t)[0]; ev = D[e].values[idx]
            for a, i in enumerate(idx):
                if a >= MIN_RUN_ERR: run[i] = rms(ev[:a])
        D["sr_" + p] = run
    for est in ("sl", "sr"):
        for p in PREDS:
            for k in KS:
                s = D[f"{est}_{p}"]; ok = s.notna(); b = k * s; e = D["e_" + p]; yh = D["yh_" + p]
                cov = (e <= b); qup = D.qmod * np.exp(yh + b)
                for cell, m in [("POOLED", ok)] + [(t, ok & (D.tag == t)) for t in cells]:
                    if m.sum() == 0: continue
                    out.append(dict(est=est, pred=p, k=k, cell=cell, n=int(m.sum()), sigma=rms(e[m]), coverage=float(cov[m].mean()), nominal=NOM[k], alarm=float((qup[m] > Q_HI).mean()),
                                    buf_ln=float(b[m].mean()), buf_km=float((D.H * b)[m].mean() / 1e3), buf_frac_band=float(b[m].mean() / W_Q), decay_factor=float(np.exp(-b[m]).mean())))
    # P0f fixed margin
    b = MARGIN_M / D.H; cov = D.y <= b; qup = D.qmod * np.exp(b)
    for cell, m in [("POOLED", np.ones(len(D), bool))] + [(t, (D.tag == t).values) for t in cells]:
        out.append(dict(est="fixed", pred="P0f", k=np.nan, cell=cell, n=int(m.sum()), sigma=rms(D.y[m]), coverage=float(cov[m].mean()), nominal=np.nan, alarm=float((qup[m] > Q_HI).mean()),
                        buf_ln=float(b[m].mean()), buf_km=MARGIN_M / 1e3, buf_frac_band=float(b[m].mean() / W_Q), decay_factor=float(np.exp(-b[m]).mean())))
    return pd.DataFrame(out), D


def lag1(D, col):
    a, b = [], []
    for (t, r), G in D.groupby(["tag", "run"]):
        if len(G) < MIN_RUN_ERR: continue
        v = G[col].values - D.loc[D.tag == t, col].mean(); a += list(v[:-1]); b += list(v[1:])
    return float(np.corrcoef(a, b)[0, 1]) if len(a) > 2 else np.nan


def periodogram(x):
    x = np.asarray(x, float) - np.mean(x); P = np.abs(np.fft.rfft(x)) ** 2; f = np.fft.rfftfreq(len(x))
    i = 1 + int(np.argmax(P[1:])); return 1.0 / f[i], P[i] / P[1:].sum(), float(P[1:][f[1:] >= 0.4].sum() / P[1:].sum())


def structure(D):
    rows = []
    for t, G in D.groupby("tag"):
        G = G.sort_values("pass_"); r = G.groupby("run").size(); big = r.idxmax()
        if r.max() < 16: continue
        R = G[G.run == big]; dl = np.degrees(np.abs(np.angle(np.exp(1j * np.diff(np.unwrap(R.lon.values))))))
        py, fy, hy = periodogram(R.y.values); pe, fe, he = periodogram(R.e_P1.values)
        rows.append(dict(tag=t, run_len=len(R), lon_period_passes=360 / np.median(dl), y_peak_period=py, y_peak_share=fy, y_share_above_0p4=hy, e_peak_period=pe, e_peak_share=fe, e_share_above_0p4=he))
    return pd.DataFrame(rows)


def verdict(D, S):
    """Pre-registered T1-T8 (primary set C15_window). PREREG values from the pre-registration file; nothing here is tuned."""
    pool = S[(S.cell == "POOLED") & (S.est == "sl")]; g = lambda p, k: pool[(pool.pred == p) & (pool.k == k)].iloc[0]
    sg = {p: g(p, 1.0).sigma for p in PREDS}; V = {}
    V["T1"] = {k: (g("P1", k).coverage, NOM[k], abs(g("P1", k).coverage - NOM[k]) <= 0.05) for k in KS}
    q = S[(S.est == "sl") & (S.pred == "P1") & (S.k == 1.64) & (S.cell != "POOLED") & (S.n >= 20)]
    V["T2"] = dict(cells_with_20=len(q), within_10pts=int((abs(q.coverage - NOM[1.64]) <= 0.10).sum()), required=12, evaluable=len(q) >= 12)
    V["T3"] = {k: (k * sg["P1"], 0.5 * W_Q, k * sg["P1"] <= 0.5 * W_Q) for k in (1.64, 2.33)}
    V["T4"] = (sg["P1"] / sg["P0s"], sg["P1"] / sg["P0s"] <= 0.5)
    V["T5"] = {p: (sg[p] / sg["P1"], sg[p] / sg["P1"] >= 0.9) for p in PREDS if p.startswith("P2")}
    V["T6"] = (sg["P3"] / sg["P1"] - 1, abs(sg["P3"] / sg["P1"] - 1) < 0.10)
    keq = (MARGIN_M / D.H) / sg["P1"]; f = S[(S.pred == "P0f") & (S.cell == "POOLED")].iloc[0]
    V["T7"] = dict(p0f_coverage=f.coverage, below_95=f.coverage < 0.95, keq_median=float(keq.median()), keq_ge_2p33=bool(keq.median() >= 2.33))
    st = structure(D); V["T8"] = dict(lag1_y=lag1(D, "y"), lag1_e=lag1(D, "e_P1"), cells_periodogram=len(st), y_peak_at_run_length=int((st.y_peak_period >= 0.5 * st.run_len).sum()), y_share_above_0p4_max=float(st.y_share_above_0p4.max()))
    return V


def unfrozen(D):
    """UNFROZEN diagnostics added after the first run: error quantiles, tail concentration, robust sigma and empirical-quantile buffer (leave-one-cell-out)."""
    e = D.e_P1.values; o = {}
    o["abs_e_quantiles_50_90_99_max"] = [float(np.quantile(np.abs(e), x)) for x in (.5, .9, .99, 1)]
    srt = np.sort(e ** 2)[::-1]; o["share_sumsq_top_1pct_of_passes"] = float(srt[:max(1, len(e) // 100)].sum() / srt.sum())
    o["share_sumsq_in_passes_with_q_above_qhi"] = float((e[D.q.values > Q_HI] ** 2).sum() / (e ** 2).sum()); o["fraction_actual_q_above_qhi"] = float((D.q > Q_HI).mean())
    rob, emp = {}, {}
    for k in KS:
        c1, c2, b1, b2 = [], [], [], []
        for t in D.tag.unique():
            oth = D.e_P1[D.tag != t].values; mine = D[D.tag == t]
            s = 1.4826 * np.median(np.abs(oth - np.median(oth))); qk = np.quantile(oth, NOM[k])
            c1 += list(mine.e_P1 <= k * s); c2 += list(mine.e_P1 <= qk); b1 += list(mine.H * k * s / 1e3); b2 += list(mine.H * qk / 1e3)
        rob[k] = (float(np.mean(c1)), float(np.mean(b1))); emp[k] = (float(np.mean(c2)), float(np.mean(b2)))
    o["robust_sigma_(coverage,mean_buffer_km)"] = rob; o["empirical_quantile_(coverage,mean_buffer_km)"] = emp
    return o


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--results-root", required=True); ap.add_argument("--selftest", action="store_true"); a = ap.parse_args()
    root = pathlib.Path(a.results_root); out = root / "venus" / "density_buffer"; out.mkdir(exist_ok=True)
    P = load(root, True); done = set(P[P.done].tag); allc = set(P.tag)
    sets = {"C15_window": (done, True), "S2_all32_window": (allc, True), "S3_C15_nowindow": (done, False)}
    for nm, (cs, w) in sets.items():
        D = build(P, cs, w); S, D = score(D); S.to_csv(out / f"scores_{nm}.csv", index=False); D.to_csv(out / f"errors_{nm}.csv", index=False)
        st = structure(D); st.to_csv(out / f"structure_{nm}.csv", index=False)
        print(f"== {nm}: cells {D.tag.nunique()}, scored passes {len(D)}, lag1 y {lag1(D, 'y'):.3f}, lag1 e(P1) {lag1(D, 'e_P1'):.3f}")
        pool = S[(S.cell == "POOLED") & (S.est.isin(["sl", "fixed"]))]
        if nm == "C15_window":
            import json; V = verdict(D, S); U = unfrozen(D); print("VERDICT", json.dumps(V, default=str, indent=1)); print("UNFROZEN", json.dumps(U, default=str, indent=1))
            (out / "verdict_and_unfrozen.json").write_text(json.dumps(dict(verdict=V, unfrozen=U), default=str, indent=1))
        print(pool[["pred", "k", "n", "sigma", "coverage", "nominal", "alarm", "buf_km", "buf_frac_band", "decay_factor"]].round(4).to_string(index=False))
    if a.selftest:
        y = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8.]); lon = np.arange(9) * 0.1; p = predict_run(y, lon)
        assert p["P1"][3] == 2 and abs(p["P3"][8] - 8.0) < 1e-9      # exact linear y(lon): P3 extrapolates y_n + m dlon = y_{n+1} (target 8 from y[7]=7, m = 10)
        y2 = y.copy(); y2[8] = 99; assert predict_run(y2, lon)["P3"][8] == p["P3"][8] and predict_run(y2, lon)["P2_3"][8] == p["P2_3"][8]   # no future leakage
        print("selftest ok")


if __name__ == "__main__": main()
