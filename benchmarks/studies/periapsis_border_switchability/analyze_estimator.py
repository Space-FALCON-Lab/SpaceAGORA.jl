"""Part 3 statistics and figure from estimator_series.csv. Usage: analyze_estimator.py OUT"""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = sys.argv[1]; S = pd.read_csv(f"{OUT}/estimator_series.csv")
def detr(t, y): return y - np.polyval(np.polyfit(t, y, 2), t)
def sine_fit(t, y, P):
    X = np.array([np.sin(2*np.pi*t/P), np.cos(2*np.pi*t/P), np.ones_like(t)]).T
    c = np.linalg.lstsq(X, y, rcond=None)[0]; return np.hypot(c[0], c[1]), np.arctan2(c[1], c[0])
rows = []
for (tag, v), g in S.groupby(["tag", "variant"]):
    t = g.time_s.values / 86400; tr = g.rp_truth.values / 1e3; tr = tr - tr[0]
    d = dict(tag=tag, variant=v, n=len(g), t_end=t[-1], truth_change=tr[-1], truth_range=tr.max()-tr.min())
    res_t = detr(t, tr + 0)
    from scipy.signal import lombscargle
    per = np.linspace(8, 30, 800); p = lombscargle(t, res_t - res_t.mean(), 2*np.pi/per); Pk = per[np.argmax(p)]
    d["truth_peak_period"] = Pk; at, ph_t = sine_fit(t, res_t, Pk); d["truth_osc_amp"] = at
    for m in ("rp_sa", "rp_sq", "rp_sa_start", "rp_hand"):
        y = g[m].values / 1e3; y = y - y[0]; err = y - tr
        d[f"{m}_rms"] = np.sqrt(np.mean(err**2)); d[f"{m}_max"] = np.abs(err).max(); d[f"{m}_final_err"] = err[-1]
        d[f"{m}_change"] = y[-1]
        ap, ph = sine_fit(t, detr(t, y), Pk); d[f"{m}_osc_amp"] = ap
        d[f"{m}_osc_phase_deg"] = np.degrees((ph - ph_t + np.pi) % (2*np.pi) - np.pi)
        d[f"{m}_corr"] = np.corrcoef(detr(t, y), res_t)[0, 1] if np.std(detr(t, y)) > 0 else np.nan
    # one-step: per-pass increments
    inc_t = np.diff(g.rp_truth.values) / 1e3; 
    for m in ("rp_1step", "rp_1step_q"):
        y = g[m].values[1:] / 1e3; d[f"{m}_inc_rms_err"] = np.sqrt(np.mean((y - g.rp_truth.values[1:]/1e3)**2))
    d["truth_inc_rms"] = np.sqrt(np.mean(inc_t**2))
    rows.append(d)
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/estimator_stats.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
cols = ["tag", "variant", "t_end", "truth_change", "truth_range", "rp_sa_rms", "rp_sa_max", "rp_sa_final_err", "rp_sq_rms", "rp_hand_rms", "rp_hand_max"]
print(R[R.variant.isin(["v5", "v8", "v9"])][cols].round(2).to_string())
print(R[R.variant.isin(["v5", "v8", "v9"])][["tag","variant","truth_peak_period","truth_osc_amp","rp_sa_osc_amp","rp_sa_osc_phase_deg","rp_sa_corr","rp_sq_osc_amp","rp_hand_osc_amp","rp_1step_inc_rms_err","truth_inc_rms"]].round(2).to_string())
print(R[R.variant=="v0"][cols].round(2).to_string())
# J2 check, v3
for tag, g in S[S.variant == "v3"].groupby("tag"):
    print(tag, "dOmega truth/pred deg", np.degrees(g.raan_truth.iloc[-1]-g.raan_truth.iloc[0]).round(3), np.degrees(g.raan_sa.iloc[-1]-g.raan_sa.iloc[0]).round(3),
          "domega", np.degrees(g.om_truth.iloc[-1]-g.om_truth.iloc[0]).round(3), np.degrees(g.om_sa.iloc[-1]-g.om_sa.iloc[0]).round(3))
# figures: 4 representative cells, drag-off Sun+Moon
sel = ["I0234_R000_W090_E1", "I0234_R000_W270_E1", "I1166_R000_W090_E2", "I0500_R270_W090_E3"]
fig, ax = plt.subplots(4, 2, figsize=(13, 12))
for k, tag in enumerate(sel):
    g = S[(S.tag == tag) & (S.variant == "v5")]; t = g.time_s / 86400
    for m, lab, c in (("rp_truth", "truth (drag off, Sun+Moon)", "k"), ("rp_sa", "single-averaged exact", "tab:red"), ("rp_sq", "single-averaged quadrupole", "tab:purple"), ("rp_hand", "double-averaged hand", "tab:green")):
        y = g[m] / 1e3; ax[k, 0].plot(t, y - y.iloc[0], c, label=lab, lw=1.4)
    ax[k, 0].set_title(tag, fontsize=9); ax[k, 0].set_ylabel("r_p - r_p0 (km)"); ax[k, 0].grid(alpha=.3)
    for m, c in (("rp_sa", "tab:red"), ("rp_sq", "tab:purple"), ("rp_hand", "tab:green")):
        y = (g[m] - g.rp_truth) / 1e3; y = y - y.iloc[0]; ax[k, 1].plot(t, y, c)
    ax[k, 1].set_ylabel("prediction error (km)"); ax[k, 1].grid(alpha=.3)
ax[0, 0].legend(fontsize=7); ax[-1, 0].set_xlabel("time (d)"); ax[-1, 1].set_xlabel("time (d)"); fig.tight_layout(); fig.savefig(f"{OUT}/fig_estimator_vs_truth.png", dpi=110)
