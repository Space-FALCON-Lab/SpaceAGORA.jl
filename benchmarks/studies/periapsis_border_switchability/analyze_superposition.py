"""Part 3 (task B1): superposition check, Python half.

Usage (numpy, pandas, scipy, matplotlib; python3 -I):
  python3 -I analyze_superposition.py select   [--results-root DIR]
  python3 -I analyze_superposition.py schedule --iter N [--results-root DIR]     # N = 0: from the passive truth; N >= 1: repair increments
                                                                                 # from the closed-loop histories of iteration N-1
  python3 -I analyze_superposition.py analyze  [--results-root DIR] [--finite-sched NAME]

`select` writes <root>/superposition/sup_cells.csv by the rule of superposition_config.toml. `schedule` writes
<root>/superposition/schedules/<tag>__<sched>__it<N>.csv (columns pass, dv_mps; signed tangential dV before pass `pass`; pass 1 is always 0
for the policies used) with the A1 library (oracle_lib). `analyze` reads the Julia outputs (run_superposition.jl) and writes tables and figures.
"""
import argparse, importlib.util, pathlib, sys, tomllib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root
spec = importlib.util.spec_from_file_location("oracle_lib", HERE / "oracle_lib.py")
lib = importlib.util.module_from_spec(spec); sys.modules["oracle_lib"] = lib; spec.loader.exec_module(lib)
CFG = tomllib.load(open(HERE / "superposition_config.toml", "rb"))
MARGIN = CFG["schedules"]["margin_km"] * 1000.0
DAY = 86400.0
SCHEDS = CFG["schedules"]["names"]
COL = {"passive": "#555555", "super": "#0072B2", "cl": "#D55E00", "final": "#009E73", "band": "#cccccc"}


def old_dir(tag): return tag[6:]


def passive_path(root, tag, inc):
    """Passive truth table with band columns. The old 50-degree grid (border/) has none; border2/ref50 holds the same passes with them."""
    return root / "border2" / "ref50" / old_dir(tag) / "passes_3b.csv" if inc == 50.0 else root / "border2" / "grid" / tag / "passes_3b.csv"


def band(df):
    return df.h_qhigh_m.values + MARGIN, df.h_qlow_m.values - MARGIN


def policy_steps(df, sched):
    """Periapsis steps (m, signed) and signed dV (m/s) of the A1 policy `sched` on the pass table `df`, per pass."""
    x = df.alt_m.values; lo, hi = band(df); L = lib.lever_kepler(df.a_m.values, df.e.values)
    if sched == "lazy_s0": c, n, steps = lib.lazy(x, lo, hi, L, s=0.0)
    elif sched == "near_passive": c, n, steps = lib.near_passive(x, lo, hi, L)
    else: raise ValueError(sched)
    return steps, steps * L, L, c, n


def write_schedule(path, dv):
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [(k + 1, v) for k, v in enumerate(dv) if v != 0.0] or [(1, 0.0)]
    pd.DataFrame(rows, columns=["pass", "dv_mps"]).to_csv(path, index=False)


def read_schedule(path, n):
    d = pd.read_csv(path); v = np.zeros(max(n, int(d["pass"].max())))
    v[d["pass"].values - 1] = d.dv_mps.values
    return v


def cmd_select(root, out):
    L = pd.read_csv(root / "oracle" / "oracle_cell_labels.csv"); c = L[L.status == "completed"]
    rows = []
    for inc in CFG["selection"]["inclinations_deg"]:
        g = c[c.inc_deg == inc].sort_values(["dp_dv", "tag"]).reset_index(drop=True)
        for role, q in (("quiet", CFG["selection"]["quantile_quiet"]), ("high_drift", CFG["selection"]["quantile_high_drift"])):
            t = g.dp_dv.quantile(q); k = (g.dp_dv - t).abs().idxmin(); r = g.loc[k]
            rows.append(dict(tag=r.tag, inc_deg=r.inc_deg, raan_deg=r.raan_deg, omega_deg=r.omega_deg, epoch=int(r.epoch), role=role,
                             quantile_target_mps=t, dp_dv=r.dp_dv, lazy_s0_dv=r.lazy_s0_dv, lazy_s0_burns=r.lazy_s0_burns,
                             near_passive_dv=r.near_passive_dv, near_passive_burns=r.near_passive_burns, a0_km=r.a0_km, e0=r.e0))
    S = pd.DataFrame(rows).sort_values("lazy_s0_dv", kind="stable").reset_index(drop=True)
    S["finite_burn_cell"] = False
    S.loc[(len(S) - 1) // 2, "finite_burn_cell"] = True        # lower median of 8 -> index 3
    out.mkdir(parents=True, exist_ok=True); S.to_csv(out / "sup_cells.csv", index=False); print(S.round(3).to_string())


def cmd_schedule(root, out, it):
    S = pd.read_csv(out / "sup_cells.csv"); rows = []
    for c in S.itertuples():
        if it == 0:
            df = pd.read_csv(passive_path(root, c.tag, c.inc_deg))
        for sched in SCHEDS:
            if it == 0:
                steps, dv, L, cost, nb = policy_steps(df, sched)
                assert dv[0] == 0.0, f"pass-1 burn needed in {c.tag} {sched}: not supported by the time-event kick"
                total = dv
            else:
                prev_name = f"{c.tag}__{sched}__it{it - 1}"
                cl = pd.read_csv(out / "runs" / prev_name / "passes_cl.csv"); prev = read_schedule(out / "schedules" / (prev_name + ".csv"), len(cl))
                steps, inc_dv, L, cost, nb = policy_steps(cl, sched)
                total = prev.copy(); n = min(len(total), len(inc_dv)); total[:n] += inc_dv[:n]
                if len(inc_dv) > len(total): total = np.concatenate([total, inc_dv[len(total):]])
                total[0] = 0.0 if total[0] == 0.0 else total[0]
            write_schedule(out / "schedules" / f"{c.tag}__{sched}__it{it}.csv", total)
            rows.append(dict(tag=c.tag, sched=sched, iter=it, n_burns=int((total != 0).sum()), dv_total_mps=float(np.abs(total).sum()), increment_dv_mps=float(cost) if it else np.nan))
    R = pd.DataFrame(rows); print(R.round(4).to_string()); R.to_csv(out / f"schedule_it{it}_summary.csv", index=False)


# ------------------------------------------------------------------------------------------------ analysis

def excursions(df):
    """Passes outside the shrunk band and outside the full band: (count, max excursion m) each; plus peak Q."""
    x = df.alt_m.values; lo, hi = band(df); flo, fhi = df.h_qhigh_m.values, df.h_qlow_m.values
    e_s = np.maximum(lo - x, x - hi); e_f = np.maximum(flo - x, x - fhi)
    return dict(n_out_shrunk=int((e_s > 1e-6).sum()), max_exc_shrunk_m=float(max(e_s.max(), 0.0)), min_margin_shrunk_m=float(-e_s.max()),
                n_out_full=int((e_f > 1e-6).sum()), max_exc_full_m=float(max(e_f.max(), 0.0)), q_max=float(df.q_w_cm2.max()), n_passes=len(df),
                n_q_above_hard=int((df.q_w_cm2 > 0.54).sum()))


def fit_moon(t_days, y, period):
    """x(t) = c0 + c1 t + c2 t^2 + A cos(w t) + B sin(w t); returns amplitude (m) and phase (deg) of A cos + B sin = R cos(w t - phi)."""
    w = 2 * np.pi / period
    M = np.column_stack([np.ones_like(t_days), t_days, t_days ** 2, np.cos(w * t_days), np.sin(w * t_days)])
    c, *_ = np.linalg.lstsq(M, y, rcond=None)
    return float(np.hypot(c[3], c[4])), float(np.degrees(np.arctan2(c[4], c[3])))


def angdiff(a, b): return (a - b + np.pi) % (2 * np.pi) - np.pi


def load_cell(root, out, c, sched, it):
    name = f"{c.tag}__{sched}__it{it}"
    d = dict(cl=pd.read_csv(out / "runs" / name / "passes_cl.csv"), sched=read_schedule(out / "schedules" / (name + ".csv"), 0))
    kp = out / "runs" / name / "kicks.csv"
    d["kicks"] = pd.read_csv(kp) if kp.exists() else None
    return d


def at_time(p, col, t_s): return np.interp(t_s, p.time_s.values, p[col].values)


def fd_ratio_function(fd, p):
    """Measured / Keplerian lever ratio from the finite-difference probes of one cell, as a function of pass time (constant outside the probes)."""
    f = fd.sort_values("t_apo_days"); f = f[f.ok.astype(str).str.lower() == "true"]
    shift_k = np.array([float(lib.shift_per_dv_km(a, e)) * 1000.0 for a, e in zip(f.a_m, f.e)])
    ratio = f.dalt_per_mps.values / shift_k
    return (lambda t_s: np.interp(t_s / DAY, f.t_apo_days.values, ratio)), f.assign(kepler_m_per_mps=shift_k, ratio=ratio)


def analyze(root, out, finite_sched):
    S = pd.read_csv(out / "sup_cells.csv"); iters = [0, 1, 2]
    FD = pd.concat([pd.read_csv(out / "runs" / f"{t}__null" / "fd_lever.csv") for t in S.tag]); FD.to_csv(out / "fd_lever_all.csv", index=False)
    null = pd.read_csv(out / "null_run_summary.csv"); null.to_csv(out / "null_check.csv", index=False)
    rows_dev, rows_exc, rows_dv, rows_cpl, rows_moon, rows_kick = [], [], [], [], [], []
    store = {}
    period = CFG["analysis"]["moon_oscillation_period_days"]
    for c in S.itertuples():
        p = pd.read_csv(passive_path(root, c.tag, c.inc_deg)); fd = FD[FD.tag == c.tag]
        rfun, fdt = fd_ratio_function(fd, p)
        for sched in SCHEDS:
            steps, dv0, L0, cost0, nb0 = policy_steps(p, sched)         # A1 estimate on the passive history (iteration 0)
            u = np.cumsum(steps)
            hist = {it: load_cell(root, out, c, sched, it) for it in iters}
            cl = hist[0]["cl"]; n = min(len(p), len(cl))
            x = p.alt_m.values; xc = cl.alt_m.values
            t_cl = cl.time_s.values
            # (1) closed loop vs superposition: pass-aligned (A1 convention), time-aligned, and measured-lever
            pred_pass = x[:n] + u[:n]
            xp_t = at_time(p, "alt_m", t_cl)
            u_cl = np.concatenate([u, np.full(max(0, len(cl) - len(u)), u[-1])])[:len(cl)]
            pred_time = xp_t + u_cl
            du = np.diff(np.concatenate([[0.0], u])); r_k = rfun(p.time_s.values)
            u_meas = np.cumsum(du * r_k)
            pred_meas = x[:n] + u_meas[:n]
            dev_pass = xc[:n] - pred_pass; dev_time = xc - pred_time; dev_meas = xc[:n] - pred_meas
            width = float(np.median((p.h_qlow_m.values - p.h_qhigh_m.values) - 2 * MARGIN))
            rows_dev.append(dict(tag=c.tag, role=c.role, inc_deg=c.inc_deg, sched=sched, n_passive=len(p), n_cl=len(cl), t_end_passive_d=p.time_s.iloc[-1] / DAY, t_end_cl_d=cl.time_s.iloc[-1] / DAY,
                                 n_bursts=int((steps != 0).sum()), u_end_m=float(u[-1]), u_range_m=float(u.max() - u.min()),
                                 dev_pass_rms_m=float(np.sqrt((dev_pass ** 2).mean())), dev_pass_max_m=float(np.abs(dev_pass).max()), dev_pass_end_m=float(dev_pass[-1]),
                                 dev_time_rms_m=float(np.sqrt((dev_time ** 2).mean())), dev_time_max_m=float(np.abs(dev_time).max()), dev_time_end_m=float(dev_time[-1]),
                                 dev_meas_rms_m=float(np.sqrt((dev_meas ** 2).mean())), dev_meas_max_m=float(np.abs(dev_meas).max()),
                                 shrunk_width_m=width, dev_time_max_over_width=float(np.abs(dev_time).max() / width)))
            store[(c.tag, sched)] = dict(p=p, hist=hist, u=u, pred_pass=pred_pass, pred_time=pred_time, dev_pass=dev_pass, dev_time=dev_time, steps=steps, dv0=dv0, L0=L0)
            # (2) excursions per iteration (and the passive history)
            for it in iters:
                e = excursions(hist[it]["cl"]); rows_exc.append(dict(tag=c.tag, sched=sched, iter=it, **e))
            rows_exc.append(dict(tag=c.tag, sched=sched, iter=-1, **excursions(p)))        # passive reference
            # (3) dV per iteration
            for it in iters:
                sc_ = hist[it]["sched"]
                rows_dv.append(dict(tag=c.tag, role=c.role, sched=sched, iter=it, dv_total_mps=float(np.abs(sc_).sum()), n_burns=int((sc_ != 0).sum()),
                                    ratio_to_A1=float(np.abs(sc_).sum() / np.abs(dv0).sum())))
            # (4) coupling: element differences, time aligned at the closed-loop end, and drift-rate changes
            k0 = hist[0]["kicks"]; t_end = cl.time_s.iloc[-1]; ts = cl.time_s.values / DAY
            da_kicks = float((k0.a_after_m - k0.a_before_m).sum())
            a_end_cl = float(cl.a_m.iloc[-1]); a_end_p = float(at_time(p, "a_m", t_end))
            def dang(col): return float(np.degrees(angdiff(cl[col].iloc[-1], np.interp(t_end, p.time_s.values, np.unwrap(p[col].values)))))
            dres = dev_time                                           # periapsis residual wrt superposition (time aligned)
            slope_res = float(np.polyfit(ts, dres, 1)[0]); slope_p = float(np.polyfit(p.time_s.values / DAY, x, 1)[0])
            rate_p = np.diff(x) / np.diff(p.time_s.values / DAY)
            om_p = np.degrees(np.unwrap(p.omega_rad.values)); om_rate_p = (om_p[-1] - om_p[0]) / (p.time_s.iloc[-1] / DAY)
            rows_cpl.append(dict(tag=c.tag, role=c.role, inc_deg=c.inc_deg, sched=sched, da_kicks_m=da_kicks, da_end_actual_m=a_end_cl - a_end_p,
                                 da_from_dynamics_m=(a_end_cl - a_end_p) - da_kicks, a_passive_decay_m=float(p.a_m.iloc[-1] - p.a_m.iloc[0]),
                                 a_cl_decay_m=float(cl.a_m.iloc[-1] - cl.a_m.iloc[0]), da_dyn_over_passive_decay=((a_end_cl - a_end_p) - da_kicks) / (p.a_m.iloc[-1] - p.a_m.iloc[0]),
                                 de_end=float(cl.e.iloc[-1] - np.interp(t_end, p.time_s.values, p.e.values)),
                                 domega_end_deg=dang("omega_rad"), draan_end_deg=dang("raan_rad"), di_end_deg=float(np.degrees(cl.i_rad.iloc[-1] - at_time(p, "i_rad", t_end))),
                                 omega_passive_total_deg=float(om_p[-1] - om_p[0]), domega_over_passive_total=dang("omega_rad") / float(om_p[-1] - om_p[0]) if abs(om_p[-1] - om_p[0]) > 1e-9 else np.nan,
                                 raan_passive_total_deg=float(np.degrees(np.unwrap(p.raan_rad.values)[-1] - np.unwrap(p.raan_rad.values)[0])),
                                 draan_over_passive_total=dang("raan_rad") / float(np.degrees(np.unwrap(p.raan_rad.values)[-1] - np.unwrap(p.raan_rad.values)[0])),
                                 omega_rate_passive_deg_d=om_rate_p, domega_rate_deg_d=dang("omega_rad") / (t_end / DAY),
                                 resid_slope_m_per_day=slope_res, passive_linear_slope_m_per_day=slope_p, passive_rms_rate_m_per_day=float(np.sqrt((rate_p ** 2).mean())),
                                 resid_slope_over_passive_rms_rate=slope_res / float(np.sqrt((rate_p ** 2).mean()))))
            # effective lever per kick (closed loop, time aligned): realized shift of the pass-to-pass residual against the commanded shift
            d_ = xc - xp_t
            dd = np.diff(np.concatenate([[0.0], d_])); du_cl = np.diff(np.concatenate([[0.0], u_cl]))
            m = np.abs(du_cl) >= CFG["analysis"]["min_kick_shift_m"]
            rat = dd[m] / du_cl[m]
            rows_kick.append(dict(tag=c.tag, sched=sched, n_kicks_used=int(m.sum()), eff_lever_ratio_median=float(np.median(rat)) if m.any() else np.nan,
                                  eff_lever_ratio_q1=float(np.quantile(rat, .25)) if m.any() else np.nan, eff_lever_ratio_q3=float(np.quantile(rat, .75)) if m.any() else np.nan,
                                  eff_lever_ratio_weighted=float(dd[m].sum() / du_cl[m].sum()) if m.any() else np.nan,
                                  fd_ratio_mean=float(fdt.ratio.mean()), fd_ratio_min=float(fdt.ratio.min()), fd_ratio_max=float(fdt.ratio.max())))
            # (5) Moon-oscillation phase of the periapsis, passive vs closed loop with the superposed offset removed
            A_p, ph_p = fit_moon(p.time_s.values / DAY, x, period); A_c, ph_c = fit_moon(ts, xc - u_cl, period)
            rows_moon.append(dict(tag=c.tag, role=c.role, epoch=c.epoch, sched=sched, amp_passive_m=A_p, phase_passive_deg=ph_p, amp_cl_minus_u_m=A_c, phase_cl_minus_u_deg=ph_c,
                                  dphase_deg=float((ph_c - ph_p + 180) % 360 - 180), amp_ratio=A_c / A_p if A_p > 0 else np.nan))
    D = pd.DataFrame(rows_dev); X = pd.DataFrame(rows_exc); V = pd.DataFrame(rows_dv); C = pd.DataFrame(rows_cpl); Mo = pd.DataFrame(rows_moon); K = pd.DataFrame(rows_kick)
    for nm, t in (("deviation", D), ("excursions", X), ("dv_iterations", V), ("coupling", C), ("moon_phase", Mo), ("effective_lever", K)):
        t.to_csv(out / f"sup_{nm}.csv", index=False); print("\n==", nm); print(t.round(4).to_string())
    # dV ratio table: final vs A1 estimate
    R = V.pivot_table(index=["tag", "role", "sched"], columns="iter", values=["dv_total_mps", "n_burns"]).reset_index()
    R.columns = ["_".join(str(x) for x in col if x != "") for col in R.columns]
    R["ratio_final_to_A1"] = R["dv_total_mps_2"] / R["dv_total_mps_0"]; R["ratio_it1_to_A1"] = R["dv_total_mps_1"] / R["dv_total_mps_0"]
    R.to_csv(out / "sup_dv_ratio.csv", index=False); print("\n== dv ratio"); print(R.round(4).to_string())
    for sched in SCHEDS:
        g = R[R.sched == sched]
        print(sched, "final/A1: median", g.ratio_final_to_A1.median(), "min", g.ratio_final_to_A1.min(), "max", g.ratio_final_to_A1.max(),
              "| it1/A1 median", g.ratio_it1_to_A1.median())
    # policy ranking: is lazy cheaper than near_passive after repair? and A1 ranking
    rk = R.pivot_table(index="tag", columns="sched", values=["dv_total_mps_0", "dv_total_mps_2"]); rk["rank_A1_lazy_cheaper"] = rk[("dv_total_mps_0", "lazy_s0")] < rk[("dv_total_mps_0", "near_passive")]
    rk["rank_final_lazy_cheaper"] = rk[("dv_total_mps_2", "lazy_s0")] < rk[("dv_total_mps_2", "near_passive")]; print(rk.round(3).to_string()); rk.to_csv(out / "sup_policy_ranking.csv")
    # cell ordering by dV: Spearman between A1 estimate and final
    from scipy import stats as st
    for sched in SCHEDS:
        g = R[R.sched == sched]; print(sched, "Spearman(A1, final) over 8 cells:", st.spearmanr(g.dv_total_mps_0, g.dv_total_mps_2)[0])
    figures_sup(out, S, store, D, V, R, C)
    finite_compare(out, S, store, finite_sched)


def figures_sup(out, S, store, D, V, R, C):
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    for part, cells in enumerate((S.iloc[:4], S.iloc[4:])):
        fig, ax = plt.subplots(2, 4, figsize=(17, 6.6), sharex="col", gridspec_kw=dict(height_ratios=[2.2, 1]))
        for j, c in enumerate(cells.itertuples()):
            d = store[(c.tag, "lazy_s0")]; p = d["p"]; cl0 = d["hist"][0]["cl"]; cl2 = d["hist"][2]["cl"]
            t = p.time_s.values / DAY; lo, hi = band(p)
            a = ax[0, j]; a.fill_between(t, lo / 1e3, hi / 1e3, color=COL["band"], alpha=.6, label="shrunk band")
            a.plot(t, p.alt_m / 1e3, color=COL["passive"], lw=1, label="passive")
            a.plot(t[:len(d["pred_pass"])], d["pred_pass"] / 1e3, color=COL["super"], lw=1.6, label="superposition (passive + u)")
            a.plot(cl0.time_s / DAY, cl0.alt_m / 1e3, color=COL["cl"], lw=1.2, ls="--", label="closed loop, A1 schedule")
            a.plot(cl2.time_s / DAY, cl2.alt_m / 1e3, color=COL["final"], lw=1.2, label="closed loop, repaired (iter 2)")
            a.set_title(f"{c.tag} ({c.role}, {c.lazy_s0_dv:.2f} m/s)", fontsize=9); a.set_ylabel("periapsis altitude (km)" if j == 0 else "")
            b = ax[1, j]; b.plot(t[:len(d["dev_pass"])], d["dev_pass"], color=COL["cl"], lw=1.2, label="pass-aligned")
            b.plot(cl0.time_s / DAY, d["dev_time"], color=COL["super"], lw=1.2, ls=":", label="time-aligned")
            b.axhline(0, color="k", lw=.5); b.set_xlabel("time (d)"); b.set_ylabel("closed loop - superposition (m)" if j == 0 else "")
            if j == 0: a.legend(fontsize=7, frameon=False, loc="best"); b.legend(fontsize=7, frameon=False)
        fig.tight_layout(); fig.savefig(out / f"fig_sup_overlay_{part + 1}.png", dpi=125); plt.close(fig)
    # dV ratio per cell
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2), sharey=True)
    for a, sched in zip(ax, SCHEDS):
        g = R[R.sched == sched].sort_values("dv_total_mps_0").reset_index(drop=True); x = np.arange(len(g))
        a.bar(x - .2, g.dv_total_mps_1 / g.dv_total_mps_0, .38, color=COL["cl"], label="after 1 repair / A1")
        a.bar(x + .2, g.ratio_final_to_A1, .38, color=COL["final"], label="after 2 repairs / A1")
        a.axhline(1, color="k", lw=.7); a.set_xticks(x); a.set_xticklabels([t[1:] for t in g.tag], rotation=70, ha="right", fontsize=8)
        a.set_title(f"{sched}: repaired dV over the superposition estimate"); a.set_ylabel("ratio" if sched == SCHEDS[0] else "")
        a.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(out / "fig_sup_dv_ratio.png", dpi=125); plt.close(fig)
    # coupling: a decay difference and excursions
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    g = C[C.sched == "lazy_s0"].sort_values("tag"); x = np.arange(len(g))
    ax[0].bar(x - .2, g.da_kicks_m / 1e3, .38, color=COL["super"], label="from the kicks (impulse)"); ax[0].bar(x + .2, g.da_end_actual_m / 1e3, .38, color=COL["cl"], label="actual (closed loop - passive)")
    ax[0].set_xticks(x); ax[0].set_xticklabels([t[1:] for t in g.tag], rotation=70, ha="right", fontsize=7); ax[0].set_ylabel("semi-major axis change at end (km)"); ax[0].legend(frameon=False, fontsize=8); ax[0].set_title("lazy s=0: a changed by the steps")
    dd = D[D.sched == "lazy_s0"].sort_values("tag")
    ax[1].bar(x - .2, dd.dev_pass_max_m / 1e3, .38, color=COL["cl"], label="pass-aligned"); ax[1].bar(x + .2, dd.dev_time_max_m / 1e3, .38, color=COL["super"], label="time-aligned")
    ax[1].axhline(1.0, color="k", lw=.7, ls="--"); ax[1].text(len(x) - .5, 1.05, "1 km margin", ha="right", fontsize=8)
    ax[1].set_xticks(x); ax[1].set_xticklabels([t[1:] for t in dd.tag], rotation=70, ha="right", fontsize=7); ax[1].set_ylabel("max |closed loop - superposition| (km)"); ax[1].legend(frameon=False, fontsize=8); ax[1].set_title("lazy s=0: deviation from superposition")
    dn = D[D.sched == "near_passive"].sort_values("tag")
    ax[2].bar(x - .2, dn.dev_pass_max_m / 1e3, .38, color=COL["cl"], label="pass-aligned"); ax[2].bar(x + .2, dn.dev_time_max_m / 1e3, .38, color=COL["super"], label="time-aligned")
    ax[2].axhline(1.0, color="k", lw=.7, ls="--"); ax[2].set_xticks(x); ax[2].set_xticklabels([t[1:] for t in dn.tag], rotation=70, ha="right", fontsize=7); ax[2].set_title("near_passive: deviation from superposition")
    fig.tight_layout(); fig.savefig(out / "fig_sup_coupling.png", dpi=125); plt.close(fig)


def finite_compare(out, S, store, sched):
    rows = []; plots = []
    for c in S.itertuples():
        fdir = out / "runs" / f"{c.tag}__{sched}__it0__finite"
        if not fdir.exists(): continue
        fb = pd.read_csv(fdir / "passes_cl.csv"); plans = pd.read_csv(fdir / "burn_plans.csv"); ev = pd.read_csv(fdir / "apsis_events.csv")
        imp = store[(c.tag, sched)]["hist"][0]["cl"]; kicks = store[(c.tag, sched)]["hist"][0]["kicks"]; p = store[(c.tag, sched)]["p"]
        n = min(len(fb), len(imp)); dxx = fb.alt_m.values[:n] - imp.alt_m.values[:n]
        apo = ev[ev.kind == "apoapsis"].time_s.values
        mid = 0.5 * (plans.start_s + plans.stop_s).values; dur = (plans.stop_s - plans.start_s).values
        off = np.array([m - apo[np.argmin(np.abs(apo - m))] for m in mid])
        # kick time of the impulsive run against the finite-burn mid time of the same firing
        km = kicks[kicks.dv_mps != 0].set_index("firing"); pl = plans.set_index("firing"); j = km.index.intersection(pl.index)
        j = j[:5]   # the timelines drift apart (period difference), so only the first five firings are compared by time
        dt = (0.5 * (pl.loc[j, "start_s"] + pl.loc[j, "stop_s"]) - km.loc[j, "t_s"]).values
        # realized shift per dV: pass-to-pass residual against passive (time aligned), finite vs impulsive
        def shift_per_dv(cl, u_sched):
            xp_t = at_time(p, "alt_m", cl.time_s.values); d_ = cl.alt_m.values - xp_t; dd = np.diff(np.concatenate([[0.0], d_]))
            return dd
        u_ = store[(c.tag, sched)]["u"]; u_ = np.concatenate([u_, np.full(max(0, n - len(u_)), u_[-1])])[:n]
        du = np.diff(np.concatenate([[0.0], u_]))
        ddf = shift_per_dv(fb, None)[:n]; ddi = shift_per_dv(imp, None)[:n]; m = np.abs(du) >= CFG["analysis"]["min_kick_shift_m"]
        rows.append(dict(tag=c.tag, n_burns=len(plans), dv_total_mps=float(plans.dv_cmd_mps.abs().sum()), n_passes_fb=len(fb), n_passes_imp=len(imp),
                         alt_diff_rms_m=float(np.sqrt((dxx ** 2).mean())), alt_diff_max_m=float(np.abs(dxx).max()), alt_diff_end_m=float(dxx[-1]),
                         burn_duration_median_s=float(np.median(dur)), burn_duration_max_s=float(dur.max()), burn_mid_minus_apoapsis_s_median=float(np.median(off)), burn_mid_minus_apoapsis_s_max_abs=float(np.abs(off).max()),
                         burn_mid_minus_impulse_time_first5_s_median=float(np.median(dt)) if len(dt) else np.nan,
                         shift_ratio_fb_over_imp_median=float(np.median(ddf[m] / ddi[m])) if m.any() else np.nan, shift_ratio_fb_over_imp_weighted=float(ddf[m].sum() / ddi[m].sum()) if m.any() else np.nan,
                         propellant_used_kg=float(plans.propellant_kg.sum()), a_end_diff_m=float(fb.a_m.iloc[-1] - imp.a_m.iloc[-1]), t_end_diff_s=float(fb.time_s.iloc[-1] - imp.time_s.iloc[-1])))
        plots.append((c, fb, imp, p, dxx))
    F = pd.DataFrame(rows); F.to_csv(out / "sup_finite_vs_impulsive.csv", index=False); print("\n== finite vs impulsive"); print(F.round(4).to_string())
    if not plots: return
    fig, ax = plt.subplots(2, len(plots), figsize=(6.2 * len(plots), 5.4), squeeze=False)
    for j, (c, fb, imp, p, dxx) in enumerate(plots):
        a = ax[0, j]; a.plot(p.time_s / DAY, p.alt_m / 1e3, color=COL["passive"], lw=1, label="passive"); a.plot(imp.time_s / DAY, imp.alt_m / 1e3, color=COL["cl"], lw=1.4, label="impulsive replay")
        a.plot(fb.time_s / DAY, fb.alt_m / 1e3, color=COL["final"], lw=1, ls="--", label="finite burn (4 N, Isp 300 s)"); a.set_ylabel("periapsis altitude (km)"); a.set_title(f"{c.tag}, {sched}, A1 schedule"); a.legend(frameon=False, fontsize=8)
        b = ax[1, j]; b.plot(fb.time_s.values[:len(dxx)] / DAY, dxx, color=COL["final"]); b.axhline(0, color="k", lw=.5); b.set_xlabel("time (d)"); b.set_ylabel("finite - impulsive (m)")
    fig.tight_layout(); fig.savefig(out / "fig_sup_finite_vs_impulsive.png", dpi=125); plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("cmd", choices=["select", "schedule", "analyze"]); add_results_root_argument(p)
    p.add_argument("--iter", type=int, default=0); p.add_argument("--finite-sched", default="lazy_s0")
    a = p.parse_args(); root = pathlib.Path(results_root(a.results_root)); out = root / "superposition"
    if a.cmd == "select": cmd_select(root, out)
    elif a.cmd == "schedule": cmd_schedule(root, out, a.iter)
    else: analyze(root, out, a.finite_sched)
