"""Static-region sweep analysis. Usage: python3 analyze_static.py OUT_DIR STEP   (STEP = metrics | race | plots)
metrics: per-cell metrics for the 288 new cells + the 50-deg reference, selection files (control_cells.csv, race_cells.csv).
race:    180-day tables, authority-down times, heat-limited altitude ramp, FD lever.
plots:   all figures (also writes the correction-cost table)."""
import sys, os, glob, tomllib
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1]; STEP = sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
S = tomllib.load(open(f"{HERE}/study_config.toml", "rb")); C = tomllib.load(open(f"{HERE}/static_config.toml", "rb"))
QLO, QHI, QHARD = S["band"]["q_low_w_cm2"], S["band"]["q_high_w_cm2"], S["band"]["q_hard_w_cm2"]
X = C["static_window"]["x_km"] * 1000.0
DAY = 86400.0; MU = 3.986004418e14
FAMILY = {23.4: "F2 (23.4)", 63.4: "F1 (63.4)", 116.6: "F1 (116.6)", 50.0: "ref (50)"}

def cell_dirs():
    for d in sorted(glob.glob(f"{OUT}/grid/I*")):
        if os.path.exists(f"{d}/passes_3b.csv"): yield os.path.basename(d), d
    for d in sorted(glob.glob(f"{OUT}/ref50/R*")):
        yield "I0500_" + os.path.basename(d), d

def unwrap_deg(x): return np.rad2deg(np.unwrap(np.asarray(x)))

def static_window(df):
    """Longest run of consecutive passes with Q in band and |alt - alt0| <= X. Returns (length_d, start_d, end_d)."""
    q = df.q_w_cm2.values; ok = (q >= QLO) & (q <= QHI) & (np.abs(df.alt_m.values - df.alt_m.values[0]) <= X)
    t = df.time_s.values / DAY; best = (0.0, np.nan, np.nan); j = 0
    while j < len(ok):
        if ok[j]:
            k = j
            while k + 1 < len(ok) and ok[k + 1]: k += 1
            if t[k] - t[j] > best[0] or np.isnan(best[1]): best = (t[k] - t[j], t[j], t[k])
            j = k + 1
        else: j += 1
    return best

def metrics(tag, d, inc):
    df = pd.read_csv(f"{d}/passes_3b.csv"); t = df.time_s.values / DAY; alt = df.alt_m.values
    w = static_window(df); q = df.q_w_cm2.values
    inb = (q >= QLO) & (q <= QHI); dt = np.append(np.diff(t), 0.0)
    fe = np.argmax(~inb) if (~inb).any() else -1
    hard = np.argmax(q > QHARD) if (q > QHARD).any() else -1
    rad = df.radius_m.values
    kozai = np.sum(((df.moon_kozai_drp_mps + df.sun_kozai_drp_mps).values * np.append(np.diff(df.time_s.values), 0.0)))
    om = unwrap_deg(df.omega_rad); raan = unwrap_deg(df.raan_rad)
    return dict(tag=tag, inc_deg=inc, raan_deg=float(tag.split("_R")[1][:3]), omega_deg=float(tag.split("_W")[1][:3]), epoch_idx=int(tag[-1]),
        family=FAMILY[inc], t_end_days=t[-1], passes=len(df), alt0_m=alt[0], alt_end_m=alt[-1], alt_min_m=alt.min(), alt_max_m=alt.max(),
        drp_total_m=alt[-1] - alt[0], max_dev_m=np.abs(alt - alt[0]).max(), pp_m=alt.max() - alt.min(),
        apo0_m=df.apo_alt_m.values[0], apo_end_m=df.apo_alt_m.values[-1], a0_m=df.a_m.values[0], a_end_m=df.a_m.values[-1],
        e0=df.e.values[0], e_end=df.e.values[-1], q_max=q.max(), q_min=q.min(), in_band_passes=int(inb.sum()), in_band_days=float(dt[inb].sum()),
        first_exit_days=t[fe] if fe >= 0 else np.nan, first_exit_edge=("" if fe < 0 else ("high_Q" if q[fe] > QHI else "low_Q")),
        exceeds_hard=bool(hard >= 0 or q.max() > QHARD), hard_first_days=t[hard] if hard >= 0 else np.nan,
        window_days=w[0], window_start_days=w[1], window_end_days=w[2],
        omega_p_moon_start=np.rad2deg(df.moon_omega_p_rad.values[0]), omega_p_sun_start=np.rad2deg(df.sun_omega_p_rad.values[0]),
        omega_p_moon_change=unwrap_deg(df.moon_omega_p_rad)[-1] - unwrap_deg(df.moon_omega_p_rad)[0],
        omega_p_sun_change=unwrap_deg(df.sun_omega_p_rad)[-1] - unwrap_deg(df.sun_omega_p_rad)[0],
        i_p_moon_start=np.rad2deg(df.moon_i_p_rad.values[0]), i_p_sun_start=np.rad2deg(df.sun_i_p_rad.values[0]),
        i_p_moon_min=np.rad2deg(df.moon_i_p_rad.min()), i_p_moon_max=np.rad2deg(df.moon_i_p_rad.max()),
        i_p_sun_min=np.rad2deg(df.sun_i_p_rad.min()), i_p_sun_max=np.rad2deg(df.sun_i_p_rad.max()),
        raan_drift_deg=raan[-1] - raan[0], omega_drift_deg=om[-1] - om[0], i_drift_deg=np.rad2deg(df.i_rad.values[-1] - df.i_rad.values[0]),
        drr_total_m=rad[-1] - rad[0], kozai_pred_drr_m=kozai, lat_amp_deg=np.rad2deg((df.lat_rad.max() - df.lat_rad.min()) / 2),
        sin2_ip_moon=np.sin(df.moon_i_p_rad.values[0]) ** 2, sin2_ip_sun=np.sin(df.sun_i_p_rad.values[0]) ** 2)

def rank_key(df, best=True):
    # best: longest window, then smallest max deviation; worst: shortest window, then largest max deviation.
    return df.sort_values(["window_days", "max_dev_m"], ascending=[False, True] if best else [True, False])

def do_metrics():
    rows = []
    for tag, d in cell_dirs():
        inc = float(int(tag[1:5])) / 10
        rows.append(metrics(tag, d, inc))
    df = pd.DataFrame(rows)
    run = pd.read_csv(f"{OUT}/static_run_summary.csv")[["tag", "retcode_3b", "cal_iterations", "wall_3b_s", "cal_wall_s"]]
    ref = pd.read_csv(os.path.join(os.environ.get("REF_DIR", C["grid"]["reference_dir"]), "grid_summary.csv"))
    ref["tag"] = "I0500_" + ref.tag; ref = ref[["tag", "retcode_3b"]]
    df = df.merge(pd.concat([run, ref]), on="tag", how="left")
    df["ended_early"] = df.retcode_3b.eq("Terminated")
    df.to_csv(f"{OUT}/static_cells.csv", index=False)
    new = df[df.inc_deg != 50.0]
    best = rank_key(new, True).head(C["control"]["n_best"]); worst = rank_key(new, False).head(C["control"]["n_worst"])
    # F2 null geometries (orbit plane near the ecliptic): every inclination-23.4 cell with RAAN 0 or 180 deg joins the control subset.
    f2n = new[(new.inc_deg == 23.4) & new.raan_deg.isin([0.0, 180.0])].assign(subset="f2_null")
    ctl = pd.concat([best.assign(subset="best"), worst.assign(subset="worst"), f2n]).drop_duplicates("tag")
    ctl[["tag", "inc_deg", "raan_deg", "omega_deg", "epoch_idx", "subset", "window_days"]].to_csv(f"{OUT}/control_cells.csv", index=False)
    fam = new.assign(f1=new.inc_deg.isin([63.4, 116.6]))
    rk = rank_key(new, True)
    pick = []
    for isf1 in (True, False):
        pick += list(rk[new.loc[rk.index, "inc_deg"].isin([63.4, 116.6]) == isf1].head(C["race"]["min_best_per_family"]).tag)
    pick += list(rk[~rk.tag.isin(pick)].head(C["race"]["n_best"] - len(pick)).tag)
    rw = list(rank_key(new, False).head(C["race"]["n_worst"]).tag)
    rc = new.set_index("tag").loc[pick + rw].reset_index()
    rc["role"] = ["best"] * len(pick) + ["worst"] * len(rw)
    rc[["tag", "inc_deg", "raan_deg", "omega_deg", "epoch_idx", "role", "window_days"]].to_csv(f"{OUT}/race_cells.csv", index=False)
    # distributions
    g = df.groupby("family")
    pd.set_option("display.width", 250)
    print(g.window_days.describe().round(2).to_string())
    print(df.groupby(["family", "omega_deg"]).window_days.agg(["mean", "median", "max"]).round(1).to_string())
    print("early-end:\n", g.ended_early.agg(["sum", "count"]).to_string(), "\nexceeds hard:\n", g.exceeds_hard.agg(["sum", "count"]).to_string())
    print("pp/max_dev/raan drift medians:\n", g[["max_dev_m", "pp_m", "raan_drift_deg", "omega_drift_deg", "lat_amp_deg", "drp_total_m"]].median().round(1).to_string())
    print("\nBEST\n", best.head(10)[["tag", "window_days", "max_dev_m", "drp_total_m"]].to_string(), "\nWORST\n", worst[["tag", "window_days", "max_dev_m", "drp_total_m"]].to_string())

# ---- correction cost ---------------------------------------------------------------------------------
def lever_mps_per_km(a, e):
    ra = a * (1 + e); va = np.sqrt(MU * (2 / ra - 1 / a))
    return 1000.0 * MU / (4 * a * a * va)       # tangential kick at apoapsis: d r_p = 4 a^2 v_a dv / mu

def correction_cost(df, horizon_days, margin_m):
    d = df[df.time_s / DAY <= horizon_days]
    lev = lever_mps_per_km(d.a_m.values, d.e.values); off = 0.0; dv = 0.0; n = 0
    for j in range(len(d)):
        lo = d.h_qhigh_m.values[j] + margin_m; hi = d.h_qlow_m.values[j] - margin_m
        a = d.alt_m.values[j] + off
        if a > hi or a < lo:
            tgt = lo if a > hi else hi
            dv += abs(tgt - a) / 1000.0 * lev[j]; off += tgt - a; n += 1
    return dv, n, float(d.time_s.values[-1] / DAY)

def do_cost(cells, race):
    m = C["cost"]["deadband_margin_km"] * 1000; rows = []
    for _, r in cells.iterrows():
        d = f"{OUT}/grid/{r.tag}" if r.inc_deg != 50.0 else f"{OUT}/ref50/{r.tag[6:]}"
        df = pd.read_csv(f"{d}/passes_3b.csv")
        dv, n, cov = correction_cost(df, 60.0, m)
        rows.append(dict(tag=r.tag, inc_deg=r.inc_deg, family=r.family, horizon_days=60.0, covered_days=cov, dv_mps=dv, n_corrections=n, ended_early=r.ended_early))
    for _, r in race.iterrows():
        df = pd.read_csv(f"{OUT}/race/{r.tag}/passes_race.csv")
        for hd in (60.0, 180.0):
            dv, n, cov = correction_cost(df, hd, m)
            rows.append(dict(tag=r.tag, inc_deg=r.inc_deg, family=FAMILY[r.inc_deg], horizon_days=hd + 0.001, covered_days=cov, dv_mps=dv, n_corrections=n, ended_early=cov < hd - 2, race=True))
    return pd.DataFrame(rows)

def band_shade(ax, df, td):
    ax.fill_between(td, df.h_qhigh_m / 1000, df.h_qlow_m / 1000, color="tab:green", alpha=0.18, lw=0, label="band 0.14-0.32 W/cm$^2$")

def race_load(race):
    return {r.tag: pd.read_csv(f"{OUT}/race/{r.tag}/passes_race.csv") for _, r in race.iterrows()}

def roll_rate(df, win_d):
    t = df.time_s.values / DAY; r = df.radius_m.values / 1000; out = np.full(len(t), np.nan)
    for j in range(len(t)):
        m = np.abs(t - t[j]) <= win_d / 2
        if m.sum() >= 4: out[j] = np.polyfit(t[m], r[m], 1)[0]
    return out

def do_race():
    race = pd.read_csv(f"{OUT}/race_cells.csv"); cells = pd.read_csv(f"{OUT}/static_cells.csv").set_index("tag")
    summ = pd.read_csv(f"{OUT}/race_run_summary.csv").set_index("tag"); thr = C["race"]["authority_rate_threshold_km_day"]; win = C["race"]["rate_window_days"]
    rows = []
    for _, r in race.iterrows():
        df = pd.read_csv(f"{OUT}/race/{r.tag}/passes_race.csv"); t = df.time_s.values / DAY
        rate = roll_rate(df, win); ok = ~np.isnan(rate)
        down = np.nan
        for j in range(len(t) - 1, -1, -1):
            if ok[j] and abs(rate[j]) >= thr: down = t[min(j + 1, len(t) - 1)] if j < len(t) - 1 else np.nan; break
        else: down = t[np.argmax(ok)] if ok.any() else np.nan
        hc = (df.h_qhigh_m + df.h_qlow_m).values / 2
        dev = df.alt_m.values - df.alt_m.values[0]; inb = (df.q_w_cm2 >= QLO) & (df.q_w_cm2 <= QHI)
        # first day the periapsis leaves the band (and stays out for good = last in-band pass)
        last_in = t[np.where(inb)[0][-1]] if inb.any() else np.nan; first_out = t[np.argmax(~inb)] if (~inb).any() else np.nan
        w = static_window(df); s = summ.loc[r.tag]
        rows.append(dict(tag=r.tag, role=r.role, inc_deg=r.inc_deg, retcode=s.retcode_race, t_end_days=t[-1], window60_days=r.window_days,
            window_race_days=w[0], first_out_of_band_days=first_out, last_in_band_days=last_in,
            alt0_km=df.alt_m.values[0] / 1000, alt_end_km=df.alt_m.values[-1] / 1000, alt_pp_km=(df.alt_m.max() - df.alt_m.min()) / 1000,
            a0_km=df.a_m.values[0] / 1000, a_end_km=df.a_m.values[-1] / 1000, e0=df.e.values[0], e_end=df.e.values[-1],
            apo0_km=df.apo_alt_m.values[0] / 1000, apo_end_km=df.apo_alt_m.values[-1] / 1000, q_end=df.q_w_cm2.values[-1],
            authority_down_days=down, max_abs_rate_km_day=np.nanmax(np.abs(rate)), rate_last30_km_day=np.nanmean(np.abs(rate[t > t[-1] - 30])),
            lat_amp_deg=np.rad2deg((df.lat_rad.max() - df.lat_rad.min()) / 2), lat_amp_first30_deg=np.rad2deg((df.lat_rad[t < 30].max() - df.lat_rad[t < 30].min()) / 2),
            band_center_move_km=(hc[-1] - hc[0]) / 1000, h_qhigh_move_km=(df.h_qhigh_m.values[-1] - df.h_qhigh_m.values[0]) / 1000,
            h_qlow_move_km=(df.h_qlow_m.values[-1] - df.h_qlow_m.values[0]) / 1000, vrel0=df.vrel_mps.values[0], vrel_end=df.vrel_mps.values[-1],
            i_p_moon_end=np.rad2deg(df.moon_i_p_rad.values[-1]), i_p_sun_end=np.rad2deg(df.sun_i_p_rad.values[-1]),
            omega_p_moon_change=unwrap_deg(df.moon_omega_p_rad)[-1] - unwrap_deg(df.moon_omega_p_rad)[0],
            omega_drift_deg=unwrap_deg(df.omega_rad)[-1] - unwrap_deg(df.omega_rad)[0], raan_drift_deg=unwrap_deg(df.raan_rad)[-1] - unwrap_deg(df.raan_rad)[0],
            lever_ok=s.get("lever_ok"), lever_a_km=s.get("lever_a_m", np.nan) / 1000, lever_mps_per_km=s.get("lever_mps_per_km", np.nan),
            lever_analytic=float(lever_mps_per_km(s.get("lever_a_m", np.nan), s.get("lever_e", np.nan))) if s.get("lever_ok") else np.nan, lever_t_apo_days=s.get("lever_t_apo_days", np.nan)))
    out = pd.DataFrame(rows); out.to_csv(f"{OUT}/race_summary.csv", index=False)
    pd.set_option("display.width", 250); print(out.round(3).T.to_string())

def do_plots():
    cells = pd.read_csv(f"{OUT}/static_cells.csv"); race = pd.read_csv(f"{OUT}/race_cells.csv")
    cost = do_cost(cells, race); cost.to_csv(f"{OUT}/correction_cost.csv", index=False)
    # 1. heat maps
    incs = [23.4, 50.0, 63.4, 116.6]; vmax = 60
    fig, ax = plt.subplots(len(incs), 3, figsize=(10, 12), constrained_layout=True)
    for i, inc in enumerate(incs):
        for e in (1, 2, 3):
            s = cells[(cells.inc_deg == inc) & (cells.epoch_idx == e)]
            rr = sorted(s.raan_deg.unique()); ww = sorted(s.omega_deg.unique())
            M = np.full((len(rr), len(ww)), np.nan)
            for _, r in s.iterrows(): M[rr.index(r.raan_deg), ww.index(r.omega_deg)] = r.window_days
            a = ax[i, e - 1]; im = a.imshow(M, origin="lower", vmin=0, vmax=vmax, cmap="viridis", aspect="auto")
            a.set_xticks(range(len(ww)), [f"{int(x)}" for x in ww]); a.set_yticks(range(len(rr)), [f"{int(x)}" for x in rr])
            for (y, x), v in np.ndenumerate(M):
                if not np.isnan(v): a.text(x, y, f"{v:.0f}", ha="center", va="center", color="w" if v < 40 else "k", fontsize=7)
            a.set_title(f"i = {inc} deg, epoch {e}", fontsize=9); a.set_xlabel("omega (deg)", fontsize=8); a.set_ylabel("RAAN (deg)", fontsize=8)
    fig.colorbar(im, ax=ax, label="static window (days)", shrink=0.6)
    fig.savefig(f"{OUT}/fig_static_window_heatmaps.png", dpi=130); plt.close(fig)
    # 2. best/worst histories (grid cells, 60 d): the 6 best + 2 worst race picks
    fig, axs = plt.subplots(2, 4, figsize=(15, 7), constrained_layout=True, sharex=True)
    for a, (_, r) in zip(axs.ravel(), race.iterrows()):
        df = pd.read_csv(f"{OUT}/grid/{r.tag}/passes_3b.csv"); td = df.time_s / DAY
        band_shade(a, df, td); a.plot(td, df.alt_m / 1000, "k.-", ms=3, lw=0.8, label="periapsis altitude")
        a.axhspan(df.alt_m.iloc[0] / 1000 - X / 1000, df.alt_m.iloc[0] / 1000 + X / 1000, fc="none", ec="tab:blue", ls="--", lw=0.8)
        a.set_title(f"{r.role}: {r.tag}\nwindow {r.window_days:.1f} d", fontsize=9); a.set_xlabel("day"); a.set_ylabel("periapsis geodetic altitude (km)")
    axs[0, 0].legend(fontsize=7, loc="best"); fig.savefig(f"{OUT}/fig_best_worst_histories.png", dpi=130); plt.close(fig)
    # 3. race plots
    rs = pd.read_csv(f"{OUT}/race_summary.csv"); D = race_load(race)
    for name, panels in (("a", [("alt", "periapsis alt (km)"), ("q", "peak Q (W/cm$^2$)"), ("rate", "d r_p/dt, 15-d fit (km/day)"), ("lat", "periapsis latitude (deg)")]),
                         ("b", [("a", "a (km)"), ("e", "e"), ("om", "omega' (deg; Moon solid, Sun dashed)"), ("ip", "i' (deg; Moon solid, Sun dashed)")])):
        fig, axs = plt.subplots(4, 8, figsize=(26, 12), constrained_layout=True, sharex=True)
        for c, (_, r) in enumerate(race.iterrows()):
            df = D[r.tag]; t = df.time_s / DAY
            for k, (key, lab) in enumerate(panels):
                a = axs[k, c]
                if key == "alt": band_shade(a, df, t); a.plot(t, df.alt_m / 1000, "k-", lw=0.9)
                elif key == "q": a.semilogy(t, df.q_w_cm2, "k-", lw=0.9); a.axhspan(QLO, QHI, color="tab:green", alpha=0.18); a.axhline(QHARD, color="r", lw=0.7, ls=":"); a.set_ylim(1e-3, 3)
                elif key == "rate": a.plot(t, roll_rate(df, C["race"]["rate_window_days"]), "k-", lw=0.9); a.axhline(0, color="gray", lw=0.5); [a.axhline(s * C["race"]["authority_rate_threshold_km_day"], color="tab:orange", lw=0.6, ls="--") for s in (1, -1)]
                elif key == "lat": a.plot(t, np.rad2deg(df.lat_rad), "k-", lw=0.9)
                elif key == "a": a.plot(t, df.a_m / 1000, "k-", lw=0.9)
                elif key == "e": a.plot(t, df.e, "k-", lw=0.9)
                elif key == "om": a.plot(t, unwrap_deg(df.moon_omega_p_rad), "tab:blue"); a.plot(t, unwrap_deg(df.sun_omega_p_rad), "tab:orange", ls="--")
                elif key == "ip": a.plot(t, np.rad2deg(df.moon_i_p_rad), "tab:blue"); a.plot(t, np.rad2deg(df.sun_i_p_rad), "tab:orange", ls="--")
                if c == 0: a.set_ylabel(lab, fontsize=8)
                if k == 0: a.set_title(f"{r.role}: {r.tag}", fontsize=8)
                if k == 3: a.set_xlabel("day")
        fig.savefig(f"{OUT}/fig_race_{name}.png", dpi=90); plt.close(fig)
    # 4. cost bars by family (60-day, all cells) + race 180-day
    c60 = cost[(cost.horizon_days == 60.0) & ~cost.get("race", pd.Series(False, index=cost.index)).fillna(False).astype(bool)]
    fams = ["F2 (23.4)", "ref (50)", "F1 (63.4)", "F1 (116.6)"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    med = [c60[c60.family == f].dv_mps.median() for f in fams]; lo = [c60[c60.family == f].dv_mps.min() for f in fams]; hi = [c60[c60.family == f].dv_mps.max() for f in fams]
    x = np.arange(len(fams)); ax[0].bar(x, med, color="tab:blue", label="median over cells"); ax[0].errorbar(x, med, yerr=[np.subtract(med, lo), np.subtract(hi, med)], fmt="none", ecolor="k", capsize=3, label="min-max")
    ax[0].set_xticks(x, fams); ax[0].set_ylabel("dV to stay in band over 60 d (m/s)"); ax[0].set_title("all cells, deadband policy"); ax[0].legend(fontsize=8)
    rc = cost[cost.get("race", pd.Series(False, index=cost.index)).fillna(False).astype(bool)]
    labs = []; 
    for j, (_, r) in enumerate(race.iterrows()):
        for hd, off, col in ((60.001, -0.2, "tab:blue"), (180.001, 0.2, "tab:red")):
            v = rc[(rc.tag == r.tag) & (rc.horizon_days == hd)].dv_mps.iloc[0]; ax[1].bar(j + off, v, 0.4, color=col, label=("60 d" if hd < 100 else "180 d") if j == 0 else None)
        labs.append(f"{r.role[0]}:{r.tag[:6]}...{r.tag[-6:]}")
    ax[1].set_xticks(range(len(race)), labs, rotation=60, fontsize=7); ax[1].set_ylabel("dV (m/s)"); ax[1].set_title("race cells"); ax[1].legend(fontsize=8)
    fig.savefig(f"{OUT}/fig_correction_dv.png", dpi=130); plt.close(fig)
    pd.set_option("display.width", 250)
    print(c60.groupby("family").dv_mps.describe().round(2).to_string()); print(rc.round(2).to_string())

{"metrics": do_metrics, "race": do_race, "plots": do_plots}[STEP]()
