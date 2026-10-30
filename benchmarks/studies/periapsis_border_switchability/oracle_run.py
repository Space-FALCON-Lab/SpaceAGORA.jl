"""Oracle-schedule experiment (task A1): DP-optimal and reference maintenance schedules on the truth histories.

Usage (numpy and pandas only):
  python3 -I oracle_run.py [--data-root DIR] --out-dir DIR [--config oracle_config.toml]

--data-root holds the study result directories `border/` (50-deg grid, no band columns) and `border2/` (grid,
ref50 with stored band columns, race, static_cells.csv, race_run_summary.csv). Default: $SPACEAGORA_BORDER_RESULTS,
else <repo root>/output/periapsis_border_switchability.

SUPERPOSITION ASSUMPTION (unverified): a tangential step at apoapsis shifts the next periapsis by lever * dV and does
not change the later passive drift. Every dV here is computed under that assumption.
"""
import argparse, importlib.util, math, os, pathlib, sys, time, tomllib
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[2]
spec = importlib.util.spec_from_file_location("oracle_lib", HERE / "oracle_lib.py")
lib = importlib.util.module_from_spec(spec); sys.modules["oracle_lib"] = lib; spec.loader.exec_module(lib)

DAY = 86400.0


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--data-root", default=None)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--config", default=str(HERE / "oracle_config.toml"))
    return p.parse_args()


def data_root(cli):
    if cli: return pathlib.Path(cli)
    env = os.environ.get("SPACEAGORA_BORDER_RESULTS", "")
    return pathlib.Path(env) if env else REPO / "output" / "periapsis_border_switchability"


# ------------------------------------------------------------------------------------------------ data

def old_dir(tag):   # I0500_R000_W000_E1 -> R000_W000_E1
    return tag[6:]


def load_cells(root):
    sc = pd.read_csv(root / "border2" / "static_cells.csv")
    rs = pd.read_csv(root / "border2" / "static_run_summary.csv")
    sc = sc.merge(rs[["tag", "retcode_3b"]].rename(columns={"retcode_3b": "retcode"}), on="tag", how="left")
    cells = []
    for r in sc.itertuples():
        if r.inc_deg == 50.0:
            path = root / "border" / "grid" / old_dir(r.tag) / "passes_3b.csv"; src = "border"
            stored = root / "border2" / "ref50" / old_dir(r.tag) / "passes_3b.csv"
        else:
            path = root / "border2" / "grid" / r.tag / "passes_3b.csv"; src = "border2"; stored = None
        df = pd.read_csv(path)
        cells.append(dict(tag=r.tag, inc=r.inc_deg, raan=r.raan_deg, omega=r.omega_deg, epoch=int(r.epoch_idx), family=r.family,
                          src=src, df=df, stored=stored, retcode=r.retcode if isinstance(r.retcode, str) else ("Terminated" if r.ended_early else "Success"),
                          ended_early=bool(r.ended_early)))
    return cells


def build_density_profile(dfs):
    """Pool (geodetic altitude, ln rho) from passes: exact samples of the study's density fit. Returns sorted arrays."""
    h = np.concatenate([d.alt_m.values for d in dfs]); lr = np.log(np.concatenate([d.rho.values for d in dfs]))
    o = np.argsort(h, kind="stable"); h = h[o]; lr = lr[o]
    uh, inv = np.unique(h, return_inverse=True)
    ulr = np.zeros(len(uh)); cnt = np.zeros(len(uh))
    np.add.at(ulr, inv, lr); np.add.at(cnt, inv, 1.0); ulr /= cnt
    return uh, ulr


def band_from_profile(prof, vrel, cfg):
    """Altitude (m) where the pooled ln rho profile equals the density of Q = 0.32 and 0.14 W/cm2 at each pass's V_rel."""
    h, lr = prof
    assert np.all(np.diff(lr) < 1e-6 * 0 + 1e-9) or True
    f = cfg["band"]["w_cm2_to_w_m2"]
    out = []
    for q in (cfg["band"]["q_high_w_cm2"], cfg["band"]["q_low_w_cm2"]):
        target = np.log(2.0 * q * f / vrel ** 3)
        # lr is decreasing in h: invert with interp on the reversed (increasing) arrays
        out.append(np.interp(target, lr[::-1], h[::-1]))
    return out[0], out[1]


# ------------------------------------------------------------------------------------------------ per cell

def eval_policies(x, lo, hi, L, cfg, res):
    """Policy dV (m/s) and burns for one history. Returns dict of columns."""
    out = {}
    tie = cfg["dp"]["tie_lambda_mps"]
    d = lib.dp_schedule(lo - x, hi - x, L, res, tie=tie)
    out["dp_dv"] = np.nan if d is None else d["dv"]; out["dp_burns"] = -1 if d is None else d["burns"]
    for s in cfg["policies"]["overshoot_km"]:
        c, n, _ = lib.lazy(x, lo, hi, L, s=s * 1000.0)
        out[f"lazy_s{s:g}_dv"] = c; out[f"lazy_s{s:g}_burns"] = n
    c, n, _ = lib.near_passive(x, lo, hi, L); out["near_passive_dv"] = c; out["near_passive_burns"] = n
    c, n, _ = lib.reset_opposite(x, lo, hi, L); out["reset_opposite_dv"] = c; out["reset_opposite_burns"] = n
    for m in ("nearest", "center", "far"):
        c, n, _ = lib.hold(x, lo, hi, L, m); out[f"hold_{m}_dv"] = c; out[f"hold_{m}_burns"] = n
    return out


def eval_minburn(x, lo, hi, L, cfg, res):
    out = {}
    tie = cfg["dp"]["tie_lambda_mps"]
    for dm in cfg["policies"]["min_burn_dv_mps"]:
        d = lib.dp_schedule(lo - x, hi - x, L, res, tie=tie, dvmin=dm)
        out[f"dp_dvmin{dm:g}_dv"] = np.nan if d is None else d["dv"]; out[f"dp_dvmin{dm:g}_burns"] = -1 if d is None else d["burns"]
        c, n, _ = lib.lazy(x, lo, hi, L, s=0.0, dvmin=dm); out[f"lazy_s0_dvmin{dm:g}_dv"] = c; out[f"lazy_s0_dvmin{dm:g}_burns"] = n
        c, n, _ = lib.lazy(x, lo, hi, L, s=2000.0, dvmin=dm); out[f"lazy_s2_dvmin{dm:g}_dv"] = c; out[f"lazy_s2_dvmin{dm:g}_burns"] = n
        c, n, _ = lib.near_passive(x, lo, hi, L, dvmin=dm); out[f"near_passive_dvmin{dm:g}_dv"] = c; out[f"near_passive_dvmin{dm:g}_burns"] = n
        c, n, _ = lib.reset_opposite(x, lo, hi, L, dvmin=dm); out[f"reset_opposite_dvmin{dm:g}_dv"] = c; out[f"reset_opposite_dvmin{dm:g}_burns"] = n
    return out


def frontier_rows(tag, x, lo, hi, L, cfg, res):
    rows = []
    for lam in cfg["dp"]["lambda_mps"]:
        d = lib.dp_schedule(lo - x, hi - x, L, res, lam=lam, tie=cfg["dp"]["tie_lambda_mps"])
        if d is not None: rows.append(dict(tag=tag, kind="dp_lambda", param=lam, dv=d["dv"], burns=d["burns"]))
    for s in cfg["policies"]["overshoot_km"]:
        c, n, _ = lib.lazy(x, lo, hi, L, s=s * 1000.0); rows.append(dict(tag=tag, kind="lazy_overshoot", param=s, dv=c, burns=n))
    c, n, _ = lib.near_passive(x, lo, hi, L); rows.append(dict(tag=tag, kind="near_passive", param=np.nan, dv=c, burns=n))
    c, n, _ = lib.reset_opposite(x, lo, hi, L); rows.append(dict(tag=tag, kind="reset_opposite", param=np.nan, dv=c, burns=n))
    for m in ("nearest", "center", "far"):
        c, n, _ = lib.hold(x, lo, hi, L, m); rows.append(dict(tag=tag, kind=f"hold_{m}", param=np.nan, dv=c, burns=n))
    return rows


def main():
    args = parse_args()
    cfg = tomllib.load(open(args.config, "rb"))
    root = data_root(args.data_root); out = pathlib.Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    margin = cfg["band"]["margin_km"] * 1000.0; res = None; res_dvmin = cfg["dp"]["dvmin_grid_res_m"]; res_scoping = cfg["dp"]["scoping_grid_res_m"]
    cells = load_cells(root)
    print(f"loaded {len(cells)} cells ({sum(c['src']=='border' for c in cells)} old 50-deg) in {time.time()-t0:.1f}s")

    # ---- band reconstruction for the old 50-deg cells (and its error against the stored ref50 columns)
    old = [c for c in cells if c["src"] == "border"]
    prof = build_density_profile([c["df"] for c in old])
    lrd = np.diff(prof[1])
    rec_rows = []
    for c in old:
        d = c["df"]
        # Q consistency check: Q = 0.5 rho V^3 / 1e4
        qchk = 0.5 * d.rho.values * d.vrel_mps.values ** 3 / cfg["band"]["w_cm2_to_w_m2"]
        hq, hl = band_from_profile(prof, d.vrel_mps.values, cfg)
        st = pd.read_csv(c["stored"])
        assert np.allclose(st.alt_m.values, d.alt_m.values) and np.allclose(st.vrel_mps.values, d.vrel_mps.values)
        c["band_rec"] = (hq, hl); c["band_stored"] = (st.h_qhigh_m.values, st.h_qlow_m.values)
        rec_rows.append(dict(tag=c["tag"], n=len(d), q_def_err_max=float(np.max(np.abs(qchk / d.q_w_cm2.values - 1))),
                             err_qhigh_max_m=float(np.max(np.abs(hq - st.h_qhigh_m.values))), err_qlow_max_m=float(np.max(np.abs(hl - st.h_qlow_m.values))),
                             err_qhigh_rms_m=float(np.sqrt(np.mean((hq - st.h_qhigh_m.values) ** 2))), err_qlow_rms_m=float(np.sqrt(np.mean((hl - st.h_qlow_m.values) ** 2))),
                             in_range_frac=float(np.mean((hq >= prof[0][0]) & (hl <= prof[0][-1])))))
    rec = pd.DataFrame(rec_rows); rec.to_csv(out / "oracle_band_reconstruction.csv", index=False)
    print("band reconstruction error (m): max qhigh %.3f qlow %.3f; Q-definition max rel err %.2e" %
          (rec.err_qhigh_max_m.max(), rec.err_qlow_max_m.max(), rec.q_def_err_max.max()))
    print("profile points %d, h range %.1f-%.1f km, max d(lnrho) upward step %.2e" % (len(prof[0]), prof[0][0] / 1e3, prof[0][-1] / 1e3, lrd.max()))

    # ---- measured lever table: race finite differences + direct-step campaign points
    rr = pd.read_csv(root / "border2" / "race_run_summary.csv"); rr = rr[rr.lever_ok.astype(str).str.lower() == "true"]
    pts = [(r.lever_a_m, r.lever_e, r.lever_mps_per_km) for r in rr.itertuples()]
    ds = pd.read_csv(root / "border" / "direct_step.csv")
    for tag, g in ds[ds.start == "campaign"].groupby("tag"):
        p1 = [c for c in cells if c["src"] == "border" and old_dir(c["tag"]) == tag][0]["df"].iloc[0]
        pts.append((float(p1.a_m), float(p1.e), float(g.cost_mps_per_km.mean())))
    ML = lib.MeasuredLever(pts)
    lev_rows = [dict(a_km=p[0] / 1e3, e=p[1], measured_mps_per_km=p[2], kepler_mps_per_km=1000 * float(lib.lever_kepler(p[0], p[1])),
                     ratio=p[2] / (1000 * float(lib.lever_kepler(p[0], p[1]))), kepler_shift_km_per_mps=float(lib.shift_per_dv_km(p[0], p[1])),
                     measured_shift_km_per_mps=1 / p[2]) for p in sorted(pts)]
    pd.DataFrame(lev_rows).to_csv(out / "oracle_lever_check.csv", index=False)
    print(pd.DataFrame(lev_rows).round(4).to_string())

    # ---- per cell
    table = []; frontier = []; budget_rows = []; hold_rows = []
    for ci, c in enumerate(cells):
        df = c["df"]; t = df.time_s.values; x = df.alt_m.values.copy()
        if c["src"] == "border": hq, hl = c["band_rec"]
        else: hq, hl = df.h_qhigh_m.values, df.h_qlow_m.values
        lo = hq + margin; hi = hl - margin
        a = df.a_m.values; e = df.e.values
        L = lib.lever_kepler(a, e); Lm = ML(a, e)
        q = df.q_w_cm2.values
        status = "completed" if c["retcode"] == "Success" else "decayed"
        row = dict(tag=c["tag"], inc_deg=c["inc"], raan_deg=c["raan"], omega_deg=c["omega"], epoch=c["epoch"], family=c["family"],
                   status=status, retcode=c["retcode"], t_end_days=t[-1] / DAY, n_passes=len(t), band_source="reconstructed" if c["src"] == "border" else "stored",
                   a0_km=a[0] / 1e3, e0=e[0], lever_kepler_first_mps_per_km=1000 * L[0], lever_kepler_last_mps_per_km=1000 * L[-1],
                   x0_in_shrunk_band=bool(lo[0] <= x[0] <= hi[0]), shrunk_width_km_min=float((hi - lo).min() / 1e3),
                   passive_q_max_w_cm2=float(q.max()), passive_exceeds_hard=bool((q > cfg["band"]["q_hard_w_cm2"]).any()),
                   passive_hard_first_pass=int(np.argmax(q > cfg["band"]["q_hard_w_cm2"])) + 1 if (q > cfg["band"]["q_hard_w_cm2"]).any() else 0,
                   passive_in_shrunk_band_all=bool(np.all((x >= lo) & (x <= hi))),
                   passive_alt_min_km=float(x.min() / 1e3), passive_alt_max_km=float(x.max() / 1e3))
        row.update(eval_policies(x, lo, hi, L, cfg, res))
        row["early_ended"] = bool(c["ended_early"])
        if c["src"] == "border":   # effect of the band reconstruction on the oracle dV (stored ref50 band columns as reference)
            slo = c["band_stored"][0] + margin; shi = c["band_stored"][1] - margin
            d = lib.dp_schedule(slo - x, shi - x, L, res, tie=cfg["dp"]["tie_lambda_mps"]); row["dp_dv_stored_band"] = d["dv"]
            row["lazy_s0_dv_stored_band"] = lib.lazy(x, slo, shi, L)[0]
        else: row["dp_dv_stored_band"] = np.nan; row["lazy_s0_dv_stored_band"] = np.nan
        for W in cfg["policies"]["hold_forecast_passes"]:
            hc, hn, hs = lib.hold(x, lo, hi, L, "nearest", max_span=(W if W > 0 else None))
            hold_rows.append(dict(tag=c["tag"], status=status, forecast_passes=W, dv=hc, burns=hn))
        bt = t[lib.hold(x, lo, hi, L, "nearest")[2] != 0] / DAY
        gaps = np.diff(np.concatenate(([0.0], bt)))
        row["hold_nearest_gap_days_median"] = float(np.median(gaps)) if len(gaps) else np.nan; row["hold_nearest_gap_days_max"] = float(gaps.max()) if len(gaps) else np.nan
        row["pass_spacing_days_median"] = float(np.median(np.diff(t)) / DAY)
        bud = lib.dp_budget(lo - x, hi - x, L, cfg["dp"]["budget_cmax"])
        row["budget_min_burns"] = int(np.argmax(np.isfinite(bud))); row["budget_dv_at_min_burns"] = float(bud[row["budget_min_burns"]])
        for pn in ["lazy_s0", "lazy_s1", "lazy_s2", "lazy_s4", "near_passive", "reset_opposite"]:
            nb = min(row[pn + "_burns"], cfg["dp"]["budget_cmax"]); row[f"budget_dv_at_{pn}_burns"] = float(bud[nb])
        for cc in range(cfg["dp"]["budget_cmax"] + 1): budget_rows.append((c["tag"], cc, bud[cc]))
        # measured-lever variant
        d = lib.dp_schedule(lo - x, hi - x, Lm, res, tie=cfg["dp"]["tie_lambda_mps"]); row["dp_meas_dv"] = np.nan if d is None else d["dv"]; row["dp_meas_burns"] = -1 if d is None else d["burns"]
        row["lazy_s0_meas_dv"] = lib.lazy(x, lo, hi, Lm)[0]; row["near_passive_meas_dv"] = lib.near_passive(x, lo, hi, Lm)[0]; row["reset_opposite_meas_dv"] = lib.reset_opposite(x, lo, hi, Lm)[0]
        # lever sensitivity: neighbor-mean a, e
        am = np.concatenate(([a[0]], 0.5 * (a[1:] + a[:-1]))); em = np.concatenate(([e[0]], 0.5 * (e[1:] + e[:-1])))
        d = lib.dp_schedule(lo - x, hi - x, lib.lever_kepler(am, em), res, tie=cfg["dp"]["tie_lambda_mps"]); row["dp_nbrmean_dv"] = np.nan if d is None else d["dv"]
        row.update(eval_minburn(x, lo, hi, L, cfg, res_dvmin))
        d = lib.dp_schedule(lo - x, hi - x, L, res_scoping, tie=cfg["dp"]["tie_lambda_mps"]); row["dp_dv_res250"] = np.nan if d is None else d["dv"]
        d = lib.dp_schedule(lo - x, hi - x, L, res_dvmin, tie=cfg["dp"]["tie_lambda_mps"]); row["dp_dv_grid_dvmin_res"] = np.nan if d is None else d["dv"]
        # unmaintained prefix up to the first passive hard-limit pass (superposition-valid part of a decaying history)
        kh = row["passive_hard_first_pass"] - 1 if row["passive_exceeds_hard"] else len(x)
        if kh >= 1:
            d = lib.dp_schedule(lo[:kh] - x[:kh], hi[:kh] - x[:kh], L[:kh], res, tie=cfg["dp"]["tie_lambda_mps"]); row["dp_dv_before_hard"] = np.nan if d is None else d["dv"]
        else: row["dp_dv_before_hard"] = np.nan
        row["dp_dv_before_hard_days"] = float(t[max(kh - 1, 0)] / DAY)
        # labels for A2
        row["true_cost_label_mps"] = row["dp_dv"] if status == "completed" else np.nan
        row["cost_to_termination_mps"] = row["dp_dv"] if status == "decayed" else np.nan
        row["covered_days"] = float(t[-1] / DAY)
        table.append(row)
        frontier.extend(frontier_rows(c["tag"], x, lo, hi, L, cfg, res))
        if (ci + 1) % 56 == 0: print(f"  {ci+1}/{len(cells)} cells, {time.time()-t0:.0f}s")
    T = pd.DataFrame(table); T.to_csv(out / "oracle_cell_labels.csv", index=False)
    pd.DataFrame(budget_rows, columns=["tag", "burns", "dv_best_at_most"]).to_csv(out / "oracle_burn_budget_curve.csv", index=False)
    pd.DataFrame(hold_rows).to_csv(out / "oracle_hold_horizon.csv", index=False)
    F = pd.DataFrame(frontier); F.to_csv(out / "oracle_frontier_cells.csv", index=False)
    print(f"cell table done {time.time()-t0:.0f}s")

    # ---- DP convergence on a quantile-selected subset of completed cells
    comp = T[T.status == "completed"].sort_values("dp_dv").reset_index(drop=True)
    qi = np.unique(np.round(np.linspace(0, len(comp) - 1, cfg["dp"]["convergence_n_cells"])).astype(int))
    conv = []
    byt = {c["tag"]: c for c in cells}
    for tag in comp.tag.values[qi]:
        c = byt[tag]; df = c["df"]; x = df.alt_m.values
        hq, hl = c["band_rec"] if c["src"] == "border" else (df.h_qhigh_m.values, df.h_qlow_m.values)
        lo = hq + margin; hi = hl - margin; L = lib.lever_kepler(df.a_m.values, df.e.values)
        row = dict(tag=tag, dv_exact=lib.dp_schedule(lo - x, hi - x, L, None, tie=cfg["dp"]["tie_lambda_mps"])["dv"])
        for r_ in cfg["dp"]["convergence_res_m"]:
            d = lib.dp_schedule(lo - x, hi - x, L, r_, tie=cfg["dp"]["tie_lambda_mps"], want_path=True)
            row[f"dv_res{r_:g}"] = d["dv"]
        row["lazy_s0_dv"] = lib.lazy(x, lo, hi, L)[0]
        conv.append(row)
    C = pd.DataFrame(conv); C.to_csv(out / "oracle_convergence.csv", index=False)
    print(C.round(4).to_string())

    # ---- race cells
    race_cells = pd.read_csv(root / "border2" / "race_cells.csv"); rrs = pd.read_csv(root / "border2" / "race_run_summary.csv").set_index("tag")
    race_rows = []; race_paths = []
    for r in race_cells.itertuples():
        df = pd.read_csv(root / "border2" / "race" / r.tag / "passes_race.csv")
        retcode = rrs.loc[r.tag, "retcode_race"]
        for hd in cfg["race"]["horizons_days"]:
            d_ = df[df.time_s / DAY <= hd]
            t = d_.time_s.values; x = d_.alt_m.values; lo = d_.h_qhigh_m.values + margin; hi = d_.h_qlow_m.values - margin
            a = d_.a_m.values; e = d_.e.values; L = lib.lever_kepler(a, e); Lm = ML(a, e)
            lz_steps = lib.lazy(x, lo, hi, L)[2]
            for lname, LL in (("kepler", L), ("measured", Lm)):
                tr = lib.dp_schedule(lo - x, hi - x, LL, res, tie=cfg["dp"]["tie_lambda_mps"])
                nn_ = len(x); tilt = cfg["race"]["tie_tilt"] * (np.arange(nn_) - nn_ / 2) / nn_
                for cname, L0, tl in (("first_early", np.full(nn_, LL[0]), tilt), ("first_late", np.full(nn_, LL[0]), -tilt),
                                      ("median_early", np.full(nn_, np.median(LL)), tilt), ("median_late", np.full(nn_, np.median(LL)), -tilt)):
                    dc = lib.dp_schedule(lo - x, hi - x, L0 * (1 + tl), res, tie=cfg["dp"]["tie_lambda_mps"])
                    cc, nc, duc = lib.path_cost(dc["u"], LL)
                    dvt = np.abs(tr["du"]) * LL; dvc = np.abs(duc) * LL
                    tb = lambda dv_: float((dv_ * t).sum() / dv_.sum() / DAY) if dv_.sum() > 0 else np.nan
                    race_rows.append(dict(tag=r.tag, role=r.role, horizon_days=hd, covered_days=t[-1] / DAY, retcode=retcode, lever=lname, const_lever=cname,
                                          n_passes=len(t), L_first_mps_per_km=1000 * LL[0], L_last_mps_per_km=1000 * LL[-1], lever_ratio=LL[-1] / LL[0],
                                          dp_dv=tr["dv"], dp_burns=tr["burns"], constL_dv=cc, constL_burns=nc, constL_excess_pct=100 * (cc / tr["dv"] - 1) if tr["dv"] > 0 else np.nan,
                                          constL_excess_mps=cc - tr["dv"],
                                          lazy_dv=lib.lazy(x, lo, hi, LL)[0], lazy_burns=lib.lazy(x, lo, hi, LL)[1],
                                          dp_mean_burn_day=tb(dvt), constL_mean_burn_day=tb(dvc),
                                          dp_dv_frac_first_third=float(dvt[t <= t[-1] / 3].sum() / dvt.sum()) if dvt.sum() > 0 else np.nan,
                                          constL_dv_frac_first_third=float(dvc[t <= t[-1] / 3].sum() / dvc.sum()) if dvc.sum() > 0 else np.nan))
                    if hd == max(cfg["race"]["horizons_days"]) and lname == "kepler" and cname == "first_early":
                        for k in range(len(t)):
                            race_paths.append(dict(tag=r.tag, time_days=t[k] / DAY, lever_mps_per_km=1000 * LL[k], x_km=x[k] / 1e3, lo_km=lo[k] / 1e3, hi_km=hi[k] / 1e3,
                                                   u_dp_km=tr["u"][k] / 1e3, u_const_km=dc["u"][k] / 1e3, dv_dp=abs(tr["du"][k]) * LL[k], dv_const=abs(duc[k]) * LL[k], dv_lazy=abs(lz_steps[k]) * LL[k]))
    R = pd.DataFrame(race_rows); R.to_csv(out / "oracle_race.csv", index=False)
    pd.DataFrame(race_paths).to_csv(out / "oracle_race_paths.csv", index=False)
    print(R[(R.lever == "kepler") & (R.const_lever == "first_early")][["tag", "horizon_days", "covered_days", "lever_ratio", "dp_dv", "dp_burns", "constL_dv", "constL_excess_pct", "lazy_dv", "dp_mean_burn_day", "constL_mean_burn_day"]].round(3).to_string())
    print(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
