"""Follow-up to venus_fixed_epoch_analysis.py (offline Python): does the closed-form Sun path survive an independent numerical integration, and does the re-evolving second-order form fix
what the frozen-element form gets wrong? Also: planner picks re-scored by the Cowell integration.

Usage: python3 -I venus_fixed_epoch_followup.py TRUTH_DIR [--jobs N]     (TRUTH_DIR holds truth.csv and analysis/ from the main run)
Reads analysis/truth_and_parameters.csv, analysis/planner_*.csv; writes analysis/followup_*.csv, followup_results.json, followup_report.txt.
"""
import json, pathlib, sys, time
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import geometry_rank_venus as V
import venus_fixed_epoch as E
import venus_fixed_epoch_cowell as C

LOG = []


def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.append(s)


def reg(y, p):
    return dict(n=len(y), R2=E.r2(y, p), MAE=E.mae(y, p), medAE=float(np.median(np.abs(y - p))), slope=float((y * p).sum() / (p * p).sum()))


def main(td, jobs):
    import multiprocessing as mp
    out = pathlib.Path(td); an = out / "analysis"; t0 = time.time(); ctx, info = V.make_ctx(); nom = E.Nominal(ctx); C._W.update(ctx=ctx, nom=nom); res = {}
    D = pd.read_csv(an / "truth_and_parameters.csv"); gr = D[D.set == "grid"].reset_index(drop=True)
    # ---- A. evolving closed form vs the ranking model (which freezes the elements) over the whole grid
    ev_floor = np.zeros(len(gr)); ev_R = np.zeros(len(gr)); ev_N = np.zeros(len(gr))
    for ep, g in gr.groupby("epoch_off"):
        p = E.closed_form_evolving(ctx, nom, ep, np.radians(g.inc.values), np.radians(g.omega0.values), np.radians(g.raan0.values)) + nom.other[None]; s = E.window(p)
        ev_floor[g.index.values] = [E.dp_floor(x, nom.w, nom.W_eff) for x in s]; ev_R[g.index.values] = np.ptp(s, axis=1); ev_N[g.index.values] = s[:, -1] - s[:, 0]
    gr["floor_ev"], gr["R_ev"], gr["N_ev"] = ev_floor, ev_R, ev_N; gr.to_csv(an / "followup_grid_with_ev.csv", index=False)
    res["ev_vs_ranking_model"] = dict(R=reg(gr.R.values, gr.R_ev.values), N=reg(gr.N.values, gr.N_ev.values), floor_mae=E.mae(gr.floor.values, gr.floor_ev.values), floor_spearman=float(pd.Series(gr.floor).corr(pd.Series(gr.floor_ev), method="spearman")),
                                      frozen_form_R=reg(gr.R.values, gr.R_cf.values))
    say("evolving form vs ranking model (which freezes the elements):", json.dumps(res["ev_vs_ranking_model"]))
    # ---- B. Cowell check on a list: polar and other-inclination sets at the headline epochs, plus the planner picks
    cands = []
    for ep in E.HEAD:
        for om in range(0, 360, 30):
            for ra in (0, 90, 180, 270): cands.append(dict(set="polar", inc=89.876, epoch_off=ep, omega=float(om), raan=float(ra)))
        for inc in (39.2, 63.4, 116.6):
            for om in (0, 90):
                for ra in (0, 180): cands.append(dict(set="other_inc", inc=inc, epoch_off=ep, omega=float(om), raan=float(ra)))
        g = gr[gr.epoch_off == ep]
        for nm, col in (("top_cf", "floor_cf"), ("top_ev", "floor_ev")):
            for r in g.sort_values(col).head(30).itertuples(): cands.append(dict(set=nm, inc=float(r.inc), epoch_off=ep, omega=float(r.omega0), raan=float(r.raan0)))
        for r in g.sample(30, random_state=1).itertuples(): cands.append(dict(set="random", inc=float(r.inc), epoch_off=ep, omega=float(r.omega0), raan=float(r.raan0)))
    pk = pd.read_csv(an / "planner_fine_pool_picks.csv")
    for r in pk[pk["rank"] < 3].itertuples(): cands.append(dict(set="fine_pick", inc=float(r.inc), epoch_off=float(r.epoch_off), omega=float(r.omega0), raan=float(r.raan0)))
    for ep in E.EPOCHS: cands.append(dict(set="vex", inc=89.876, epoch_off=ep, omega=75.505, raan=104.115))
    for i, c in enumerate(cands): c["id"] = "C%04d" % i
    say(len(cands), "Cowell candidates"); t1 = time.time()
    with mp.get_context("fork").Pool(jobs) as pool: rows = pool.map(C.one, cands, chunksize=2)
    R = pd.DataFrame(rows)
    # ev path and floor for the same candidates
    ev_f, ev_r = [], []
    for c in cands:
        a = [np.array([x]) for x in (np.radians(c["inc"]), np.radians(c["omega"]), np.radians(c["raan"]))]; p = E.closed_form_evolving(ctx, nom, c["epoch_off"], *a)[0] + nom.other; s = E.window(p); ev_f.append(E.dp_floor(s, nom.w, nom.W_eff)); ev_r.append(float(np.ptp(s)))
    R["floor_ev"], R["R_ev"] = ev_f, ev_r
    # ranking-model truth lookup for grid-members
    key = lambda d: d.inc.round(3).astype(str) + "_" + d.epoch_off.astype(str) + "_" + d.omega0.astype(str) + "_" + d.raan0.astype(str)
    lk = pd.Series(gr.floor.values, index=key(gr).values); R["floor_model"] = key(R).map(lk).values
    R.to_csv(an / "followup_cowell.csv", index=False); say("Cowell integrations: %.0f s wall" % (time.time() - t1))
    ok = R[R.set != "vex"]
    res["cowell_vs_forms"] = dict(frozen_form_R=reg(ok.R_num.values, ok.R_cf.values), evolving_form_R=reg(ok.R_num.values, ok.R_ev.values), frozen_form_floor_MAE=E.mae(ok.floor_num.values, ok.floor_cf.values), evolving_form_floor_MAE=E.mae(ok.floor_num.values, ok.floor_ev.values),
                                  frozen_sun_path_rms_median=float(ok.sun_path_rms_err.median()), n=len(ok))
    m = ok.dropna(subset=["floor_model"]); res["cowell_vs_ranking_model_floor"] = dict(n=len(m), MAE=E.mae(m.floor_num.values, m.floor_model.values), spearman=float(m.floor_num.corr(m.floor_model, method="spearman")))
    say(json.dumps(res["cowell_vs_forms"], indent=1)); say(json.dumps(res["cowell_vs_ranking_model_floor"]))
    by = []
    for st, g in ok.groupby("set"):
        by.append(dict(set=st, n=len(g), floor_num_median=g.floor_num.median(), floor_model_median=g.floor_model.median(), floor_cf_median=g.floor_cf.median(), floor_ev_median=g.floor_ev.median(), floor_num_min=g.floor_num.min(), frac_num_lt0p1=(g.floor_num < 0.1).mean(), frac_num_lt1=(g.floor_num < 1).mean(),
                       MAE_cf=E.mae(g.floor_num.values, g.floor_cf.values), MAE_ev=E.mae(g.floor_num.values, g.floor_ev.values)))
    B = pd.DataFrame(by); B.to_csv(an / "followup_cowell_by_set.csv", index=False); say(B.round(3).to_string())
    # ---- C. per-epoch planner re-scoring in the Cowell world (headline epochs)
    rows = []
    for ep in E.HEAD:
        g = R[(R.epoch_off == ep)]; vx = R[(R.set == "vex") & (R.epoch_off == ep)].iloc[0]; cf_pick = gr[gr.epoch_off == ep].sort_values(["floor_cf", "R_cf"]).iloc[0]; ev_pick = gr[gr.epoch_off == ep].sort_values(["floor_ev", "R_ev"]).iloc[0]
        def cw(p):
            h = R[(R.epoch_off == ep) & (R.inc.round(3) == round(p.inc, 3)) & (R.omega0 == p.omega0) & (R.raan0 == p.raan0)]; return h.floor_num.iloc[0] if len(h) else np.nan
        pool = g[g.set.isin(["top_cf", "top_ev", "random", "polar", "other_inc"])]
        rows.append(dict(epoch_off=ep, vex_cowell=vx.floor_num, cf_pick_model=cf_pick.floor, cf_pick_cowell=cw(cf_pick), ev_pick_model=ev_pick.floor, ev_pick_cowell=cw(ev_pick), cowell_best_in_evaluated=pool.floor_num.min(), random_cowell_median=g[g.set == "random"].floor_num.median(),
                         fine_pick1_cowell=g[g.set == "fine_pick"].sort_values("floor_cf").floor_num.iloc[0] if (g.set == "fine_pick").any() else np.nan, cf_pick_geom="i%.1f/w%.0f/O%.0f" % (cf_pick.inc, cf_pick.omega0, cf_pick.raan0), ev_pick_geom="i%.1f/w%.0f/O%.0f" % (ev_pick.inc, ev_pick.omega0, ev_pick.raan0)))
    PL = pd.DataFrame(rows); PL.to_csv(an / "followup_planner_cowell.csv", index=False); say(PL.round(3).to_string())
    # ---- D. sensitivity of the pick to errors in the initial angles (ranking model, closed form; one full-model eval per point on the best polar pick per epoch)
    sens = []
    for ep in E.HEAD:
        b = gr[gr.epoch_off == ep].sort_values(["floor_cf", "R_cf"]).iloc[0]
        for dw in (-10, -5, -2, 0, 2, 5, 10):
            for dO in (-10, -5, 0, 5, 10):
                a = [np.array([x]) for x in (np.radians(b.inc), np.radians(b.omega0 + dw), np.radians(b.raan0 + dO))]; p = E.closed_form(ctx, nom, ep, *a)["path3"][0] + nom.other; s = E.window(p)
                sens.append(dict(epoch_off=ep, base="i%.1f/w%.0f/O%.0f" % (b.inc, b.omega0, b.raan0), d_omega0=dw, d_raan0=dO, floor_cf=E.dp_floor(s, nom.w, nom.W_eff)))
    S = pd.DataFrame(sens); S.to_csv(an / "followup_pick_angle_sensitivity.csv", index=False)
    piv = S.groupby(["epoch_off", "d_omega0"]).floor_cf.max().unstack(); say("worst floor over Omega0 error of +-10 deg, by omega0 error (closed form), rows epoch:\n" + piv.round(2).to_string())
    piv2 = S.groupby(["epoch_off", "d_raan0"]).floor_cf.max().unstack(); say("worst floor over omega0 error of +-10 deg, by Omega0 error:\n" + piv2.round(2).to_string())
    (an / "followup_results.json").write_text(json.dumps(res, indent=1, default=float)); (an / "followup_report.txt").write_text("\n".join(LOG)); say("done %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    a = sys.argv[1:]; main(a[0], int(a[a.index("--jobs") + 1]) if "--jobs" in a else 6)
