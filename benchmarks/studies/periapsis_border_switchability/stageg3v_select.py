"""Venus gate V-G(V): geometry selection BY RULE (plan R8.3). Python only, no flights, no GRAM. New file; imports the frozen Venus ranking and fixed-epoch code unchanged.
Epoch: offset 0 = the Venus Express anchor epoch fixed by the Venus ranking (2014-05-19T04:00:00; the same epoch for all three geometries, one GRAM epoch).
Rule (all constants of the ranking model; nothing new is sourced here):
  VG1   fine pool (i 20..160 every 5 deg, omega0 and Omega0 every 5 deg, 150,336 candidates; the planner's pool of venus_fixed_epoch_analysis.py): top 10 by closed-form dV floor (ties: smaller range),
        each verified with the full ranking model (deterministic floor and the 12-run jitter ensemble); VG1 = lowest verified ensemble floor (ties: smaller sd). Read from planner_fine_pool_picks.csv (same function).
  VG2   the same rule on the pool restricted to i in {63.4, 116.6} (omega0, Omega0 every 5 deg, 10,368 candidates).
  VGref the Venus Express geometry (i 89.876, omega0 75.505, Omega0 104.115) at the same epoch.
Writes geometry_gate_table_venus.json (+ .sha256) to the task directory. Usage: python3 -I stageg3v_select.py [--jobs N]"""
import hashlib, json, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import geometry_rank as G
import geometry_rank_venus as V
import venus_fixed_epoch as E

EP = 0.0
BASE_UTC = "2014-05-19T04:00:00"


def verify(ctx, c):
    r, n, S, d = G.evaluate(ctx, c, None, None); fl = G.jitter_floors(ctx, c, None, None)
    return dict(inc=c["inc"], omega0=c["omega"], raan0=c["raan"], floor_det=r["floor_mps"], floor_ens=float(fl.mean()), floor_ens_sd=float(fl.std()), R=r["R"], N=r["N"], regime=r["regime"],
                q_req=r["q_req"], W_eff=r["W_eff"], n_passes=r["n_passes"], days=r["days"], days_at_qcap=r["days_at_qcap"], floor_timesafe_det=r["floor_timesafe_mps"], floor_driftfree=r["floor_driftfree_mps"])


def pool_top10(ctx, nom, incs):
    oms = np.arange(0.0, 360.0, 5.0); Oms = np.arange(0.0, 360.0, 5.0)
    II, WW, OO = np.meshgrid(np.asarray(incs, float), oms, Oms, indexing="ij"); II, WW, OO = II.ravel(), WW.ravel(), OO.ravel()
    Rcf = np.empty(len(II))
    for s0 in range(0, len(II), 20000):
        sl = slice(s0, s0 + 20000); F, _ = E.features(ctx, nom, EP, np.radians(II[sl]), np.radians(WW[sl]), np.radians(OO[sl])); Rcf[sl] = F["R_cf"]
    cand = np.argsort(Rcf)[:3000]
    F, _ = E.features(ctx, nom, EP, np.radians(II[cand]), np.radians(WW[cand]), np.radians(OO[cand])); fl = np.array([E.dp_floor(s, nom.w, nom.W_eff) for s in F["_slow"]])
    order = np.lexsort((F["R_cf"], fl))[:10]
    return [dict(rank=rk, inc=float(II[cand[j]]), omega0=float(WW[cand[j]]), raan0=float(OO[cand[j]]), floor_cf=float(fl[j]), R_cf=float(F["R_cf"][j])) for rk, j in enumerate(order)], len(II)


def choose(ctx, top):
    for t in top: t.update(verify(ctx, dict(id="V", inc=t["inc"], epoch_off=EP, omega=t["omega0"], raan=t["raan0"])))
    best = sorted(top, key=lambda t: (round(t["floor_ens"], 6), t["floor_ens_sd"]))[0]; return best, top


def main():
    ctx, info = V.make_ctx(); nom = E.Nominal(ctx)
    pk = pd.read_csv(pathlib.Path(V.notes_dir()) / "venus_fixed_epoch_data/truth/analysis/planner_fine_pool_picks.csv"); pk = pk[pk.epoch_off == EP].sort_values("rank")
    top1 = [dict(rank=int(r.rank), inc=float(r.inc), omega0=float(r.omega0), raan0=float(r.raan0), floor_cf=float(r.floor_cf), R_cf=float(r.R_cf)) for r in pk.itertuples()]
    # recompute the closed-form top 10 for the whole fine pool as a check against the stored picks (150,336 candidates)
    incs = np.arange(20.0, 160.1, 5.0); chk, npool = pool_top10(ctx, nom, incs)
    same = all((a["inc"], a["omega0"], a["raan0"]) == (b["inc"], b["omega0"], b["raan0"]) for a, b in zip(top1, chk)); print("fine pool recomputed (%d candidates); top 10 equal to the stored picks: %s" % (npool, same), flush=True)
    vg1, t1 = choose(ctx, chk)
    top2, n2 = pool_top10(ctx, nom, [63.4, 116.6]); vg2, t2 = choose(ctx, top2); print("critical pool %d candidates" % n2, flush=True)
    ref = verify(ctx, dict(id="V", inc=89.876, epoch_off=EP, omega=75.505, raan=104.115)); ref.update(rank=-1, floor_cf=float("nan"), R_cf=float("nan"))
    # existing-table cross-check of VG2: best 63.4 / 116.6 candidate of the 7,296-candidate fine ranking (any epoch) and of the fixed-epoch truth grid at epoch 0
    rf = pd.read_csv(pathlib.Path(V.notes_dir()) / "venus_rank_data/rank_fine/geometry_rank.csv"); rc = rf[rf.inc_deg.round(1).isin([63.4, 116.6])].sort_values("floor_mps").iloc[0]
    tr = pd.read_csv(pathlib.Path(V.notes_dir()) / "venus_fixed_epoch_data/truth/truth.csv"); tc = tr[(tr.set == "grid") & (tr.epoch_off == EP) & tr.inc.round(1).isin([63.4, 116.6])].sort_values("floor").iloc[0]
    geoms = []
    for role, g in (("VG1", vg1), ("VG2", vg2), ("VGref", ref)):
        geoms.append(dict(id=role, role=role, in_default_set=True, inclination_deg=g["inc"], omega0_deg=g["omega0"], raan0_deg=g["raan0"], epoch_offset_days=EP, epoch_utc=BASE_UTC,
                          floor_ens_mps=round(g["floor_ens"], 4), floor_ens_sd_mps=round(g["floor_ens_sd"], 4), floor_det_mps=round(g["floor_det"], 4), regime=g["regime"], drift_range_ln=round(g["R"], 4),
                          net_drift_ln=round(g["N"], 4), q_req_w_cm2=round(g["q_req"], 4), w_eff_ln=round(g["W_eff"], 4), n_passes_model=int(g["n_passes"]), days_model=round(g["days"], 3),
                          days_at_qcap=round(g["days_at_qcap"], 3), floor_timesafe_det_mps=round(g["floor_timesafe_det"], 4)))
    T = dict(description="Venus gate V-G(V) geometries (plan R8.3), selected by rule from the Venus ranking model (constant H = 3.86 km; the mean-field replacement is in stageg3v_predictions.json).",
             selection_rule=dict(epoch="offset 0 = the Venus Express anchor epoch of the Venus ranking, " + BASE_UTC + " (all three geometries)",
                                 VG1="fine pool (i 20 to 160 every 5 deg, omega0 and Omega0 every 5 deg, %d candidates): top 10 by closed-form floor (ties: smaller range), each verified with the full ranking model; lowest 12-run jitter-ensemble floor (ties: smaller sd)" % npool,
                                 VG2="same rule on the pool i in {63.4, 116.6} (%d candidates)" % n2, VGref="the Venus Express geometry (i 89.876, omega0 75.505, Omega0 104.115) at the same epoch"),
             start=dict(ra0_m=float(ctx.cfg["orbit"]["ra0_m"]), ra_target_m=float(ctx.cfg["orbit"]["ra_target_m"]), h_ref_m=float(ctx.cfg["orbit"]["h_ref_m"]), true_anomaly_deg=180.0, d_ref_days=float(ctx.cfg["band"]["d_ref_days"]),
                        note="ra0 SOURCED vex_orvm_f.toml; ra_target and d_ref DERIVED by geometry_rank_venus.derive_campaign (flown VEX passes 2950 to 3005); start at apoapsis, epoch_utc is the starting-apoapsis time (ASSUMED as at Mars)"),
             geometries=geoms, verified_top10=dict(VG1=t1, VG2=t2),
             cross_checks=dict(fine_pool_top10_equals_stored_picks=bool(same), rank_fine_best_critical=dict(id=str(rc["id"]), inc=float(rc.inc_deg), omega0=float(rc.omega0_deg), raan0=float(rc.raan0_deg), epoch_off_d=float(rc.epoch_off_d), floor_det=float(rc.floor_mps)),
                               fixed_epoch_grid_best_critical_epoch0=dict(inc=float(tc.inc), omega0=float(tc.omega0), raan0=float(tc.raan0), floor_det=float(tc.floor))))
    p = pathlib.Path(V.notes_dir()) / "geometry_gate_table_venus.json"; raw = json.dumps(T, indent=1).encode(); p.write_bytes(raw); h = hashlib.sha256(raw).hexdigest()
    (p.parent / "geometry_gate_table_venus.json.sha256").write_text(h + "  geometry_gate_table_venus.json\n"); print(json.dumps([{k: g[k] for k in ("id", "inclination_deg", "omega0_deg", "raan0_deg", "floor_ens_mps", "floor_ens_sd_mps", "floor_det_mps", "regime", "drift_range_ln", "n_passes_model", "days_model")} for g in geoms], indent=1))
    print("cross-check", json.dumps(T["cross_checks"])); print("sha256", h)


if __name__ == "__main__":
    main()
