"""Revision 20d (no final periapsis maneuver; time only from the campaign and the orbit period). Analysis only: Python, no Julia, no flights.

wait_rule(): the arrival-time options of Revision 20d item 2 as a pure function. The run stops at the first apoapsis at or below the target (error err_km <= 0); waiting k more final orbits (each a further
drag pass) is admissible while err_km - k * last_drop_km >= -tol_km. Arrival at the closing apoapsis passage k = t_close + k * T_final.
Commands (python3 -I planner_rev20d.py CMD): selftest | baselines | replay | replay_tables | curves | corridor | p20_10 | all (all = everything; `replay` takes about 25 min).   Output: tasks/.../rev20_data/rev20d/*.json and *.md (tables for results_rev20d.md).
Labels. ASSUMED: tol_km = 25 (Rev 14/16b), next drops = the last flown drop, K_w = 3 final-orbit periods (user's; 1, 3, 5, 10 swept), the period of a waited orbit = the closing orbit's (the apoapsis changes by <= tol).
SOURCED: corridor [0.05, 0.25] W/cm2 (planner_config [common.band]), Q_s never-exceed 0.65 W/cm2.
Replay: planner_rev20_replay.replay(trace=True) with GUARD_V2 unchanged; same seeds, ic and t* as its evaluate/rev17 commands (checked against the stored files by command `replay`)."""
import argparse, glob, json, math, os, pathlib, sys
import numpy as np, pandas as pd

HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_analysis as _AN
TD = pathlib.Path(PLB.task_dir()); OUT = TD / "rev20_data" / "rev20d"
ROOT = pathlib.Path(os.environ.get("SPACEAGORA_PLANNER_ROOT", str(pathlib.Path.home() / "Documents" / "SpaceAGORA-planner-data" / "planner")))
TOL_KM = 25.0            # ASSUMED (Revision 16b): final apoapsis error tolerance
K_W = 3.0                # ASSUMED (user): window half-width in final-orbit periods
K_SWEEP = (1.0, 3.0, 5.0, 10.0)
K_CAP = 50               # numerical cap on waited orbits
Q_BAND = (0.05, 0.25)    # SOURCED planner_config [common.band] q_low, q_high (W/cm2)
QS_HARD = 0.65           # SOURCED planner_config q_never_exceed_qs_w_cm2
OFFS = (0, 1, 3, 5, 10, 20, 40, 80)       # |t* - reference| in final orbits (numerical grid)


def wait_rule(err_km, t_close_d, last_drop_km, T_final_d, t_star_d, K_w=K_W, tol_km=TOL_KM, k_cap=K_CAP):
    """Arrival options and window verdict. err_km: final apoapsis error at the first apoapsis at or below the target (<= 0); t_close_d: its arrival time; last_drop_km: last flown apoapsis drop (assumed for
    every further drop); T_final_d: final-orbit period; t_star_d: target time; K_w: window half-width in periods. Waiting is admissible only if the first closing is itself inside the tolerance
    (apo_ok); k_max = floor((err + tol) / drop) (0 when the drop is not positive or not finite). Returns dict(apo_ok, k_max, times_d, met_natural, met_wait, k_used, dt_best_d, t_arr_d):
    met_natural = |t_close - t*| <= K_w T; met_wait = some k in 0..k_max does; k_used = smallest such k (else the k nearest to t*); dt_best_d = signed t_arr - t* for that k."""
    W = K_w * T_final_d; apo_ok = bool(-tol_km <= err_km <= 0.0)
    drop = last_drop_km if (last_drop_km is not None and np.isfinite(last_drop_km)) else 0.0
    k_max = min(int(math.floor((err_km + tol_km) / drop + 1e-12)), k_cap) if (apo_ok and drop > 0) else 0
    times = [t_close_d + k * T_final_d for k in range(k_max + 1)]; dts = [t - t_star_d for t in times]
    ok = [k for k, d in enumerate(dts) if abs(d) <= W + 1e-12]; k_used = ok[0] if ok else int(np.argmin(np.abs(dts)))
    return dict(apo_ok=apo_ok, k_max=k_max, times_d=times, met_natural=abs(dts[0]) <= W + 1e-12, met_wait=bool(ok), k_used=k_used, dt_best_d=dts[k_used], t_arr_d=times[k_used])


def selftest():
    T = 0.142
    r = wait_rule(-5.0, 60.0, 10.0, T, 60.0 + 2 * T)
    assert r["k_max"] == 2 and r["met_wait"] and r["k_used"] == 0                      # 2T late is inside K_w = 3 already
    r = wait_rule(-5.0, 60.0, 10.0, T, 60.0 + 5 * T); assert not r["met_natural"] and r["met_wait"] and r["k_used"] == 2 and abs(r["dt_best_d"] + 3 * T) < 1e-12   # 5T away: k = 2 brings it to 3T
    r = wait_rule(-5.0, 60.0, 10.0, T, 60.0 + 5.1 * T); assert not r["met_wait"] and r["k_used"] == 2                  # still 3.1T away at k = 2
    r = wait_rule(-5.0, 60.0, 40.0, T, 60.0 + 5 * T); assert r["k_max"] == 0 and not r["met_wait"]                      # drop 40 > margin 20: no wait
    r = wait_rule(-30.0, 60.0, 5.0, T, 60.0 + 5 * T); assert r["k_max"] == 0 and not r["apo_ok"] and not r["met_wait"]   # already outside the tolerance
    r = wait_rule(-0.0, 60.0, 25.0, T, 60.0 + 4 * T); assert r["k_max"] == 1 and r["met_wait"]                          # margin 25 / drop 25: exactly one more orbit
    r = wait_rule(-5.0, 60.0, 10.0, T, 60.0 - 4 * T); assert not r["met_wait"]                                          # waiting never moves earlier
    r = wait_rule(-5.0, 60.0, 0.0, T, 60.0 + 9 * T); assert r["k_max"] == 0                                             # no positive drop: no wait
    r = wait_rule(-5.0, 60.0, 10.0, T, 60.0 + 3 * T, K_w=1.0); assert not r["met_natural"] and r["met_wait"] and r["k_used"] == 2
    print("planner_rev20d selftest ok: wait rule (k_max, window by natural and waited arrival, edge cases)")


def sha16(p): return __import__("hashlib").sha256(pathlib.Path(p).read_bytes()).hexdigest()[:16]
def md(head, rows): return ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)] + ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
def f2(x, p=2): return "n/a" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{p}f}"
def med(v):
    v = [x for x in v if x is not None and np.isfinite(x)]; return float(np.median(v)) if v else float("nan")
def jdump(o, name): OUT.mkdir(parents=True, exist_ok=True); (OUT / name).write_text(json.dumps(o, indent=1, default=lambda x: None if isinstance(x, float) and not math.isfinite(x) else float(x)))


# ---------------------------------------------------------------------------------------------- 2. baselines and flown arms (WP3 flights)

def base_flights():
    import planner_rev20_baselines as B
    fl = B.load_flights()
    for f in fl:
        fin = f["ra"][np.isfinite(f["ra"])]; f["last_drop_km"] = float((fin[-1] - f["ra_f"]) / 1e3) if len(fin) else float("nan"); f["t_c_d"] = f["t_c"] / 86400.0; f["T_d"] = f["T_f"] / 86400.0
    return fl, B


def score_flight(f, ts, K=K_W):
    w = wait_rule(f["err_km"], f["t_c_d"], f["last_drop_km"], f["T_d"], ts, K)
    s3 = bool(f["relaxed"] and w["met_wait"]); s3n = bool(f["relaxed"] and w["met_natural"])      # reached + no hard + window (apoapsis tolerance not required)
    s4 = bool(f["success"] and w["met_wait"]); s4n = bool(f["success"] and w["met_natural"])      # ... + final apoapsis in [-25, 0] km (WP3 mission success)
    return dict(w=w, s3=s3, s3n=s3n, s4=s4, s4n=s4n)


def cell_stats(fs, ts_fn, K=K_W):
    """fs: flights of one cell; ts_fn(f) -> t* (days) list; counts over (flight, t*) pairs, fractions and median dV of successes."""
    P = [(f, score_flight(f, t, K)) for f in fs for t in ts_fn(f)]; n = len(P)
    g = lambda key: sum(r[key] for _, r in P)
    d3 = [f["dv"] for f, r in P if r["s3"]]; d4 = [f["dv"] for f, r in P if r["s4"]]
    return dict(pairs=n, flights=len(fs), s3_nat=g("s3n"), s3=g("s3"), s4_nat=g("s4n"), s4=g("s4"), dv_s3=med(d3), dv_s4=med(d4), dv_all_reached=med([f["dv"] for f in fs if f["relaxed"]]))


def baselines():
    fl, B = base_flights(); ARMS = B.ARMS; LAB = B.LABEL; out = dict(cells={}, flights={})
    ic_nom = {g: json.load(open(TD / "rev20_data" / "ic" / f"curves_{x}.json"))["edges"]["ref"]["nominal_d"] for g, x in (("ODY", "ODY"), ("G1", "G1"), ("GC", "G0975"))}   # DERIVED (WP1, reference start)
    nom_sets = dict(ic_ref=ic_nom, rev19a=dict(B.NOMINAL_D))
    cells = [(g, a) for g in ARMS for a in ARMS[g] if any(f["geom"] == g and f["arm"] == a for f in fl)]
    tmed = {}; Tmed = {}
    for g, a in cells:
        S = [f for f in fl if f["geom"] == g and f["arm"] == a and f["success"]]; R = [f for f in fl if f["geom"] == g and f["arm"] == a and f["relaxed"]]
        tmed[(g, a)] = med([f["t_c_d"] for f in S]); Tmed[(g, a)] = med([f["T_d"] for f in S]) if S else med([f["T_d"] for f in R])
    kdist = {}
    for g, a in cells:
        fs = [f for f in fl if f["geom"] == g and f["arm"] == a]; key = f"{g}|{a}"
        t0 = tmed[(g, a)]
        out["cells"][key] = dict(flights=len(fs), reached_nohard=sum(f["relaxed"] for f in fs), strict=sum(f["success"] for f in fs), t_cell_median_d=t0, T_final_d=Tmed[(g, a)],
                                 at_median={f"K{int(K)}": cell_stats(fs, lambda f: [t0], K) for K in K_SWEEP},
                                 sweep={f"{o}": {f"K{int(K)}": cell_stats(fs, lambda f, o=o: ([t0] if o == 0 else [t0 - o * Tmed[(g, a)], t0 + o * Tmed[(g, a)]]), K) for K in (K_W,)} for o in OFFS},
                                 nominal={ns: {f"K{int(K)}": cell_stats(fs, lambda f: [nm[g]], K) for K in K_SWEEP} for ns, nm in nom_sets.items() if g in nm})
        km = [score_flight(f, t0)["w"]["k_max"] for f in fs if f["success"]]
        out["cells"][key]["kmax_of_strict"] = dict(n=len(km), median=med(km), p90=float(np.percentile(km, 90)) if km else float("nan"), max=max(km) if km else None, frac_ge1=float(np.mean([k >= 1 for k in km])) if km else float("nan"),
                                                   drop_median_km=med([f["last_drop_km"] for f in fs if f["success"]]))
    out["nominal_sets"] = nom_sets
    # flight-level window-by-wait changes at the registered t* (P20-10 evidence)
    jdump(out, "baselines.json")
    L = []
    hd = ["geometry", "arm", "flights", "reach+no hard", "+apo tol", "t* = cell median (d)", "S3 first / wait (K=3)", "S4 first / wait (K=3)", "S3 wait K=1/5/10", "S4 wait K=1/5/10", "median aero dV of S4 wait (m/s)", "median aero dV of S3 wait"]
    rows = []
    for key, c in out["cells"].items():
        g, a = key.split("|"); m = c["at_median"]; k3 = m["K3"]
        rows.append([g, LAB.get(a, a), c["flights"], c["reached_nohard"], c["strict"], f2(c["t_cell_median_d"], 1), f"{k3['s3_nat']} / {k3['s3']}", f"{k3['s4_nat']} / {k3['s4']}", "/".join(str(m[f"K{k}"]["s3"]) for k in (1, 5, 10)), "/".join(str(m[f"K{k}"]["s4"]) for k in (1, 5, 10)), f2(k3["dv_s4"]), f2(k3["dv_s3"])])
    L += ["### Table B1. t* = each cell's own median arrival (of the strict successes), aerobraking dV only", ""] + md(hd, rows) + [""]
    for ns in ("ic_ref", "rev19a"):
        rows = []
        for key, c in out["cells"].items():
            g, a = key.split("|")
            if ns not in c["nominal"]: continue
            m = c["nominal"][ns]; k3 = m["K3"]
            rows.append([g, LAB.get(a, a), c["flights"], f2(nom_sets[ns][g], 2), f"{k3['s3_nat']} / {k3['s3']}", f"{k3['s4_nat']} / {k3['s4']}", "/".join(str(m[f"K{k}"]["s3"]) for k in (1, 5, 10)), "/".join(str(m[f"K{k}"]["s4"]) for k in (1, 5, 10)), f2(k3["dv_s4"]), f2(k3["dv_s3"])])
        L += [f"### Table B2 ({ns}). t* fixed at the geometry's model-world nominal arrival ({'WP1 ic_plan reference start; GC uncalibrated (kappa not measured)' if ns == 'ic_ref' else 'Revision 19a values; no GC value'})", ""] + md(["geometry", "arm", "flights", "t* (d)", "S3 first / wait (K=3)", "S4 first / wait (K=3)", "S3 wait K=1/5/10", "S4 wait K=1/5/10", "median aero dV S4 wait", "median aero dV S3 wait"], rows) + [""]
    hd = ["geometry", "arm"] + [f"|dt|={o} orb" for o in OFFS]
    rows = []
    for key, c in out["cells"].items():
        g, a = key.split("|"); rows.append([g, LAB.get(a, a)] + [f"{c['sweep'][str(o)]['K3']['s4']}/{c['sweep'][str(o)]['K3']['pairs']} ({f2(c['sweep'][str(o)]['K3']['dv_s4'])})" for o in OFFS])
    L += ["### Table B3. S4 successes with the wait rule (K=3) / (flight x t*) pairs, and median aero dV of those successes (m/s), versus |t* - cell median| in final orbits (both signs pooled; the cell median is the strict-success median arrival)", ""] + md(hd, rows) + [""]
    rows = []
    for key, c in out["cells"].items():
        g, a = key.split("|"); rows.append([g, LAB.get(a, a)] + [f"{c['sweep'][str(o)]['K3']['s3']}/{c['sweep'][str(o)]['K3']['pairs']} ({f2(c['sweep'][str(o)]['K3']['dv_s3'])})" for o in OFFS])
    L += ["### Table B4. Same as B3 for S3 (reached + no hard + window; apoapsis tolerance not required)", ""] + md(hd, rows) + [""]
    rows = []
    for key, c in out["cells"].items():
        g, a = key.split("|"); k = c["kmax_of_strict"]; rows.append([g, LAB.get(a, a), k["n"], f2(k["drop_median_km"], 1), f2(k["median"], 0), f2(k["p90"], 0), k["max"], f2(k["frac_ge1"], 2)])
    L += ["### Table B5. Admissible extra final orbits of the strict successes (k_max; last flown drop, 25 km tolerance)", ""] + md(["geometry", "arm", "n", "median last drop (km)", "median k_max", "p90", "max", "fraction k_max >= 1"], rows) + [""]
    (OUT / "tables_baselines.md").write_text("\n".join(L) + "\n"); print("\n".join(L))


# ---------------------------------------------------------------------------------------------- 3. replay (guard v2) with the wait rule

def replay_rows(R, ids, geom, guard, ts, ic, noise):
    out = []
    for rid in ids:
        F = R.get_flight(rid); m = R.split_id(rid)
        r = R.replay(F, geom, "OLRT3", ic=ic, noise=noise, guard=guard, tstar_d=ts, h_des_m=R.H_DES_M, seed=m["seed"] * 100 + (m["level"] or 0), trace=True)
        row = dict(run_id=rid, reached=bool(r["reached"]), hard=bool(r["hard_violation"]), tstar_d=ts, dv_campaign=r["dv_campaign"], dv_pre=r.get("dv_pre", 0.0))
        if r["reached"]:
            tr = r["trace"]; ld = (tr[-2][3] - tr[-1][3]) / 1e3 if len(tr) >= 2 else float("nan")
            row.update(err_km=r["final_err_km"], t_arr_d=r["t_arr_d"], T_d=r["T_final_s"] / 86400.0, last_drop_km=ld, fb_unused=r["fb"], w_max_d=r["w_max_d"])
        out.append(row)
    return out


def score_replay(rows, K=K_W):
    S = []
    for r in rows:
        if not r["reached"]: S.append(dict(first=False, wait=False, k_max=0, apo=False)); continue
        w = wait_rule(r["err_km"], r["t_arr_d"], r["last_drop_km"], r["T_d"], r["tstar_d"], K)
        S.append(dict(first=w["met_natural"], wait=w["met_wait"], k_max=w["k_max"], apo=w["apo_ok"], k_used=w["k_used"]))
    return S


def cell_summary(rows):
    out = dict(n=len(rows), reached=sum(r["reached"] for r in rows), hard=sum(r["hard"] for r in rows), apo_ok=0, kmax_median=float("nan"), by_K={})
    S3 = score_replay(rows); out["apo_ok"] = sum(s["apo"] for s in S3); out["kmax_median"] = med([s["k_max"] for s in S3 if True]); out["frac_kmax_ge1"] = float(np.mean([s["k_max"] >= 1 for s in S3]))
    for K in K_SWEEP:
        S = score_replay(rows, K)
        ok_first = [r["reached"] and not r["hard"] and r["err_km"] >= -TOL_KM and s["first"] for r, s in zip(rows, S)]; ok_wait = [r["reached"] and not r["hard"] and s["apo"] and s["wait"] for r, s in zip(rows, S)]
        dv_f = [r["dv_campaign"] + r["dv_pre"] for r, ok in zip(rows, ok_first) if ok]; dv_w = [r["dv_campaign"] + r["dv_pre"] for r, ok in zip(rows, ok_wait) if ok]
        out["by_K"][f"K{int(K)}"] = dict(window_first=sum(s["first"] for s in S), window_wait=sum(s["wait"] for s in S), success_first=int(sum(ok_first)), success_wait=int(sum(ok_wait)), dv_success_first=med(dv_f), dv_success_wait=med(dv_w))
    out["med_campaign_dv"] = med([r["dv_campaign"] for r in rows if r["reached"]]); out["med_total_dv"] = med([r["dv_campaign"] + r["dv_pre"] for r in rows if r["reached"]])
    return out


def paired_extra(base, guard):
    d = [g["dv_campaign"] - b["dv_campaign"] for b, g in zip(base, guard) if b["reached"] and g["reached"]]
    return dict(n=len(d), median=med(d), q25=float(np.percentile(d, 25)) if d else float("nan"), q75=float(np.percentile(d, 75)) if d else float("nan"))


def replay_cmd():
    import planner_rev20_replay as R
    gd = {**R.GUARD_DEFAULT, **R.GUARD_V2}; res = dict(held={}, rev17={}, check={})
    maxdiff = 0.0
    for noise in ("pass", "time", "state"):
        for sc in ("ref", "S-fixed"):
            stored = {"pass": f"evaluate_held_ref+S-fixed_v2_pass.json", "time": "evaluate_held_ref+S-fixed_v2_time.json", "state": "evaluate_held_ref+S-fixed_v2_state.json"}[noise]
            SJ = json.load(open(R.OUT / stored))["cells"]
            for g in ("ODY", "G1"):
                ids = R.held_paths(g); ic = R.get_ic(g, sc); ts = R.tstar_for(g, ic)
                b = replay_rows(R, ids, g, None, ts, ic, noise); gg = replay_rows(R, ids, g, gd, ts, ic, noise)
                for tag, rows, key in (("baseline", b, "rows_baseline"), ("guard", gg, "rows_guard")):
                    for r, s in zip(rows, SJ[f"{g}|{sc}"][key]):
                        if r["reached"]: maxdiff = max(maxdiff, abs(r["t_arr_d"] - s["t_arr_d"]), abs(r["dv_campaign"] - s["dv_campaign"]))
                res["held"][f"{noise}|{g}|{sc}"] = dict(tstar_d=ts, baseline=cell_summary(b), guard=cell_summary(gg), paired=paired_extra(b, gg), rows_baseline=b, rows_guard=gg)
                print(noise, sc, g, "done", flush=True)
    res["check"]["held_max_abs_diff_vs_stored_t_arr_d_or_dv"] = maxdiff
    for noise in ("pass", "time", "state"):
        for sc in ("ref", "S-fixed"):
            ic = R.get_ic("ODY", sc); ts = R.tstar_for("ODY", ic); allb, allg = [], []
            cells = {}
            for lvl in range(5):
                for src in ("OL", "OLR", "OLRT3"):
                    ids = R.rev17_paths(lvl, src)
                    if not ids: continue
                    b = replay_rows(R, ids, "ODY", None, ts, ic, noise); gg = replay_rows(R, ids, "ODY", gd, ts, ic, noise); allb += b; allg += gg
                    cells[f"L{lvl}|{src}"] = dict(n=len(ids), baseline=cell_summary(b), guard=cell_summary(gg), paired=paired_extra(b, gg))
            res["rev17"][f"{noise}|{sc}"] = dict(tstar_d=ts, cells=cells, pooled=dict(n=len(allb), baseline=cell_summary(allb), guard=cell_summary(allg), paired=paired_extra(allb, allg)))
            print("rev17", noise, sc, len(allb), flush=True)
    jdump(res, "replay.json")


def replay_tables():
    J = json.load(open(OUT / "replay.json")); L = []
    rows = []
    for key, c in J["held"].items():
        noise, g, sc = key.split("|"); b, gd_, p = c["baseline"], c["guard"], c["paired"]
        rows.append([noise, g, sc, f"{b['by_K']['K3']['window_first']} / {b['by_K']['K3']['window_wait']}", f"{gd_['by_K']['K3']['window_first']} / {gd_['by_K']['K3']['window_wait']}", f"{gd_['by_K']['K3']['success_first']} / {gd_['by_K']['K3']['success_wait']}",
                     f"{b['hard']} / {gd_['hard']}", f2(gd_["by_K"]["K3"]["dv_success_wait"]), f2(b["by_K"]["K3"]["dv_success_wait"]), f2(p["median"]), f"{f2(gd_['kmax_median'], 0)} ({f2(gd_['frac_kmax_ge1'], 2)})"])
    L += ["### Table R1. Held paths (16 per cell): window first-passage / with wait rule (K=3), OLRT3 alone vs OLRT3 + guard v2; aerobraking dV (campaign + pre-burn, no FB)", ""] + md(
        ["indexing", "geom", "start", "window OLRT3 alone: first / wait", "window guard v2: first / wait", "success guard (reach+no hard+apo+window): first / wait", "hard alone / guard", "median aero dV of guard successes (wait)", "median aero dV of alone successes (wait)", "paired extra campaign dV (median)", "guard k_max median (frac >= 1)"], rows) + [""]
    rows = []
    for key, c in J["held"].items():
        noise, g, sc = key.split("|"); gd_ = c["guard"]["by_K"]; rows.append([noise, g, sc] + [f"{gd_[f'K{int(K)}']['window_first']} / {gd_[f'K{int(K)}']['window_wait']}" for K in K_SWEEP])
    L += ["### Table R2. Guard v2 window count first-passage / wait, by K_w (final orbits), held paths", ""] + md(["indexing", "geom", "start"] + [f"K={int(K)}" for K in K_SWEEP], rows) + [""]
    rows = []
    for key, c in J["rev17"].items():
        noise, sc = key.split("|"); p = c["pooled"]; b, gd_ = p["baseline"], p["guard"]
        rows.append([noise, sc, p["n"], f"{b['by_K']['K3']['window_first']} / {b['by_K']['K3']['window_wait']}", f"{gd_['by_K']['K3']['window_first']} / {gd_['by_K']['K3']['window_wait']}", f"{gd_['by_K']['K3']['success_first']} / {gd_['by_K']['K3']['success_wait']}", f"{b['hard']} / {gd_['hard']}", f2(gd_["by_K"]["K3"]["dv_success_wait"]), f2(p["paired"]["median"])])
    L += ["### Table R3. Rev 17 paths pooled (OL, OLR, OLRT3 source arms, L0-L4): window first-passage / wait (K=3), alone vs guard", ""] + md(["indexing", "start", "paths", "window alone: first / wait", "window guard: first / wait", "guard success: first / wait", "hard alone / guard", "median aero dV of guard successes (wait)", "paired extra campaign dV (median)"], rows) + [""]
    adds = [(k, c["guard"]["by_K"]["K3"]["window_wait"] - c["guard"]["by_K"]["K3"]["window_first"], c["baseline"]["by_K"]["K3"]["window_wait"] - c["baseline"]["by_K"]["K3"]["window_first"]) for k, c in J["held"].items()]
    L += [f"Maximum window gain from the wait rule over the 12 held cells (guard / alone): {max(a[1] for a in adds)} / {max(a[2] for a in adds)}. Max |held replay - stored| in t_arr_d or campaign dV: {J['check']['held_max_abs_diff_vs_stored_t_arr_d_or_dv']:.2e}.", ""]
    (OUT / "tables_replay.md").write_text("\n".join(L) + "\n"); print("\n".join(L))


# ---------------------------------------------------------------------------------------------- 4. WP1 curves: aerobraking dV (campaign + pre-burn) versus |dt| in final orbits

def curves():
    out = {}; L = []; XS = (0, 1, 3, 5, 10, 20, 40, 80)
    for g, tag in (("ODY", "ODY"), ("G1", "G1"), ("GC", "G0975")):
        J = json.load(open(TD / "rev20_data" / "ic" / f"curves_{tag}.json")); T = J["T_final_s"] / 86400.0; out[g] = {}
        for sc in ("ref", "S-fixed", "S-free"):
            nom = J["edges"][sc]["nominal_d"]
            rows = sorted([r for r in J["rows"] if r["scenario"] == sc and r["s_d"] == 0 and r["fb_term"] is False and r["feasible"]], key=lambda r: r["t_star_d"])
            x = np.array([(r["t_star_d"] - nom) / T for r in rows]); d = np.array([r["campaign_dv_mps"] + r["dv_pre_mps"] for r in rows])
            # the 'nominal' row (t* = nominal) is on the grid; D0(x) interpolates linearly between grid points (numerical, ASSUMED)
            D0 = lambda xx: float(np.interp(xx, x, d)); dnom = D0(0.0); res = {}
            xmin, xmax = float(x.min()), float(x.max())
            for K in K_SWEEP:
                row = {}
                for o in XS:
                    for sg in ((1,) if o == 0 else (-1, 1)):
                        xx = sg * o
                        if xx < xmin - 1e-9 or xx > xmax + 1e-9: row[f"{sg * o:+d}"] = None; continue
                        e = sg * max(o - K, 0.0); row[f"{sg * o:+d}"] = D0(e) - dnom if (xmin - 1e-9 <= e <= xmax + 1e-9) else None
                res[f"K{int(K)}"] = row
            # check: the stored s = 3-period rows (0.4236 d ~ 3 T) against the shifted s = 0 curve at the grid t*
            wr = sorted([r for r in J["rows"] if r["scenario"] == sc and abs(r["s_d"] - J["W_max_d"]) < 1e-9 and r["fb_term"] is False and r["feasible"]], key=lambda r: r["t_star_d"]); chk = []
            for r in wr:
                xx = (r["t_star_d"] - nom) / T
                e = math.copysign(max(abs(xx) - 3.0, 0.0), xx)
                if xmin <= e <= xmax: chk.append(abs((r["campaign_dv_mps"] + r["dv_pre_mps"]) - D0(e)))
            out[g][sc] = dict(nominal_d=nom, T_final_d=T, dv_nominal_mps=dnom, x_min_orbits=xmin, x_max_orbits=xmax, excess=res, check_K3_max_abs_diff_mps=max(chk) if chk else None, check_K3_median_abs_diff_mps=med(chk) if chk else None, n_grid=len(x))
    for K in K_SWEEP:
        rows = []
        for g in out:
            for sc, o_ in out[g].items():
                e = o_["excess"][f"K{int(K)}"]; rows.append([g, sc, f2(o_["nominal_d"], 2), f2(o_["dv_nominal_mps"], 3)] + [("-" if e.get(k) is None else f2(e[k], 3)) for k in ("-80", "-40", "-20", "-10", "-5", "-3", "-1", "+0", "+1", "+3", "+5", "+10", "+20", "+40", "+80")])
        L += [f"### Table C{int(K)}. Excess aerobraking dV (campaign + pre-burn) over the nominal-arrival plan (m/s) versus t* - nominal in final orbits (negative = earlier), K_w = {int(K)} orbits of slack; '-' = infeasible (earlier than the fastest feasible campaign) or beyond the computed range", ""] + md(
            ["geom", "scenario", "nominal (d)", "dV at nominal (m/s)", "-80", "-40", "-20", "-10", "-5", "-3", "-1", "0", "+1", "+3", "+5", "+10", "+20", "+40", "+80"], rows) + [""]
    rows = [[g, sc, f2(o_["x_min_orbits"], 1), f2(o_["x_max_orbits"], 1), o_["n_grid"], f2(o_["check_K3_max_abs_diff_mps"], 3), f2(o_["check_K3_median_abs_diff_mps"], 3)] for g in out for sc, o_ in out[g].items()]
    L += ["### Table C0. Range of the curves (orbits relative to nominal) and consistency check of the K=3 shifted curve against WP1's stored s = 3-period solves", ""] + md(["geom", "scenario", "earliest feasible t* (orbits)", "latest computed t* (orbits)", "grid points", "max |shifted - stored s=W| (m/s)", "median"], rows) + [""]
    jdump(out, "curves.json"); (OUT / "tables_curves.md").write_text("\n".join(L) + "\n"); print("\n".join(L))


# ---------------------------------------------------------------------------------------------- 5. AC5 corridor check

def flight_dirs():
    """(group, geom, arm, level, run_id, dir) for 366 indexed + trx_pull4 dev OLRT3 + rev17 + rev20olt (+ rev20val flights that have a summary row)."""
    import planner_rev18_flights as FL
    TOK = FL.TOK; seen = set(); out = []
    for r in json.load(open(TD / "dv_correction" / "flights.json")):
        if r.get("status") == "not_corrected": continue
        d = FL.run_dir(r["root"], r["phase"], r["run_id"])
        if d is None: continue
        out.append((r["phase"], TOK[str(r["geometry"])], r["arm"], "", r["run_id"], d)); seen.add(r["run_id"])
    for d in sorted((ROOT.parent / "trx_pull4" / "planner" / "devR").glob("*_PL_OLRT3_*")):
        if d.name not in seen and (d / "passes.csv").exists():
            gid = d.name.split("_PL_")[0]; out.append(("devR_trx", {"ODY": "ODY", "G5863": "G1"}[gid], "OLRT3", "", d.name, d)); seen.add(d.name)
    for blk, pat in (("rev17", "summary_rev17*.csv"), ("rev20olt", "summary_rev20olt*.csv"), ("rev20val", "summary_rev20val*.csv")):
        done = set()
        for f in glob.glob(str(ROOT / pat)):
            s = pd.read_csv(f); done |= set(s.run_id[s.status == "ok"])
        for d in sorted((ROOT / blk).glob("*_PL_*")):
            if d.name in done and (d / "passes.csv").exists():
                gid, rest = d.name.split("_PL_"); arm = rest.split("_pert_")[0]; lv = d.name.split("_")[-1] if d.name.split("_")[-1].startswith("L") else ""
                out.append((blk, {"ODY": "ODY", "G5863": "G1"}.get(gid, gid), arm, lv, d.name, d))
    return out


def corridor():
    fl = flight_dirs(); rows = []
    for ph, g, arm, lv, rid, d in fl:
        P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra = P.ra_m.values.astype(float); hit = np.where(np.isfinite(ra) & (ra <= _AN.RA_TARGET))[0]
        reached = len(hit) > 0; i = int(hit[0]) if reached else len(P) - 1
        qs = 0.5 * P.rho_peak_kg_m3.values * P.vrel_mps.values ** 3 / 1e4; q = float(P.q_w_cm2.iloc[i])
        rows.append(dict(group=ph, geom=g, arm=arm, level=lv, run_id=rid, reached=reached, q_last=q, qs_last=float(qs[i]), qs_max_any=float(np.nanmax(qs[:i + 1])), alt_km=float(P.alt_m.iloc[i]) / 1e3,
                         in_band=bool(Q_BAND[0] <= q <= Q_BAND[1]), below=bool(q < Q_BAND[0]), above=bool(q > Q_BAND[1]), qs_over=bool(qs[i] > QS_HARD)))
    df = pd.DataFrame(rows); df.to_csv(OUT / "corridor_rows.csv", index=False)
    R = df[df.reached]; g = R.groupby(["group", "geom", "arm", "level"])
    T = g.agg(n=("run_id", "size"), in_band=("in_band", "sum"), below=("below", "sum"), above=("above", "sum"), qs_over=("qs_over", "sum"), q_med=("q_last", "median"), q_min=("q_last", "min"), q_max=("q_last", "max"), qs_med=("qs_last", "median"), qs_max=("qs_last", "max"),
              alt_med=("alt_km", "median"), alt_min=("alt_km", "min"), alt_max=("alt_km", "max")).reset_index()
    nr = df[~df.reached].groupby(["group", "geom", "arm", "level"]).size().rename("not_reached").reset_index(); T = T.merge(nr, how="left", on=["group", "geom", "arm", "level"]).fillna({"not_reached": 0})
    tot = dict(flights=len(df), reached=int(len(R)), not_reached=int((~df.reached).sum()), in_band=int(R.in_band.sum()), below=int(R.below.sum()), above=int(R.above.sum()), qs_over=int(R.qs_over.sum()), q_med=float(R.q_last.median()), q_p05=float(R.q_last.quantile(.05)), q_p95=float(R.q_last.quantile(.95)),
               alt_med=float(R.alt_km.median()), alt_min=float(R.alt_km.min()), alt_max=float(R.alt_km.max()), by_group={k: dict(n=int(len(v)), in_band=int(v.in_band.sum()), below=int(v.below.sum()), above=int(v.above.sum()), qs_over=int(v.qs_over.sum())) for k, v in R.groupby("group")},
               by_geom={k: dict(n=int(len(v)), in_band=int(v.in_band.sum()), below=int(v.below.sum()), above=int(v.above.sum()), qs_over=int(v.qs_over.sum()), q_med=float(v.q_last.median()), alt_med=float(v.alt_km.median())) for k, v in R.groupby("geom")})
    jdump(dict(total=tot, cells=T.to_dict("records")), "corridor.json")
    L = ["### Table Q0. Totals", "", f"Flights {tot['flights']} (reached {tot['reached']}, not reached {tot['not_reached']}); last-pass mean q in [0.05, 0.25]: {tot['in_band']}, below {tot['below']}, above {tot['above']}; last-pass Q_s > 0.65: {tot['qs_over']}; median q {tot['q_med']:.3f} (5-95%: {tot['q_p05']:.3f}-{tot['q_p95']:.3f}); final periapsis altitude median {tot['alt_med']:.1f} km (range {tot['alt_min']:.1f}-{tot['alt_max']:.1f}).", ""]
    rows = [[k, v["n"], v["in_band"], v["below"], v["above"], v["qs_over"]] for k, v in tot["by_group"].items()] + [[f"geom {k}", v["n"], v["in_band"], v["below"], v["above"], v["qs_over"]] for k, v in tot["by_geom"].items()]
    L += ["### Table Q1. By group / geometry (reached flights)", ""] + md(["group", "reached", "q in band", "q below 0.05", "q above 0.25", "last-pass Q_s > 0.65"], rows) + [""]
    rows = [[r.group, r.geom, r.arm, r.level or "-", r.n, int(r.not_reached), int(r.in_band), int(r.below), int(r.above), int(r.qs_over), f"{r.q_med:.3f} [{r.q_min:.3f}, {r.q_max:.3f}]", f"{r.qs_med:.3f} (max {r.qs_max:.3f})", f"{r.alt_med:.1f} [{r.alt_min:.1f}, {r.alt_max:.1f}]"] for r in T.itertuples()]
    L += ["### Table Q2. Per cell (group, geometry, arm, level)", ""] + md(["group", "geom", "arm", "level", "reached", "not reached", "in band", "below", "above", "Qs>0.65", "last-pass q: median [min, max] (W/cm2)", "last-pass Q_s: median (max)", "final periapsis alt (km): median [min, max]"], rows) + [""]
    (OUT / "tables_corridor.md").write_text("\n".join(L) + "\n"); print("\n".join(L))


# ---------------------------------------------------------------------------------------------- P20-10 evidence

def p20_10():
    """Largest window gain of the wait rule over a scan of t* (flights: cell median +- 200 orbits; replay: nominal +- 100 orbits), per cell, and the k_max distribution. Written to p20_10_scan.json."""
    fl, B = base_flights(); out = dict(flights=[], replay={})
    for g in B.ARMS:
        for a in B.ARMS[g]:
            fs = [f for f in fl if f["geom"] == g and f["arm"] == a]
            if not fs: continue
            S = [f for f in fs if f["success"]]; t0 = float(np.median([f["t_c_d"] for f in S])); T = float(np.median([f["T_d"] for f in S])); mx3 = mx4 = 0; arg = None
            for o in range(-200, 201):
                r = [score_flight(f, t0 + o * T) for f in fs]; d3 = sum(x["s3"] for x in r) - sum(x["s3n"] for x in r); d4 = sum(x["s4"] for x in r) - sum(x["s4n"] for x in r)
                if d3 > mx3 or d4 > mx4: arg = o
                mx3 = max(mx3, d3); mx4 = max(mx4, d4)
            out["flights"].append(dict(geom=g, arm=a, n=len(fs), max_add_s3=mx3, max_add_s4=mx4, at_offset_orbits=arg))
    J = json.load(open(OUT / "replay.json")); ks = []
    for k, c in J["held"].items():
        mx = {}
        for tag in ("rows_baseline", "rows_guard"):
            rows = [r for r in c[tag] if r["reached"]]; m = 0
            for r in rows: ks.append(wait_rule(r["err_km"], r["t_arr_d"], r["last_drop_km"], r["T_d"], r["tstar_d"])["k_max"])
            for o in range(-100, 101):
                w = [wait_rule(r["err_km"], r["t_arr_d"], r["last_drop_km"], r["T_d"], r["tstar_d"] + o * r["T_d"]) for r in rows]
                m = max(m, sum(x["met_wait"] for x in w) - sum(x["met_natural"] and x["apo_ok"] for x in w))
            mx[tag] = m
        out["replay"][k] = mx
    out["replay_kmax"] = dict(n=len(ks), median=float(np.median(ks)), max=int(max(ks)), frac_ge1=float(np.mean(np.array(ks) >= 1)))
    jdump(out, "p20_10_scan.json"); print(json.dumps(out, indent=1))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["selftest", "baselines", "replay", "replay_tables", "curves", "corridor", "p20_10", "all"]); a = ap.parse_args(); OUT.mkdir(parents=True, exist_ok=True)
    if a.cmd == "selftest": return selftest()
    if a.cmd in ("baselines", "all"): baselines()
    if a.cmd in ("replay", "all"): replay_cmd()
    if a.cmd in ("replay_tables", "all"): replay_tables()
    if a.cmd in ("curves", "all"): curves()
    if a.cmd in ("corridor", "all"): corridor()
    if a.cmd in ("p20_10", "all"): p20_10()


if __name__ == "__main__":
    main()
