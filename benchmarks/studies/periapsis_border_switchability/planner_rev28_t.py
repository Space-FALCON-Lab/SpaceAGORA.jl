"""Revision 28A (rev28t): default t* (OLRT3-W3T, fast end + 6 d) against the nominal t* (OLRT3-W3E), both EI 150, aim 12.5 km. Analysis only; Python, read-only on flight data.
Per flight (planner_rev26_val2.load): success, hard, window against the flight's OWN t* (W3T: summary pos_tstar_d; W3E: logged aim minus the offset), corridor, aerobraking dV, arrival minus t*, position raw and with L1, mission success.
Pairs W3T minus W3E on the same seed: W3E from rev28t (G0975M s13-16; ODY and G1230 s17-20) or the earlier blocks (rev27ei ODY s13-16, rev27low G1038 and G1230 s13-16).
Statistics as results_rev27_sensitivity.md: median of the paired differences, bootstrap 95% CI (10,000 resamples, seed 20261030), exact two-sided sign test, verdict by the sensitivity rule (one-pair scale 1.6 m/s).
Usage: python3 -I planner_rev28_t.py --out-dir DIR"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev26_val2 as V
import planner_rev20_score as SC
import planner_rev21_position as PS
NB = 10000; SEED = 20261030; S_PAIR = 1.6


def sign_p(d):
    d = [x for x in d if x != 0]; n = len(d); k = sum(1 for x in d if x > 0)
    if n == 0: return 1.0
    kk = min(k, n - k); return min(1.0, 2 * sum(math.comb(n, i) for i in range(kk + 1)) / 2 ** n)


def paired(d):
    d = np.asarray(d, float); rng = np.random.default_rng(SEED); bs = np.median(rng.choice(d, size=(NB, len(d)), replace=True), axis=1); lo, hi = np.percentile(bs, [2.5, 97.5]); p = sign_p(d)
    sgn = np.sign(np.median(d)); agree = float(np.mean(np.sign(d) == sgn)); excl = lo > 0 or hi < 0
    strong = excl and (p < 0.05 or (len(d) <= 4 and agree == 1.0 and (lo if sgn > 0 else -hi) > S_PAIR))
    return dict(n=len(d), med=float(np.median(d)), mean=float(np.mean(d)), ci=[float(lo), float(hi)], p=p, pos=int((d > 0).sum()), neg=int((d < 0).sum()), verdict="still holds" if strong else ("holds only in sign" if agree >= 2 / 3 else "not resolved"), d=[float(x) for x in d])


def nominal_for_tstar(tok, tstar_d, ic=None):
    """Nominal plan for a default t*: the table's zero-burn level path shifted by a constant c in ln q so that the model arrival (kappa of the flights, ASSUMED same as the guard's) equals t*; its (a, e) history drives the J2 track of the target direction (DERIVED, ASSUMED construction).
    The logged position columns instead evaluate the UNSHIFTED nominal plan (arrival at the nominal t*) at the shifted t*."""
    RP = PS._rp(); tab = PS.get_tab(tok); x0 = math.log(RP.QFIRST_CAL.get(tok, PS.O.Q_FIRST)) + RP.ic_shift(PS.headline_ic(tok)); path = np.array(tab["path"], float); rt = PS.target_radius(tok); kap = PS.kappa_for(tok)
    f = lambda c: PS.arrival_history(tab, x0 + path - path[0] + c, x0, kap, rt)["days"] - tstar_d
    lo, hi = -1.5, 1.5
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0: lo = mid      # arrival too late -> more heating
        else: hi = mid
    c = 0.5 * (lo + hi); h = PS.arrival_history(tab, x0 + path - path[0] + c, x0, kap, rt); i0, om0, Om0 = PS.geometry_elements(tok)
    tr = PS.J2Track(h["t"], h["a"], h["e"], np.full(len(h["t"]), i0), om0, Om0, i0)
    return dict(track=tr, days=h["days"], hist=h, tok=tok, N=h["N"], shift=c)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    B = {b: V.load(b, None) for b in ("rev28t", "rev27ei", "rev27low")}; rows = []
    for b, F in B.items():
        for k, f in F.items():
            p = f["pos"]; rows.append(dict(block=b, arm=k[0], geom=k[1], level=k[2], seed=k[3], ok=f["ok"], hard=f["hard"], window=f["window"], corr=f["corr"], err=f["err"], dv=f["dv_sum"], burns=f["burns"], t_arr=f["t_arr_d"], tstar=f["tstar"], dt_d=f["t_arr_d"] - f["tstar"], dt_w=f["dt_nom_w"], wmax=f["wmax_d"],
                                           raw=None if not p else p["raw"], l1=None if not p else p["l1"], all=f["all"], stop=f["stop"], log_pos=float(f["summ"].get("pos_err_deg", float("nan"))), log_l1=float(f["summ"].get("l1_err_deg", float("nan"))), fast=float(f["summ"].get("tstar_fast_d", float("nan"))), slow=float(f["summ"].get("tstar_slow_d", float("nan")))))
    for r_ in rows:       # position against a nominal rebuilt for the flight's own default t* (W3T only)
        if r_["arm"] == "OLRT3W3T" and r_["raw"] is not None:
            f_ = B[r_["block"]][(r_["arm"], r_["geom"], r_["level"], r_["seed"])]; nm = nominal_for_tstar(f_["tok"], r_["tstar"]); q_ = V.position(f_["dir"], f_["tok"], r_["tstar"], f_["rt"], None, nm)
            r_["raw_def"] = q_["raw"]; r_["l1_def"] = q_["l1"]; r_["shift_lnq"] = nm["shift"]; r_["all_def"] = bool(f_["ok"] and not f_["hard"] and f_["window"] and f_["corr"] and q_["l1"] < 1.0)
        else: r_["raw_def"] = r_["raw"]; r_["l1_def"] = r_["l1"]; r_["shift_lnq"] = None; r_["all_def"] = r_["all"]
    df = pd.DataFrame(rows); (out / "rev28t_rows.json").write_text(json.dumps(SC.nan_none(rows), indent=1))
    v = df[(df.block == "rev28t") & df.log_pos.notna() & df.raw.notna()]
    print("cross-check W3T (mine vs logged): n", len(v), "max |raw|", float((v.raw - v.log_pos).abs().max()), "max |L1|", float((v.l1 - v.log_l1).abs().max()))
    # pairs
    def pick(geom, seed, arm_w3e):
        for blk in ("rev28t", "rev27ei", "rev27low"):
            x = df[(df.block == blk) & (df.arm == "OLRT3W3E") & (df.geom == geom) & (df.seed == seed) & (df.level == "L0")]
            if len(x): return x.iloc[0]
        return None
    P = {}
    for geom in ("G0975M", "G1038", "ODY", "G1230"):
        pr = []
        for s in range(13, 21):
            t = df[(df.block == "rev28t") & (df.arm == "OLRT3W3T") & (df.geom == geom) & (df.seed == s) & (df.level == "L0")]
            if not len(t): continue
            e = pick(geom, s, None)
            if e is None: continue
            pr.append((int(s), t.iloc[0], e))
        P[geom] = pr
    R = {}
    for geom, pr in P.items():
        d = [e.dv - t.dv for s, t, e in pr]       # saving = W3E minus W3T (positive = default t* cheaper)
        R[geom] = dict(seeds=[s for s, _, _ in pr], src=[e.block for _, _, e in pr], saving=paired(d), dv_T=[float(t.dv) for _, t, _ in pr], dv_E=[float(e.dv) for _, _, e in pr], diff_of_medians=float(np.median([e.dv for _, _, e in pr]) - np.median([t.dv for _, t, _ in pr])),
                       T=dict(ok=int(sum(t.ok for _, t, _ in pr)), hard=int(sum(t.hard for _, t, _ in pr)), win=int(sum(t.window for _, t, _ in pr)), corr=int(sum(t['corr'] for _, t, _ in pr)), l1=int(sum((t.l1 is not None and t.l1 < 1) for _, t, _ in pr)), allr=int(sum(t['all'] for _, t, _ in pr)), n=len(pr)),
                       E=dict(ok=int(sum(e.ok for _, _, e in pr)), hard=int(sum(e.hard for _, _, e in pr)), win=int(sum(e.window for _, _, e in pr)), corr=int(sum(e['corr'] for _, _, e in pr)), l1=int(sum((e.l1 is not None and e.l1 < 1) for _, _, e in pr)), allr=int(sum(e['all'] for _, _, e in pr)), n=len(pr)))
    # whole-arm cells (all W3T seeds, all W3E seeds)
    cells = {}
    for geom in ("G0975M", "G1038", "ODY", "G1230"):
        for arm in ("OLRT3W3T", "OLRT3W3E"):
            x = df[(df.arm == arm) & (df.geom == geom) & (df.level == "L0")]
            if not len(x): continue
            cells[f"{geom} {arm}"] = dict(n=len(x), blocks=sorted(set(x.block)), ok=int(x.ok.sum()), hard=int(x.hard.sum()), win=int(x.window.sum()), corr=int(x['corr'].sum()), dv_med=float(x.dv.median()), dv=[float(z) for z in x.dv], dt_d_med=float(x.dt_d.median()), dt_d=[float(z) for z in x.dt_d], dt_w_med=float(x.dt_w.median()), tstar_med=float(x.tstar.median()), arr_med=float(x.t_arr.median()),
                                         raw_med=float(x.raw.median()), l1_med=float(x.l1.median()), l1_lt1=int((x.l1 < 1).sum()), raw_lt1=int((x.raw < 1).sum()), allr=int(x['all'].sum()), fast=float(x.fast.median()), burns_med=float(x.burns.median()), err=[float(z) for z in x.err])
    # natural (unguarded) arrivals for the margin discussion
    S = SC.load_summary(V.ROOT, "summary_rev21c_*.csv"); nat = {}
    for blk, arm, geom in (("rev21c", "OLRT3", "ODY"), ("rev17", "OLRT3", "ODY")):
        Sx = SC.load_summary(V.ROOT, f"summary_{blk}_*.csv"); F = V.load(blk, None); x = [f["t_arr_d"] for k, f in F.items() if k[0] == "OLRT3" and k[1] == geom and f["reached"] and (k[2] in ("L0",))]
        nat[f"{blk} {geom} L0"] = x
    # regression
    c = SC.compare_passes(V.ROOT / "rev28t" / "ODY_PL_OLRT3_pert_ODY_s14_L4", V.ROOT / "rev17" / "ODY_PL_OLRT3_pert_ODY_s14_L4"); reg = dict(equal=c["equal"], ncols=c["ncols"], rows=(c["rows_a"], c["rows_b"]), first=c["first"])
    (out / "rev28t_summary.json").write_text(json.dumps(SC.nan_none(dict(pairs=R, cells=cells, natural=nat, regression=reg)), indent=1))
    for g, r in R.items(): print(g, r["saving"]["med"], r["saving"]["ci"], r["saving"]["p"], r["saving"]["verdict"], "T", r["T"], "E", r["E"])
    for k, c_ in cells.items(): print(k, {a: (round(b, 2) if isinstance(b, float) else b) for a, b in c_.items() if a not in ("dv", "dt_d", "err")})
    print(nat, reg)


if __name__ == "__main__":
    main()
