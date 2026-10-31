"""Revision 29 (rev29aim): 5 km terminal aim (OLRT3-W3E5) against 12.5 km (OLRT3-W3E), EI 150, nominal t*. Analysis only; Python, read-only on flight data.
Pairs by seed: Odyssey W3E5 s13-16 against rev27ei W3E, s17-20 against rev28t W3E; G1038 W3E5 s13-16 against rev27low W3E, s17-20 against rev29aim W3E (12.5 km).
Per flight: success, low-side miss (final error < -25 km), final error, window (own t*), hard, passes, terminal-approach passes (passes from the first apoapsis within 40 km of the target to the stop), aerobraking dV, position raw and with L1, mission success.
Statistics as results_rev27_sensitivity.md (paired median, bootstrap 95% CI of the median with 10,000 resamples seed 20261030, exact two-sided sign test, 1.6 m/s sensitivity rule). Usage: python3 -I planner_rev29_aim.py --out-dir DIR"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev26_val2 as V
import planner_rev28_t as T28
import planner_rev20_score as SC
import planner_rev21_position as PS


def term_passes(d, rt):
    P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra = P.ra_m.values.astype(float); i = np.where(np.isfinite(ra) & (ra - rt <= 40e3))[0]
    return int(len(P) - i[0]) if len(i) else 0


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    B = {b: V.load(b, None) for b in ("rev29aim", "rev27ei", "rev28t", "rev27low")}; rows = []
    for b, F in B.items():
        for k, f in F.items():
            if k[2] != "L0" or k[0] not in ("OLRT3W3E5", "OLRT3W3E"): continue
            if b == "rev28t" and k[1] != "ODY": continue
            p = f["pos"]; rows.append(dict(block=b, arm=k[0], geom=k[1], seed=k[3], ok=f["ok"], hard=f["hard"], window=f["window"], corr=f["corr"], err=f["err"], dv=f["dv_sum"], passes=int(f["passes"]), term=term_passes(f["dir"], f["rt"]), burns=f["burns"], dt_w=f["dt_nom_w"], raw=None if not p else p["raw"], l1=None if not p else p["l1"], all=f["all"], low=bool(f["err"] < -25.0)))
    df = pd.DataFrame(rows); (out / "rev29aim_rows.json").write_text(json.dumps(SC.nan_none(rows), indent=1))
    R = {}
    for geom in ("ODY", "G1038"):
        n5 = df[(df.arm == "OLRT3W3E5") & (df.geom == geom)].sort_values("seed")
        pr = []
        for _, r in n5.iterrows():
            c = df[(df.arm == "OLRT3W3E") & (df.geom == geom) & (df.seed == r.seed)]
            c = c[c.block.isin(("rev27ei", "rev28t") if geom == "ODY" else ("rev27low", "rev29aim"))]
            if len(c): pr.append((r, c.iloc[0]))
        d = [r.dv - c.dv for r, c in pr]; P5 = T28.paired(d)
        e5 = n5.err.values; ec = np.array([c.err for _, c in pr])
        R[geom] = dict(n5=len(n5), n_pairs=len(pr), dv=P5, src=[c.block for _, c in pr], dp=[int(r.passes - c.passes) for r, c in pr], dterm=[int(r.term - c.term) for r, c in pr],
                       W5=dict(low=int(n5.low.sum()), err_med=float(np.median(e5)), err_worst=float(e5.min()), err_best=float(e5.max()), win=int(n5.window.sum()), ok=int(n5.ok.sum()), hard=int(n5.hard.sum()), corr=int(n5['corr'].sum()), l1=int((n5.l1 < 1).sum()), all=int(n5['all'].sum()), passes_med=float(n5.passes.median()), term_med=float(n5.term.median()), dv_med=float(n5.dv.median()), l1_med=float(n5.l1.median()), raw_med=float(n5.raw.median())),
                       WE=dict(n=len(pr), low=int(sum(c.low for _, c in pr)), err_med=float(np.median(ec)), err_worst=float(ec.min()), win=int(sum(c.window for _, c in pr)), ok=int(sum(c.ok for _, c in pr)), hard=int(sum(c.hard for _, c in pr)), l1=int(sum(c.l1 < 1 for _, c in pr)), all=int(sum(c['all'] for _, c in pr)), passes_med=float(np.median([c.passes for _, c in pr])), term_med=float(np.median([c.term for _, c in pr])), dv_med=float(np.median([c.dv for _, c in pr]))),
                       W5_pairs=dict(low=int(sum(r.low for r, _ in pr)), err_med=float(np.median([r.err for r, _ in pr])), win=int(sum(r.window for r, _ in pr))))
    c = SC.compare_passes(V.ROOT / "rev29aim" / "ODY_PL_OLRT3_pert_ODY_s14_L4", V.ROOT / "rev17" / "ODY_PL_OLRT3_pert_ODY_s14_L4"); reg = dict(equal=c["equal"], ncols=c["ncols"], rows=(c["rows_a"], c["rows_b"]), first=c["first"])
    (out / "rev29aim_summary.json").write_text(json.dumps(SC.nan_none(dict(cells=R, regression=reg)), indent=1))
    for g, r in R.items(): print(g, "dV pair med %.2f CI %s p %.3f %s" % (r["dv"]["med"], np.round(r["dv"]["ci"], 2), r["dv"]["p"], r["dv"]["verdict"]), "d:", np.round(r["dv"]["d"], 2).tolist(), "\n  W5", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r["W5"].items()}, "\n  WE", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r["WE"].items()}, "\n  dp", r["dp"], "dterm", r["dterm"])
    print(reg)
    print(df[(df.arm == "OLRT3W3E5")][["geom", "seed", "ok", "err", "window", "dv", "passes", "term", "l1", "all"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
