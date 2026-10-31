"""Revision 31.4 (rev31g1230): OLRT3-W3E7 (v3, 7.5 km aim, EI 150) at G1230 L0 seeds 13-20 against rev30seed W3E5 (5 km, same seeds) and the 12.5 km flights (rev27low s13-16, rev28t W3E s17-20). Analysis only; Python, read-only on flight data.
Last-pass q against 0.05, low-side misses (< -25 km), final error, window (own t*), position raw / L1, ALL with the corridor not counted / counted, paired dV (W3E7 minus W3E5). Usage: python3 -I planner_rev31_g1230.py --out-dir DIR"""
import argparse, json, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev26_val2 as V
import planner_rev28_t as T28
import planner_rev30_seed as S30
import planner_rev31_pos as P31
import planner_rev20_score as SC


def sumry(fs):
    s = P31.stats(fs); q = [float(f["q_last"]) for f in fs]; e = [float(f["err"]) for f in fs]
    s.update(q=q, q_below=int(sum(x < 0.05 for x in q)), q_above=int(sum(x > 0.25 for x in q)), q_min=min(q), q_max=max(q), err=e, err_med=float(np.median(e)), err_worst=min(e), low=int(sum(x < -25.0 for x in e)),
             all_corr=int(sum(bool(f["ok"] and not f["hard"] and f["window"] and f["pos"]["l1"] < 1.0 and f["corr"]) for f in fs)), dt=[float(f["dt_nom_w"]) for f in fs], passes=[int(f["n"]) for f in fs]); s["wilson_corr"] = S30.wilson(s["all_corr"], s["n"]); return s


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    B = {b: V.load(b, None) for b in ("rev31g1230", "rev30seed", "rev27low", "rev28t")}
    g = lambda blk, arm, seeds: [f for k, f in sorted(B[blk].items(), key=lambda t: t[0][3]) if k[0] == arm and k[1] == "G1230" and k[2] == "L0" and k[3] in seeds]
    f7 = g("rev31g1230", "OLRT3W3E7", range(13, 21)); f5 = g("rev30seed", "OLRT3W3E5", range(13, 21)); f125 = g("rev27low", "OLRT3W3E", range(13, 17)) + g("rev28t", "OLRT3W3E", range(17, 21))
    R = dict(w7=sumry(f7), w5=sumry(f5), w125=sumry(f125)); R["dv_pair_7_5"] = T28.paired([x["dv_sum"] - y["dv_sum"] for x, y in zip(f7, f5)]); R["dv_pair_7_125"] = T28.paired([x["dv_sum"] - y["dv_sum"] for x, y in zip(f7, f125)])
    c = SC.compare_passes(V.ROOT / "rev31g1230" / "ODY_PL_OLRT3_pert_ODY_s14_L4", V.ROOT / "rev17" / "ODY_PL_OLRT3_pert_ODY_s14_L4"); reg = dict(equal=c["equal"], ncols=c["ncols"], rows=(c["rows_a"], c["rows_b"]), first=c["first"])
    S = pd.read_csv(next(V.ROOT.glob("summary_rev31g1230_*.csv"))); rr = S[S.arm == "OLRT3"].iloc[0]; reg["ei"] = float(rr.ei_km); reg["stop"] = str(rr.stop_reason)
    (out / "rev31g1230.json").write_text(json.dumps(SC.nan_none(dict(cells=R, regression=reg)), indent=1, default=float))
    for k in ("w7", "w5", "w125"):
        s = R[k]; print(k, f"n {s['n']} ok {s['ok']} low {s['low']} err med {s['err_med']:.1f} worst {s['err_worst']:.1f} win {s['win']} q<0.05 {s['q_below']} qrange {s['q_min']:.3f}-{s['q_max']:.3f} raw {s['raw_med']:.2f} l1 {s['l1_med']:.2f} l11 {s['l11']} ALLnc {s['all']} ALLc {s['all_corr']} dV {s['dv_med']:.2f}")
    for k in ("dv_pair_7_5", "dv_pair_7_125"): p = R[k]; print(k, np.round(p["d"], 2).tolist(), round(p["med"], 2), np.round(p["ci"], 2).tolist(), p["verdict"], p["p"])
    print("w7 q", np.round(R["w7"]["q"], 3).tolist(), "l1", np.round(R["w7"]["l1"], 2).tolist(), "dt", np.round(R["w7"]["dt"], 2).tolist(), "err", np.round(R["w7"]["err"], 1).tolist()); print(reg)


if __name__ == "__main__":
    main()
