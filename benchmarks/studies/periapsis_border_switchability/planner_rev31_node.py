"""Revision 31.2 (rev31node): OLRT3-W4LN5 (v4 + L1 + node check v2 at 1.0 deg, window-clamped) against OLRT3-W4L5 at G1 and G0975M, L0, seeds 13-20 (W4L5 s13-16 from rev31pos, s17-20 from rev31node). Analysis only; Python, read-only on flight data.
Per geometry: within 1 deg with L1, cross-track (node) and along-track errors, windows, paired dV, ALL; engagement of the node check: first differing decision pass against the paired W4L5 flight, differing ttg_engaged / ttg_solves counts, logged ttg_tt_s against the nominal t*.
Usage: python3 -I planner_rev31_node.py --out-dir DIR"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev26_val2 as V
import planner_rev28_t as T28
import planner_rev31_pos as P31
import planner_rev30_seed as S30
import planner_rev20_score as SC


def first_diff(a, b, f="decisions", col="delta_v_mps"):
    A = pd.read_csv(a / f"{f}.csv"); B = pd.read_csv(b / f"{f}.csv"); n = min(len(A), len(B)); x = A[col].values[:n].astype(float); y = B[col].values[:n].astype(float); i = np.where(~((x == y) | (np.isnan(x) & np.isnan(y))))[0]
    return (int(A.target_pass.values[i[0]]) if f == "decisions" and len(i) else (int(i[0]) + 1 if len(i) else None)), len(A), len(B)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    B = {b: V.load(b, None) for b in ("rev31node", "rev31pos")}; R = {}
    for geom in ("G1", "G0975M"):
        N = {k[3]: f for k, f in B["rev31node"].items() if k[0] == "OLRT3W4LN5" and k[1] == geom and k[2] == "L0"}
        L = {k[3]: f for k, f in B["rev31pos"].items() if k[0] == "OLRT3W4L5" and k[1] == geom and k[2] == "L0"}; L.update({k[3]: f for k, f in B["rev31node"].items() if k[0] == "OLRT3W4L5" and k[1] == geom and k[2] == "L0"})
        seeds = sorted(set(N) & set(L)); fn = [N[s] for s in seeds]; fl = [L[s] for s in seeds]; sn = P31.stats(fn); sl = P31.stats(fl); d = [x["dv_sum"] - y["dv_sum"] for x, y in zip(fn, fl)]
        eng = []
        for s in seeds:
            fd, na, nb = first_diff(N[s]["dir"], L[s]["dir"]); rn = pd.read_csv(N[s]["dir"] / "replan.csv"); rl = pd.read_csv(L[s]["dir"] / "replan.csv")
            tt = rn.ttg_tt_s.dropna(); eng.append(dict(seed=s, first_dv_diff_pass=fd, identical=bool(fd is None and na == nb), ttg_tt_min_d=float(tt.min() / 86400), ttg_tt_max_d=float(tt.max() / 86400), solves_N=int(rn.ttg_solves.fillna(0).sum()), solves_L=int(rl.ttg_solves.fillna(0).sum()), engaged_N=int(rn.ttg_engaged.fillna(0).sum()), engaged_L=int(rl.ttg_engaged.fillna(0).sum()),
                                      acted_N=float(V.SC.pd.read_csv(next(V.ROOT.glob("summary_rev31node_*.csv"))).set_index("run_id").loc[N[s]["dir"].name].pos_stage_acted), win_N=bool(N[s]["window"]), win_L=bool(L[s]["window"]), l1_N=N[s]["pos"]["l1"], l1_L=L[s]["pos"]["l1"], dv_N=N[s]["dv_sum"], dv_L=L[s]["dv_sum"], dt_N=N[s]["dt_nom_w"], dt_L=L[s]["dt_nom_w"]))
        R[geom] = dict(seeds=seeds, node=sn, base=sl, dv_pair=T28.paired(d), eng=eng, tstar=float(fn[0]["tstar"]), within_change=int(sn["l11"] - sl["l11"]), windows_lost=int(sl["win"] - sn["win"]))
    c = SC.compare_passes(V.ROOT / "rev31node" / "ODY_PL_OLRT3_pert_ODY_s14_L4", V.ROOT / "rev17" / "ODY_PL_OLRT3_pert_ODY_s14_L4"); reg = dict(equal=c["equal"], ncols=c["ncols"], rows=(c["rows_a"], c["rows_b"]), first=c["first"])
    S = pd.read_csv(next(V.ROOT.glob("summary_rev31node_*.csv"))); rr = S[S.arm == "OLRT3"].iloc[0]; reg["ei"] = float(rr.ei_km); reg["stop"] = str(rr.stop_reason)
    (out / "rev31node.json").write_text(json.dumps(SC.nan_none(dict(cells=R, regression=reg)), indent=1, default=float))
    for g, r in R.items():
        for tag in ("node", "base"):
            s = r[tag]; print(g, tag, f"n {s['n']} ok {s['ok']} win {s['win']} corr {s['corr']} raw {s['raw_med']:.2f} l1 {s['l1_med']:.2f} l11 {s['l11']} raw1 {s['raw1']} ALL {s['all']} cross {s['cross_med']:.2f} alongafter {s['along_after_med']:.2f} share {s['cross_share']:.2f} dV {s['dv_med']:.2f}")
        print("  dV pair", np.round(r["dv_pair"]["d"], 2).tolist(), round(r["dv_pair"]["med"], 2), np.round(r["dv_pair"]["ci"], 2).tolist(), r["dv_pair"]["verdict"], "within change", r["within_change"], "windows lost", r["windows_lost"])
        for e in r["eng"]: print("   ", e["seed"], e["first_dv_diff_pass"], "ident", e["identical"], "tt", round(e["ttg_tt_min_d"], 3), round(e["ttg_tt_max_d"], 3), "solves", e["solves_N"], e["solves_L"], "eng", e["engaged_N"], e["engaged_L"], "win", e["win_N"], e["win_L"], "l1", round(e["l1_N"], 2), round(e["l1_L"], 2), "dV", round(e["dv_N"], 2), round(e["dv_L"], 2))
    print(reg)


if __name__ == "__main__":
    main()
