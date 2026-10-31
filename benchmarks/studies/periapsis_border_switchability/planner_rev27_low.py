"""rev27low supporting analysis (Python only, read-only on flight data): (1) regression flight against the Rev 17 flight, (2) the dV anatomy per low-altitude cell: table drift, tube, model-world arrival range against the guard t*, flown level, burns.
Usage: python3 -I planner_rev27_low.py --out-dir DIR"""
import argparse, json, math, os, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_rev21_position as PS
import planner_rev20_score as SC
ROOT = pathlib.Path(os.environ.get("SPACEAGORA_PLANNER_ROOT", str(pathlib.Path.home() / "Documents" / "SpaceAGORA-planner-data" / "planner")))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True); R = {}
    d1 = ROOT / "rev27low" / "ODY_PL_OLRT3_pert_ODY_s14_L4"; d0 = ROOT / "rev17" / "ODY_PL_OLRT3_pert_ODY_s14_L4"; c = SC.compare_passes(d1, d0)
    S = SC.load_summary(ROOT, "summary_rev27low_*.csv"); S17 = SC.load_summary(ROOT, "summary_rev17_*.csv"); f1 = SC.flight(d1, S.get(d1.name)); f0 = SC.flight(d0, S17.get(d0.name))
    R["regression"] = dict(equal=c["equal"], ncols=c["ncols"], rows=(c["rows_a"], c["rows_b"]), first=c["first"], stop=(f1["stop"], f0["stop"]), same_stop=f1["stop"] == f0["stop"], ei=(float(S[d1.name]["ei_km"]) if "ei_km" in S[d1.name] else None))
    RP = PS._rp(); sw, P = RP.arm_def("OLR"); X0 = math.log(0.212); R["cells"] = {}
    for tok in ("G0486", "G1230", "G1038"):
        tab = PS.get_tab(tok); inp = PLB.ol_inputs(sw, tab, P, X0); lo = np.array(inp["lo"]); hi = np.array(inp["hi"]); path = np.array(tab["path"], float); tr = PS.target_radius(tok); e_t = (tr - 3.504e6) / (tr + 3.504e6)
        def arr(x, kap): return PS.arrival_history(tab, np.asarray(x, float), X0, kap, tr)["days"]
        x_nat = X0 + path - path[0]; n = len(x_nat); xs = lambda v: np.concatenate(([X0], np.asarray(v, float)))[:n]
        K = 0.925      # the measured drop-law ratio of the block (results_rev27low_droplaw)
        cell = dict(n_table=int(tab["n"]), q_req=float(tab["q_req"]), W=float(tab["W"]), drift_ln=float(path[-1] - path[0]), tube_ln=float(np.median(hi - lo)), ln_q_cap=float(math.log(P["q_cap_mean_w_cm2"])), ln_q_req=float(math.log(tab["q_req"])),
                    e_target=e_t, arr_nat=arr(x_nat, 1.0), arr_hi=arr(xs(hi), 1.0), arr_lo=arr(xs(lo), 1.0), arr_req=arr(np.full(3 * n, math.log(tab["q_req"])), 1.0), arr_hi_k=arr(xs(hi), K), arr_lo_k=arr(xs(lo), K), arr_req_k=arr(np.full(3 * n, math.log(tab["q_req"])), K), tstar=PS.nominal(tok)["days"])
        fl = []
        for dd in sorted((ROOT / "rev27low").glob(f"{tok}_PL_OLRT3-W3E_pert_ODY_s*_L0")):
            Pp = pd.read_csv(dd / "passes.csv").sort_values("pass"); D = pd.read_csv(dd / "decisions.csv"); dv = D.delta_v_mps.fillna(0).values; nz = dv[np.abs(dv) > 0]; lq = np.log(Pp.q_w_cm2.values.astype(float))
            fl.append(dict(run=dd.name, passes=len(Pp), days=float(Pp.time_s.iloc[-1] / 86400), burns=len(nz), sum_abs=float(np.abs(nz).sum()), net=float(abs(nz.sum())), up=float(nz[nz > 0].sum()), down=float(nz[nz < 0].sum()), med_lnq=float(np.median(lq)), frac_above_cap=float(np.mean(lq > math.log(P["q_cap_mean_w_cm2"]))), med_q=float(np.exp(np.median(lq))), lnq_first=float(lq[0]), lnq_last=float(lq[-1]),
                           lnq_early=float(np.median(lq[:60])), lnq_late=float(np.median(lq[-60:]))))
        cell["flights"] = fl; R["cells"][tok] = cell
    (out / "rev27low_anatomy.json").write_text(json.dumps(R, indent=1, default=float))
    print(R["regression"])
    for t, c in R["cells"].items():
        print(t, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in c.items() if k != "flights"})
        for f in c["flights"]: print("   ", f["run"][-6:], f["passes"], round(f["days"], 1), f["burns"], "sum|dv| %.2f net %.2f up %.2f down %.2f" % (f["sum_abs"], f["net"], f["up"], f["down"]), "med q %.3f lnq early %.2f late %.2f above cap %.2f" % (f["med_q"], f["lnq_early"], f["lnq_late"], f["frac_above_cap"]))


if __name__ == "__main__":
    main()
