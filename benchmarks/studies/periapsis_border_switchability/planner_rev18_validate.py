"""Revision 18: in-model arrival time (planner_rev18_offline.arrival) driven by the FLOWN ln q path of each flight, against the flown arrival time (analysis only). Usage: python3 -I planner_rev18_validate.py --flights rev18_flights.json --out F"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev18_offline as O
import planner_rev18_flights as FL


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--flights", required=True); ap.add_argument("--flights-src", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
    src = {r["run_id"]: r for r in json.load(open(a.flights_src)) if r.get("status") != "not_corrected"}; J = json.load(open(a.flights))["flights"]; tabs = {g: O.load_tab(g) for g in ("ODY", "G1")}; rows = []
    for f in J:
        if f["geom"] not in tabs: continue
        r = src[f["run_id"]]; d = FL.run_dir(r["root"], r["phase"], r["run_id"]); P = pd.read_csv(d / "passes.csv").sort_values("pass"); n = f["n"]
        x = np.log(P.q_w_cm2.values[:n].astype(float)); m = O.arrival(tabs[f["geom"]], x, math.log(float(P.q_w_cm2.values[0])))
        rows.append(dict(cell=f["arm"].replace("-", "") + "@" + f["geom"], flown_d=f["t_arr_d"], model_d=m["days"], flown_n=n, model_n=m["N"]))
    D = pd.DataFrame(rows); D["ratio"] = D.flown_d / D.model_d; g = D.groupby("cell").agg(n=("ratio", "size"), ratio_med=("ratio", "median"), ratio_min=("ratio", "min"), ratio_max=("ratio", "max"), flown_med=("flown_d", "median"), model_med=("model_d", "median"))
    out = dict(cells=g.reset_index().to_dict("records"), overall_ratio_median=float(D.ratio.median()), overall_ratio_iqr=[float(D.ratio.quantile(.25)), float(D.ratio.quantile(.75))], n=int(len(D)))
    json.dump(out, open(a.out, "w"), indent=1, default=float); print(g.round(3).to_string()); print(out["overall_ratio_median"], out["overall_ratio_iqr"], out["n"])


if __name__ == "__main__":
    main()
