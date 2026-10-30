"""Revision 21, WP21-A (analysis only, Python, no Julia): guard v3 work on DEV paths. Imports planner_rev20_replay (guard v2 = GUARD_V2, frozen sha in results_rev20_guard_v2.md).

score(overrides, indexings): the 27 dev paths (devT OLRT + trx_pull4 flown OLRT3, ODY 11 and G1 16) under the given indexings; per indexing the pooled window count, hard violations, apoapsis successes, the
paired median extra campaign dV over OLRT3 alone, and the median arrival error against t* in units of W_max per geometry. `bias_probe`: prediction of the CHOSEN plan against the realized arrival.
Held and Rev 17 paths are never opened here. Usage: python3 -I planner_rev21_guard.py tune NAME  (reads configs from a python expression on argv[2])
"""
import json, math, pathlib, sys, time
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_rev20_replay as R

R.DEV_EXT = True
OUT = R.OUT
BASE = {}

# Guard v3 (Revision 21): FROZEN after the dev-only tuning (v3_stage*.jsonl, 33 configurations x 3 indexings on 27 dev paths): v2 plus a half as wide engagement band (s_frac 0.5 -> 0.25) and the terminal aim emulated in
# the arrival prediction (term_model). Everything else is GUARD_V2 (round-1 parameters, ODY aim offset +0.15 d, deviation threshold 0.05 m/s). See results_rev21_guard.md.
GUARD_V3 = dict(R.GUARD_V2, s_frac=0.25, term_model=True)


def baseline(nm, g):
    if (nm, g) not in BASE:
        BASE[(nm, g)] = R.dev_cell(g, noise=nm)
    return BASE[(nm, g)]


def score(ov, indexings=("pass", "time", "state"), geoms=("ODY", "G1"), base_guard=None):
    gdp = {**(R.GUARD_V2 if base_guard is None else base_guard), **ov}; row = dict(guard=ov, per={})
    for nm in indexings:
        win = hard = apo = 0; extra = []; z = {}; bias = {}
        for g in geoms:
            res = R.dev_cell(g, guard=gdp, noise=nm); b = baseline(nm, g); sm = R.summarize(res)
            win += sm["window_ok"]; hard += sm["hard"]; apo += sm["apo_ok"]
            extra += [t["dv_campaign"] - x["dv_campaign"] for x, t in zip(b, res) if x["reached"] and t["reached"]]
            z[g] = float(np.median([r["dt_d"] / r["w_max_d"] for r in res if r["reached"]]))
        row["per"][nm] = dict(window=win, hard=hard, apo=apo, med_extra=float(np.median(extra)), med_dt_over_W=z)
    row["worst_window"] = min(v["window"] for v in row["per"].values()); row["hard_total"] = sum(v["hard"] for v in row["per"].values())
    row["mean_extra"] = float(np.mean([v["med_extra"] for v in row["per"].values()]))
    return row


def tune(configs, name, indexings=("pass", "time", "state")):
    OUT.mkdir(parents=True, exist_ok=True)
    for ov in configs:
        t0 = time.time(); row = score(ov, indexings); row["wall_s"] = time.time() - t0
        with open(OUT / f"{name}.jsonl", "a") as fh:
            fh.write(json.dumps(row, default=float) + "\n")
        print("v3", ov, "| window %s | worst %d | hard %d | extra dV %s | dt/W ODY,G1 (pass) %.2f %.2f | %.0fs" % ("/".join(str(row["per"][n]["window"]) for n in indexings), row["worst_window"], row["hard_total"],
              "/".join("%.2f" % row["per"][n]["med_extra"] for n in indexings), row["per"][indexings[0]]["med_dt_over_W"]["ODY"], row["per"][indexings[0]]["med_dt_over_W"]["G1"], row["wall_s"]), flush=True)


def bias_probe(ov, geom="ODY"):
    """Median over dev paths of (prediction of the plan the guard CHOSE) minus (realized arrival) at n-30, n-15, n-8 passes before the end (days), pass indexing."""
    gdp = {**R.GUARD_V2, **ov}; ts = R.tstar_for(geom); out = []
    for rid in R.dev_paths_ext(geom):
        F = R.get_flight(rid); r = R.replay(F, geom, "OLRT3", guard=gdp, tstar_d=ts); n = r["n"]; act = r["t_arr_d"]
        pl = {int(a): b for a, b, c in r["info"]["pred_log"]}; lam = {int(l[0]): l[4] for l in r["lam_log"] if len(l) > 4}
        out.append([lam.get(j, pl[j]) - act for j in (max(x for x in pl if x <= n - k) for k in (30, 15, 8))])
    return [float(x) for x in np.median(out, 0)]


if __name__ == "__main__":
    if sys.argv[1] == "tune":
        tune(eval(sys.argv[3]), sys.argv[2])
