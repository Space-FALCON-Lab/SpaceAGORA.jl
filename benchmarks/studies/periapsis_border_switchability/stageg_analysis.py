"""Stage G analysis: asymmetric adaptive band versus the fixed 1000 m margin and versus Stage F adaptive, on the held-out cells of the Stage F at 20x20 tree.
Plan: tasks/2026-10-15_periapsis-region-border-targeting/stageg_plan.md (frozen before any run). Reuses stagef20_analysis.py (run_metrics, summaries, boot_ci) unchanged.
Usage: python3 -I stageg_analysis.py --results-root DIR [--ctrl adaptive_asym --phase f2g] [--tag g]   (--ctrl adaptive --phase f2 --tag plumbing reproduces the Stage F rows: a plumbing check)"""
import argparse, json, pathlib, sys
import numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import stagef20_analysis as sf

AN, X = sf.AN, sf.AN["win_cost_tolerance"]
MARS_DV_FACTOR = 2.0   # G3 (plan section 4)


def cell_table(root, pl, ctrls):
    """One row per held-out cell with run metrics of every arm in ctrls = {name: (ctrl_dir_prefix, summary_phase)}; None if any run is missing."""
    cells = sf.cells_of(root, pl); rows = []
    S = {k: sf.summaries(root, ph) for k, (c, ph) in ctrls.items()}
    for t in cells[cells.split == "heldout"].tag:
        d = {}
        for k, (c, ph) in ctrls.items():
            f = sf.sfdir(root) / pl / "runs" / f"{c}_pert_{t}" / "passes.csv"
            sr = S[k][(S[k].planet == pl) & (S[k].run_id == f"{c}_pert_{t}")] if len(S[k]) else []
            if not f.exists() or len(sr) == 0: d = None; break
            d[k] = sf.run_metrics(pd.read_csv(f), sr.iloc[0].to_dict())
        if d and all(d.values()): rows.append(dict(tag=t, **{f"{k}_{m}": v for k, dd in d.items() for m, v in dd.items()}))
    return pd.DataFrame(rows)


def arm_summary(R, a):
    tot = lambda k: float(R[f"{a}_{k}"].sum())
    bad = int(((R[f"{a}_retcode"] != "Success") | R[f"{a}_stopped"]).sum())
    return dict(passes=tot("n_pass"), low=int(tot("low")), soft=int(tot("soft")), hard=int(tot("hard")), frac_low=tot("low") / tot("n_pass"), frac_soft=tot("soft") / tot("n_pass"),
                frac_hard=tot("hard") / tot("n_pass"), dv_mean=float(R[f"{a}_dv"].mean()), burns_mean=float(R[f"{a}_burns"].mean()), decay_median_m_per_day=float(R[f"{a}_decay_m_per_day"].median()),
                t_end_days_mean=float(R[f"{a}_t_end_days"].mean()), bad_runs=bad, collapsed=int(R[f"{a}_collapsed"].sum()))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--results-root", required=True); ap.add_argument("--ctrl", default="adaptive_asym"); ap.add_argument("--phase", default="f2g"); ap.add_argument("--tag", default="g")
    a = ap.parse_args(); root = a.results_root; rng = np.random.default_rng(AN["bootstrap_seed"]); res = {}
    for pl in sf.SF["planets"]:
        if not (sf.sfdir(root) / pl).exists(): continue
        R = cell_table(root, pl, {"fixed": ("fixed", "f1"), "sf": ("adaptive", "f2"), "asym": (a.ctrl, a.phase)})
        if len(R) == 0: continue
        R.to_csv(sf.sfdir(root) / "analysis" / f"{a.tag}_cells_{pl}.csv", index=False)
        S = {k: arm_summary(R, k) for k in ("fixed", "sf", "asym")}
        dsoft = (R.asym_soft / R.asym_n_pass - R.fixed_soft / R.fixed_n_pass).values; lo, hi = sf.boot_ci(dsoft, rng, AN["bootstrap_draws"])
        f, g = S["fixed"], S["asym"]
        # Stage F outcome classes of asym against fixed (same formulas as stagef20_analysis.phase_f2)
        w1 = bool(g["frac_soft"] < f["frac_soft"] and hi < 0 and g["frac_hard"] <= f["frac_hard"])
        w2 = bool(g["decay_median_m_per_day"] >= (1 - X) * f["decay_median_m_per_day"] and g["dv_mean"] <= (1 + X) * f["dv_mean"] and g["bad_runs"] <= f["bad_runs"])
        loss = bool(lo > 0 or g["frac_hard"] > f["frac_hard"])
        cost_win = bool((not loss) and g["decay_median_m_per_day"] > (1 + X) * f["decay_median_m_per_day"] and g["dv_mean"] <= (1 + X) * f["dv_mean"])
        outcome = "LOSS" if loss else "WIN" if (w1 and w2) else "SAFETY_TRADE" if w1 else "COST_WIN" if cost_win else "TIE"
        v = dict(G1_hard_le_fixed=bool(g["hard"] <= f["hard"]), G2_soft_lt_fixed=bool(g["frac_soft"] < f["frac_soft"]), G3_dv_le_2x_fixed=bool(g["dv_mean"] <= MARS_DV_FACTOR * f["dv_mean"]),
                 G5_bad_runs_le_fixed=bool(g["bad_runs"] <= f["bad_runs"]),
                 G6_no_loss=bool(not loss), G6_dv_le_1p1x_fixed=bool(g["dv_mean"] <= (1 + X) * f["dv_mean"]))
        v["Mars_verdict_G1_G2_G3_G5"] = bool(v["G1_hard_le_fixed"] and v["G2_soft_lt_fixed"] and v["G3_dv_le_2x_fixed"] and v["G5_bad_runs_le_fixed"])
        v["EarthVenus_verdict_G6"] = bool(v["G6_no_loss"] and v["G6_dv_le_1p1x_fixed"])
        res[pl] = dict(n_cells=len(R), arms=S, dsoft_asym_minus_fixed_mean=float(np.mean(dsoft)), dsoft_ci90=[lo, hi], W1=w1, W2=w2, outcome_vs_fixed=outcome,
                       low_q_fraction_cost=dict(fixed=f["frac_low"], stagef_adaptive=S["sf"]["frac_low"], asym=g["frac_low"]), criteria=v)
    (sf.sfdir(root) / "analysis" / f"{a.tag}_results.json").write_text(json.dumps(res, indent=1, default=str)); print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__": main()
