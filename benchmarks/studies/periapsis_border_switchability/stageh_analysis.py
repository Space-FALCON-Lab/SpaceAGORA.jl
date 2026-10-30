"""Stage H analysis: asymmetric band with a partial low-heat relaxation e_lo_frac * |e_lo|. Plan: tasks/2026-10-15_periapsis-region-border-targeting/stageh_plan.md (frozen before any run).
Reuses stagef20_analysis.py (run_metrics, summaries, boot_ci) unchanged and the Stage G criterion formulas (G1-G6).
  --mode dev  --frac F          Mars development cells, flown at relaxation F, against the fixed 1000 m arm on the same cells; prints the pre-registered selection criteria D1-D4
  --mode held --frac F          held-out cells of every planet with runs, criteria G1-G6 against fixed, with Stage F adaptive and Stage G reported alongside
Usage: python3 -I stageh_analysis.py --results-root DIR --mode dev|held --frac 0.25 [--planets mars]"""
import argparse, json, pathlib, sys
import numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import stagef20_analysis as sf

AN, X = sf.AN, sf.AN["win_cost_tolerance"]
DV_FACTOR = 2.0   # G3 / D3, unchanged from Stage G


def cell_table(root, pl, split, ctrls):
    """One row per cell of the split with run metrics of every arm in ctrls = {name: (ctrl_dir_prefix, summary_phase)}; cells where any arm is missing are dropped."""
    cells = sf.cells_of(root, pl); rows = []
    S = {k: sf.summaries(root, ph) for k, (c, ph) in ctrls.items()}
    for t in cells[cells.split == split].tag:
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


def mars_criteria(f, g):
    c = dict(G1_hard_le_fixed=bool(g["hard"] <= f["hard"]), G2_soft_lt_fixed=bool(g["frac_soft"] < f["frac_soft"]), G3_dv_le_2x_fixed=bool(g["dv_mean"] <= DV_FACTOR * f["dv_mean"]),
             G5_bad_runs_le_fixed=bool(g["bad_runs"] <= f["bad_runs"]))
    c["all"] = all(c.values()); return c


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--results-root", required=True); ap.add_argument("--mode", required=True, choices=["dev", "held"]); ap.add_argument("--frac", type=float, required=True)
    ap.add_argument("--planets", default="venus,mars,earth"); a = ap.parse_args(); root = a.results_root; code = f"{round(100 * a.frac):03d}"; ctrl = f"adaptive_h{code}"; res = {}
    rng = np.random.default_rng(AN["bootstrap_seed"])
    for pl in a.planets.split(","):
        if not (sf.sfdir(root) / pl).exists(): continue
        if a.mode == "dev":
            R = cell_table(root, pl, "dev", {"fixed": ("fixed", "f1"), "h": (ctrl, f"hdev{code}")})
            if len(R) == 0: continue
            f, g = arm_summary(R, "fixed"), arm_summary(R, "h"); c = mars_criteria(f, g)
            res[pl] = dict(frac=a.frac, n_cells=len(R), arms=dict(fixed=f, h=g), criteria={k.replace("G", "D"): v for k, v in c.items()})
            continue
        arms = {"fixed": ("fixed", "f1"), "sf": ("adaptive", "f2"), "g": ("adaptive_asym", "f2g"), "h": (ctrl, f"hheld{code}")}
        R = cell_table(root, pl, "heldout", arms)
        if len(R) == 0: continue
        R.to_csv(sf.sfdir(root) / "analysis" / f"h{code}_cells_{pl}.csv", index=False)
        S = {k: arm_summary(R, k) for k in arms}
        dsoft = (R.h_soft / R.h_n_pass - R.fixed_soft / R.fixed_n_pass).values; lo, hi = sf.boot_ci(dsoft, rng, AN["bootstrap_draws"]); f, g = S["fixed"], S["h"]
        w1 = bool(g["frac_soft"] < f["frac_soft"] and hi < 0 and g["frac_hard"] <= f["frac_hard"])
        w2 = bool(g["decay_median_m_per_day"] >= (1 - X) * f["decay_median_m_per_day"] and g["dv_mean"] <= (1 + X) * f["dv_mean"] and g["bad_runs"] <= f["bad_runs"])
        loss = bool(lo > 0 or g["frac_hard"] > f["frac_hard"])
        cost_win = bool((not loss) and g["decay_median_m_per_day"] > (1 + X) * f["decay_median_m_per_day"] and g["dv_mean"] <= (1 + X) * f["dv_mean"])
        outcome = "LOSS" if loss else "WIN" if (w1 and w2) else "SAFETY_TRADE" if w1 else "COST_WIN" if cost_win else "TIE"
        c = mars_criteria(f, g); c["G6_no_loss"] = bool(not loss); c["G6_dv_le_1p1x_fixed"] = bool(g["dv_mean"] <= (1 + X) * f["dv_mean"]); c["G6_all"] = bool(c["G6_no_loss"] and c["G6_dv_le_1p1x_fixed"])
        res[pl] = dict(frac=a.frac, n_cells=len(R), arms=S, dsoft_h_minus_fixed_mean=float(np.mean(dsoft)), dsoft_ci90=[lo, hi], W1=w1, W2=w2, outcome_vs_fixed=outcome, criteria=c,
                       decay_vs_fixed_pct=100 * (g["decay_median_m_per_day"] / f["decay_median_m_per_day"] - 1))
    (sf.sfdir(root) / "analysis").mkdir(parents=True, exist_ok=True)
    (sf.sfdir(root) / "analysis" / f"h_{a.mode}{code}_results.json").write_text(json.dumps(res, indent=1, default=str)); print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__": main()
