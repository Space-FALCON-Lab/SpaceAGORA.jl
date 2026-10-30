"""POST HOC re-score under Revision 13 (user, 2026-10-27): hard violations are failures; extra dV spent to avoid violating a heating constraint is fine. No new flights. Original scoring is kept (results.md).
HARD VIOLATION (only existing frozen quantities, no new threshold): a flight is a hard failure when ANY pass has pass-maximum Q_s on V_rel above the never-exceed limit 0.65 W/cm2 (planner_config q_never_exceed_qs_w_cm2,
SOURCED qualification limit, used as never-exceed by assumption U12) -- equivalently the flight stopped on the safety_q stop (which fires at the closing apoapsis of such a pass). The altitude stop (3 H below the start, OL only)
is NOT a heating quantity and is counted separately (none occurred). SOFT: passes with Q_s above 0.54 (the C2 budget E x n_seeds = 8, unchanged): reported as before, never a failure by itself in the re-score.
Re-scored rules: C2' = no hard failure in the cell (soft count reported beside the original C2); C1' = original C1 computed on the successful (non-hard) flights of the cell only; D1' = ratios among jointly successful seeds < 1
(upper 90% bound) AND C1' AND C2'; D6' (comparator vs reference): the arm with fewer hard failures is preferred and the dV difference is the price of that; with equal failure counts the paired median dV difference decides as in D6.
SURVIVOR BIAS: dV is compared only among successful flights, and a flight that fails is often one that stayed low-dV (never corrected a rising heat rate), so the successful subset of a failing arm is biased toward
cases where it needed fewer corrections or got benign atmosphere; its dV is therefore optimistic relative to its true cost of staying safe, and a paired dV difference against an arm with no failures understates the price of avoiding violations.
Usage: python3 -I planner_rescore13.py --planner-root R/planner --stageo-root S/stageo --out-json F --out-md F"""
import argparse, json, pathlib, sys
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_analysis as A
import planner_lib as PLB
import planner_results as PR


def hard(r):
    return bool(r["v065"] > 0 or "safety_q" in str(r["stop_reason"]))


def alt_stop(r):
    return "safety_alt" in str(r["stop_reason"])


def okset(runs):
    return {k: r for k, r in runs.items() if not hard(r)}


def fmt(b, nd=3):
    return PR.fmt(b, nd)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--planner-root", required=True); ap.add_argument("--stageo-root"); ap.add_argument("--out-json", required=True); ap.add_argument("--out-md", required=True); a = ap.parse_args()
    cfg = PLB.load_cfg(); td = A.task_dir(); roles = A.gate_roles(td); proot = pathlib.Path(a.planner_root); PR.link_sens(proot)
    preds = json.loads((td / cfg["planner"]["predictions"]["path"]).read_text()); tabs = json.loads((td / cfg["planner"]["ground_table"]["path"]).read_text())
    runs = A.collect([proot] + ([pathlib.Path(a.stageo_root)] if a.stageo_root else []), roles, preds, tabs, ("held", "sens")); rng = A.rng_of(); E = cfg["planner"]["E_budget"]
    ok = {k: okset(v) for k, v in runs.items()}; out = {"definition": __doc__.split("Usage:")[0]}
    # per cell table
    cells = {}
    for k, v in sorted(runs.items()):
        s = ok[k]; cells[k] = dict(flights=len(v), hard_fail=len(v) - len(s), safety_q_stops=sum("safety_q" in r["stop_reason"] for r in v.values()), alt_stops=sum(alt_stop(r) for r in v.values()), soft_gt054_passes=int(sum(r["v054"] for r in v.values())),
                                   hard_passes_gt065=int(sum(r["v065"] for r in v.values())), dv_success_median=A.boot([r["dv"] for r in s.values()], np.median, rng), burns_success_median=float(np.median([r["burns"] for r in s.values()])) if s else None,
                                   days_success_median=float(np.median([r["days"] for r in s.values()])) if s else None, dv_all_median=float(np.median([r["dv"] for r in v.values()])))
    out["cells"] = cells
    body, fo = runs.get("BODY@ODY", {}), runs.get("FO@ODY", {}); bases = {n: ok[f"{n}@ODY"] for n in ("BODY", "FO", "N2I") if f"{n}@ODY" in ok}
    # C1', C2', D1'
    resc = {}
    for key in sorted(runs):
        arm = key.split("@")[0]
        if arm in ("FO", "N2I"):
            continue
        s = ok[key]; c1o = A.c1(runs[key], body, fo, rng) if body and fo else None; c1n = A.c1(s, ok["BODY@ODY"], ok["FO@ODY"], rng) if s else None
        c2o = A.c2(runs[key]); resc[key] = dict(C1_original=c1o["pass"] if c1o else None, C1_rescored=c1n["pass"] if c1n else None, C1_rescored_detail=dict(days=c1n["days"], passes=c1n["passes"], n=c1n["n_seeds"]) if c1n else None,
                                               C2_original=c2o["pass"], C2_soft_count=c2o["count_gt054"], C2_soft_budget=E * len(runs[key]), C2_rescored_no_hard_failure=bool(cells[key]["hard_fail"] == 0), hard_fail=cells[key]["hard_fail"])
    out["C1_C2"] = resc
    d1 = {}
    for tok, name in (("G1", "D1a"), ("ODY", "D1b")):
        key = f"L3@{tok}"
        if key in ok and ok[key] and bases:
            ratios = {}
            for n, X in bases.items():
                _, r = A.pair(ok[key], X, lambda x, y: x["dv"] / y["dv"]); ratios[n] = A.boot(r, np.median, rng)
            d1[name] = dict(original=None, ratios_success_only=ratios, C1_rescored=resc[key]["C1_rescored"], C2_rescored=resc[key]["C2_rescored_no_hard_failure"],
                            **{"pass_rescored": bool(all(v["hi"] < 1 for v in ratios.values()) and resc[key]["C1_rescored"] and resc[key]["C2_rescored_no_hard_failure"])})
    out["D1"] = d1

    def d6p(comp, ref):
        C, R = runs.get(comp), runs.get(ref)
        if not C or not R:
            return None
        fc, fr = cells[comp]["hard_fail"], cells[ref]["hard_fail"]
        s, d = A.pair(ok[comp], ok[ref], lambda x, y: x["dv"] - y["dv"]); b = A.boot(d, np.median, rng)
        if fc > fr:
            v = f"{ref.split('@')[0]} preferred: fewer hard failures ({fr} vs {fc}); dV difference is the price of avoiding them"
        elif fc < fr:
            v = f"{comp.split('@')[0]} preferred: fewer hard failures ({fc} vs {fr}); its dV difference is the price"
        else:
            v = "equal hard failures: dV decides: " + ("comparator costs less (SUPPORTS comparator)" if b["hi"] < 0 else "reference costs less" if b["lo"] > 0 else "NULL")
        return dict(comparator=comp, reference=ref, hard_fail_comparator=fc, hard_fail_reference=fr, n_joint_success=len(s), dv_diff_comparator_minus_reference=b, verdict=v)
    out["D6"] = [x for x in [d6p(f"OL@{t}", f"L3@{t}") for t in ("ODY", "G1")] + [d6p(f"OLR@{t}", f"L3@{t}") for t in ("ODY", "G1")] + [d6p(f"OL@{t}", f"OLR@{t}") for t in ("ODY", "G1")] if x]
    out["price_of_avoiding"] = {t: {"OLR_minus_OL_dv_joint_success": A.boot(A.pair(ok[f"OLR@{t}"], ok[f"OL@{t}"], lambda x, y: x["dv"] - y["dv"])[1], np.median, rng), "OLR_minus_OL_burns": A.boot(A.pair(ok[f"OLR@{t}"], ok[f"OL@{t}"], lambda x, y: x["burns"] - y["burns"])[1], np.median, rng),
                                    "OL_hard_failures": cells[f"OL@{t}"]["hard_fail"], "OLR_hard_failures": cells[f"OLR@{t}"]["hard_fail"]} for t in ("ODY", "G1") if f"OL@{t}" in ok and f"OLR@{t}" in ok}
    pathlib.Path(a.out_json).write_text(json.dumps(out, indent=1, default=float))
    L = ["# Results re-score under Revision 13 (POST HOC, by user instruction, 2026-10-27)\n", "No new flights. The original scoring is in results.md and is unchanged; this file re-scores with hard violations counted as failures.\n", "## Definitions\n", out["definition"].replace("\n", " ").strip() + "\n"]
    L.append("## Failures and successful-flight dV per cell (original = all flights, re-score = successful flights only)\n\n" + PR.table([[k, c["flights"], c["hard_fail"], c["safety_q_stops"], c["alt_stops"], c["soft_gt054_passes"], c["hard_passes_gt065"], f"{c['dv_all_median']:.3f}", fmt(c["dv_success_median"]),
                                                                                                                       None if c["burns_success_median"] is None else f"{c['burns_success_median']:.0f}", None if c["days_success_median"] is None else f"{c['days_success_median']:.1f}"] for k, c in cells.items()],
                                                                                                                     ["arm@geometry", "flights", "hard failures", "safety_q stops", "altitude stops", "soft passes >0.54", "passes >0.65", "dV median, all (original)", "dV median, successful [90%]", "burns (successful)", "days (successful)"]))
    L.append("\n## C1 and C2, original and re-scored\n\n" + PR.table([[k, v["C1_original"], v["C1_rescored"], v["C2_original"], f"{v['C2_soft_count']} of {v['C2_soft_budget']:.0f}", v["C2_rescored_no_hard_failure"], v["hard_fail"]] for k, v in resc.items()],
                                                                       ["arm@geometry", "C1 original", "C1 re-scored (successful flights)", "C2 original", "soft passes >0.54 (budget)", "C2 re-scored (no hard failure)", "hard failures"]))
    L.append("\n## D1a and D1b, re-scored (original: both FAIL, see results.md)\n\n")
    for n, v in d1.items():
        L.append(f"{n}: pass (re-scored) = {v['pass_rescored']} (C1' {v['C1_rescored']}, C2' {v['C2_rescored']}); median dV ratio L3 / X among successful seeds: " + ", ".join(f"{x} {fmt(r)}" for x, r in v["ratios_success_only"].items()) + "\n\n")
    L.append("## D6 and OL versus OL-R, re-scored (original verdicts in results.md)\n\n" + PR.table([[f"{d['comparator']} vs {d['reference']}", d["hard_fail_comparator"], d["hard_fail_reference"], d["n_joint_success"], fmt(d["dv_diff_comparator_minus_reference"]), d["verdict"]] for d in out["D6"]],
                                                                                                                 ["pair", "comparator hard failures", "reference hard failures", "jointly successful seeds", "dV diff comparator - reference [90%]", "re-scored verdict"]))
    L.append("\n## Price of avoiding hard violations (OL-R over OL)\n\n" + PR.table([[t, v["OL_hard_failures"], v["OLR_hard_failures"], fmt(v["OLR_minus_OL_dv_joint_success"]), fmt(v["OLR_minus_OL_burns"], 1)] for t, v in out["price_of_avoiding"].items()],
                                                                                   ["geometry", "OL hard failures (of 16)", "OL-R hard failures (of 16)", "extra dV of OL-R over OL, jointly successful seeds, m/s [90%]", "extra burns [90%]"]) +
             "\nThe extra dV is the price of reacting to the measured heating and is reported as such, not as a loss. Survivor bias: OL's failing flights are excluded from this difference, and they are likely those that most needed correction, so the true price is larger than this paired number.\n")
    L.append("\n## Gates\n\nMars and Venus gate rows are not affected: the only gate flight stopped by the Q_s limit (Mars G2 zero-burn s102, safety_q at pass 27) was already scored INSUFFICIENT and stays so; no gate verdict depends on counting 0.65 exceedances differently.\n")
    pathlib.Path(a.out_md).write_text("\n".join(L)); print("wrote", a.out_md)


if __name__ == "__main__":
    main()
