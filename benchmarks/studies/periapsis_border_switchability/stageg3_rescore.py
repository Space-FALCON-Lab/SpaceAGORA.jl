"""Post hoc re-scoring of the Mars gate (stageg3_data) and the Venus gate (stageg3v_data) under plan Revision 10 (user, 2026-10-27): a flight that reaches the terminal condition (end radius) with no constraint violation
(no pass above the V_rel never-exceed cap, no safety stop) is NOT insufficient because of its length. It is scored on the passes it has; metrics that need a window (15-pass smoothed correlation, drift-range ratio, floor on the
smoothed path) are 'not scored' when the flight is shorter than one smoothing window (then the flight is judged on terminal condition and no violation). Safety-stopped or violating flights keep their original treatment.
THIS IS A POST HOC CHANGE made after the Venus verdict; the frozen verdicts are reported beside the re-scored ones. NO frozen file is edited and no flight is flown: stageg3_analysis.py and stageg3v_analysis.py are IMPORTED
(criteria, burn removal, model fed with the flown epochs); only the minimum-pass threshold is lowered on a COPY of the config. Usage: python3 -I stageg3_rescore.py [--out FILE]"""
import argparse, copy, glob, json, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import stageg3_analysis as A3
import stageg3v_analysis as SV
import geometry_rank as GR

TD = A3.SOA.task_dir()
NEVER_EXCEED = 0.65      # ASSUMED (U19), the safety-stop value of both gates


def early_success(row, run_dir):
    """Revision 10 test on one flight: end radius reached, no safety stop, no pass-maximum Q_s on V_rel above the never-exceed cap."""
    P = pd.read_csv(run_dir / "passes.csv").dropna(subset=["ra_m"]); qr = A3.qs_rel(P.rho_peak_kg_m3, P.vrel_mps)
    ok = row.stop_reason == "end_radius" and not (isinstance(row.safety_reason, str) and row.safety_reason) and not (qr > NEVER_EXCEED).any()
    return bool(ok), float(qr.max()), len(P)


def rescore_gate(name, geoms, cfg_for_score, S, runs, ctx, field_of, ratio, cr, outcome_fn, roles, orig_path):
    orig = json.load(open(orig_path)); window = int(cr["smooth_passes"]); cfgs = copy.deepcopy(cfg_for_score); cfgs["stageg3"]["criteria"]["min_passes_scored"] = window
    ok = S[S.status == "ok"]; per = {}; recs = []
    for g in geoms:
        per[g["id"]] = []
        for r in ok[ok.geometry == g["id"]].itertuples():
            rd = runs / r.run_id; o = next(f for f in orig["flights"] if f["run"] == r.run_id)
            if not o.get("insufficient"):
                rec = dict(o, rescored=False, treatment="unchanged (scored in the frozen analysis)"); per[g["id"]].append(rec); recs.append(rec); continue
            es, qmax, npass = early_success(r, rd)
            if not es:
                rec = dict(o, rescored=False, treatment="unchanged (safety stop or violation: still INSUFFICIENT)", stop_reason=r.stop_reason, qs_rel_max=qmax); per[g["id"]].append(None); recs.append(rec); continue
            if npass >= window:
                rec = A3.score_flight(ctx, cfgs, g, field_of(g), rd, ratio); rec["treatment"] = "early success: scored on %d available passes (Revision 10)" % npass
            else:
                rec = dict(run=r.run_id, n_flown=npass, insufficient=False, corr=None, R_ratio=None, F_flown=None, F_pred=None, a_pass=None, b_pass=None, c_pass=None, **{"pass": True},
                           treatment="early success, %d passes < one smoothing window (%d): window metrics NOT SCORED; judged on terminal condition and no violation" % (npass, window))
            rec.update(geometry=g["id"], role=g["role"], keeper=r.keeper, stop_reason=r.stop_reason, qs_rel_max=qmax, rescored=True, original="INSUFFICIENT (%s)" % o.get("reason"))
            per[g["id"]].append(None if rec.get("insufficient") else rec); recs.append(rec)
    ver = {g["id"]: A3.verdict(per[g["id"]], cr["flights_rule"]) for g in geoms}
    return dict(gate=name, original_verdict={k: v["passed"] for k, v in orig["verdict"].items()}, original_outcome=orig["outcome"], rescored_verdict=ver, flights=recs, **outcome_fn(ver))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default=str(TD / "stageg3_rescore.json")); a = ap.parse_args()
    res = dict(note="POST HOC (plan Revision 10, user 2026-10-27). Frozen verdicts kept beside the re-scored ones; no frozen file edited; no flight flown.")
    # Mars gate
    cfg3 = A3.load3(); ctx = A3.make_ctx(cfg3); root = TD / "stageg3_data" / "stageg3"; S = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(root / "summary_gate_mars_*.csv")))], ignore_index=True)
    mf = lambda g: A3.meanfield_for(cfg3, ctx, None, g)[0]
    def mars_out(ver):
        r = {g["role"]: ver[g["id"]] for g in cfg3["stageg3"]["geometry"]}; return dict(rescored_outcome=A3.outcome(r["G1"], r["G2"]), gworst_informative=r["Gworst"])
    res["mars"] = rescore_gate("Mars gate (G1, G2, Gworst)", cfg3["stageg3"]["geometry"], cfg3, S, root / "runs", ctx, mf, float(cfg3["stageg3"]["burn_response_ratio"]), cfg3["stageg3"]["criteria"], mars_out, None, root / "analysis" / "gate_verdict.json")
    # Venus gate
    cfgv = SV.loadv(); ctxv, _ = SV.make_ctx(cfgv); rv = TD / "stageg3v_data"; Sv = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(str(rv / "summary_gate_venus_*.csv")))], ignore_index=True)
    F, _ = SV.load_field(cfgv)
    def ven_out(ver): return dict(rescored_outcome=SV.outcome(ver["VG1"], ver["VGref"]), vg2_diagnostic=ver["VG2"])
    res["venus"] = rescore_gate("Venus gate (VG1, VG2, VGref)", cfgv["stageg3v"]["geometry"], cfgv, Sv, rv / "runs", ctxv, lambda g: SV.field_of(F, g, "full"), float(cfgv["stageg3v"]["burn_response_ratio"]), cfgv["stageg3v"]["criteria"], ven_out, None, rv / "analysis" / "gate_verdict.json")
    pathlib.Path(a.out).write_text(json.dumps(res, indent=1, default=float))
    for k in ("mars", "venus"):
        r = res[k]; print(r["gate"], "| original", r["original_verdict"], r["original_outcome"], "| rescored", r["rescored_verdict"], r["rescored_outcome"])
        for f in r["flights"]:
            if f.get("rescored") or f.get("insufficient"): print("  ", f.get("run"), f.get("treatment"), {x: (round(f[x], 3) if isinstance(f.get(x), float) else f.get(x)) for x in ("corr", "R_pred", "R_flown", "R_ratio", "F_pred", "F_flown", "F_tol", "regime_pred", "regime_flown", "a_pass", "b_pass", "c_pass", "pass")})


if __name__ == "__main__": main()
