"""Analysis only (no flights): the final apoapsis error of every held and sensitivity flight, and what a terminal-accuracy-first scoring would change. Definitions (read from the code): planner_analysis.run_metrics.
  target = ra_target_m = 7,408,000 m (stageo_config / planner_config [stageo] ra_target_m, SOURCED Table III end of main phase 2).
  value  = ra_m of the last logged pass of the flight: passes.csv column ra_m of pass k is the apoapsis radius at the apoapsis that CLOSES pass k (PassBook.book_apoapsis!, stageo_common.jl). The flight stops at the first closing
           apoapsis with radius <= target, so for a flight that reaches the terminal condition the value is the radius at that stop apoapsis (the first one at or below the target).
  error  = ra_final - target (signed, m). Negative = below (inside) the target, i.e. overshoot past it; positive = still above it (flight stopped early). Reached flights therefore always have error <= 0, bounded below by
           about -(apoapsis drop of the last pass).
  stops:  safety_q: the Q_s test runs at the closing apoapsis AFTER the pass has been logged, so the value is that apoapsis (above target, error > 0); safety_alt (OL only, none occurred): at a periapsis event, so the value is the
           previous closing apoapsis; horizon: the last pass can be open (ra NaN), the code uses the last finite ra_m (none of the flights hit the horizon).
Usage: python3 -I planner_terminal_apoapsis.py --planner-root R/planner --stageo-root S/stageo --out-md F --out-json F"""
import argparse, json, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_analysis as A
import planner_lib as PLB
import planner_results as PR

TARGET = A.RA_TARGET; TOLS_KM = (10.0, 25.0, 50.0, 100.0)    # ASSUMED tolerances: no apoapsis tolerance is defined in any frozen file (Stage O plan: "main-phase scope stops by rule, so this is not a metric")


def per_run(runs_dir, pre):
    rows = {}
    return rows


def collect_terminal(roots, roles):
    import glob
    out = {}
    for root in roots:
        root = pathlib.Path(root)
        for f in sorted(glob.glob(str(root / "summary_held_*.csv")) + glob.glob(str(root / "summary_sens_*.csv"))):
            S = pd.read_csv(f); S = S[S.status == "ok"]
            for r in S.itertuples():
                arm = str(r.arm).replace("-", ""); tok = roles.get(getattr(r, "geometry", "ODY"), "ODY"); d = A._run_dir(root, r.run_id) or (root / "sens" / r.run_id)
                if not (d / "passes.csv").exists():
                    continue
                P = pd.read_csv(d / "passes.csv").sort_values("pass"); ra = P.ra_m.values.astype(float); fin = ra[np.isfinite(ra)]
                seed = int(str(r.run_id).split("_pert_")[1].split("_")[1].lstrip("s")); reached = bool(len(fin) and fin[-1] <= TARGET) or str(r.stop_reason) == "end_radius"
                qs = 0.5 * P.rho_peak_kg_m3.values * P.vrel_mps.values ** 3 / 1e4; hardf = bool(np.nanmax(qs) > 0.65 or "safety_q" in str(r.stop_reason))
                drop = float(fin[-2] - fin[-1]) if len(fin) > 1 else float("nan"); D = d / "decisions.csv"
                Dd = pd.read_csv(D) if D.exists() and D.stat().st_size > 2 else pd.DataFrame(dict(target_pass=[], delta_v_mps=[]))
                lastb = int(Dd[Dd.delta_v_mps.abs() >= 1e-3].target_pass.max()) if (Dd.delta_v_mps.abs() >= 1e-3).any() else 0
                out.setdefault(f"{arm}@{tok}", {})[seed] = dict(err_km=(fin[-1] - TARGET) / 1e3, ra_final_km=fin[-1] / 1e3, reached=reached, passes=len(P), days=float(r.t_end_days), stop=str(r.stop_reason), hard=hardf, last_drop_km=drop / 1e3,
                                                                frac=(TARGET - fin[-1]) / drop if drop > 0 else float("nan"), last_burn_pass=lastb, q_last10=float(np.mean(P.q_w_cm2.values[-10:])), alt_last_km=float(P.alt_m.values[-1] / 1e3), dv=float(r.total_dv_mps))
    return out


def med(x):
    x = [v for v in x if np.isfinite(v)]; return float(np.median(x)) if x else float("nan")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--planner-root", required=True); ap.add_argument("--stageo-root"); ap.add_argument("--out-md", required=True); ap.add_argument("--out-json", required=True); a = ap.parse_args()
    td = A.task_dir(); roles = A.gate_roles(td); proot = pathlib.Path(a.planner_root); T = collect_terminal([proot] + ([pathlib.Path(a.stageo_root)] if a.stageo_root else []), roles); rng = A.rng_of()
    pathlib.Path(a.out_json).write_text(json.dumps(T, indent=1, default=float))
    L = ["# Terminal apoapsis: definition, per-cell results and what a terminal-accuracy-first scoring would change\n", "Analysis only: no flights, no frozen file touched. Source: passes.csv, decisions.csv, summaries of the held and sensitivity flights (planner_terminal_apoapsis.py, JSON beside this file in the results root).\n",
         "## 1. Definition (from the code)\n\n" + __doc__.split("Usage:")[0].split("Definitions (read from the code):")[1].strip() + "\n\n**Tolerance.** No final-apoapsis tolerance is defined in any frozen file. The Stage O plan says the main-phase scope 'stops by rule, so this is not a metric' (odyssey_comparison_plan.md). The tolerances used below (10, 25, 50, 100 km) are ASSUMED and only illustrate the effect.\n"]
    rows = []
    for k, v in sorted(T.items()):
        s = {x: r for x, r in v.items() if not r["hard"]}; h = {x: r for x, r in v.items() if r["hard"]}
        e = [r["err_km"] for r in s.values()]; b = A.boot(e, np.median, rng)
        rows.append([k, len(v), f"{np.mean([r['reached'] for r in v.values()]):.2f}", len(h), f"{med(e):.1f}", f"{min(e):.1f} to {max(e):.1f}" if e else "n/a", f"{med([r['last_drop_km'] for r in s.values()]):.1f}", f"{med([r['frac'] for r in s.values()]):.2f}",
                     f"{med([r['passes'] for r in s.values()]):.0f}", f"{med([r['days'] for r in s.values()]):.1f}", ", ".join(sorted({r['stop'] for r in v.values()})), (", ".join(f"s{x}:{r['err_km']:+.0f}" for x, r in sorted(h.items())) or "-")])
    L.append("## 2. Per cell (successful = no hard violation; error in km, negative = below the target)\n\n" + PR.table(rows, ["arm@geometry", "flights", "fraction reached", "hard failures", "median error, successful", "range, successful", "median last-pass apoapsis drop", "median overshoot as fraction of last drop", "median passes", "median days", "stop reasons", "hard-failed flights (seed: error km)"]))
    L.append("\nFull per-seed errors (km), successful flights first, hard-failed flagged with *:\n")
    for k, v in sorted(T.items()):
        L.append(f"- {k}: " + ", ".join(f"s{x}:{r['err_km']:+.1f}{'*' if r['hard'] else ''}" for x, r in sorted(v.items())))
    # mechanism tests
    def cmp(key):
        return {x: r for x, r in T.get(key, {}).items() if not r["hard"]}
    L.append("\n## 3. Why OL-R's median error at Odyssey (-40.5 km) is worse than OL's (-14.9 km)\n")
    ol, olr = cmp("OL@ODY"), cmp("OLR@ODY"); both = sorted(set(ol) & set(olr))
    f = lambda d, k: [d[s][k] for s in both]
    L.append(f"Numbers on the {len(both)} seeds where both are successful (the -14.9 km in the original table is the median over all 16 OL flights, 4 of which stopped early with positive errors): median error OL {med(f(ol,'err_km')):.1f} km, OL-R {med(f(olr,'err_km')):.1f} km; median last-pass apoapsis drop OL {med(f(ol,'last_drop_km')):.1f} km, OL-R {med(f(olr,'last_drop_km')):.1f} km; median overshoot fraction OL {med(f(ol,'frac')):.2f}, OL-R {med(f(olr,'frac')):.2f}; mean heat rate over the last 10 passes OL {med(f(ol,'q_last10')):.3f}, OL-R {med(f(olr,'q_last10')):.3f} W/cm2; pass of the last burn OL {med(f(ol,'last_burn_pass')):.0f}, OL-R {med(f(olr,'last_burn_pass')):.0f} of about {med(f(olr,'passes')):.0f}.\n")
    L.append("""
Facts from the code (tested by reading it): (a) the DP objective contains only the heat-rate tube and the burn cost (lever weights, optional burn penalty); there is NO terminal-apoapsis term, no target r_a in the planner, and no end-game (final raise and trim are out of scope, plan R1.6 U6). (b) A tangential apoapsis burn does not change the apoapsis radius at the burn; it moves the next periapsis, hence the heat rate and so the decay of the apoapsis over later passes. The planner can therefore only change when the apoapsis reaches 7,408 km (time and passes), not the value at the stop. (c) The flight stops at the first closing apoapsis at or below the target, so the error of a reached flight is the part of the last pass's drop that overshoots the target; it lies between minus the last-pass drop and zero. The error is therefore set by the size of the last-pass drop (driven by the heat rate at the end) and by where the target falls inside that last step (arbitrary phase), not by a planning decision. The 'stop short of the target' case happens only on a hard stop (safety_q), where the error is positive and large.

Tested here (numbers above): (i) the gap in the headline medians is mostly a median artifact. The -14.9 km quoted for OL is the median over all 16 OL flights, four of which stopped on safety_q with POSITIVE errors (+36, +464, +3,627, +4,539 km); on the 12 seeds where both arms are successful the medians are close (OL -32.0 km, OL-R -36.5 km; paired seeds, see the numbers above), and OL-R is the arm that reaches the terminal condition in 16 of 16 flights. (ii) The error of a reached flight is the fraction of the last pass's apoapsis drop that overshoots the target (overshoot fraction 0.35 for OL and 0.49 for OL-R on common seeds, a phase quantity), times the last drop (OL 84 km, OL-R 77 km): so the remaining 4 km difference comes from the overshoot phase, not from the drop size. Hypothesis (not tested): any systematic difference in the heat level near the end (the logs give the last-10-pass mean heat rate above) would act only through the drop size. The planner does not stop short of the target by choice and the error is not the error at a last decision: it is measured at the closing apoapsis of the last flown pass, after all decisions.
""")
    # terminal-first scoring
    L.append("## 4. Does the planner treat the terminal apoapsis as the goal, and what would terminal-first scoring change\n\nAs built, the terminal apoapsis is the stopping condition of the flight, not an objective: the planner (and B-ODY, L0 to L3, OL, OL-R) minimizes burn dV subject to the heat-rate tube; the campaign ends when the apoapsis first reaches the target. Duration is a constraint through q_req (the tube's lower edge, D_ref) and is scored by C1. So 'did it reach the terminal condition without a hard violation' is already implicit in the stop reason; the size of the overshoot is not controlled.\n")
    L.append("Terminal-accuracy-first scoring (ASSUMED tolerance tol; success = reached, |error| <= tol, no hard violation; then dV among successes; then duration). Success counts out of 16 (8 for sensitivity cells):\n")
    rows2 = []
    for k, v in sorted(T.items()):
        cnt = [sum(1 for r in v.values() if r["reached"] and not r["hard"] and abs(r["err_km"]) <= t) for t in TOLS_KM]; s = [r for r in v.values() if r["reached"] and not r["hard"]]
        rows2.append([k, len(v), sum(r["hard"] for r in v.values()), sum(not r["reached"] for r in v.values())] + cnt + [f"{med([r['dv'] for r in s]):.3f}", f"{med([r['days'] for r in s]):.1f}"])
    L.append(PR.table(rows2, ["arm@geometry", "flights", "hard failures", "not reached"] + [f"success at {int(t)} km" for t in TOLS_KM] + ["median dV of reached, no hard violation", "median days"]))
    # D1 / D6 under terminal-first success at each ASSUMED tolerance (paired seeds where both arms succeed)
    def okc(key, tol):
        return {x: r for x, r in T.get(key, {}).items() if r["reached"] and not r["hard"] and abs(r["err_km"]) <= tol}
    rows3 = []
    for tol in TOLS_KM:
        for tok in ("ODY", "G1"):
            l3 = okc(f"L3@{tok}", tol); rr = []
            for X in ("BODY", "FO", "N2I"):
                Xs = okc(f"{X}@ODY", tol); sd = sorted(set(l3) & set(Xs))
                rr.append(f"{X} {np.median([l3[x]['dv'] / Xs[x]['dv'] for x in sd]):.3f} (n={len(sd)})" if sd else f"{X} n/a (n=0)")
            ol, olr = okc(f"OL@{tok}", tol), okc(f"OLR@{tok}", tol); sd = sorted(set(ol) & set(olr))
            rows3.append([int(tol), tok, ", ".join(rr), f"{len(okc(f'OL@{tok}', tol))} / {len(okc(f'OLR@{tok}', tol))} / {len(okc(f'L3@{tok}', tol))}", (f"{np.median([olr[x]['dv'] - ol[x]['dv'] for x in sd]):+.3f} (n={len(sd)})" if sd else "n/a")])
    L.append("\nTerminal-first success at each ASSUMED tolerance, then dV among jointly successful seeds: median dV ratio L3 / X (n = seeds where both succeed); successes OL / OL-R / L3 out of 16; median dV of OL-R minus OL on common successes (m/s):\n\n" + PR.table(rows3, ["tolerance km (ASSUMED)", "geometry", "L3 / B-ODY, L3 / F-O, L3 / N2-I", "successes OL / OL-R / L3", "OL-R minus OL dV"]))
    L.append("""
What it changes (plain): a reached flight overshoots by less than one last-pass drop (Odyssey cells: 55 to 87 km at the median; G1 cells: 31 to 37 km), so a tolerance of 100 km is met by every reached flight except two F-O flights, 50 km separates the Odyssey arms (F-O 6 of 16, L3 9, L0 10, L1 12, OL-R 12, N2-I 13, L2 13, B-ODY 15 successes), and 10 or 25 km mostly scores the phase of the last drop, not the planning. Under a 50 km tolerance L3 at Odyssey would be a success in 9 of 16 flights, so the D1b dV ratios would rest on those 9 (table above); the direction (ratios below 1) is unchanged at every tolerance where pairs exist, but the sample shrinks and the conclusion becomes sensitive to an ASSUMED number. D6 and OL against OL-R: OL fails to reach the terminal condition (safety_q stops) in 4 of 16 flights at Odyssey and 1 at G1, OL-R in none, so terminal-first scoring (like the Revision 13 hard-violation rule) prefers OL-R; at 50 km OL-R succeeds in 12 of 16 against OL's 10 at Odyssey. The extra dV of OL-R over OL is again the price of reaching the terminal condition. A genuinely terminal-first controller (a terminal-apoapsis term or end-game trim in the DP) would be a design change and is not evaluated here. Note: the 'pass of the last burn' above is the target pass of the last decision, which can be the decision made at the closing orbit that the flight ends before executing (hence it can exceed the pass count).
""")
    pathlib.Path(a.out_md).write_text("\n".join(L)); print("wrote", a.out_md)


if __name__ == "__main__":
    main()
