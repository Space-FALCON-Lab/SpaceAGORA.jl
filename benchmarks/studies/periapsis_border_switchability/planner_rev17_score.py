"""Revision 17b scoring (analysis only). Flights: summary_rev17_*.csv + <root>/rev17/<run_id>. Order (user): mission success (reached, final apoapsis error in [-25, 0] km ASSUMED, no hard violation: pass-max Q_s on V_rel > 0.65 or a
safety_q stop; soft 0.54 reported only), then dV among successes (dV to the last counted pass = ve ln(m0/m_n), and net |sum of signed logged dV| up to pass n), plus table-error diagnostics: realized offset (median over passes 5..n of
measured.csv y = ln(rho_mean / baseline fit)), drift update (replan.csv density_offset_yhat, drift_bias_per_pass, nu_cum), D8 deviation (planner_analysis.replan_metrics on rows up to pass n), hard and soft counts. Paired by seed.
Usage: python3 -I planner_rev17_score.py --planner-root R/planner --out F.md [--json F]"""
import argparse, glob, json, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_analysis as A
TARGET = A.RA_TARGET; VE = 220.0 * 9.80665; ARMS = ["OL", "OLR", "OLRT3"]; LEVELS = ["L0", "L1", "L2", "L3", "L4"]
LABEL = {"L0": "nominal", "L1": "noise 1.5", "L2": "offset +0.22", "L3": "offset -0.22", "L4": "offset +0.22 and noise 1.5"}


def flight(d, srow):
    P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra = P.ra_m.values.astype(float); hit = np.where(np.isfinite(ra) & (ra <= TARGET))[0]
    n = int(hit[0] + 1) if len(hit) else len(P); Pn = P.iloc[:n]; fin = ra[:n][np.isfinite(ra[:n])]; m = P.mass_kg.values.astype(float)
    qs = 0.5 * Pn.rho_peak_kg_m3.values * Pn.vrel_mps.values ** 3 / 1e4; hard = bool(np.nanmax(qs) > 0.65 or "safety_q" in str(srow.stop_reason)); reached = bool(len(hit))
    err = (fin[-1] - TARGET) / 1e3 if len(fin) else float("nan"); D = pd.read_csv(d / "decisions.csv"); Dn = D[D.target_pass <= n]; dv = Dn.delta_v_mps.fillna(0).values
    M = pd.read_csv(d / "measured.csv"); y = M.y.values.astype(float); R = pd.read_csv(d / "replan.csv"); Rn = R[R["pass"] <= n]
    out = dict(n=n, reached=reached, hard=hard, err=err, ok=bool(reached and not hard and -25.0 <= err <= 0.0), dv_sum=VE * np.log(m[0] / m[n - 1]), dv_net=abs(float(dv.sum())), dv_check=abs(float(np.abs(dv).sum()) - VE * np.log(m[0] / m[n - 1])),
               burns=int((np.abs(dv) > 0).sum()), soft=int((qs > 0.54).sum()), qmax=float(np.nanmax(qs)), days=float(P.time_s.iloc[n - 1] / 86400), stop=str(srow.stop_reason),
               offset=float(np.nanmedian(y[4:n])) if n > 5 else float("nan"), offset_first5=float(np.nanmean(y[:5])), y_sd=float(np.nanstd(y[:n])))
    if "density_offset_yhat" in Rn:
        yh = Rn.density_offset_yhat.values.astype(float); out.update(yhat_late=float(np.nanmedian(yh[20:])) if len(yh) > 21 else float("nan"), bias_abs_max=float(np.nanmax(np.abs(Rn.drift_bias_per_pass.values.astype(float)))), nu_cum_final=float(Rn.nu_cum.iloc[-1]))
    if "apriori_dv_mps" in Rn and Rn.apriori_dv_mps.notna().any():
        out.update({"d8_" + k: v for k, v in A.replan_metrics(Rn).items() if k in ("n_diverge", "dv_diff_abs_sum", "n_burn_replan", "n_burn_apriori")})
    return out


def med(x):
    x = [v for v in x if np.isfinite(v)]; return float(np.median(x)) if x else float("nan")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--planner-root", required=True); ap.add_argument("--out", required=True); ap.add_argument("--json"); a = ap.parse_args(); root = pathlib.Path(a.planner_root); F = {}
    for f in sorted(glob.glob(str(root / "summary_rev17_*.csv"))):
        S = pd.read_csv(f); S = S[S.status == "ok"]
        for r in S.itertuples():
            lv = r.disp; arm = r.arm; seed = int(r.run_id.split("_pert_")[1].split("_")[1].lstrip("s")); F[(arm, lv, seed)] = flight(root / "rev17" / r.run_id, r)
    L = ["# Results: Revision 17 (table error, open loop versus replanner), Odyssey, held seeds 13-16", "", "Pre-registered in the plan (Revisions 17, 17a, 17b). Scored in the user's order. Levels: " + "; ".join(f"{k} = {v}" for k, v in LABEL.items()) + ". Success = reached, final apoapsis error in [-25, 0] km (ASSUMED), no hard violation (Q_s on V_rel > 0.65 or safety stop). dV = dV to the last counted pass (sum |dv| and net |sum dv|). 4 seeds per cell: preliminary; the smallest sign-test p is 0.125.", "",
         "## 1. Success, hard violations, dV among successes", "", "| arm | level | flights | success | hard | not reached | passes Q_s>0.54 (median) | median sum dV of successes (m/s) | median net dV | median error (km) | median burns | median days |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for arm in ARMS:
        for lv in LEVELS:
            v = {s: F[(arm, lv, s)] for s in (13, 14, 15, 16) if (arm, lv, s) in F}
            if not v: continue
            ok = [r for r in v.values() if r["ok"]]
            L.append(f"| {arm} | {lv} | {len(v)} | {len(ok)}/{len(v)} | {sum(r['hard'] for r in v.values())} | {sum(not r['reached'] for r in v.values())} | {med([r['soft'] for r in v.values()]):.0f} | {med([r['dv_sum'] for r in ok]):.2f} | {med([r['dv_net'] for r in ok]):.2f} | {med([r['err'] for r in v.values()]):.1f} | {med([r['burns'] for r in v.values()]):.0f} | {med([r['days'] for r in v.values()]):.1f} |")
    L += ["", "Per flight (arm level seed: OK/FAIL [HARD], error km, sum dV, net dV, burns, soft passes, stop):"]
    L += [f"- {arm} {lv} s{s}: {'OK' if r['ok'] else 'FAIL'}{' HARD' if r['hard'] else ''}, {r['err']:+.1f}, {r['dv_sum']:.2f}, {r['dv_net']:.2f}, {r['burns']}, {r['soft']}, {r['stop']}" for (arm, lv, s), r in sorted(F.items(), key=lambda kv: (ARMS.index(kv[0][0]), LEVELS.index(kv[0][1]), kv[0][2]))]
    L += ["", "## 2. Table-error diagnostics", "", "Realized offset = median over passes 5..n of y = ln(rho_mean / baseline fit) (measured.csv). Expected: L0 about 0, L2/L4 about +0.22, L3 about -0.22 (the probe values at the first 40 periapsis points). Drift update columns are for the replanning arms (median of the estimator density_offset_yhat after pass 20; largest |drift bias| per pass; final cumulative level shift nu_cum in ln q). D8: orbits whose burn/no-burn decision differs from the a priori schedule, and sum |executed - a priori| (m/s).", "",
          "| arm | level | realized offset (median of flights) | y sd | yhat after pass 20 | max drift bias per pass | nu_cum final | D8 orbits differing | D8 sum abs dv diff (m/s) |", "|---|---|---|---|---|---|---|---|---|"]
    for arm in ARMS:
        for lv in LEVELS:
            v = [F[(arm, lv, s)] for s in (13, 14, 15, 16) if (arm, lv, s) in F]
            if v: L.append(f"| {arm} | {lv} | {med([r['offset'] for r in v]):+.3f} | {med([r['y_sd'] for r in v]):.3f} | {med([r.get('yhat_late', float('nan')) for r in v]):+.3f} | {med([r.get('bias_abs_max', float('nan')) for r in v]):.4f} | {med([r.get('nu_cum_final', float('nan')) for r in v]):+.2f} | {med([r.get('d8_n_diverge', float('nan')) for r in v]):.0f} | {med([r.get('d8_dv_diff_abs_sum', float('nan')) for r in v]):.2f} |")
    L += ["", "## 3. Paired comparisons by seed (dV on seeds where both members succeed; n shown)", "", "| level | comparison | n paired | median dV difference (m/s) | success A / B |", "|---|---|---|---|---|"]
    for lv in LEVELS:
        for A_, B_ in (("OLR", "OL"), ("OLRT3", "OLR"), ("OLRT3", "OL")):
            sd = [s for s in (13, 14, 15, 16) if (A_, lv, s) in F and (B_, lv, s) in F and F[(A_, lv, s)]["ok"] and F[(B_, lv, s)]["ok"]]
            sa = sum(F[(A_, lv, s)]["ok"] for s in (13, 14, 15, 16) if (A_, lv, s) in F); sb = sum(F[(B_, lv, s)]["ok"] for s in (13, 14, 15, 16) if (B_, lv, s) in F)
            L.append(f"| {lv} | {A_} - {B_} | {len(sd)} | {med([F[(A_, lv, s)]['dv_sum'] - F[(B_, lv, s)]['dv_sum'] for s in sd]):+.2f} | {sa} / {sb} |")
    L += ["", "Level effect within an arm (level minus L0, same seed, both succeed):", "", "| arm | level | n paired | median sum dV difference (m/s) | median net dV difference | burn difference |", "|---|---|---|---|---|---|"]
    for arm in ARMS:
        for lv in LEVELS[1:]:
            sd = [s for s in (13, 14, 15, 16) if (arm, lv, s) in F and (arm, "L0", s) in F and F[(arm, lv, s)]["ok"] and F[(arm, "L0", s)]["ok"]]
            L.append(f"| {arm} | {lv} | {len(sd)} | {med([F[(arm, lv, s)]['dv_sum'] - F[(arm, 'L0', s)]['dv_sum'] for s in sd]):+.2f} | {med([F[(arm, lv, s)]['dv_net'] - F[(arm, 'L0', s)]['dv_net'] for s in sd]):+.2f} | {med([F[(arm, lv, s)]['burns'] - F[(arm, 'L0', s)]['burns'] for s in sd]):+.0f} |")
    L += ["", "Survivor-bias note: dV among successes is conditioned on success. Where an arm fails more often at a level, its surviving flights are the seeds that were easy, so its dV is biased low and a comparison with an arm that succeeded on all seeds is not like for like; compare only the paired rows (n shown) and read the success counts first. With n of 4 or fewer no ranking is claimed.", ""]
    # predictions
    def cnt(arm, lv, key): return sum(F[(arm, lv, s)][key] for s in (13, 14, 15, 16) if (arm, lv, s) in F)
    def nok(arm, lv): return cnt(arm, lv, "ok")
    P1 = all(nok("OLRT3", lv) >= 3 for lv in LEVELS); P2 = all(nok("OLRT3", lv) > nok("OL", lv) for lv in LEVELS)
    P3 = {lv: (med([F[(arm, lv, s)]["offset"] for arm in ("OL",) for s in (13, 14, 15, 16) if (arm, lv, s) in F]), med([F[(arm, lv, s)].get("yhat_late", np.nan) for arm in ("OLR", "OLRT3") for s in (13, 14, 15, 16) if (arm, lv, s) in F])) for lv in LEVELS}
    P4a = all(cnt(arm, "L0", "hard") == 0 for arm in ("OLR", "OLRT3")) and all(cnt(arm, lv, "hard") <= 1 for arm in ("OLR", "OLRT3") for lv in LEVELS)
    P4b = (any(cnt("OL", lv, "hard") >= 1 for lv in ("L0", "L2", "L4"))) and cnt("OL", "L2", "hard") >= cnt("OL", "L3", "hard") and cnt("OL", "L4", "hard") >= cnt("OL", "L3", "hard")
    P6 = {arm: nok(arm, "L0") - nok(arm, "L1") for arm in ("OLR", "OLRT3")}
    L += ["## 4. Pre-registered predictions against the result (mechanical readings; refutation rules in results_rev17_design.md and plan Revision 17b)", "",
          f"- P1 (OLRT3 success at least 3/4 at every level): {'held' if P1 else 'REFUTED'}; successes by level: " + ", ".join(f"{lv} {nok('OLRT3', lv)}/4" for lv in LEVELS) + ".",
          f"- P2 (OLRT3 succeeds more often than OL at every level; confounded by the terminal aim): {'held' if P2 else 'REFUTED'}; OL: " + ", ".join(f"{lv} {nok('OL', lv)}/4" for lv in LEVELS) + ".",
          "- P3 (realized offset): OL median y by level " + ", ".join(f"{lv} {P3[lv][0]:+.3f}" for lv in LEVELS) + "; replanner yhat after pass 20 " + ", ".join(f"{lv} {P3[lv][1]:+.3f}" for lv in LEVELS) + ". Reading: the realized offset y is the same for all arms (it is the measured mean-field density against the baseline fit); the test is whether the replanner's estimate tracks it (yhat close to y) while the open loop schedule does not use it.",
          f"- P4 (rewritten in Revision 17b), part a (OLR/OLRT3: 0 hard at L0, at most 1 of 4 at any level): {'held' if P4a else 'REFUTED'}; hard flights OLR " + ", ".join(f"{lv} {cnt('OLR', lv, 'hard')}" for lv in LEVELS) + "; OLRT3 " + ", ".join(f"{lv} {cnt('OLRT3', lv, 'hard')}" for lv in LEVELS) + f". Part b (OL has a hard flight at L0, L2 or L4 and at least as many at L2 and L4 as at L3): {'held' if P4b else 'REFUTED'}; OL " + ", ".join(f"{lv} {cnt('OL', lv, 'hard')}" for lv in LEVELS) + ".",
          f"- P6 (success of OLR/OLRT3 falls by at most one flight from L0 to L1): OLR {P6['OLR']}, OLRT3 {P6['OLRT3']} flights lost: {'held' if max(P6.values()) <= 1 else 'REFUTED'}.",
          "- P5 (OLRT3 median dV within 0.8 m/s of its L0 median at every level): " + ", ".join(f"{lv} {med([F[('OLRT3', lv, s)]['dv_sum'] for s in (13, 14, 15, 16) if ('OLRT3', lv, s) in F and F[('OLRT3', lv, s)]['ok']]):.2f}" for lv in LEVELS) + " m/s (read against the 0.8 band).", ""]
    # L0 reproduction of the held flights
    rep = []
    for arm in ("OL", "OLR"):
        for s in (13, 14, 15, 16):
            h = root / "runs" / f"ODY_PL_{ {'OL': 'OL', 'OLR': 'OLR'}[arm] }_pert_ODY_s{s}"; n_ = root / "rev17" / f"ODY_PL_{arm}_pert_ODY_s{s}_L0"
            if (h / "passes.csv").exists() and (n_ / "passes.csv").exists():
                a_ = pd.read_csv(h / "passes.csv"); b_ = pd.read_csv(n_ / "passes.csv"); sh = [c for c in a_.columns if c in b_.columns]
                same = a_.shape[0] == b_.shape[0] and all((a_[c].astype(str) == b_[c].astype(str)).all() for c in sh); rep.append((arm, s, same, a_.shape[0], b_.shape[0]))
    L += ["## 5. L0 check against the existing held flights (same arm, seed, geometry; informal)", ""] + [f"- {arm} s{s}: {'bit-identical on all shared passes.csv columns' if same else 'DIFFERENT'} (passes {na} held, {nb} now)" for arm, s, same, na, nb in rep]
    txt = "\n".join(L); print(txt); pathlib.Path(a.out).write_text(txt)
    if a.json: pathlib.Path(a.json).write_text(json.dumps({f"{k[0]}|{k[1]}|{k[2]}": v for k, v in F.items()}, indent=1, default=float))


if __name__ == "__main__":
    main()
