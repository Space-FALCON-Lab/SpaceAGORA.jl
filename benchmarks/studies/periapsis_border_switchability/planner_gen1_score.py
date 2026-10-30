"""Revision 16b phase 1 scoring (one geometry, seeds 13-16). Success = reached, final apoapsis error in [-25, 0] km (ASSUMED tolerance), no hard violation (Q_s on V_rel > 0.65 W/cm2).
dV among successes: net = |sum of signed executed dV|, sum = sum |dV| (replan.csv executed_dv_mps; cross-checked against the summary total).
Usage: python3 -I planner_gen1_score.py <planner_root> [out.md]"""
import glob, pathlib, sys
import numpy as np, pandas as pd
TARGET = 7.408e6; VE = 220.0 * 9.80665   # Isp 220 s: DERIVED from the flights (ve = 2157.463 m/s reproduces the summary dV of flights with no post-final-pass burn)
root = pathlib.Path(sys.argv[1]); rows = []
for f in sorted(glob.glob(str(root / "summary_gen1_*.csv"))):
    S = pd.read_csv(f); S = S[S.status == "ok"]
    for r in S.itertuples():
        d = root / "gen1" / r.run_id; P = pd.read_csv(d / "passes.csv").sort_values("pass"); ra = P.ra_m.values.astype(float); fin = ra[np.isfinite(ra)]
        qs = 0.5 * P.rho_peak_kg_m3.values * P.vrel_mps.values ** 3 / 1e4; hard = bool(np.nanmax(qs) > 0.65 or "safety_q" in str(r.stop_reason))
        D = pd.read_csv(d / "decisions.csv"); n = len(P); dv = D[D.target_pass <= n].delta_v_mps.fillna(0).values; err = (fin[-1] - TARGET) / 1e3
        R = pd.read_csv(d / "replan.csv"); m = P.mass_kg.values; dv_thru = VE * np.log(m[0] / m[-1])
        rows.append(dict(arm=r.arm, seed=int(r.run_id.split("_s")[-1]), err=err, hard=hard, reached=bool(fin[-1] <= TARGET), net=abs(dv.sum()), sm=dv_thru, flown=float(r.total_dv_mps), chk=abs(np.abs(dv).sum() - dv_thru),
                         burns=int(r.burns), days=float(r.t_end_days), term=int(R.reason.astype(str).str.startswith("terminal").sum())))
T = pd.DataFrame(rows); T["ok"] = T.reached & ~T.hard & (T.err >= -25) & (T.err <= 0)
out = ["| arm | flights | success | hard | not reached | median net dV (m/s) | median sum dV (m/s) | median final error (km) | median burns | median days |", "|---|---|---|---|---|---|---|---|---|---|"]
for a in ["BODYT", "N2IT", "OLRT3", "OLRT5", "OLRT35"]:
    g = T[T.arm == a]; s = g[g.ok]
    if len(g): out.append(f"| {a} | {len(g)} | {len(s)}/{len(g)} | {int(g.hard.sum())} | {int((~g.reached).sum())} | {s.net.median() if len(s) else float('nan'):.2f} | {s.sm.median() if len(s) else float('nan'):.2f} | {g.err.median():.1f} | {g.burns.median():.0f} | {g.days.median():.1f} |")
out.append("\nPer flight (arm seed: ok, error km, net dV, sum dV, burns, terminal burns):")
out += [f"- {r.arm} s{r.seed}: {'OK' if r.ok else 'FAIL'}{' HARD' if r.hard else ''}, {r.err:+.1f}, {r.net:.2f}, {r.sm:.2f} (summary total incl. post-target burn {r.flown:.2f}), {r.burns}, {r.term}" for r in T.sort_values(["arm", "seed"]).itertuples()]
out.append(f"\nSum dV = dV delivered through the closing apoapsis that first reached the target (rocket equation on the mass at the last pass). Net dV = |sum of signed commanded dV| for decisions up to that pass. The summary total_dv_mps adds a burn commanded at that final apoapsis and cut off when the run stops (see Revision 17 note); cross-check max |sum|commanded| - delivered| = {T.chk.max():.1e} m/s.")
txt = "\n".join(out); print(txt)
if len(sys.argv) > 2: pathlib.Path(sys.argv[2]).write_text(txt)
