"""Mission-success scoring (Revision 14, user): success = reaches the target apoapsis within the stated tolerance (ASSUMED 25 km below the target) AND no hard heating violation (Q_s on V_rel > 0.65); a failed mission has no
dV credit; then minimum dV among successes; duration is not a criterion (reported). Reads the summaries of the given phases and prints / writes a markdown table. Error = last closing apoapsis - 7,408 km (m), negative below.
Usage: python3 -I planner_terminal_score.py --planner-root R/planner --phases devT [--stageo-root S] [--tols 10,25,50] [--out F] [--json F]"""
import argparse, glob, json, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_analysis as A

TARGET = A.RA_TARGET


def collect(roots, roles, phases):
    out = {}
    for root in roots:
        root = pathlib.Path(root)
        for ph in phases:
            for f in sorted(glob.glob(str(root / f"summary_{ph}_*.csv"))):
                S = pd.read_csv(f); S = S[S.status == "ok"]
                for r in S.itertuples():
                    d = root / ph / r.run_id
                    if not (d / "passes.csv").exists():
                        d = A._run_dir(root, r.run_id)
                    if d is None or not (d / "passes.csv").exists():
                        continue
                    P = pd.read_csv(d / "passes.csv").sort_values("pass"); ra = P.ra_m.values.astype(float); fin = ra[np.isfinite(ra)]
                    qs = 0.5 * P.rho_peak_kg_m3.values * P.vrel_mps.values ** 3 / 1e4; hard = bool(np.nanmax(qs) > 0.65 or "safety_q" in str(r.stop_reason))
                    arm = str(r.arm).replace("-", ""); tok = roles.get(getattr(r, "geometry", "ODY"), "ODY"); seed = int(str(r.run_id).split("_pert_")[1].split("_")[1].lstrip("s"))
                    out.setdefault(f"{arm}@{tok}", {})[seed] = dict(err_km=(fin[-1] - TARGET) / 1e3, reached=bool(fin[-1] <= TARGET), hard=hard, dv=float(r.total_dv_mps), days=float(r.t_end_days), passes=len(P), stop=str(r.stop_reason),
                                                                    n_burns=int(r.burns), terminal_burns=0)
                    D = d / "replan.csv"
                    if D.exists() and D.stat().st_size > 2:
                        R = pd.read_csv(D); out[f"{arm}@{tok}"][seed]["terminal_burns"] = int(R.reason.astype(str).str.startswith("terminal").sum()); out[f"{arm}@{tok}"][seed]["terminal_dv"] = float(R[R.reason.astype(str).str.startswith("terminal")].executed_dv_mps.abs().sum())
    return out


def success(r, tol):
    return bool(r["reached"] and not r["hard"] and -tol <= r["err_km"] <= 0.0)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--planner-root", required=True); ap.add_argument("--phases", required=True); ap.add_argument("--stageo-root"); ap.add_argument("--tols", default="10,25,50"); ap.add_argument("--out"); ap.add_argument("--json"); a = ap.parse_args()
    roles = A.gate_roles(A.task_dir()); tols = [float(x) for x in a.tols.split(",")]; T = collect([pathlib.Path(a.planner_root)] + ([pathlib.Path(a.stageo_root)] if a.stageo_root else []), roles, a.phases.split(","))
    rows = ["| arm@geometry | flights | hard failures | not reached | " + " | ".join(f"success {int(t)} km" for t in tols) + " | median error (km) | median dV of successes at 25 km | median dV all | median days | terminal-aim burns (median) |", "|" + "---|" * (8 + len(tols))]
    for k, v in sorted(T.items()):
        s25 = [r for r in v.values() if success(r, 25.0)]
        rows.append(f"| {k} | {len(v)} | {sum(r['hard'] for r in v.values())} | {sum(not r['reached'] for r in v.values())} | " + " | ".join(f"{sum(success(r, t) for r in v.values())}" for t in tols) +
                    f" | {np.median([r['err_km'] for r in v.values()]):.1f} | {(np.median([r['dv'] for r in s25]) if s25 else float('nan')):.3f} | {np.median([r['dv'] for r in v.values()]):.3f} | {np.median([r['days'] for r in v.values()]):.1f} | {np.median([r['terminal_burns'] for r in v.values()]):.0f} |")
    txt = "\n".join(rows) + "\n\nPer-seed errors (km; * = hard failure):\n" + "\n".join(f"- {k}: " + ", ".join(f"s{x}:{r['err_km']:+.1f}{'*' if r['hard'] else ''}" for x, r in sorted(v.items())) for k, v in sorted(T.items()))
    print(txt)
    if a.out: pathlib.Path(a.out).write_text(txt)
    if a.json: pathlib.Path(a.json).write_text(json.dumps(T, indent=1, default=float))


if __name__ == "__main__":
    main()
