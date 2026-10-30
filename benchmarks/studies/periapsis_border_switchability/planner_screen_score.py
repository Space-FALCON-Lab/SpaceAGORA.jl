"""Revision 15 screen / full-dev scoring (user order): mission success (reached, error in [-25, 0] km, no hard violation), then dV among successes (also net |sum dv| and sum |dv|, final error), then burns and sign
reversals. Baseline = the unmodified OLRT dev flights (devT). Usage: python3 -I planner_screen_score.py --planner-root R/planner [--extra-root R2/planner] --seeds 1,2,3,5,7 --geometry ODY [--out F]"""
import argparse, glob, json, pathlib, sys
import numpy as np, pandas as pd
TARGET = 7.408e6


def flight(d):
    P = pd.read_csv(d / "passes.csv").sort_values("pass"); R = pd.read_csv(d / "replan.csv"); ex = R.executed_dv_mps.values; b = np.abs(ex) >= 1e-3; sg = np.sign(ex[b]); ra = P.ra_m.dropna().values
    qs = 0.5 * P.rho_peak_kg_m3.values * P.vrel_mps.values ** 3 / 1e4; err = (ra[-1] - TARGET) / 1e3
    hard = bool(np.nanmax(qs) > 0.65); reached = bool(ra[-1] <= TARGET)
    return dict(err=err, reached=reached, hard=hard, ok=bool(reached and not hard and -25 <= err <= 0), sumdv=float(np.abs(ex).sum()), netdv=float(abs(ex.sum())), burns=int(b.sum()), rev=int((np.diff(sg) != 0).sum()) if b.sum() > 1 else 0, days=float(P.time_s.iloc[-1] / 86400), passes=len(P))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--planner-root", required=True); ap.add_argument("--extra-root", action="append", default=[]); ap.add_argument("--seeds", default="1,2,3,5,7"); ap.add_argument("--geometry", default="ODY"); ap.add_argument("--out"); a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]; roots = [pathlib.Path(a.planner_root)] + [pathlib.Path(x) for x in a.extra_root]; g = "ODY" if a.geometry == "ODY" else "G5863"; res = {}
    for arm in ("OLRT", "OLRT1", "OLRT2", "OLRT3", "OLRT4", "OLRT5", "OLRT35"):
        for s in seeds:
            for root in roots:
                for sub in ("devT", "devR"):
                    d = root / sub / f"{g}_PL_{arm}_pert_ODY_s{s}"
                    if (d / "passes.csv").exists():
                        res.setdefault(arm, {})[s] = flight(d)
    rows = ["| variant | flights | successes | mean dV, successes | mean dV, all | net |sum dv| (mean) | burns (mean) | reversals (mean) | final error median (km) | days (mean) | seeds with hard failure |", "|" + "---|" * 11]
    for arm, v in res.items():
        ok = [r for r in v.values() if r["ok"]]; f = lambda k, S=v.values(): np.mean([r[k] for r in S])
        rows.append(f"| {arm} | {len(v)} | {len(ok)} | {np.mean([r['sumdv'] for r in ok]) if ok else float('nan'):.2f} | {f('sumdv'):.2f} | {f('netdv'):.2f} | {f('burns'):.1f} | {f('rev'):.1f} | {np.median([r['err'] for r in v.values()]):.1f} | {f('days'):.1f} | {[s for s, r in v.items() if r['hard']]} |")
    rows.append("\nPer seed (sum dV, burns, reversals, error km):\n")
    for arm, v in res.items():
        rows.append(f"- {arm}: " + "; ".join(f"s{s}: {r['sumdv']:.2f}, {r['burns']}, {r['rev']}, {r['err']:+.1f}{'' if r['ok'] else ' FAIL'}" for s, r in sorted(v.items())))
    txt = "\n".join(rows); print(txt)
    if a.out: pathlib.Path(a.out).write_text(txt)


if __name__ == "__main__":
    main()
