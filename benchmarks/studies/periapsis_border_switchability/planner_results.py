"""Assemble results.md from the held and sensitivity flights (planner_analysis.report plus tables). Reports the pre-registered criteria plainly; no interpretation beyond them. Python only.
Usage: python3 -I planner_results.py --planner-root R/planner --stageo-root S/stageo [--dev-diagnosis F] [--calibration F] --out results.md [--json F]
sens-phase runs live in <planner-root>/sens; they are linked into runs/ (symlinks) so the shared loader finds them. The drag ledger is the analytic one (geometry_rank channels at the calibrated context): the dV-equivalent
of each natural-drift channel (positive = pushes the periapsis down, needs burns up), per geometry; the measured split of the flown paths is not computed here."""
import argparse, glob, json, os, pathlib, sys
import numpy as np, pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_analysis as A
import planner_lib as PLB


def link_sens(proot):
    runs = proot / "runs"; runs.mkdir(exist_ok=True)
    for d in glob.glob(str(proot / "sens" / "*")):
        t = runs / pathlib.Path(d).name
        if not t.exists():
            os.symlink(d, t)


def fmt(b, nd=3):
    if b is None or not isinstance(b, dict) or "point" not in b:
        return "n/a"
    f = lambda x: "nan" if x is None or not np.isfinite(x) else f"{x:.{nd}f}"
    return f"{f(b['point'])} [{f(b['lo'])}, {f(b['hi'])}]"


def med(v, key):
    x = [r[key] for r in v.values() if np.isfinite(r.get(key, np.nan))]; return float(np.median(x)) if x else float("nan")


def drag_ledger(cfg):
    import planner_predictions as PP
    GR, G3A = PP._imports(); cfg3 = G3A.load3(); ctx = PP.make_ctx_cal(cfg3, cfg); geoms = PP.geometries(cfg, cfg3); out = {}
    for tok in ("ODY", "G1"):
        g = geoms[tok]; field, _ = G3A.meanfield_for(cfg3, ctx, None, g["g"]); b = ctx.cfg["band"]; Hk = None if field is None else field["H"]
        q = GR.solve_q_req(ctx, g["cand"], b["d_ref_days"], Hk); S = GR.schedule(ctx, g["cand"], q, Hk); d = GR.channels(ctx, g["cand"], S, field)
        n = S["n"]; w = np.zeros(n); w[1:] = S["H"][1:] / (ctx.ratio * np.array([GR.lever(ctx, S["a"][k - 1], S["e"][k - 1]) for k in range(1, n)]))
        led = {}
        for ch, inc in d.items():
            c = np.asarray(inc) * w[1:len(inc) + 1]; led[ch] = dict(down_dv=float(c[c > 0].sum()), up_dv=float(-c[c < 0].sum()), net_dv=float(c.sum()))
        led["total_net"] = float(sum(v["net_dv"] for v in led.values())); out[tok] = led
    return out


def table(rows, head):
    s = "| " + " | ".join(head) + " |\n|" + "---|" * len(head) + "\n"
    for r in rows:
        s += "| " + " | ".join(str(x) for x in r) + " |\n"
    return s


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--planner-root", required=True); ap.add_argument("--stageo-root"); ap.add_argument("--dev-diagnosis"); ap.add_argument("--calibration"); ap.add_argument("--out", required=True)
    ap.add_argument("--json"); a = ap.parse_args()
    cfg = PLB.load_cfg(); td = A.task_dir(); roles = A.gate_roles(td); proot = pathlib.Path(a.planner_root); link_sens(proot)
    preds = json.loads((td / cfg["planner"]["predictions"]["path"]).read_text()); tabs = json.loads((td / cfg["planner"]["ground_table"]["path"]).read_text())
    roots = [proot] + ([pathlib.Path(a.stageo_root)] if a.stageo_root else [])
    runs = A.collect(roots, roles, preds, tabs, ("held", "sens")); rep = A.report(runs, preds)
    # sensitivity arms against their penalty-0.01 counterparts
    sens = {}
    for z, base in (("L3Z", "L3"), ("OLRZ", "OLR")):
        for tok in ("ODY", "G1"):
            Z, Bq = runs.get(f"{z}@{tok}"), runs.get(f"{base}@{tok}")
            if Z and Bq:
                rng = A.rng_of(); s = {}
                for k, f in (("dv", lambda x, y: x["dv"] - y["dv"]), ("burns", lambda x, y: x["burns"] - y["burns"]), ("v032", lambda x, y: x["v032"] - y["v032"])):
                    seeds, dd = A.pair(Z, Bq, f); s[k] = A.boot(dd, np.median, rng); s["n"] = len(seeds)
                sens[f"{z}@{tok} - {base}@{tok}"] = s
    led = drag_ledger(cfg)
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(dict(report=rep, sensitivity=sens, drag_ledger=led), indent=1, default=float))
    L = ["# Results: planner arms, held seeds 13 to 28 (paired), sensitivity arms\n",
         "Pre-registered criteria only, reported plainly. Every ASSUMED value is a named field of planner_config.toml; calibrated values (sigma_p 0.0452, planning D_ref 70.4 d, q_cap 0.2503) are MEASURED / DERIVED from dev seeds 1 to 8 (calibration_freeze.md). The never-exceed 0.65 is a SOURCED value used as a terminator by assumption (U12); the burn penalty 0.01 m/s is ASSUMED. Intervals are 90% percentile bootstraps over seeds; medians unless stated. One Julia flight per seed and arm; walls are contended.\n"]
    L.append("## Flights per arm\n\n" + table([[k, len(v), sum(not r['reached'] for r in v.values()), ", ".join(sorted({r['stop_reason'] for r in v.values()}))] for k, v in sorted(runs.items())], ["arm@geometry", "flights", "not reached", "stop reasons"]))
    L.append("## dV, time, violations (medians over seeds)\n\n" + table([[k, f"{med(v,'dv'):.3f}", f"{med(v,'burns'):.0f}", f"{med(v,'days'):.1f}", f"{med(v,'n_passes'):.0f}", sum(r['v032'] for r in v.values()), sum(r['v054'] for r in v.values()), sum(r['v065'] for r in v.values())] for k, v in sorted(runs.items())],
                                                                      ["arm@geometry", "dV m/s", "burns", "days", "passes", "passes Q_s>0.32 (total)", ">0.54", ">0.65"]))
    for name in ("D1a", "D1b"):
        r = rep.get(name)
        L.append(f"## {name}\n\n" + (f"PASS = {r['pass']} (C1 {r['c1']}, C2 {r['c2']}). Median dV ratio L3 / X, upper bound must be < 1:\n\n" + table([[x, fmt(v)] for x, v in r["ratios"].items()], ["X", "ratio [90%]"]) if r else "not computable (arms missing)\n"))
    at = rep.get("attribution", {}).get("components", {})
    L.append("## D2 and D3 (attribution, paired means, m/s; positive = saving)\n\n")
    for X, c in at.items():
        if "C" in c:
            L.append(f"First rung {X} (n = {c['n_seeds']} seeds; the ladder telescopes exactly per seed):\n\n" + table([[k, fmt(c[k])] for k in ("C", "N", "P_ff", "P_t", "P_g", "P_g_main", "P_g_interaction", "P", "total")], ["component", "mean [90%]"]) + f"\nP share of total: {c['P_share_of_total']:.3f}\n")
        else:
            L.append(f"First rung {X}: {c}\n")
    L.append(f"D2: {json.dumps(rep.get('D2'), default=float)}\n\nD3: components above.\n")
    L.append("## D5\n\n" + (f"B-ODY at G1 vs Odyssey: diff {fmt(rep['D5']['diff_median'])}, ratio {fmt(rep['D5']['ratio_median'])}\n" if "D5" in rep else "n/a\n"))
    for tok in ("ODY", "G1"):
        for comp in ("OL", "OLR"):
            r = rep.get(f"D6:{comp}@{tok}")
            if r:
                L.append(f"## D6 {comp} - L3 at {tok}\n\nverdict {r['verdict']}; dV diff {fmt(r['dv'])}; passes>0.32 diff {fmt(r['v032'])}; >0.65 diff {fmt(r['v065'])}\n")
    L.append("## D7 and the OL versus OL-R table (OL minus OL-R, paired medians)\n")
    for tok in ("ODY", "G1"):
        r = rep.get(f"D7:{tok}")
        if r:
            L.append(f"\nGeometry {tok}:\n\n" + table([[v["label"], v["n_seeds"], fmt(v)] for k, v in r.items()], ["metric", "seeds", "OL - OL-R [90%]"]))
            ol, olr = runs.get(f"OL@{tok}", {}), runs.get(f"OLR@{tok}", {})
            L.append("\n" + table([[nm, f"{med(v,'track_rms_m'):.0f}", f"{med(v,'track_max_m'):.0f}", f"{med(v,'in_tube'):.3f}", f"{med(v,'final_apo_err_m'):.0f}", sum(r_['v032'] for r_ in v.values()), sum(r_['v065'] for r_ in v.values()), f"{med(v,'dv'):.3f}", f"{med(v,'burns'):.0f}"] for nm, v in (("OL", ol), ("OL-R", olr))],
                                  ["arm", "track rms m", "track max m", "in tube", "final apo err m", "Q_s>0.32 (total)", ">0.65", "dV m/s", "burns"]))
    L.append("\n## D8 replan versus a priori (medians over seeds, per arm and geometry)\n")
    for k, r in sorted(rep.items()):
        if k.startswith("D8:") and "n_diverge" in r:
            L.append(f"\n{k} ({r['n_flights']} flights, {r['n_no_divergence']} never diverged):\n\n" + table([[m, fmt(r[m], 4)] for m in A.D8_METRICS], ["metric", "median [90%]"]))
        elif k.startswith("D8:"):
            L.append(f"\n{k}: {r}\n")
    L.append("\n## Criteria C1, C2, C4, C4abs, C5 to C7 (per arm and geometry)\n\n")
    rows = []
    for k in sorted(runs):
        if k.split("@")[0] in ("FO", "N2I"):
            continue
        c1, c2, c4, c4a, c5, c6, c7 = (rep.get(f"{n}:{k}") for n in ("C1", "C2", "C4", "C4abs", "C5", "C6", "C7"))
        rows.append([k, c1 and c1.get("pass"), c2 and f"{c2['pass']} (n>0.54: {c2['count_gt054']}, exp {c2['expected']:.1f})", c4 and f"{c4['pass']} (med {c4['median_dv']:.2f}, upper {c4['upper']:.2f})",
                     c4a and f"{c4a['pass']} (diff {c4a['difference']:+.2f})", c5 and c5["share_pass"], c6 and c6.get("status"), c7 and f"{c7['pass']} (q99 {c7['q99']:.3f})"])
    L.append(table(rows, ["arm@geometry", "C1", "C2", "C4 (upper bound only unless floors attached)", "C4abs (reported)", "C5 share passing", "C6", "C7"]))
    L.append("\n## Sensitivity: penalty 0 minus penalty 0.01 (paired medians)\n\n" + (table([[k, v["n"], fmt(v["dv"]), fmt(v["burns"], 1), fmt(v["v032"], 1)] for k, v in sens.items()], ["pair", "seeds", "dV diff m/s", "burns diff", "passes>0.32 diff"]) if sens else "arms not flown or incomplete\n"))
    L.append("\n## Drag ledger (analytic, dV-equivalent m/s per main phase, model path at the calibrated context)\n\nDrag cannot be switched off: it is the mission. 'down' = burns up the channel would need, 'up' = burns down, net = signed. The drag descent is given to every keeper for free; it is not a saving relative to the baselines.\n")
    for tok, l in led.items():
        L.append(f"\n{tok}:\n\n" + table([[ch, f"{v['down_dv']:.3f}", f"{v['up_dv']:.3f}", f"{v['net_dv']:.3f}"] for ch, v in l.items() if ch != "total_net"], ["channel", "down", "up", "net"]) + f"\ntotal net {l['total_net']:.3f}\n")
    if a.calibration and pathlib.Path(a.calibration).exists():
        cal = json.loads(pathlib.Path(a.calibration).read_text()); L.append(f"\n## Calibration (dev)\n\nsigma_p measured {cal['sigma_p_measured_pooled']:.4f}; D_ref median {cal['d_ref_median']:.2f} d; per-arm e_p sd: " + ", ".join(f"{k} {v['sd']:.3f}" for k, v in cal["e_p"].items()) + "\n")
    if a.dev_diagnosis and pathlib.Path(a.dev_diagnosis).exists():
        dd = json.loads(pathlib.Path(a.dev_diagnosis).read_text()); m = dd["medians"]; t = dd["totals"]
        L.append("\n## POST HOC diagnosis: why L3 cost more than L1 at Odyssey in dev (not pre-registered; flagged)\n\n"
                 f"Dev, seeds 1 to 8, from replan.csv (planner_dev_diagnosis.py). A burn is 'forced' when the one-pass prediction z0 lies outside the full tube (then L1 burns too) and 'anticipatory' when z0 is inside it (the burn comes from the planned horizon: reserve-shrunk tube for later steps, drift feed-forward from the ground table, early-timing weights).\n\n"
                 + table([[arm, f"{t[arm+'_sum']['burns']:.0f}", f"{t[arm+'_sum']['dv']:.2f}", f"{t[arm+'_sum']['forced_burns']:.0f}", f"{t[arm+'_sum']['forced_dv']:.2f}", f"{t[arm+'_sum']['antic_burns']:.0f}", f"{t[arm+'_sum']['antic_dv']:.2f}", f"{t[arm+'_sum']['down_dv']:.2f}", f"{t[arm+'_sum']['sign_reversals']:.0f}"] for arm in ("L1", "L3")],
                         ["arm", "burns (8 seeds)", "dV", "forced burns", "forced dV", "anticipatory burns", "anticipatory dV", "burn-down dV (adds heat)", "sign reversals"])
                 + "\nReading (hypothesis consistent with these counts, not tested by a separate experiment): L1 acts only when the next pass leaves the tube (all its burns are forced, 16 sign reversals in 8 flights). L3 has the same forced burns plus anticipatory ones; its extra dV is mostly anticipatory burns, over half of them burn-down (heat-adding) burns, and the sign reversals are concentrated in two seeds (2 and 5), i.e. the planner steps toward the lower edge of the reserve-shrunk tube and then reverses when the measured next pass differs from the table's slow path. The planner design was not changed.\n")
    pathlib.Path(a.out).write_text("\n".join(L))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
