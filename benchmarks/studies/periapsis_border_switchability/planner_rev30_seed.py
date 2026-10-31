"""Revision 30.3 (rev30seed) and the one-page status table. Analysis only; Python, read-only on flight data.
rev30seed: OLRT3-W3L5 (G0528, v3 + L1), OLRT3-W3E5 (G1230, v3), OLRT3-W4L5 (ODY, v4 + L1), 5 km aim, EI 150, seeds 13-20, plus the regression flight. G0528 history-term shift (rev26val2 EI 150 pooled, rev23geo EI 2000).
Status table rows: per class/geometry the most recent consistent setting; ALL requirements per flight = success [-25, 0] km, no hard, window (own t*), position < 1.0 deg with L1 (the corridor is REPORTED, Revision 20d item 3; all_corr also requires it) (offline for arms without L1). Wilson 95% intervals.
Usage: python3 -I planner_rev30_seed.py --out-dir DIR"""
import argparse, json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev26_val2 as V
import planner_rev28_t as T28
import planner_rev30_g0528 as G
import planner_rev21_position as PS
import planner_rev20_score as SC


def wilson(k, n, z=1.96):
    if n == 0: return (float("nan"), float("nan"))
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)); return ((c - h) / d, (c + h) / d)


def boot_unpaired(a, b, nb=10000, seed=20261030):
    rng = np.random.default_rng(seed); a = np.asarray(a, float); b = np.asarray(b, float)
    bs = np.median(rng.choice(a, (nb, len(a))), axis=1) - np.median(rng.choice(b, (nb, len(b))), axis=1); return float(np.median(a) - np.median(b)), [float(x) for x in np.percentile(bs, [2.5, 97.5])]


def cell(fs):
    n = len(fs); l1 = [f["pos"]["l1"] if f["pos"] else float("inf") for f in fs]; raw = [f["pos"]["raw"] if f["pos"] else float("inf") for f in fs]
    return dict(n=n, ok=int(sum(f["ok"] for f in fs)), hard=int(sum(f["hard"] for f in fs)), win=int(sum(f["window"] for f in fs)), corr=int(sum(f["corr"] for f in fs)), dv_med=float(np.median([f["dv_sum"] for f in fs])), dv=[float(f["dv_sum"]) for f in fs],
                raw_med=float(np.median(raw)), l1_med=float(np.median(l1)), raw1=int(sum(x < 1 for x in raw)), l11=int(sum(x < 1 for x in l1)), all=int(sum(bool(f["ok"] and not f["hard"] and f["window"] and f["pos"] and f["pos"]["l1"] < 1.0) for f in fs)), all_corr=int(sum(f["all"] for f in fs)),
                wilson=wilson(int(sum(bool(f["ok"] and not f["hard"] and f["window"] and f["pos"] and f["pos"]["l1"] < 1.0) for f in fs)), n))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    B = {b: V.load(b, None) for b in ("rev30seed", "rev29aim", "rev28t", "rev27low", "rev27ei", "rev26val2", "rev23geo", "rev21v", "rev20val")}
    sel = lambda blk, arm, geom, lvl="L0", seeds=None: [f for k, f in sorted(B[blk].items(), key=lambda t: t[0][3]) if k[0] == arm and k[1] == geom and k[2] == lvl and (seeds is None or k[3] in seeds)]
    rows = {}
    spec = [
        ("Odyssey 93.5 deg, 7,408 km", "rev30seed", "OLRT3W4L5", "ODY", "v4 + L1", 150, 5.0, "nominal"),
        ("Odyssey 93.5 deg, 7,408 km (v3 alternative)", "rev29aim", "OLRT3W3E5", "ODY", "v3", 150, 5.0, "nominal"),
        ("G0528 63.4 deg, 7,408 km", "rev30seed", "OLRT3W3L5", "G0528", "v3 + L1", 150, 5.0, "nominal"),
        ("G1230 116.6 deg, 1,000 km circle", "rev30seed", "OLRT3W3E5", "G1230", "v3", 150, 5.0, "nominal"),
        ("G1038 93.5 deg, 500 km circle", "rev29aim", "OLRT3W3E5", "G1038", "v3", 150, 5.0, "nominal"),
        ("G0486 63.4 deg, 2,000 km circle", "rev27low", "OLRT3W3E", "G0486", "v3", 150, 12.5, "nominal"),
        ("G1 116.6 deg, 7,408 km", "rev26val2", "OLRT3W4L", "G1", "v4 + L1", 150, 12.5, "nominal"),
        ("G0975M 93.5 deg, omega0 30, 7,408 km (mean field)", "rev28t", "OLRT3W3E", "G0975M", "v3", 150, 12.5, "nominal"),
        ("G0975M, default t*", "rev28t", "OLRT3W3T", "G0975M", "v3", 150, 12.5, "fast end + 6 d"),
        ("G0063M 25.2 deg, 7,408 km (mean field)", "rev23geo", "OLRT3W3", "G0063M", "v3", 2000, 12.5, "nominal"),
        ("G0890 93.5 deg, omega0 90, 7,408 km", "rev23geo", "OLRT3W3", "G0890", "v3", 2000, 12.5, "nominal"),
    ]
    for lab, blk, arm, geom, guard, ei, aim, tst in spec:
        fs = sel(blk, arm, geom); c = cell(fs); c.update(label=lab, block=blk, arm=arm, guard=guard, ei=ei, aim=aim, tstar_rule=tst); rows[lab] = c
    # older / alternative settings for the footnote
    alt = {}
    for lab, blk, arm, geom, lv in (("Odyssey v3 S-fixed EI 2000 (rev21v, L0-L4 s13-20)", "rev21v", "OLRT3W3", "ODY", None), ("Odyssey S-free EI 2000 (rev21v)", "rev21v", "OLRT3WF3", "ODY", "L0"), ("Odyssey v4+L1 EI 150, 12.5 km (rev26val2 L0, L2)", "rev26val2", "OLRT3W4L", "ODY", None), ("G0528 v3 EI 2000 (rev23geo)", "rev23geo", "OLRT3W3", "G0528", "L0"),
                                  ("G1230 v3 EI 150, 12.5 km (rev27low s13-16)", "rev27low", "OLRT3W3E", "G1230", "L0"), ("G1038 v3 EI 150, 12.5 km (rev27low)", "rev27low", "OLRT3W3E", "G1038", "L0"), ("G0975 preview EI 2000 (rev21v, L0 + L2p)", "rev21v", "OLRT3W3", "G0975", None), ("G0063 preview EI 2000 (rev21v)", "rev21v", "OLRT3W3", "G0063", "L0"), ("G1 v3 EI 2000 (rev21v)", "rev21v", "OLRT3W3", "G1", "L0")):
        fs = [f for k, f in B[blk].items() if k[0] == arm and k[1] == geom and (lv is None or k[2] == lv) and (blk != "rev21v" or geom != "ODY" or True)]
        alt[lab] = cell(fs)
    # P30-3 cells with the registered measure
    P = {}
    for lab, geom, arm, need, use_l1 in (("G0528", "G0528", "OLRT3W3L5", 5, True), ("G1230", "G1230", "OLRT3W3E5", 6, False), ("Odyssey", "ODY", "OLRT3W4L5", 5, True)):
        fs = sel("rev30seed", arm, geom); c = cell(fs); P[lab] = dict(c, need=need, measure="L1" if use_l1 else "raw", within=c["l11"] if use_l1 else c["raw1"], within_other=c["raw1"] if use_l1 else c["l11"])
    # G0528 pooling and history shift
    def decomp(f):
        q = PS.flight_position(f["dir"], "G0528", f["tstar"]); return dict(c=q["c_hist_deg"], a=q["a_timing_deg"], b=q["b_pred_deg"], raw=f["pos"]["raw"], l1=f["pos"]["l1"], dt=q["dt_d"], dv=f["dv_sum"])
    g150a = [decomp(f) for f in sel("rev26val2", "OLRT3W3L", "G0528", seeds=range(13, 17))]; g150b = [decomp(f) for f in sel("rev30seed", "OLRT3W3L5", "G0528")]; g2000 = [decomp(f) for f in sel("rev23geo", "OLRT3W3", "G0528")]
    pool = g150a + g150b; cs150 = [d["c"] for d in pool]; cs2000 = [d["c"] for d in g2000]; sh, ci = boot_unpaired(cs150, cs2000)
    shb = boot_unpaired([d["c"] for d in g150b], cs2000)
    l1_150 = sum(d["l1"] < 1 for d in pool); l1_2000 = sum(d["l1"] < 1 for d in g2000)
    G0 = dict(pool_n=len(pool), pool_aim12=len(g150a), pool_aim5=len(g150b), c150=cs150, c2000=cs2000, c150_med=float(np.median(cs150)), c2000_med=float(np.median(cs2000)), shift=sh, ci=ci, shift_aim5_only=shb[0], ci_aim5_only=shb[1], l1_150=l1_150, l1_2000=l1_2000,
              fisher_l1=G.fisher(l1_2000, len(g2000) - l1_2000, l1_150, len(pool) - l1_150), a150_med=float(np.median([d["a"] for d in pool])), a2000_med=float(np.median([d["a"] for d in g2000])), b150_med=float(np.median([d["b"] for d in pool])), b2000_med=float(np.median([d["b"] for d in g2000])),
              l1_150_med=float(np.median([d["l1"] for d in pool])), l1_2000_med=float(np.median([d["l1"] for d in g2000])), c150_aim5=[d["c"] for d in g150b])
    c = SC.compare_passes(V.ROOT / "rev30seed" / "ODY_PL_OLRT3_pert_ODY_s14_L4", V.ROOT / "rev17" / "ODY_PL_OLRT3_pert_ODY_s14_L4"); reg = dict(equal=c["equal"], ncols=c["ncols"], rows=(c["rows_a"], c["rows_b"]), first=c["first"])
    S = pd.read_csv(next(V.ROOT.glob("summary_rev30seed_*.csv"))); regrow = S[S.arm == "OLRT3"].iloc[0]; reg["ei"] = float(regrow.ei_km); reg["stop"] = str(regrow.stop_reason)
    (out / "rev30seed.json").write_text(json.dumps(SC.nan_none(dict(rows=rows, alt=alt, P30_3=P, g0528=G0, regression=reg)), indent=1, default=float))
    for lab, r in rows.items(): print(f"{lab:55s} n={r['n']} ok={r['ok']} hard={r['hard']} win={r['win']} corr={r['corr']} l1<1={r['l11']} raw<1={r['raw1']} ALL={r['all']} W=[{r['wilson'][0]:.2f},{r['wilson'][1]:.2f}] dV={r['dv_med']:.2f}")
    for lab, r in alt.items(): print("ALT", f"{lab:60s} n={r['n']} ALL={r['all']} l1<1={r['l11']} dV={r['dv_med']:.2f}")
    for lab, r in P.items(): print("P30-3", lab, r["measure"], r["within"], "/", r["n"], "need", r["need"], "other", r["within_other"], "ok", r["ok"], "win", r["win"], "dv", round(r["dv_med"], 2), "raw/l1 med", round(r["raw_med"], 2), round(r["l1_med"], 2), "all", r["all"])
    print(G0); print(reg)


if __name__ == "__main__":
    main()
