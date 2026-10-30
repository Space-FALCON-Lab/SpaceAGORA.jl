"""Venus generality, part (a): frozen estimators v3 and v5 against the 20x20 engine truth, scored against the pre-registered expectations V1-V6
(tasks/2026-10-15_periapsis-region-border-targeting/venus_generality_plan.md). Definitions are those of Stage C (analyze_mars_generality.py) and the v5 note:
bins by truth pass time (right-closed), RMS of predicted minus truth per pass index, completed = the 60 d run ended with Success.
Usage: python3 -I venus_analyze.py [--results-root DIR]  -> <root>/venus/analysis/{bin_errors.csv, cell_errors.csv, verdict.csv, classes.csv}
"""
import argparse, math, pathlib, sys, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from results_root import add_results_root_argument, results_root

VC = tomllib.load(open(HERE / "venus_generality_config.toml", "rb"))
BINS = VC["analysis"]["lookahead_bin_edges_days"]; Q_HARD = VC["band"]["q_hard_w_cm2"]; DAY = 86400.0
SPIN = abs(2.99e-7)       # SOURCED: Venus() struct spin magnitude (rad/s), used only for the plan's V4 rationale number
EST = {"v3": "pred_v3", "v5": "pred_v5"}; EST_SC = {"v3": "pred_v3_selfcal", "v5": "pred_v5_selfcal"}
LAB = [f"{BINS[i]:g}-{BINS[i + 1]:g}" for i in range(len(BINS) - 1)]


def rms(x): x = np.asarray(x, float); return float(np.sqrt(np.mean(x ** 2))) if len(x) else float("nan")


def load(root):
    v = root / "venus"; cells = pd.read_csv(v / "cells.csv"); rs = pd.read_csv(v / "truth_run_summary.csv").set_index("tag")
    out = []
    for c in cells.itertuples():
        tr = pd.read_csv(v / "truth" / c.tag / "passes_3b.csv")
        status = "completed" if rs.loc[c.tag, "retcode"] == "Success" else "decayed"
        out.append(dict(tag=c.tag, inc=c.inc, raan=c.raan, omega=c.omega, status=status, truth=tr,
                        true_class="decayed" if status == "decayed" else ("breach" if tr.q_w_cm2.max() > Q_HARD else "safe")))
    return out


def errors(root, cells, folder):
    recs = []
    for c in cells:
        f = root / "venus" / folder / f"{c['tag']}.csv"
        if not f.exists(): continue
        pr = pd.read_csv(f); tr = c["truth"]; n = min(len(pr), len(tr))
        if n < 1: continue
        t = tr.time_s.values[:n] / DAY
        for k in range(n):
            recs.append(dict(tag=c["tag"], status=c["status"], k=k + 1, t_days=t[k], err_alt_m=pr.alt_m.values[k] - tr.alt_m.values[k],
                             err_q_rel=pr.q_w_cm2.values[k] / tr.q_w_cm2.values[k] - 1, dt_s=pr.t_s.values[k] - tr.time_s.values[k]))
    D = pd.DataFrame(recs); D["bin"] = pd.cut(D.t_days, BINS, right=True, include_lowest=True, labels=LAB); return D


def bin_table(D, band_km, H_km):
    g = D.groupby("bin", observed=False)
    t = g.agg(n=("err_alt_m", "size"), n_cells=("tag", "nunique"), rms_alt_km=("err_alt_m", lambda v: rms(v) / 1e3), max_alt_km=("err_alt_m", lambda v: np.max(np.abs(v)) / 1e3 if len(v) else np.nan),
              rms_q_rel=("err_q_rel", rms), rms_dt_s=("dt_s", rms), max_dt_s=("dt_s", lambda v: np.max(np.abs(v)) if len(v) else np.nan)).reset_index()
    t["rms_alt_over_band"] = t.rms_alt_km / band_km; t["rms_alt_over_H"] = t.rms_alt_km / H_km
    return t


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap); a = ap.parse_args(); root = pathlib.Path(results_root(a.results_root))
    out = root / "venus" / "analysis"; out.mkdir(exist_ok=True); pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    cells = load(root); H_km = float(pd.read_csv(root / "venus" / "planet_constants.csv").scale_height_m.iloc[0]) / 1e3
    band_km = float(np.median([(c["truth"].h_qlow_m.iloc[0] - c["truth"].h_qhigh_m.iloc[0]) / 1e3 for c in cells]))
    ncomp = sum(c["status"] == "completed" for c in cells); primary = "completed" if ncomp >= 8 else "all"
    print(f"cells {len(cells)}, completed {ncomp}, primary set = {primary}; band width at pass 1 (median) {band_km:.4f} km")
    rows, ctab, cm_rows, cell_rows = [], [], [], []
    for seedname, folders in (("truth_seed", EST), ("selfcal_seed", EST_SC)):
        for est, folder in folders.items():
            D = errors(root, cells, folder)
            for setname, Dsub in (("completed", D[D.status == "completed"]), ("all", D)):
                T = bin_table(Dsub, band_km, H_km); T.insert(0, "set", setname); T.insert(0, "seed", seedname); T.insert(0, "estimator", est); rows.append(T)
    B = pd.concat(rows); B.to_csv(out / "bin_errors.csv", index=False)
    # pass-1 offsets and per-cell errors, truth seed
    P1 = {}
    for est, folder in EST.items():
        off = []
        for c in cells:
            f = root / "venus" / folder / f"{c['tag']}.csv"; pr = pd.read_csv(f)
            e = pr.alt_m.values[:len(c["truth"])] - c["truth"].alt_m.values[:len(pr)]
            off.append(dict(estimator=est, tag=c["tag"], inc=c["inc"], raan=c["raan"], omega=c["omega"], status=c["status"], pass1_offset_m=float(pr.alt_m.iloc[0] - c["truth"].alt_m.iloc[0]),
                            rms_alt_km=rms(e) / 1e3 if len(e) else np.nan, n_truth=len(c["truth"]), n_pred=len(pr)))
        P1[est] = pd.DataFrame(off); cell_rows.append(P1[est])
    pd.concat(cell_rows).to_csv(out / "cell_errors.csv", index=False)
    # 3-class table
    cl = []
    for est in EST:
        pc = pd.read_csv(root / "venus" / f"pred_cells_{est}.csv").set_index("tag")
        for c in cells:
            pr = pd.read_csv(root / "venus" / EST[est] / f"{c['tag']}.csv"); dec = bool(pc.loc[c["tag"], "decayed_pred"])
            pcl = "decayed" if dec else ("breach" if pr.q_w_cm2.max() > Q_HARD else "safe")
            cl.append(dict(estimator=est, tag=c["tag"], true_class=c["true_class"], pred_class=pcl, t_end_true_days=c["truth"].time_s.iloc[-1] / DAY,
                           t_decay_pred_days=pc.loc[c["tag"], "t_decay_pred_days"], n_truth=len(c["truth"]), n_pred=len(pr)))
    CL = pd.DataFrame(cl); CL.to_csv(out / "classes.csv", index=False)
    # verdict rows (truth seed, primary set)
    V = []
    def add(i, stmt, value, thr, ok, note=""): V.append(dict(id=i, statement=stmt, value=value, threshold=thr, passed=bool(ok), note=note))
    Bp = B[(B.seed == "truth_seed") & (B.set == primary)]
    for est in EST:
        be = Bp[Bp.estimator == est].set_index("bin")
        p1 = P1[est]; p1 = p1[p1.status == "completed"] if primary == "completed" else p1
        r1 = rms(p1.pass1_offset_m) / 1e3
        add(f"V1 {est}", "pass-1 altitude offset RMS (km)", r1, "<= 0.25", r1 <= 0.25, f"{len(p1)} cells")
        for b in LAB:
            r = be.loc[b]
            add(f"V2 {est}", f"RMS altitude error {b} d (km)", r.rms_alt_km, "<= 1.0", r.rms_alt_km <= 1.0, f"n={int(r.n)} passes, {int(r.n_cells)} cells; /band {r.rms_alt_over_band:.3f}; band-scaled threshold 0.555 km: {'PASS' if r.rms_alt_km <= 0.555 else 'FAIL'}")
            add(f"V3 {est}", f"RMS relative Q error {b} d", r.rms_q_rel, "<= 0.10", r.rms_q_rel <= 0.10)
        if est == "v5":
            for b in LAB[:3]:
                r = be.loc[b]; add("V4 v5", f"pooled pass-time RMS {b} d (s)", r.rms_dt_s, "<= 1000", r.rms_dt_s <= 1000.0)
    b3 = Bp[Bp.estimator == "v3"].set_index("bin"); b5 = Bp[Bp.estimator == "v5"].set_index("bin")
    for b in LAB:
        add("V5", f"v5 altitude RMS <= 1.05 x v3, {b} d", b5.loc[b].rms_alt_km, f"<= {1.05 * b3.loc[b].rms_alt_km:.4f} (v3 {b3.loc[b].rms_alt_km:.4f})", b5.loc[b].rms_alt_km <= 1.05 * b3.loc[b].rms_alt_km)
        add("V5", f"v5 pass-time RMS <= 1.05 x v3, {b} d", b5.loc[b].rms_dt_s, f"<= {1.05 * b3.loc[b].rms_dt_s:.1f} (v3 {b3.loc[b].rms_dt_s:.1f})", b5.loc[b].rms_dt_s <= 1.05 * b3.loc[b].rms_dt_s)
    for est in EST:
        g = CL[CL.estimator == est]; cm = pd.crosstab(g.true_class, g.pred_class).reindex(index=["safe", "breach", "decayed"], columns=["safe", "breach", "decayed"], fill_value=0); cm.to_csv(out / f"confusion_{est}.csv")
        acc = np.trace(cm.values) / cm.values.sum(); nd = cm.loc["decayed"].sum(); rec = cm.loc["decayed", "decayed"] / nd if nd else np.nan
        nc = cm.loc[["safe", "breach"]].values.sum(); fd = cm.loc[["safe", "breach"], "decayed"].sum() / nc if nc else np.nan
        add(f"V6a {est}", "3-class accuracy", acc, ">= 0.95", acc >= 0.95); add(f"V6b {est}", "decay recall", rec, ">= 0.95", nd > 0 and rec >= 0.95, f"{cm.loc['decayed','decayed']}/{nd}")
        add(f"V6c {est}", "false-decay rate", fd, "<= 0.05", fd <= 0.05, f"{cm.loc[['safe','breach'],'decayed'].sum()}/{nc}")
    VT = pd.DataFrame(V); VT.to_csv(out / "verdict.csv", index=False)
    print(B[(B.seed == "truth_seed")].round(4).to_string(index=False)); print(VT.to_string(index=False))
    print("V4 rationale number: 0.3 / (4 w_rot) = %.3g s" % (0.3 / (4 * SPIN)))
    for est in EST:
        g = CL[CL.estimator == est]; d = g[(g.true_class == "decayed") & (g.pred_class == "decayed")]
        print(est, "decayed cells predicted decayed: n=%d, decay-time error days: mean %.3f rms %.3f max abs %.3f; pass count diff (pred - truth) median %.1f" % (
            len(d), (d.t_decay_pred_days - d.t_end_true_days).mean(), rms(d.t_decay_pred_days - d.t_end_true_days), (d.t_decay_pred_days - d.t_end_true_days).abs().max(), (g.n_pred - g.n_truth).median()))


if __name__ == "__main__":
    main()
