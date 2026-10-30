"""Revision 21 flight-validation tables (analysis only: Python, reads flight directories, no Julia, writes no flight data).
Scores rev21v (OLRT3-W3 = guard v3 + S-fixed, OLRT3-WF3 = guard v3 + S-free), rev21c (OLRT3 and baselines at G0975/G0063/ODY), gen1 (Revision 16b generic flights at G0975) and the Revision 20 v2 flights (rev20val OLRT3-W) with the
per-flight metrics of planner_rev20_score.flight, and builds the tables of results_rev21_val.md: v3 against v2 paired by (geometry, level, seed), the terminal anatomy of one flight, S-free against S-fixed with the first-pass heat rate against
exp(x0), G0975/G0063 cells against OLRT3 and the baselines, ODY L0/L2 over seeds 13-20, and the acceptance table. All numbers are MEASURED (flight logs) or DERIVED (formulas named in the note); thresholds are ASSUMED as registered.
Usage: python3 -I planner_rev21_val.py --out-dir DIR     (DIR gets val_flights.json and val_tables.md)
Env: SPACEAGORA_PLANNER_ROOT (flight root, default ~/Documents/SpaceAGORA-planner-data/planner)."""
import argparse, json, math, os, pathlib, sys
import numpy as np, pandas as pd

HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev20_score as SC
import planner_rev21_position as PS

ROOT = pathlib.Path(os.environ.get("SPACEAGORA_PLANNER_ROOT", str(pathlib.Path.home() / "Documents" / "SpaceAGORA-planner-data" / "planner")))
OFF = {"ODY": 0.15, "G1": 0.0}      # ASSUMED aim offsets (planner_config [planner.guard] and guard_sets.v3), days


def load(block, pattern):
    """All finished flights of a block: dict (armkey, geom, level, seed) -> flight metrics (+ t* fields, dirs)."""
    S = SC.load_summary(ROOT, pattern); D = SC.scan(ROOT, block); out = {}
    for k, d in D.items():
        if not SC.finished(d): continue
        f = SC.flight(d, S.get(d.name), dv_pre_fallback=0.0); f["dir"] = str(d); f["geom"] = k[1]
        sr = S.get(d.name) or {}; f["ic_scenario"] = str(sr.get("ic_scenario", "n/a")); f["ic_x0"] = float(sr.get("ic_x0", float("nan"))) if sr.get("ic_x0") is not None else float("nan")
        f["stop_reason"] = str(sr.get("stop_reason", "n/a")); f["total_dv_summary"] = float(sr.get("total_dv_mps", float("nan")))
        out[k] = f
    return out


def tstar_fields(F, nominal_by_geom):
    """Per flight: aimed t* (logged ttg aim) and nominal t* (aim minus offset); unguarded flights use the geometry nominal. Adds dt at both, window first passage and with the wait rule."""
    for k, f in F.items():
        g = k[1]; aim = f["tt_aim_d"]
        if np.isfinite(aim): nom = aim - OFF.get(g, 0.0)
        else: nom = nominal_by_geom[g]; aim = nom + OFF.get(g, 0.0)
        f["tstar_nom_d"] = nom; f["tstar_aim_d"] = aim
        if f["reached"]:
            f["dt_nom_w"] = (f["t_arr_d"] - nom) / f["wmax_d"]; f["dt_aim_w"] = (f["t_arr_d"] - aim) / f["wmax_d"]; f["window"] = bool(abs(f["dt_nom_w"]) <= 1.0)
            w = SC.D20.wait_rule(f["err"], f["t_arr_d"], f["last_drop_km"], f["T_final_s"] / 86400.0, nom, SC.W_MAX_PERIODS); f["window_wait"] = bool(w["met_wait"])
        else: f["dt_nom_w"] = f["dt_aim_w"] = float("nan"); f["window"] = f["window_wait"] = False
        f["mission"] = bool(f["ok"] and f["window"])
    return F


def q1_ratio(f):
    """First-pass mean-field heat rate (passes.csv q_w_cm2, pass 1) over exp(ic_x0)."""
    P = pd.read_csv(pathlib.Path(f["dir"]) / "passes.csv").sort_values("pass"); return float(P.q_w_cm2.iloc[0]) / math.exp(f["ic_x0"]) if np.isfinite(f["ic_x0"]) else float("nan")


def q(x, p): return float(np.percentile(x, p)) if len(x) else float("nan")


def stats(x):
    x = [v for v in x if np.isfinite(v)]; return dict(n=len(x), med=q(x, 50), q25=q(x, 25), q75=q(x, 75), mn=min(x) if x else float("nan"), mx=max(x) if x else float("nan"))


def fmt(s, p=2): return f"{s['med']:.{p}f} [{s['q25']:.{p}f}, {s['q75']:.{p}f}] ({s['mn']:.{p}f} to {s['mx']:.{p}f})"


def tab(head, rows): return ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)] + ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    V3 = load("rev21v", "summary_rev21v_*.csv"); C = load("rev21c", "summary_rev21c_*.csv"); G1 = load("gen1", "summary_gen1_*.csv"); V2 = load("rev20val", "summary_rev20val*.csv"); R17 = load("rev17", "summary_rev17_*.csv")
    nom = {"ODY": PS.nominal("ODY")["days"], "G1": PS.nominal("G1")["days"], "G0975": PS.nominal("G0975")["days"], "G0063": PS.nominal("G0063")["days"]}
    for F in (V3, C, G1, V2, R17): tstar_fields(F, nom)
    pos = json.load(open(out / "score_rev21.json"))["position"]["flights"] if (out / "score_rev21.json").exists() else {}
    L = []; J = {}
    # ---------- (1) v3 against v2
    L += ["### Table 1a. Per pair: OLRT3-W3 (v3, rev21v) minus OLRT3-W (v2, rev20val), same geometry, level and seed", ""]
    rows = []; pairs = []
    for g, lvls in (("ODY", ["L0", "L1", "L2", "L3", "L4"]), ("G1", ["L0"])):
        for lv in lvls:
            for s in (13, 14, 15, 16):
                f3 = V3.get(("OLRT3W3", g, lv, s)); f2 = V2.get(("OLRT3W", g, lv, s))
                if not f3 or not f2: continue
                pairs.append((g, lv, s, f3, f2))
                rows.append([g, lv, s, "OK" if f2["ok"] else "FAIL", "OK" if f3["ok"] else "FAIL", f"{f2['dv_sum']:.2f}", f"{f3['dv_sum']:.2f}", f"{f3['dv_sum'] - f2['dv_sum']:+.2f}", f2["burns"], f3["burns"], f2["reversals"], f3["reversals"],
                             f"{f2['dt_nom_w']:+.2f}", f"{f3['dt_nom_w']:+.2f}", f"{f2['dt_aim_w']:+.2f}", f"{f3['dt_aim_w']:+.2f}", f"{f2['err']:.1f}", f"{f3['err']:.1f}"])
    L += tab(["geom", "level", "seed", "v2", "v3", "v2 aero dV (m/s)", "v3 aero dV", "v3 - v2", "v2 burns", "v3 burns", "v2 reversals", "v3 reversals", "v2 dt/W_max (nominal t*)", "v3 dt/W_max (nominal)", "v2 dt/W_max (aimed t*)", "v3 dt/W_max (aimed)", "v2 final error (km)", "v3 final error (km)"], rows) + [""]
    groups = {"ODY (20 pairs)": [p for p in pairs if p[0] == "ODY"], "G1 L0 (4 pairs)": [p for p in pairs if p[0] == "G1"], "all (24 pairs)": pairs}
    rows = []; J["v3_v2"] = {}
    for nm, ps in groups.items():
        f3s = [p[3] for p in ps]; f2s = [p[4] for p in ps]
        d = stats([p[3]["dv_sum"] - p[4]["dv_sum"] for p in ps]); db = stats([p[3]["burns"] - p[4]["burns"] for p in ps]); dr = stats([p[3]["reversals"] - p[4]["reversals"] for p in ps])
        r2n, r3n = stats([f["dt_nom_w"] for f in f2s]), stats([f["dt_nom_w"] for f in f3s]); r2a, r3a = stats([f["dt_aim_w"] for f in f2s]), stats([f["dt_aim_w"] for f in f3s])
        J["v3_v2"][nm] = dict(d_dv=d, d_burns=db, d_rev=dr, v2_dt_nom=r2n, v3_dt_nom=r3n, v2_dt_aim=r2a, v3_dt_aim=r3a, v3_gt_v2_dv=int(sum(1 for p in ps if p[3]["dv_sum"] > p[4]["dv_sum"])),
                              v2=dict(ok=sum(f["ok"] for f in f2s), hard=sum(f["hard"] for f in f2s), win=sum(f["window"] for f in f2s), wait=sum(f["window_wait"] for f in f2s), abs_aim=stats([abs(f["dt_aim_w"]) for f in f2s]), abs_nom=stats([abs(f["dt_nom_w"]) for f in f2s]), dv=stats([f["dv_sum"] for f in f2s]), burns=stats([f["burns"] for f in f2s]), rev=stats([f["reversals"] for f in f2s])),
                              v3=dict(ok=sum(f["ok"] for f in f3s), hard=sum(f["hard"] for f in f3s), win=sum(f["window"] for f in f3s), wait=sum(f["window_wait"] for f in f3s), abs_aim=stats([abs(f["dt_aim_w"]) for f in f3s]), abs_nom=stats([abs(f["dt_nom_w"]) for f in f3s]), dv=stats([f["dv_sum"] for f in f3s]), burns=stats([f["burns"] for f in f3s]), rev=stats([f["reversals"] for f in f3s])))
        rows.append([nm, "v2", f"{sum(f['ok'] for f in f2s)}/{len(ps)}", sum(f["hard"] for f in f2s), f"{sum(f['window'] for f in f2s)}/{len(ps)}", fmt(stats([f["dv_sum"] for f in f2s])), fmt(stats([f["burns"] for f in f2s]), 0), fmt(stats([f["reversals"] for f in f2s]), 0), fmt(r2n), fmt(r2a)])
        rows.append([nm, "v3", f"{sum(f['ok'] for f in f3s)}/{len(ps)}", sum(f["hard"] for f in f3s), f"{sum(f['window'] for f in f3s)}/{len(ps)}", fmt(stats([f["dv_sum"] for f in f3s])), fmt(stats([f["burns"] for f in f3s]), 0), fmt(stats([f["reversals"] for f in f3s]), 0), fmt(r3n), fmt(r3a)])
        rows.append([nm, "v3 - v2 (paired)", "", "", "", fmt(d), fmt(db, 0), fmt(dr, 0), "", ""])
    L += ["### Table 1b. v3 against v2, medians [IQR] (min to max)", ""] + tab(["group", "guard", "success", "hard", "window (first)", "aero dV (m/s)", "burns", "reversals", "dt/W_max at the nominal t*", "dt/W_max at the aimed t*"], rows) + [""]
    # paired extra dV over OLRT3 (comparator: rev17 for ODY, rev20val OLRT3 for G1)
    def extra(Fw, comp, g, arm):
        r = []
        for (a_, g_, lv, s), f in Fw.items():
            if a_ != arm or g_ != g or s not in (13, 14, 15, 16): continue
            c = comp.get(("OLRT3", g, lv, s))
            if c: r.append(f["dv_sum"] - c["dv_sum"])
        return r
    ex = {}
    ex["ODY v2"] = extra(V2, R17, "ODY", "OLRT3W"); ex["ODY v3"] = extra(V3, R17, "ODY", "OLRT3W3"); ex["G1 v2"] = extra(V2, V2, "G1", "OLRT3W"); ex["G1 v3"] = extra(V3, V2, "G1", "OLRT3W3")
    J["extra_over_olrt3"] = {k: stats(v) for k, v in ex.items()}
    L += ["### Table 1c. Extra aerobraking dV over OLRT3 on the same seeds/levels (comparators: Rev 17 OLRT3 at Odyssey, rev20val OLRT3 at G1; pairs with a comparator)", ""] + tab(["cell", "pairs", "median extra dV (m/s)", "IQR", "min to max"], [[k, s["n"], f"{s['med']:.2f}", f"{s['q25']:.2f} to {s['q75']:.2f}", f"{s['mn']:.2f} to {s['mx']:.2f}"] for k, s in J["extra_over_olrt3"].items()]) + [""]
    # ---------- (2) anatomy of the failure
    f18 = V3[("OLRT3W3", "ODY", "L0", 18)]; P = pd.read_csv(pathlib.Path(f18["dir"]) / "passes.csv"); D = pd.read_csv(pathlib.Path(f18["dir"]) / "decisions.csv")
    ra = P.ra_m.values.astype(float); qq = P.q_w_cm2.values.astype(float); drop = np.concatenate(([np.nan], -np.diff(ra))) / 1e3; c = drop / qq; n = len(P)
    rows = [[int(P["pass"][i]), f"{(ra[i] - SC.TARGET) / 1e3:+.3f}", f"{drop[i]:.2f}", f"{qq[i]:.4f}", f"{c[i]:.1f}", D[D.target_pass == int(P['pass'][i])].reason.iloc[0] if (D.target_pass == int(P['pass'][i])).any() else ""] for i in range(n - 8, n)]
    c_hist = c[n - 7:n - 1]; pred_q = math.exp(float(D.iloc[-2].z0) + float(D.iloc[-2].nu0)); need = (ra[n - 2] - SC.TARGET) / 1e3 + 12.5
    anat = dict(final_err_km=(ra[-1] - SC.TARGET) / 1e3, tol_km=25.0, miss_by_m=(abs((ra[-1] - SC.TARGET) / 1e3) - 25.0) * 1e3, last_drop_km=float(drop[-1]), needed_drop_km=float(need), c_median_last6=float(np.median(c_hist)), c_last=float(c[-1]), c_prior_max=float(np.nanmax(c[n - 30:n - 1])),
                c_cv_30=float(np.nanstd(c[n - 30:n - 1]) / np.nanmean(c[n - 30:n - 1])), pred_q=pred_q, q_last=float(qq[-1]), pred_drop_km=float(np.median(c_hist) * pred_q), q_ratio=float(qq[-1] / pred_q), c_ratio=float(c[-1] / np.median(c_hist)), z_c=float((c[-1] - np.nanmean(c[n - 30:n - 1])) / np.nanstd(c[n - 30:n - 1])))
    J["anatomy_s18"] = anat
    L += ["### Table 2. ODY L0 s18 (OLRT3-W3): the last passes", ""] + tab(["pass", "apoapsis minus target (km), after the pass", "drop (km)", "pass q (W/cm2)", "drop / q (km per W/cm2)", "decision before the pass"], rows) + [""]
    # ---------- (3) S-free against S-fixed
    rows = []; J["sfree"] = {}
    for g, tok in (("ODY", "ODY"), ("G1", "G1")):
        for s in (13, 14, 15, 16):
            fw = V3.get(("OLRT3WF3", g, "L0", s)); fx = V3.get(("OLRT3W3", g, "L0", s))
            if not fw or not fx: continue
            J["sfree"][f"{g}_s{s}"] = dict(q1_free=q1_ratio(fw), q1_fixed=q1_ratio(fx), dv_free=fw["dv_sum"], dv_fixed=fx["dv_sum"], pre_free=fw["dv_pre"], pre_fixed=fx["dv_pre"], ok=(fw["ok"], fx["ok"]), hard=(fw["hard"], fx["hard"]), win=(fw["window"], fx["window"]), x0=fw["ic_x0"])
            rows.append([g, s, f"{fw['ic_x0']:.4f}", f"{math.exp(fw['ic_x0']):.4f}", f"{q1_ratio(fw):.4f}", f"{q1_ratio(fx):.4f}", f"{fw['dv_sum']:.2f}", f"{fx['dv_sum']:.2f}", f"{fw['dv_sum'] - fx['dv_sum']:+.2f}", f"{fw['dv_pre']:.3f}", f"{fx['dv_pre']:.3f}", "OK" if fw["ok"] else "FAIL", "OK" if fx["ok"] else "FAIL", int(fw["hard"]), "yes" if fw["window"] else "no", "yes" if fx["window"] else "no", f"{fw['dt_nom_w']:+.2f}", f"{fx['dt_nom_w']:+.2f}"])
    L += ["### Table 3. OLRT3-WF3 (S-free) against OLRT3-W3 (S-fixed), L0, same seeds", ""] + tab(["geom", "seed", "ic x0 (ln q)", "target exp(x0) (W/cm2)", "WF3 first-pass q / exp(x0)", "W3 first-pass q / exp(x0)", "WF3 aero dV (m/s)", "W3 aero dV", "WF3 - W3", "WF3 pre-burn", "W3 pre-burn", "WF3", "W3", "WF3 hard", "WF3 window", "W3 window", "WF3 dt/W_max", "W3 dt/W_max"], rows) + [""]
    cal = pd.read_csv(ROOT / "summary_calx0_mars_calx0.csv"); J["calx0"] = {r.tag: dict(q_first=float(r.q_first_w_cm2), q_target=float(r.q_target_w_cm2), rel=float(r.q_first_w_cm2 / r.q_target_w_cm2 - 1)) for r in cal.itertuples()}
    # ---------- (4) generic geometries
    nomg = nom
    rows = []; J["generic"] = {}
    def cell(F, arm, g, lv, label, src):
        fs = [f for k, f in F.items() if k[0] == arm and k[1] == g and (k[2] == lv or (lv == "L0" and k[2] == ""))]
        if not fs: return None
        t_nat = stats([f["t_arr_d"] for f in fs]); r = dict(label=label, src=src, n=len(fs), ok=sum(f["ok"] for f in fs), hard=sum(f["hard"] for f in fs), win=sum(f["window"] for f in fs), wait=sum(f["window_wait"] for f in fs), dv=stats([f["dv_sum"] for f in fs]), burns=stats([f["burns"] for f in fs]), t_arr=t_nat, dt=stats([f["dt_nom_w"] for f in fs]), tstar=stats([f["tstar_nom_d"] for f in fs]), pre=stats([f["dv_pre"] for f in fs]))
        return r
    spec = [("G0975", "L0", [("OLRT3W3", V3, "OLRT3-W3 (v3)", "rev21v"), ("OLRT3", C, "OLRT3", "rev21c"), ("OLRT3", G1, "OLRT3", "gen1 (rev16b)"), ("OLRT35", G1, "OLRT35", "gen1"), ("OLRT5", G1, "OLRT5", "gen1"), ("BODYT", G1, "B-ODY-T", "gen1"), ("N2IT", G1, "N2-I-T", "gen1")]),
            ("G0975", "L2p", [("OLRT3W3", V3, "OLRT3-W3 (v3)", "rev21v"), ("OLRT3", C, "OLRT3", "rev21c"), ("BODYT", C, "B-ODY-T", "rev21c"), ("N2IT", C, "N2-I-T", "rev21c")]),
            ("G0063", "L0", [("OLRT3W3", V3, "OLRT3-W3 (v3)", "rev21v"), ("OLRT3", C, "OLRT3", "rev21c"), ("BODYT", C, "B-ODY-T", "rev21c"), ("N2IT", C, "N2-I-T", "rev21c")])]
    for g, lv, arms in spec:
        for arm, F, label, src in arms:
            r = cell(F, arm, g, lv, label, src)
            if not r: continue
            J["generic"][f"{g}_{lv}_{label}_{src}"] = r
            rows.append([g, lv, label, src, r["n"], r["ok"], r["hard"], f"{r['win']}/{r['n']}", f"{r['wait']}/{r['n']}", fmt(r["dv"]), fmt(r["burns"], 0), f"{r['t_arr']['med']:.1f}", f"{r['tstar']['med']:.2f}", f"{r['dt']['med']:+.1f}", f"{r['pre']['med']:.3f}"])
    L += ["### Table 4. Generic geometries: OLRT3-W3 (v3) against OLRT3 and the baselines (window: first passage, W_max = 3 final orbits, t* = the preview-table nominal of the geometry)", ""] + tab(["geom", "level", "arm", "block", "n", "success", "hard", "window (first)", "window (wait)", "aero dV median [IQR] (min to max)", "burns", "median arrival (d)", "t* (d)", "median dt/W_max", "pre-burn (m/s)"], rows) + [""]
    pr = []
    for g, lv in (("G0975", "L0"), ("G0975", "L2p"), ("G0063", "L0")):
        w3 = [f for k, f in V3.items() if k[0] == "OLRT3W3" and k[1] == g and k[2] == lv]; cc = [f for k, f in C.items() if k[0] == "OLRT3" and k[1] == g and k[2] == lv]
        for s_ in (13, 14, 15, 16):
            a_ = next((f for k, f in V3.items() if k[0] == "OLRT3W3" and k[1] == g and k[2] == lv and k[3] == s_), None); b_ = next((f for k, f in C.items() if k[0] == "OLRT3" and k[1] == g and k[2] == lv and k[3] == s_), None)
            if a_ and b_: pr.append([g, lv, s_, f"{a_['dv_sum']:.2f}", f"{b_['dv_sum']:.2f}", f"{a_['dv_sum'] - b_['dv_sum']:+.2f}", a_["burns"], b_["burns"], f"{a_['t_arr_d']:.2f}", f"{b_['t_arr_d']:.2f}", f"{a_['tstar_nom_d'] - b_['t_arr_d']:.1f}", f"{(a_['dv_sum'] - b_['dv_sum']) / (a_['tstar_nom_d'] - b_['t_arr_d']):.3f}"])
    L += ["### Table 4b. Pairs W3 minus OLRT3 at the generic cells, with the delay the guard has to add (t* minus the natural arrival of the paired OLRT3 flight)", ""] + tab(["geom", "level", "seed", "W3 dV", "OLRT3 dV", "extra dV (m/s)", "W3 burns", "OLRT3 burns", "W3 arrival (d)", "OLRT3 arrival (d)", "delay to add (d)", "extra dV per day of delay (m/s/d)"], pr) + [""]
    # ---------- (5) ODY L0/L2 over seeds 13-20
    rows = []; J["p21_6"] = {}
    for lv in ("L0", "L2"):
        v3 = [V3[("OLRT3W3", "ODY", lv, s)] for s in range(13, 21) if ("OLRT3W3", "ODY", lv, s) in V3]
        mix = [V2[("OLRT3W", "ODY", lv, s)] for s in range(13, 17) if ("OLRT3W", "ODY", lv, s) in V2] + [V3[("OLRT3W3", "ODY", lv, s)] for s in range(17, 21) if ("OLRT3W3", "ODY", lv, s) in V3]
        for nm, fs in ((f"v3 on seeds 13-20 (all rev21v)", v3), ("v2 on 13-16 (rev20val) + v3 on 17-20", mix)):
            J["p21_6"][f"{lv} {nm}"] = dict(n=len(fs), ok=sum(f["ok"] for f in fs), hard=sum(f["hard"] for f in fs), win=sum(f["window"] for f in fs), wait=sum(f["window_wait"] for f in fs), mission=sum(f["mission"] for f in fs))
            rows.append([lv, nm, len(fs), sum(f["ok"] for f in fs), sum(f["hard"] for f in fs), sum(f["window"] for f in fs), sum(f["window_wait"] for f in fs), sum(f["mission"] for f in fs), fmt(stats([f["dv_sum"] for f in fs])), fmt(stats([f["dt_nom_w"] for f in fs]))])
    L += ["### Table 5. Odyssey L0 and L2, seeds 13-20 (seeds 17-20 flew v3; in rev21v seeds 13-16 flew v3 as well)", ""] + tab(["level", "set", "flights", "success", "hard", "window (first)", "window (wait)", "mission (success and window)", "aero dV median [IQR] (min to max)", "dt/W_max (nominal t*)"], rows) + [""]
    # ---------- (6) acceptance table
    ac = []; P9 = lambda blk, key: pos.get(blk, {}).get(key, {}).get("err_deg")
    def poscell(blk, arm_label, g, lv, F, armkey):
        e = []
        for (a_, g_, l_, s_), f in F.items():
            if a_ == armkey and g_ == g and l_ == lv:
                v = P9(blk, f"{armk}@{g}@{lv}@s{s_}")
                if v is not None: e.append(v)
        return e
    cells = [("OLRT3W3", "OLRT3-W3", "ODY", lv, V3, "rev21v") for lv in ("L0", "L1", "L2", "L3", "L4")] + [("OLRT3W3", "OLRT3-W3", "G1", "L0", V3, "rev21v"), ("OLRT3WF3", "OLRT3-WF3", "ODY", "L0", V3, "rev21v"), ("OLRT3WF3", "OLRT3-WF3", "G1", "L0", V3, "rev21v"),
                                                                                                         ("OLRT3W3", "OLRT3-W3", "G0975", "L0", V3, "rev21v"), ("OLRT3W3", "OLRT3-W3", "G0975", "L2p", V3, "rev21v"), ("OLRT3W3", "OLRT3-W3", "G0063", "L0", V3, "rev21v")]
    comp_for = lambda g, lv, s: (R17.get(("OLRT3", g, lv, s)) if g == "ODY" and s <= 16 else (C.get(("OLRT3", g, lv, s)) if g in ("ODY", "G0975", "G0063") else V2.get(("OLRT3", g, lv, s))))
    J["ac"] = {}
    for armk, lab, g, lv, F, blk in cells:
        fs = {k[3]: f for k, f in F.items() if k[0] == armk and k[1] == g and k[2] == lv}; n = len(fs)
        cmp_ = [comp_for(g, lv, s) for s in sorted(fs)]; cmp_ok = [c for c in cmp_ if c]
        hard = sum(f["hard"] for f in fs.values()); hard_c = sum(c["hard"] for c in cmp_ok); apo = sum(f["apo_ok"] for f in fs.values()); apo_c = sum(c["apo_ok"] for c in cmp_ok)
        need = math.ceil((14 if lv == "L0" else 12) / 16 * n); win = sum(f["window"] for f in fs.values()); wt = sum(f["window_wait"] for f in fs.values())
        ac2 = "PASS" if hard == 0 or (lv != "L0" and hard <= hard_c) else "FAIL"; ac3 = "PASS" if apo >= apo_c - n / 16 else "FAIL"; ac3 = ac3 if len(cmp_ok) == n else ac3 + f" (comparator {len(cmp_ok)}/{n})"
        ac4 = "PASS" if win >= need else "FAIL"; corr = sum(1 for f in fs.values() if f["reached"] and f["corridor"] and f["qs_last_ok"]); reach = sum(1 for f in fs.values() if f["reached"])
        e = poscell(blk, lab, g, lv, F, armk); w1 = sum(1 for x in e if x <= 1.0)
        J["ac"][f"{lab}@{g}@{lv}"] = dict(n=n, hard=hard, hard_c=hard_c, apo=apo, apo_c=apo_c, win=win, need=need, wait=wt, corr=corr, reach=reach, pos_med=q(e, 50), pos_p90=q(e, 90), pos_w1=w1, pos_n=len(e), AC2=ac2, AC3=ac3, AC4=ac4)
        ac.append([f"{lab} {g} {lv}", n, f"{hard} (comparator {hard_c})", ac2, f"{apo} vs {apo_c}", ac3, f"{win} (need {need}); wait {wt}", ac4, f"{corr}/{reach}", f"{q(e, 50):.2f} / {q(e, 90):.2f}; {w1}/{len(e)} within 1.0 deg"])
    L += ["### Table 6. Acceptance criteria per cell (rev21v; 4-flight and 8-flight scaling of the 14/16 and 12/16 thresholds is ASSUMED)", ""] + tab(["cell", "flights", "AC2 hard (OLRT3 comparator hard)", "AC2", "AC3 apoapsis success vs comparator", "AC3", "AC4 window (first passage)", "AC4", "AC5 corridor and Q_s", "AC9 position median / p90 (deg)"], ac) + [""]
    # baselines window (AC6): all rev21c + gen1 baseline flights
    bl = [f for F in (C, G1) for k, f in F.items() if k[0] in ("BODYT", "N2IT")]; J["ac6"] = dict(n=len(bl), win=sum(f["window"] for f in bl), wait=sum(f["window_wait"] for f in bl), ok=sum(f["ok"] for f in bl))
    ol3 = [f for F in (C,) for k, f in F.items() if k[0] == "OLRT3"]; J["ol3_alone"] = dict(n=len(ol3), win=sum(f["window"] for f in ol3), ok=sum(f["ok"] for f in ol3))
    J["w3_all"] = dict(n=len(V3), ok=sum(f["ok"] for f in V3.values()), hard=sum(f["hard"] for f in V3.values()), win=sum(f["window"] for f in V3.values()), wait=sum(f["window_wait"] for f in V3.values()))
    (out / "val_tables.md").write_text("\n".join(L) + "\n"); (out / "val_flights.json").write_text(json.dumps(SC.nan_none(J), indent=1)); print("ok", len(V3), "rev21v flights")


if __name__ == "__main__":
    main()
