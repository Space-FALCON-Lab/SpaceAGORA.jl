"""Revision 23b preparation (analysis only: Python, no Julia, no flights; imports are read-only, no existing file is edited).

  tables   planner_ground_tables_rev23.json (+ .sha256): ground tables WITH the GRAM mean-field term for G0528, G0890 (cut-256 samples) and G0975M, G0063M (fine samples, found by elements), same schema as
           planner_ground_tables_generic.json (built by planner_predictions.ground_table at the calibrated planner context, as planner_generic_tables.py does).
  ic       rev20_data/ic/headline_rev23.json (+ .sha256) = headline_rev21.json (entries unchanged) + the select_ic plans of the new tokens (3 scenarios x 2 horizons x fb_term False/True, kappa 1 ASSUMED).
  check    G0975 / G0063 consequence check: DP with mean-field vs preview tables; guard-v3 replay (planner_rev21_guard.GUARD_V3) of the existing paths with each table; predicted guard t* for G0528 / G0890.
Tokens (table keys): G0528, G0890, G0975M, G0063M.  select_ic is run on the new tables through planner_rev23_geometry.patched (the module-global swap), so planner_rev20_ic is not modified.
Usage: python3 -I planner_rev23b_prep.py tables | ic | check"""
import hashlib, json, math, pathlib, sys
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_rev23_geometry as X
import planner_rev20_ic as R

TD = pathlib.Path(PLB.task_dir())
TAB_FILE = TD / "planner_ground_tables_rev23.json"
IC_FILE = TD / "rev20_data" / "ic" / "headline_rev23.json"
NEW = {"G0528": dict(src="cut", gid="G0528", role="GE"), "G0890": dict(src="cut", gid="G0890", role="GF"),
       "G0975M": dict(src="fine", base="G0975", role="GC"), "G0063M": dict(src="fine", base="G0063", role="GA")}
KAPPA = 1.0                                   # ASSUMED (no flight calibration; as KAPPA_G0975 in planner_rev20_ic)


def cand_field(tok):
    G = X.init(); o = NEW[tok]
    if o["src"] == "cut":
        c = next(c for c in G["grid"] if c["id"] == o["gid"]); return c, G["mf"]["cut"][o["gid"]], "analytic_data/meanfield (cut-256 id %s)" % o["gid"]
    g = G["cfg"]["planner"]["geometry"][o["base"]]; c = dict(id=o["base"], inc=float(g["inc_deg"]), epoch_off=float(g["epoch_off_d"]), omega=float(g["omega_deg"]), raan=float(g["raan_deg"]))
    fid = G["fine_key"][X.key(c)]; return c, X.load_fine()[fid], "analytic_data/meanfield_fine (fine id %s)" % fid


def build_tables():
    G = X.init(); out = {}
    for tok, o in NEW.items():
        c, field, src = cand_field(tok); tab, m = G["PP"].ground_table(G["ctx"], c, field, G["cfg"])
        chk = X.table_from(X.build(c, field)); assert np.array_equal(np.asarray(tab["path"]), chk["path"]) and tab["n"] == chk["n"], tok
        tab.update(token=tok, role=o["role"], field_mode="meanfield", field_source=src, inc_deg=c["inc"], omega_deg=c["omega"], raan_deg=c["raan"], epoch_off_d=c["epoch_off"]); out[tok] = tab
    return out


def cmd_tables():
    T = build_tables()
    doc = dict(schema="planner_ground_tables/1", note="Revision 23b; WITH the GRAM mean-field term (field_mode meanfield); G0975M / G0063M = the generic geometries rebuilt with the field; kappa not stored (1 ASSUMED)", geometries=T)
    TAB_FILE.write_text(json.dumps(doc, indent=0, default=float)); h = PLB.sha256_file(TAB_FILE)
    (TD / (TAB_FILE.name + ".sha256")).write_text(h + "  " + TAB_FILE.name + "\n"); print("tables", {k: (v["n"], round(v["q_req"], 5), round(v["W"], 4)) for k, v in T.items()}, "sha256", h)


def load_new():
    G = json.load(open(TAB_FILE))["geometries"]; return {k: {a: (np.array(b, float) if isinstance(b, list) else b) for a, b in v.items()} for k, v in G.items()}


def old_tab(base):
    return {a: (np.array(b, float) if isinstance(b, list) else b) for a, b in json.load(open(TD / "planner_ground_tables_generic.json"))["geometries"][base].items()}


def plan(tab, sc, hz, fbt, tok):
    with X.patched(tab):
        try:
            c = R.Case("X", sc, 250e3, fbt, KAPPA, hz, 0.0, 0.0); ic = R.select_ic(tok, sc, fb_term=fbt, kappa=KAPPA, horizon=hz, case=c)
        except TypeError:                                   # Case has no self-consistent arrival horizon (the plan arrives beyond the table): reported as infeasible, as select_ic does for no plan
            return dict(scenario=sc, geom=tok, feasible=False, tube_binding=True, t_star_d=None, s_d=None, fb_term=bool(fbt), horizon=hz, h_des_m=250e3, delta_ln=0.0, why="no plan inside the table (arrival beyond table + tail)")
    ic["geom"] = tok; return ic


def slim(d):
    return {k: v for k, v in d.items() if k != "labels"}


def cmd_ic():
    H = json.load(open(TD / "rev20_data" / "ic" / "headline_rev21.json")); T = load_new(); new = []
    for tok, tab in T.items():
        for hz in ("arrival", "full"):
            for sc in ("ref", "S-fixed", "S-free"):
                for fbt in (False, True):
                    new.append(slim(plan(tab, sc, hz, fbt, tok)))
    H["plans"] = H["plans"] + new
    H["rev23_note"] = "Revision 23b: headline_rev21.json (entries unchanged) + select_ic plans of G0528, G0890, G0975M, G0063M on planner_ground_tables_rev23.json (mean field; kappa 1 ASSUMED; planner_rev23b_prep.py)"
    IC_FILE.write_text(json.dumps(H, indent=1, default=float)); h = PLB.sha256_file(IC_FILE); (IC_FILE.parent / (IC_FILE.name + ".sha256")).write_text(h + "  " + IC_FILE.name + "\n")
    print("plans", len(H["plans"]), "sha256", h)
    print("infeasible entries:", [(p["geom"], p["scenario"], p["horizon"], p["fb_term"]) for p in new if not p["feasible"]])
    for p in new:
        if p["scenario"] == "S-fixed" and p["horizon"] == "arrival" and not p["fb_term"] and p["feasible"]:
            print(p["geom"], "S-fixed x0 %.4f nu_pre %.4f pre %.4f campaign %.4f N %d t_arr %.2f feasible %s" % (p["x0"], p["nu_pre"], p["dv_pre_mps"], p["campaign_dv_mps"], p["N"], p["t_arr_d"], p["feasible"]))


# ----------------------------------------------------------------------------------------------------- replay helpers

def rp_register(RP, tok, tab):
    """Register a table under a replay token (derived fields as planner_rev20_replay.load_tab)."""
    t = dict(tab); t["token"] = tok; ra = t["ra"]; n = t["n"]
    t["d0"] = np.concatenate((np.diff(-ra), [float(ra[-2] - ra[-1])])); t["_ra_inc"] = ra[::-1]; t["_d0_inc"] = t["d0"][::-1]; t["_k_inc"] = np.arange(n, 0, -1, dtype=float)[::-1] * 1.0; t["_kk"] = np.arange(1, n + 1, dtype=float)
    RP._TABS[tok] = t; RP.KAPPA[tok] = KAPPA; RP.QFIRST_CAL.setdefault(tok, RP.QFIRST_CAL.get("G0975", 0.211971)); return t


def run_paths(RP, ids, tok, tab, guard, tstar):
    out = []
    for rid in ids:
        F = RP.extract(rid, tab, KAPPA); m = RP.split_id(rid)
        r = RP.replay(F, tok, "OLRT3", ic=None, noise=RP.NOISE, guard=guard, tstar_d=tstar, seed=m["seed"] * 100 + (m["level"] or 0)); r["tstar_d"] = tstar
        r["dt_d"] = (r["t_arr_d"] - tstar) if r["reached"] else float("nan"); r["ok_win"] = bool(r["reached"] and abs(r["dt_d"]) <= r["w_max_d"]); r["success"] = bool(r["reached"] and r["ok_apo"] and not r["hard_violation"] and r["ok_win"])
        out.append(r)
    return out


def cmd_check():
    import planner_rev20_replay as RP
    import planner_rev21_guard as G3
    out = dict(dp={}, replay={}, tstar={}); T = load_new()
    for tok, base in (("G0975M", "G0975"), ("G0063M", "G0063")):
        old = old_tab(base); new = T[tok]; row = {}
        for sc in ("ref", "S-fixed", "S-free"):
            a = plan(old, sc, "arrival", False, base); b = plan(new, sc, "arrival", False, tok)
            row[sc] = dict(preview=dict(aero=a["campaign_dv_mps"] + a["dv_pre_mps"], N=a["N"], t_arr_d=a["t_arr_d"], x0=a["x0"]), mean_field=dict(aero=b["campaign_dv_mps"] + b["dv_pre_mps"], N=b["N"], t_arr_d=b["t_arr_d"], x0=b["x0"]))
        row["q_req"] = dict(preview=float(old["q_req"]), mean_field=float(new["q_req"])); row["n"] = dict(preview=int(old["n"]), mean_field=int(new["n"])); out["dp"][tok] = row
    # replay: G0975 and G0063 paths with each table (guard v3, reference start, noise indexings)
    import glob
    def run_dir_ext(run_id):
        for ph in ("dev", "devT", "devR", "gen1", "sens", "rev17", "runs", "rev21c", "rev21v"):
            c = RP.ROOT / ph / run_id
            if (c / "passes.csv").exists() and all((c / f).exists() for f in ("decisions.csv", "measured.csv", "replan.csv")): return c
        return None
    RP.run_dir = run_dir_ext
    def ids_for(pre, arms):
        ids = []
        for ph in ("gen1", "rev21c", "rev21v"):
            for d in sorted((RP.ROOT / ph).glob(pre + "_PL_*")) if (RP.ROOT / ph).exists() else []:
                nm = d.name; arm = nm.split("_PL_")[1].split("_pert_")[0]; lv = nm.split("_L")[1] if "_L" in nm else ""
                if arm in arms and lv in ("", "0") and run_dir_ext(nm) is not None and nm not in ids: ids.append(nm)
        return ids
    RP.TOKEN.update({"G0063": "G0063"})
    for tok, base in (("G0975M", "G0975"), ("G0063M", "G0063")):
        ids = ids_for(base, ("OLRT3", "OLRT35", "OLRT5")); res = {}
        for noise in ("pass", "state"):
            RP.NOISE = noise; res[noise] = {}
            for ver, tab in (("preview", old_tab(base)), ("mean_field", T[tok])):
                vt = base + "_" + ver
                tab = rp_register(RP, vt, tab); RP.QFIRST_CAL[vt] = RP.QFIRST_CAL.get("G0975", 0.211971)
                ts = RP.tstar_for(vt); b = run_paths(RP, ids, vt, tab, None, ts); g = run_paths(RP, ids, vt, tab, G3.GUARD_V3, ts)
                dd = [x["dv_campaign"] - y["dv_campaign"] for y, x in zip(b, g) if y["reached"] and x["reached"]]
                res[noise][ver] = dict(tstar_d=ts, guard=[dict(id=i, ok_win=x["ok_win"], hard=bool(x["hard_violation"]), reached=bool(x["reached"]), dt_d=x["dt_d"], dv=x["dv_campaign"]) for i, x in zip(ids, g)],
                                       base=[dict(id=i, ok_win=x["ok_win"], hard=bool(x["hard_violation"]), reached=bool(x["reached"]), dt_d=x["dt_d"], dv=x["dv_campaign"]) for i, x in zip(ids, b)],
                                       extra_dv_median=float(np.median(dd)) if dd else float("nan"))
        flips = {}
        for noise in res:
            p, m = res[noise]["preview"], res[noise]["mean_field"]
            fl = []
            for kind in ("guard", "base"):
                for a, b_ in zip(p[kind], m[kind]):
                    if a["ok_win"] != b_["ok_win"] or a["hard"] != b_["hard"]: fl.append((kind, a["id"], a["ok_win"], b_["ok_win"], a["hard"], b_["hard"]))
            flips[noise] = dict(n_changed=len(fl), changed=fl, extra_dv_median_preview=p["extra_dv_median"], extra_dv_median_mean_field=m["extra_dv_median"], extra_dv_shift=m["extra_dv_median"] - p["extra_dv_median"],
                                tstar_preview=p["tstar_d"], tstar_mean_field=m["tstar_d"], window_ok_guard=(sum(x["ok_win"] for x in p["guard"]), sum(x["ok_win"] for x in m["guard"])), n=len(ids))
        out["replay"][tok] = dict(ids=ids, results=res, comparison=flips)
        for noise, f in flips.items(): print(tok, noise, "n", f["n"], "changed", f["n_changed"], "extra dV shift %.3f" % f["extra_dv_shift"], "window ok guard (preview, M)", f["window_ok_guard"], "t* %.2f -> %.2f" % (f["tstar_preview"], f["tstar_mean_field"]))
    # predicted guard t* for G0528, G0890 (kappa 1 ASSUMED, S-fixed IC of headline_rev23.json, aim offset 0.0 d: not in guard t_offset_d)
    H = json.load(open(IC_FILE))
    for tok in ("G0528", "G0890"):
        rp_register(RP, tok, T[tok]); pl = [p for p in H["plans"] if p["geom"] == tok and p["scenario"] == "S-fixed" and p["horizon"] == "arrival" and not p["fb_term"]][0]
        ts = RP.tstar_for(tok, pl); off = G3.GUARD_V3["t_offset_d"].get(tok, 0.0); sw, P = RP.arm_def("OLRT3")
        out["tstar"][tok] = dict(guard_tstar_d=ts, aim_offset_d=off, aim_d=ts + off, ic_plan_t_arr_d=pl["t_arr_d"], N=pl["N"], x0=pl["x0"], nu_pre=pl["nu_pre"], dv_pre_mps=pl["dv_pre_mps"], campaign_dv_mps=pl["campaign_dv_mps"], kappa=KAPPA,
                                 w_max_d_nominal=3.0 * 2 * math.pi * math.sqrt((0.5 * (RP.RA_T + RP.GUARD_W_RP_M)) ** 3 / RP.MU) / 86400.0)
        print(tok, out["tstar"][tok])
    p = TD / "rev23_data" / "prep_check.json"; p.write_text(json.dumps(out, indent=1, default=float)); print("sha256", PLB.sha256_file(p)[:16])


if __name__ == "__main__":
    {"tables": cmd_tables, "ic": cmd_ic, "check": cmd_check}[sys.argv[1]]()
