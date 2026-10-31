"""Revision 31.3 (analysis only: Python, no Julia, no flights): the unflown 20 deg and 80 deg inclination classes, node free.

  select   rank EVERY geometry of the class in the fine mean-field grid (12 omega0 x 8 Omega0 x 4 epochs = 384 per class; neither class is in the 1,536 grid, both are in the 8,064 fine grid and its GRAM samples) by the DP
           aerobraking dV (S-fixed, free arrival) at two targets: the Odyssey radius 7,408 km (mean-field table from the fine sample) and 2,000 km altitude (radius 5,396,190 m; PREVIEW table: no mean field along the longer track).
  request  machine-readable GRAM sampling request (tracks CSV + JSON, format of rev25_data/sampling_request) for the 2,000 km tracks of the selected geometries.
  build    after the sampling: mean-field tables for the four cells, S-fixed IC plans, guard t* (kappa 1 ASSUMED) into NEW pinned files planner_ground_tables_rev31.json, rev20_data/ic/headline_rev31.json, rev25_data/../rev31_data/guard_tstar_rev31.json.
Cell keys (token = key): <fine id>_A2000 (2,000 km) and <fine id>_AODY (Odyssey radius 7,408 km). Tie rule for the choice: S-fixed dV within 0.005 m/s are ties, broken by the dV of the plan hardened for a +0.2231 ln density offset, then the id.
Usage: python3 -I planner_rev31_classes.py select [--jobs 4] | request | build --samples DIR | build --preview"""
import argparse, hashlib, json, math, pathlib, sys
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_rev23_geometry as X
import planner_rev20_ic as R
import planner_select_geometry as SG

TD = SG.TD; OUT = TD / "rev31_data"
CLASSES = {"EQ20": 20.0, "HIGH80": 80.0}      # ASSUMED classes (Revision 24)
ALT_LOW_KM = 2000.0                            # the Revision 31.3 low target
RA_LOW = SG.radius_of(ALT_LOW_KM); RA_ODY = 7408000.0
KAPPA = 1.0; HORIZON_D = SG.HORIZON_D
TIE = 0.005


def fine_class(inc):
    G = X.init(); return [c for c in G["fine"] if c["inc"] == inc]


def evaluate(c):
    G = X.init(); fld = X.load_fine().get(c["id"]); out = dict(id=c["id"], inc=c["inc"], omega0=c["omega"], raan0=c["raan"], epoch=c["epoch_off"], field=fld is not None)
    for lv, ra, prev in (("ODY", None, False), ("LOW", RA_LOW, True)):
        tab, info = SG.candidate_table(c, ra, field=fld, preview=prev); rt = RA_ODY if ra is None else ra; d = {}
        for sc in ("S-fixed", "S-free"):
            r = X.dp(tab, sc, kappa=KAPPA, ra_t=rt, horizon_d=HORIZON_D); d[sc] = {k: r.get(k) for k in ("feasible", "why", "aero", "N", "days", "x0", "nu_pre", "rp_final")}
        d["table"] = info["field_mode"]; d["n"] = info["n"]; out[lv] = d
        if d["S-fixed"]["feasible"] and d["S-fixed"]["aero"] < 0.3:                        # robustness only for the competitive ones
            h = X.dp(tab, "S-fixed", kappa=KAPPA, harden=(0.0, X.OFFSET_LN), ra_t=rt, horizon_d=HORIZON_D); d["harden_plus"] = h["aero"] if h.get("feasible") else None
    return out


def cmd_select(jobs):
    X.init(); X.load_fine(); T = []
    for name, inc in CLASSES.items(): T += fine_class(inc)
    rows = X.pmap(evaluate, T, jobs, "rev31"); X.write_json(OUT / "selection_raw.json", rows); return summarize(rows)


def summarize(rows=None):
    rows = rows or json.load(open(OUT / "selection_raw.json")); S = dict(coverage={}, best={}, top={})
    G = X.init(); grid_incs = sorted(set(c["inc"] for c in G["grid"])); S["coverage"] = dict(grid_1536_inclinations=grid_incs, in_1536={n: (inc in grid_incs) for n, inc in CLASSES.items()},
                                                                                   fine_candidates={n: sum(1 for r in rows if r["inc"] == inc) for n, inc in CLASSES.items()}, fine_with_field={n: sum(1 for r in rows if r["inc"] == inc and r["field"]) for n, inc in CLASSES.items()})
    for name, inc in CLASSES.items():
        g = [r for r in rows if r["inc"] == inc]
        for lv in ("ODY", "LOW"):
            ok = [r for r in g if r[lv]["S-fixed"]["feasible"]]
            key = lambda r: (round(r[lv]["S-fixed"]["aero"] / TIE), r[lv].get("harden_plus") if r[lv].get("harden_plus") is not None else 9.0, r["id"])
            ok.sort(key=key); S["top"]["%s_%s" % (name, lv)] = [dict(id=r["id"], omega0=r["omega0"], raan0=r["raan0"], epoch=r["epoch"], aero_fixed=r[lv]["S-fixed"]["aero"], aero_free=r[lv]["S-free"]["aero"], harden_plus=r[lv].get("harden_plus"), N=r[lv]["S-fixed"]["N"], days=r[lv]["S-fixed"]["days"], table=r[lv]["table"]) for r in ok[:6]]
            S["best"]["%s_%s" % (name, lv)] = S["top"]["%s_%s" % (name, lv)][0]; S["best"]["%s_%s" % (name, lv)]["n_feasible"] = len(ok); S["best"]["%s_%s" % (name, lv)]["n"] = len(g)
            S["best"]["%s_%s" % (name, lv)]["median"] = float(np.median([r[lv]["S-fixed"]["aero"] for r in ok])); S["best"]["%s_%s" % (name, lv)]["n_le_0p1"] = int(sum(r[lv]["S-fixed"]["aero"] <= 0.1 for r in ok))
    X.write_json(OUT / "selection.json", S); return S


def cells(S=None):
    """The four cells: key -> dict(cand, ra_target_m, alt_km, class, level). Geometry = the best of the class at THAT target (ODY: with field; LOW: preview ranking)."""
    S = S or json.load(open(OUT / "selection.json")); G = X.init(); out = {}
    for name in CLASSES:
        for lv, ra, alt, tag in (("ODY", RA_ODY, 4011.81, "AODY"), ("LOW", RA_LOW, ALT_LOW_KM, "A2000")):
            b = S["best"]["%s_%s" % (name, lv)]; c = next(c for c in G["fine"] if c["id"] == b["id"]); out["%s_%s" % (b["id"], tag)] = dict(cand=c, ra_target_m=ra, alt_km=alt, cls=name, level=lv, aero_fixed=b["aero_fixed"])
    return out


def cmd_request(out_dir=None):
    G = X.init(); GR = G["GR"]; ctx = G["ctx"]; cfg = ctx.cfg; out = pathlib.Path(out_dir or OUT / "sampling_request"); out.mkdir(parents=True, exist_ok=True); rows = []; man = dict(candidates={}); offs = cfg["sampling"]["alt_offsets_m"]; hc = cfg["orbit"]["h_ref_m"]
    for key, o in cells().items():
        if o["level"] != "LOW": continue
        c = o["cand"]; q = GR.solve_q_req(ctx, c, cfg["band"]["d_ref_days"])
        with X.orbit_target(ctx, o["ra_target_m"], SG.MAX_PASSES): S = GR.schedule(ctx, c, q)
        S0 = GR.schedule(ctx, c, q); lon, lat, _ = GR.lon_lat(ctx, c, S); ep = (pd.Timestamp(cfg["anchor"]["epoch_utc"]) + pd.Timedelta(days=c["epoch_off"])).isoformat()
        for k in range(S["n"]): rows.append((key, ep, c["epoch_off"], c["omega"], c["raan"], c["inc"], k, S["t"][k], lat[k], lon[k], hc, 0, 0))
        man["candidates"][key] = dict(token=key, class_tag=o["cls"], fine_id=c["id"], inc_deg=c["inc"], omega0_deg=c["omega"], raan0_deg=c["raan"], epoch_off_d=c["epoch_off"], epoch_utc=ep, target_alt_km=o["alt_km"], ra_target_m=o["ra_target_m"], passes=int(S["n"]),
                                      passes_default_radius=int(S0["n"]), extra_passes=int(S["n"] - S0["n"]), calls=int(S["n"] * len(offs)), calls_extension_only=int((S["n"] - S0["n"]) * len(offs)), q_req=float(q))
    df = pd.DataFrame(rows, columns="cand epoch_utc epoch_off_d omega0_deg raan0_deg inc_deg k t_s lat_gd_rad lon_rad h_center_m in_cut256 in_cut384".split()); csv = out / "tracks_rev31low.csv"; df.to_csv(csv, index=False)
    man.update(schema="gram_sampling_request/1", tracks_csv=csv.name, tracks_sha256=hashlib.sha256(csv.read_bytes()).hexdigest(), alt_offsets_m=list(offs), h_center_m=hc, calls_total=int(len(df) * len(offs)),
               calls_extension_only_total=int(sum(v["calls_extension_only"] for v in man["candidates"].values())), est_call_seconds=float(len(df) * len(offs) * 0.157e-3),
               julia="julia --project=. benchmarks/studies/periapsis_border_switchability/meanfield_sample.jl --phase=sample --out=<DIR> --tracks=<this csv> --cut=none   (one Julia job per host, route :process, MemoryMax cap)",
               consumer="planner_rev31_classes.py build --samples <DIR> reads <DIR>/mf_*.csv; the Odyssey-radius cells use the existing fine-grid samples (meanfield_fine) and need no new sampling",
               note="every pass of the longer track is sampled (3 altitudes around h_center_m); est_call_seconds is 0.157 ms per call (SOURCED probe), not a performance number")
    (out / "sampling_request.json").write_text(json.dumps(man, indent=1, default=float)); return man


def cmd_build(samples=None, out_dir=None):
    import planner_rev20_replay as RP, planner_rev21_guard as G3, planner_rev18_offline as O
    out = pathlib.Path(out_dir) if out_dir else TD; (out / "rev20_data" / "ic").mkdir(parents=True, exist_ok=True); (out / "rev31_data").mkdir(parents=True, exist_ok=True)
    F = SG.lowalt_fields(samples) if samples else {}; tabs = {}; plans = []; tstar = {}
    for key, o in cells().items():
        c = o["cand"]; ra_t = o["ra_target_m"]
        if o["level"] == "LOW":
            field = F.get(key) if samples else None
            if samples and field is None: raise KeyError("no samples for %s" % key)
            tab, info = SG.candidate_table(c, ra_t, field=field, preview=(field is None))
        else:
            tab, info = SG.candidate_table(c, None)             # fine-grid mean-field sample of the geometry (field_mode meanfield)
        t = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in tab.items()}
        t.update(token=key, role=o["cls"], field_mode=info["field_mode"], field_source=info.get("field_source"), inc_deg=c["inc"], omega_deg=c["omega"], raan_deg=c["raan"], epoch_off_d=c["epoch_off"], ra_target_m=ra_t, target_alt_km=o["alt_km"], fine_id=c["id"]); tabs[key] = t
        for hz in ("arrival", "full"):
            for sc in ("ref", "S-fixed", "S-free"):
                for fbt in (False, True):
                    try:
                        with X.patched(tab, ra_t=ra_t, horizon_d=HORIZON_D):
                            cs = R.Case("X", sc, 250e3, fbt, KAPPA, hz, 0.0, 0.0); ic = R.select_ic(key, sc, fb_term=fbt, kappa=KAPPA, horizon=hz, case=cs)
                    except TypeError:
                        ic = dict(scenario=sc, geom=key, feasible=False, tube_binding=True, t_star_d=None, s_d=None, fb_term=bool(fbt), horizon=hz, h_des_m=250e3, delta_ln=0.0, why="no plan inside the table (arrival beyond table + tail)")
                    ic = {k: v for k, v in ic.items() if k != "labels"}; ic["geom"] = key; ic["ra_target_m"] = ra_t; plans.append(ic)
        arr = {k: (np.array(v, float) if isinstance(v, list) else v) for k, v in t.items()}; old = (RP.RA_T, O.RA_T); RP.RA_T = O.RA_T = ra_t
        try:
            RP.KAPPA[key] = KAPPA; RP.QFIRST_CAL.setdefault(key, RP.QFIRST_CAL.get("G0975", 0.211971)); a2 = dict(arr); n = a2["n"]; ra = a2["ra"]
            a2["d0"] = np.concatenate((np.diff(-ra), [float(ra[-2] - ra[-1])])); a2["_ra_inc"] = ra[::-1]; a2["_d0_inc"] = a2["d0"][::-1]; a2["_k_inc"] = np.arange(n, 0, -1, dtype=float)[::-1] * 1.0; a2["_kk"] = np.arange(1, n + 1, dtype=float); RP._TABS[key] = a2
            pl = [p for p in plans if p["geom"] == key and p["scenario"] == "S-fixed" and p["horizon"] == "arrival" and not p["fb_term"]][0]; ts = RP.tstar_for(key, pl) if pl.get("feasible") else float("nan")
        finally:
            RP.RA_T, O.RA_T = old
        tstar[key] = dict(guard_tstar_d=ts, aim_offset_d=G3.GUARD_V3["t_offset_d"].get(key, 0.0), ic_plan_t_arr_d=pl.get("t_arr_d"), N=pl.get("N"), x0=pl.get("x0"), nu_pre=pl.get("nu_pre"), dv_pre_mps=pl.get("dv_pre_mps"), campaign_dv_mps=pl.get("campaign_dv_mps"),
                          kappa=KAPPA, ra_target_m=ra_t, field_mode=info["field_mode"])
    tf = out / "planner_ground_tables_rev31.json"; hf = out / "rev20_data" / "ic" / "headline_rev31.json"; gf = out / "rev31_data" / "guard_tstar_rev31.json"
    tf.write_text(json.dumps(dict(schema="planner_ground_tables/1", note="Revision 31.3 tables for the 20 and 80 deg classes; token = key; ra_target_m per table; kappa 1 ASSUMED; %s" % ("2,000 km cells: PREVIEW (no samples given)" if not samples else "2,000 km cells: mean field from the rev31 samples"), geometries=tabs), indent=0, default=float))
    hf.write_text(json.dumps(dict(labels="Revision 31.3 IC plans (planner_rev20_ic.Case, kappa 1 ASSUMED, horizon 130 d); headline layout; extra IC file", plans=plans), indent=1, default=float)); gf.write_text(json.dumps(tstar, indent=1, default=float)); sh = {}
    for f in (tf, hf):
        h = PLB.sha256_file(f); pathlib.Path(str(f) + ".sha256").write_text(h + "  " + f.name + "\n"); sh[f.name] = h
    sh[gf.name] = PLB.sha256_file(gf); print(json.dumps(sh, indent=1)); return dict(files=sh, tstar=tstar)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["select", "summarize", "request", "build"]); ap.add_argument("--jobs", type=int, default=4); ap.add_argument("--samples"); ap.add_argument("--preview", action="store_true"); ap.add_argument("--out-dir"); a = ap.parse_args()
    if a.cmd == "select": print(json.dumps(cmd_select(a.jobs)["best"], indent=1, default=float))
    elif a.cmd == "summarize": print(json.dumps(summarize()["best"], indent=1, default=float))
    elif a.cmd == "request": print(json.dumps(cmd_request(a.out_dir), indent=1)[:1200])
    else:
        if not a.preview and not a.samples: ap.error("build needs --preview or --samples DIR")
        cmd_build(None if a.preview else a.samples, a.out_dir)
