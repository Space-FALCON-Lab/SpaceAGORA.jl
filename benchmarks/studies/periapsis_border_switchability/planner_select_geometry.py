"""Revision 25 (WP25-P): geometry selection and the low-altitude block preparation (analysis only: Python, no Julia, no flights). Every new behavior is default-off; with the defaults nothing else changes.

  select_geometry(...)   rank candidate geometries by the Revision 20 DP aerobraking dV to a given closing apoapsis radius, with an OPTIONAL plane term (target inclination, optional node) added to the net cost.
                         Tables are built per candidate from the sampled GRAM mean field (cut-256 or fine samples, found by the candidate's elements) or, where no sample exists, from the preview model (flagged).
  sampling_request()     machine-readable GRAM call set (tracks CSV for meanfield_sample.jl + JSON manifest) along the longer tracks of the block geometries to their target radii.
  build_lowalt()         after the sampling: mean-field tables to each target radius, S-fixed IC plans and guard t* into NEW pinned files (`--preview` runs the same code on preview tables, used now).
  selftest               default-off identities and the build path on preview data.
Named parameters (ASSUMED unless noted): ROBUST (density offset, L1 reserve scale, jitter: SOURCED Revision 17a / geometry_rank_config), HORIZON_D 130 (SOURCED planner horizon), KAPPA 1.0 (ASSUMED), BLOCK (the Revision 25 cells).
Usage: python3 -I planner_select_geometry.py selftest | request | build [--preview | --samples DIR] [--out-dir DIR]"""
import argparse, hashlib, json, math, os, pathlib, sys
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_lib as PLB
import planner_analysis as AN
import planner_rev23_geometry as X
import planner_rev20_ic as R
import planner_rev22_rendezvous as R22

TD = pathlib.Path(PLB.task_dir())
R_EQ = 3396190.0                                  # SOURCED geometry_rank_config.toml [mars] r_eq_m
KAPPA = 1.0                                       # ASSUMED (no flight calibration at low apoapsis; measured by the block)
HORIZON_D = 130.0                                 # SOURCED planner_config.toml horizon_days
MAX_PASSES = 5000                                 # ASSUMED guard of the extended schedule loop
ROBUST = dict(offset_ln=X.OFFSET_LN, l1_scale=X.L1_SCALE, reserve_ln=X.RESERVE0, jitter_frac=X.JIT_FRAC, jitter_n=X.JIT_N)      # SOURCED Revision 17a / geometry_rank_config.toml
# The Revision 25 block cells (WP24 targets; geometry elements are the candidates' own, asserted against rev23_data/cut256_table.json by WP25-J)
BLOCK = {"G0486": dict(alt_km=2000.0, inc_deg=63.4, omega0=0.0, raan0=270.0, epoch_off=90.0, cell="GG", target_id="CRIT63_2000km"),
         "G1230": dict(alt_km=1000.0, inc_deg=116.6, omega0=270.0, raan0=270.0, epoch_off=0.0, cell="GH", target_id="CRIT117_1000km"),
         "G1038": dict(alt_km=500.0, inc_deg=93.5, omega0=270.0, raan0=270.0, epoch_off=180.0, cell="GI", target_id="ODY_500km")}
OUT_DIR = TD / "rev25_data"


def radius_of(alt_km):
    return R_EQ + float(alt_km) * 1e3


def exact_plane_increment(di_deg, ra_m, rp_m):
    """Plane term of an apoapsis-node burn into a circle of radius ra_m from the final ellipse (r_p rp_m): sqrt((vc - va)^2 + 4 vc va sin^2(di/2)) - (vc - va), m/s. DERIVED (Revision 24)."""
    mu = R22.MU
    vc = math.sqrt(mu / ra_m); va = math.sqrt(mu * (2 / ra_m - 2 / (ra_m + rp_m)))
    return math.sqrt((vc - va) ** 2 + 4 * vc * va * math.sin(math.radians(di_deg) / 2) ** 2) - (vc - va)


# ----------------------------------------------------------------------------------------------------- tables per candidate

def find_field(cand):
    """Mean-field sample of a candidate by its elements: the cut-256 sample (id of the 1,536 grid), else the fine sample (id found by elements). Returns (field, source) or (None, None)."""
    G = X.init(); k = X.key(cand)
    if cand.get("id") == "ODY" and "ODY" in G["mf"]["cut"]: return G["mf"]["cut"]["ODY"], "cut:ODY"        # the Odyssey anchor is not on the grid
    for c in G["grid"]:
        if X.key(c) == k and c["id"] in G["mf"]["cut"]: return G["mf"]["cut"][c["id"]], "cut:" + c["id"]
    fid = G["fine_key"].get(k)
    if fid is not None:
        f = X.load_fine().get(fid)
        if f is not None: return f, "fine:" + fid
    return None, None


def candidate_table(cand, ra_target_m=None, field="auto", preview=False, samples=None):
    """Ground table of a candidate to closing apoapsis radius ra_target_m (None = the default RA_TARGET, i.e. the existing table). field: 'auto' = look the sample up; a dict = that field; None = none.
    Same corridor: q_req and tube W of the default-radius campaign (as Revision 24). Returns (table dict with arrays, info)."""
    G = X.init(); ctx = G["ctx"]; GR = G["GR"]; b = ctx.cfg["band"]; src = "given"
    if field == "auto":
        field, src = find_field(cand)
    use_field = (field is not None) and not preview
    H = field["H"] if use_field else None
    q_req = GR.solve_q_req(ctx, cand, b["d_ref_days"], H)
    with X.orbit_target(ctx, ra_target_m, MAX_PASSES if ra_target_m is not None else None):
        S = GR.schedule(ctx, cand, q_req, H); d = GR.channels(ctx, cand, S, field if use_field else None)
    tab = X.table_from(dict(q_req=q_req, S=S, d=d))
    return tab, dict(field_mode="meanfield" if use_field else "preview_no_field", field_source=src if use_field else None, n=tab["n"], q_req=q_req, ra_target_m=float(AN.RA_TARGET if ra_target_m is None else ra_target_m))


def plane_cost(cand, plane, ra_m, rp_m):
    """Optional plane term for a target plane dict(inc_deg, node_deg=None): node None = node free (inclination difference only); a node needs the candidate's arrival node (not available here: ValueError)."""
    if plane is None: return 0.0
    if plane.get("node_deg") is not None: raise ValueError("node-fixed plane term needs Omega0-matched samples (Revision 23 section 5); not implemented")
    return exact_plane_increment(abs(cand["inc"] - plane["inc_deg"]), ra_m, rp_m)


def select_geometry(candidates, scenario="S-fixed", ra_target_m=None, plane=None, preview=False, kappa=KAPPA, horizon_d=HORIZON_D, robust=False, tstar_d=None):
    """Ranked list of dict(candidate, aero, plane, net, x0, nu_pre, t_arr_d, N, feasible, field_mode, ...). plane=None: no plane term (default-off). robust=True adds the hardened +/- offset and L1-reserve dV (ROBUST).
    Decision variables returned per candidate: the geometry (the candidate), x0, nu_pre (S-fixed pre-burn), t* (the plan's arrival, or the window plan when tstar_d is given). Objective: DP aerobraking dV (+ plane term)."""
    G = X.init(); ra_t = float(AN.RA_TARGET if ra_target_m is None else ra_target_m); rows = []
    for cand in candidates:
        tab, info = candidate_table(cand, None if ra_t == float(AN.RA_TARGET) else ra_t, preview=preview)
        r = X.dp(tab, scenario, kappa=kappa, ra_t=ra_t, horizon_d=horizon_d, tstar=tstar_d)
        row = dict(id=cand.get("id"), inc=cand["inc"], omega0=cand["omega"], raan0=cand["raan"], epoch_off=cand["epoch_off"], feasible=bool(r.get("feasible")), why=r.get("why"), field_mode=info["field_mode"], ra_target_m=ra_t)
        if r.get("feasible"):
            pl = plane_cost(cand, plane, ra_t, r["rp_final"]); row.update(aero=r["aero"], plane=pl, net=r["aero"] + pl, x0=r["x0"], nu_pre=r["nu_pre"], dv_pre=r["pre"], t_arr_d=r["days"], N=r["N"], burns=r["burns"], rp_final=r["rp_final"])
            if robust:
                ro = {}
                for nm, kw in (("harden_plus", dict(harden=(0.0, ROBUST["offset_ln"]))), ("harden_minus", dict(harden=(ROBUST["offset_ln"], 0.0))), ("reserve_l1", dict(reserve=ROBUST["reserve_ln"] * ROBUST["l1_scale"]))):
                    q = X.dp(tab, scenario, kappa=kappa, ra_t=ra_t, horizon_d=horizon_d, **kw); ro[nm] = q["aero"] if q.get("feasible") else None
                row["robust"] = ro
        rows.append(row)
    return sorted(rows, key=lambda x: (not x["feasible"], x.get("net", float("inf")), str(x["id"])))


# ----------------------------------------------------------------------------------------------------- sampling request

def block_candidates():
    G = X.init(); out = {}
    for tok, o in BLOCK.items():
        c = next(c for c in G["grid"] if c["id"] == tok); assert abs(c["inc"] - o["inc_deg"]) < 1e-9 and c["omega"] == o["omega0"] and c["raan"] == o["raan0"] and c["epoch_off"] == o["epoch_off"], tok
        out[tok] = c
    return out


def sampling_request(out_dir=None):
    """Tracks CSV (header as geometry_rank.write_tracks: cand epoch_utc epoch_off_d omega0_deg raan0_deg inc_deg k t_s lat_gd_rad lon_rad h_center_m in_cut256 in_cut384) of the three block geometries along the
    campaign to their target radii (cand id = '<token>_A<alt_km 4 digits>'), the same predicted-track convention as the existing samples (fallback scale height), and a JSON manifest with the call counts."""
    G = X.init(); GR = G["GR"]; ctx = G["ctx"]; cfg = ctx.cfg; out = pathlib.Path(out_dir or OUT_DIR / "sampling_request"); out.mkdir(parents=True, exist_ok=True); rows = []; man = dict(candidates={})
    offs = cfg["sampling"]["alt_offsets_m"]; hc = cfg["orbit"]["h_ref_m"]
    for tok, c in block_candidates().items():
        ra_t = radius_of(BLOCK[tok]["alt_km"]); cid = "%s_A%04d" % (tok, int(BLOCK[tok]["alt_km"]))
        q = GR.solve_q_req(ctx, c, cfg["band"]["d_ref_days"])
        with X.orbit_target(ctx, ra_t, MAX_PASSES):
            S = GR.schedule(ctx, c, q)
        S0 = GR.schedule(ctx, c, q)                                   # the default-radius campaign: its track is the prefix of the longer one
        lon, lat, _ = GR.lon_lat(ctx, c, S); ep = (pd.Timestamp(cfg["anchor"]["epoch_utc"]) + pd.Timedelta(days=c["epoch_off"])).isoformat()
        assert np.array_equal(S["t"][:S0["n"] - 1], S0["t"][:S0["n"] - 1])
        for k in range(S["n"]):
            rows.append((cid, ep, c["epoch_off"], c["omega"], c["raan"], c["inc"], k, S["t"][k], lat[k], lon[k], hc, 0, 0))
        man["candidates"][cid] = dict(token=tok, cell=BLOCK[tok]["cell"], inc_deg=c["inc"], omega0_deg=c["omega"], raan0_deg=c["raan"], epoch_off_d=c["epoch_off"], epoch_utc=ep, target_alt_km=BLOCK[tok]["alt_km"], ra_target_m=ra_t,
                                      passes=int(S["n"]), passes_default_radius=int(S0["n"]), extra_passes=int(S["n"] - S0["n"]), calls=int(S["n"] * len(offs)), calls_extension_only=int((S["n"] - S0["n"]) * len(offs)), q_req=float(q))
    df = pd.DataFrame(rows, columns="cand epoch_utc epoch_off_d omega0_deg raan0_deg inc_deg k t_s lat_gd_rad lon_rad h_center_m in_cut256 in_cut384".split()); csv = out / "tracks_rev25low.csv"; df.to_csv(csv, index=False)
    man.update(schema="gram_sampling_request/1", tracks_csv=csv.name, tracks_sha256=hashlib.sha256(csv.read_bytes()).hexdigest(), alt_offsets_m=list(offs), h_center_m=hc, calls_total=int(len(df) * len(offs)),
               calls_extension_only_total=int(sum(v["calls_extension_only"] for v in man["candidates"].values())), est_call_seconds=float(len(df) * len(offs) * 0.157e-3),
               julia="julia --project=. benchmarks/studies/periapsis_border_switchability/meanfield_sample.jl --phase=sample --out=<DIR> --tracks=<this csv> --cut=none   (one Julia job per host, route :process, MemoryMax cap)",
               consumer="planner_select_geometry.py build --samples <DIR> reads <DIR>/mf_*.csv (cand, k, alt_idx, alt_m, lat_rad, lon_rad, t_s, rho_kg_m3); expected files: mf_<epochtag>_<chunk>.csv per epoch",
               note="every pass of the longer track is sampled (3 altitudes around h_center_m); the first passes duplicate the existing default-radius samples; est_call_seconds is 0.157 ms per call (SOURCED probe), not a performance number")
    (out / "sampling_request.json").write_text(json.dumps(man, indent=1, default=float)); return man


# ----------------------------------------------------------------------------------------------------- post-sampling build

def lowalt_fields(samples):
    if samples is None: return {}
    G = X.init()
    return G["GR"].load_meanfield(samples, G["gcfg"])


def build_lowalt(samples=None, out_dir=None, tag="preview"):
    """Tables to each block target radius (mean field from `samples`, else preview), S-fixed (and all) IC plans, guard t* (kappa 1) into NEW files: planner_ground_tables_rev25<tag>.json, rev20_data/ic/headline_rev25<tag>.json,
    rev25_data/guard_tstar_rev25<tag>.json (+ .sha256 for the first two). tag 'preview' -> files named ..._rev25.json (the names WP25-J pins), 'mf' -> ..._rev25mf.json."""
    import planner_rev20_replay as RP, planner_rev21_guard as G3, planner_rev18_offline as O
    out = pathlib.Path(out_dir) if out_dir else TD; (out / "rev20_data" / "ic").mkdir(parents=True, exist_ok=True); (out / "rev25_data").mkdir(parents=True, exist_ok=True)
    suffix = "" if tag == "preview" else tag; F = lowalt_fields(samples); cands = block_candidates(); tabs = {}; plans = []; tstar = {}
    for tok, c in cands.items():
        o = BLOCK[tok]; ra_t = radius_of(o["alt_km"]); cid = "%s_A%04d" % (tok, int(o["alt_km"])); field = F.get(cid) if samples is not None else None
        if samples is not None and field is None: raise KeyError("no samples for %s in %s" % (cid, samples))
        tab, info = candidate_table(c, ra_t, field=field, preview=(field is None))
        t = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in tab.items()}
        t.update(token=tok, role=o["cell"], field_mode=info["field_mode"], field_source=(None if field is None else "samples:" + cid), inc_deg=c["inc"], omega_deg=c["omega"], raan_deg=c["raan"], epoch_off_d=c["epoch_off"],
                 ra_target_m=ra_t, target_alt_km=o["alt_km"], target_id=o["target_id"]); tabs[tok] = t
        for hz in ("arrival", "full"):
            for sc in ("ref", "S-fixed", "S-free"):
                for fbt in (False, True):
                    try:
                        with X.patched(tab, ra_t=ra_t, horizon_d=HORIZON_D):
                            cs = R.Case("X", sc, 250e3, fbt, KAPPA, hz, 0.0, 0.0); ic = R.select_ic(tok, sc, fb_term=fbt, kappa=KAPPA, horizon=hz, case=cs)
                    except TypeError:
                        ic = dict(scenario=sc, geom=tok, feasible=False, tube_binding=True, t_star_d=None, s_d=None, fb_term=bool(fbt), horizon=hz, h_des_m=250e3, delta_ln=0.0, why="no plan inside the table (arrival beyond table + tail)")
                    ic = {k: v for k, v in ic.items() if k != "labels"}; ic["geom"] = tok; ic["ra_target_m"] = ra_t; plans.append(ic)
        # guard t* (planner_ttg_setup definition: nominal closing-apoapsis time of the S-fixed plan; kappa 1; aim offset 0)
        arr = {k: (np.array(v, float) if isinstance(v, list) else v) for k, v in t.items()}
        old = (RP.RA_T, O.RA_T); RP.RA_T = O.RA_T = ra_t
        try:
            RP.KAPPA[tok] = KAPPA; RP.QFIRST_CAL.setdefault(tok, RP.QFIRST_CAL.get("G0975", 0.211971)); a2 = dict(arr); n = a2["n"]; ra = a2["ra"]
            a2["d0"] = np.concatenate((np.diff(-ra), [float(ra[-2] - ra[-1])])); a2["_ra_inc"] = ra[::-1]; a2["_d0_inc"] = a2["d0"][::-1]; a2["_k_inc"] = np.arange(n, 0, -1, dtype=float)[::-1] * 1.0; a2["_kk"] = np.arange(1, n + 1, dtype=float); RP._TABS[tok] = a2
            pl = [p for p in plans if p["geom"] == tok and p["scenario"] == "S-fixed" and p["horizon"] == "arrival" and not p["fb_term"]][0]
            ts = RP.tstar_for(tok, pl) if pl.get("feasible") else float("nan")
        finally:
            RP.RA_T, O.RA_T = old
        tstar[tok] = dict(guard_tstar_d=ts, aim_offset_d=G3.GUARD_V3["t_offset_d"].get(tok, 0.0), ic_plan_t_arr_d=pl.get("t_arr_d"), N=pl.get("N"), x0=pl.get("x0"), nu_pre=pl.get("nu_pre"), dv_pre_mps=pl.get("dv_pre_mps"),
                          campaign_dv_mps=pl.get("campaign_dv_mps"), kappa=KAPPA, ra_target_m=ra_t, field_mode=info["field_mode"])
    tf = out / ("planner_ground_tables_rev25%s.json" % suffix); hf = out / "rev20_data" / "ic" / ("headline_rev25%s.json" % suffix); gf = out / "rev25_data" / ("guard_tstar_rev25%s.json" % suffix)
    tf.write_text(json.dumps(dict(schema="planner_ground_tables/1", note="Revision 25 low-altitude block tables (%s); token = key; ra_target_m per table; kappa not stored (1 ASSUMED)" % ("preview, no mean field" if samples is None else "mean field from the rev25 samples"), geometries=tabs), indent=0, default=float))
    hf.write_text(json.dumps(dict(labels="Revision 25 IC plans (planner_rev20_ic.Case on the rev25 tables, kappa 1 ASSUMED, horizon 130 d); headline layout; extra IC file, rev23 stays pinned", plans=plans, rev25_note="geom = table key = geometry token; ra_target_m in every entry"), indent=1, default=float))
    gf.write_text(json.dumps(tstar, indent=1, default=float)); sh = {}
    for f in (tf, hf):
        h = PLB.sha256_file(f); pathlib.Path(str(f) + ".sha256").write_text(h + "  " + f.name + "\n"); sh[f.name] = h
    sh[gf.name] = PLB.sha256_file(gf); print(json.dumps(sh, indent=1)); return dict(files=sh, tstar=tstar)


# ----------------------------------------------------------------------------------------------------- selftest

def selftest():
    import tempfile
    G = X.init()
    # default-off identities: with the default target and no plane the ranking equals the Revision 23 DP of the same tables
    ody = G["GR"].anchor_cand(G["ctx"].cfg); r = select_geometry([ody])[0]; t0 = X.table_from(X.build(ody, G["mf"]["cut"]["ODY"])); a = X.dp(t0, "S-fixed")
    assert r["feasible"] and r["plane"] == 0.0 and abs(r["aero"] - a["aero"]) < 1e-12 and r["ra_target_m"] == AN.RA_TARGET and abs(r["aero"] - 1.159413233391487) < 1e-9, (r, a)
    # plane term: zero at equal inclination, positive and monotone otherwise; node-fixed refused
    assert exact_plane_increment(0.0, 7.408e6, 3.504e6) == 0.0 and exact_plane_increment(2.0, 7.408e6, 3.504e6) > exact_plane_increment(1.0, 7.408e6, 3.504e6) > 0
    rp = select_geometry([ody], plane=dict(inc_deg=63.4))[0]; assert rp["plane"] > 700 and rp["net"] > rp["aero"]
    try:
        select_geometry([ody], plane=dict(inc_deg=63.4, node_deg=10.0)); raise SystemExit("node-fixed must be refused")
    except ValueError:
        pass
    # candidate lookup by elements: ODY not on the grid -> cut anchor sample; a fine-only geometry (80 deg) found in the fine sample
    c80 = next(c for c in G["fine"] if c["inc"] == 80.0); f, src = find_field(c80); assert f is not None and src.startswith("fine:")
    # sampling request: prefix property, counts, machine-readable manifest
    with tempfile.TemporaryDirectory() as d:
        m = sampling_request(d); assert set(m["candidates"]) == {"G0486_A2000", "G1230_A1000", "G1038_A0500"} and m["calls_total"] == sum(v["calls"] for v in m["candidates"].values())
        assert all(v["extra_passes"] > 0 for v in m["candidates"].values()); df = pd.read_csv(pathlib.Path(d) / "tracks_rev25low.csv"); assert len(df) * 3 == m["calls_total"] and list(df.columns)[:2] == ["cand", "epoch_utc"]
        # build on preview data into a scratch dir (the test of the post-sampling path)
        res = build_lowalt(None, d, "preview"); T = json.load(open(pathlib.Path(d) / "planner_ground_tables_rev25.json"))["geometries"]
        assert set(T) == {"G0486", "G1230", "G1038"} and all(T[k]["token"] == k and T[k]["field_mode"] == "preview_no_field" and abs(T[k]["ra"][-1] - T[k]["ra_target_m"]) < 6e4 + 5e4 for k in T) and all(T[k]["ra"][-1] < radius_of(BLOCK[k]["alt_km"]) + 1e5 for k in T)
        H = json.load(open(pathlib.Path(d) / "rev20_data" / "ic" / "headline_rev25.json"))["plans"]; assert len(H) == 36 and all(p["ra_target_m"] > 3.8e6 for p in H)
        assert all(math.isfinite(v["guard_tstar_d"]) for v in res["tstar"].values()), res["tstar"]
    assert R.RA_T == AN.RA_TARGET and R.HORIZON_D == 130.0 and G["ctx"].cfg["orbit"]["ra_target_m"] == 7408000.0
    print("planner_select_geometry selftest ok")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["selftest", "request", "build"]); ap.add_argument("--preview", action="store_true"); ap.add_argument("--samples"); ap.add_argument("--out-dir"); a = ap.parse_args()
    if a.cmd == "selftest": selftest()
    elif a.cmd == "request": print(json.dumps(sampling_request(a.out_dir), indent=1)[:1500])
    else:
        if not a.preview and not a.samples: ap.error("build needs --preview or --samples DIR")
        build_lowalt(None if a.preview else a.samples, a.out_dir, "preview" if a.preview else "mf")
