"""Venus gate V-G(V) analysis (plan R8.3, R8.4; criteria R4.3 and U19; heat rate on V_rel, R4.4, R5.1). Python only, no Julia. The Venus counterpart of stageg3_analysis.py, which is IMPORTED unchanged (A3: criteria, verdict,
burn removal, model fed with the flown pass epochs, U15 rotation identity); only the Venus context (geometry_rank_venus.make_ctx), the Venus config/table loader, the mean-field loader, the R8.3 outcome rule and the Venus
integrity checks are new. The Venus ranking and fixed-epoch files are imported unchanged too.
  selftest                                  config and table, U15 identity at Venus, mean-field loader, criteria, outcomes, burn removal and epoch feeding on the Venus model
  meanfield [--mf DIR]                      H(h) and the profile shape from the GRAM samples; floors of VG1, VG2, VGref with the constant H = 3.86 km and with the sampled H(k) (and H(k) with ln rho)
  predict   [--mf DIR] [--out F]            pre-flight predictions (stageg3v_predictions.json): drift path, range, floor, regime per geometry, deterministic and jitter ensemble, with the sampled field
  gate      --results-root R [--mf DIR]     criteria (a) to (c) per flight with the model fed the flown pass epochs, geometry verdicts, the R8.3 outcome
  integrity --results-root R                route, failures, epoch, NaN, burns, stop reasons, host stamps, scoring table (Q_s on V_rel and, for comparison, inertial)
Heat cap: NO sourced Venus value exists (config [stageg3v.heat_cap]); the scoring thresholds 0.32, 0.54, 0.65 and the ranking cap 0.234 are ASSUMED (Odyssey numbers carried).
Usage: python3 -I stageg3v_analysis.py selftest | meanfield | predict | gate --results-root R | integrity --results-root R"""
import argparse, glob, hashlib, json, math, pathlib, sys, tempfile, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import geometry_rank as GR
import geometry_rank_venus as V
import stageg3_analysis as A3

CFGV_PATH = HERE / "stageg3v_config.toml"


def loadv(path=CFGV_PATH):
    """Config with the geometry list from the Venus table (sha256 verified against [stageg3v.table].sha256; a mismatch is an ERROR). cfg["stageg3v"]["geometry"] rows: id, role, inc_deg, omega_deg, raan_deg, epoch_off_d, epoch_utc, field_id."""
    cfg = tomllib.loads(pathlib.Path(path).read_text()); t = cfg["stageg3v"]["table"]; raw = (V.notes_dir() / t["path"]).read_bytes(); h = hashlib.sha256(raw).hexdigest()
    if h != t["sha256"]: raise ValueError(f"geometry table sha256 mismatch: {t['path']} is {h}, config has {t['sha256']}")
    T = json.loads(raw); st = cfg["stageg3v"]["start"]; ts = T["start"]
    if (ts["ra0_m"], ts["ra_target_m"], ts["true_anomaly_deg"]) != (st["ra0_m"], st["ra_target_m"], st["true_anomaly_deg"]): raise ValueError("table start differs from [stageg3v.start]")
    cfg["stageg3v"]["geometry"] = [dict(id=g["id"], role=g["role"], inc_deg=g["inclination_deg"], omega_deg=g["omega0_deg"], raan_deg=g["raan0_deg"], epoch_off_d=g["epoch_offset_days"], epoch_utc=g["epoch_utc"], field_id=g["id"])
                                   for g in T["geometries"] if g["in_default_set"]]
    cfg["stageg3v"]["table_hash"] = h; cfg["stageg3"] = cfg["stageg3v"]      # A3.score_flight reads cfg["stageg3"]["criteria"]
    return cfg


def make_ctx(cfgv):
    ctx, info = V.make_ctx(); ctx.ratio = float(cfgv["stageg3v"]["burn_response_ratio"])
    assert abs(ctx.cfg["orbit"]["ra_target_m"] - cfgv["stageg3v"]["start"]["ra_target_m"]) < 1e-3, "ranking campaign differs from the table start"
    return ctx, info


def load_field(cfgv, mf_dir=None):
    """GRAM mean field per geometry id: dict(H per pass, lnrho0 per pass, k) from mf_*.csv (geometry_rank.load_meanfield with the Venus altitude offsets of [stageg3v.sampling])."""
    d = pathlib.Path(mf_dir) if mf_dir else V.notes_dir() / cfgv["stageg3v"]["criteria"]["meanfield_dir"]
    cfg = dict(GR.load_cfg()); cfg["sampling"] = dict(cfg["sampling"], alt_offsets_m=list(cfgv["stageg3v"]["sampling"]["alt_offsets_m"]))
    return GR.load_meanfield(d, cfg), str(d)


def field_of(F, g, mode="full"):
    """Field dict of one geometry: mode 'full' (H and ln rho), 'H' (sampled H(k), field channel zero), 'none' (constant H)."""
    if mode == "none" or g["field_id"] not in F: return None
    f = F[g["field_id"]]; return f if mode == "full" else dict(H=f["H"], lnrho0=np.zeros(len(f["lnrho0"])), k=f["k"])


# ---- predictions ------------------------------------------------------------------------------------------------------------------------------

def predict_geometry(ctx, g, field):
    """Deterministic and jitter-ensemble prediction of one geometry (ranking model): drift path (slow ln q), range, floor, regime."""
    c = A3.cand_of(g); r, n, S, d = GR.evaluate(ctx, c, field, ctx.s_t_ody); path = GR.path_of(d); slow, fast = GR.split_slow_fast(path, ctx.cfg["band"]["fast_window_passes"])
    jc = ctx.cfg["floor"]; ens = []
    for j in np.linspace(-jc["phase_jitter_frac"], jc["phase_jitter_frac"], int(jc["phase_jitter_n"])):
        rj, _, _, dj = GR.evaluate(ctx, c, field, ctx.s_t_ody, j); sj = GR.split_slow_fast(GR.path_of(dj), ctx.cfg["band"]["fast_window_passes"])[0]
        ens.append(dict(jitter=float(j), floor_mps=rj["floor_mps"], R_ln=rj["R"], regime=rj["regime"], n_passes=rj["n_passes"], slow_path_ln=[float(x) for x in sj]))
    fl = np.array([e["floor_mps"] for e in ens]); Rs = np.array([e["R_ln"] for e in ens])
    return dict(id=g["id"], role=g["role"], inc_deg=g["inc_deg"], omega0_deg=g["omega_deg"], raan0_deg=g["raan_deg"], epoch_utc=g["epoch_utc"], field_mode=r["field_mode"],
                deterministic=dict(floor_mps=r["floor_mps"], R_ln=r["R"], N_ln=r["N"], regime=r["regime"], q_req_w_cm2=r["q_req"], W_eff=r["W_eff"], W=r["W"], fast_w=r["fast_w"], n_passes=r["n_passes"], days=r["days"],
                                   floor_driftfree_mps=r["floor_driftfree_mps"], slow_path_ln=[float(x) for x in slow], t_s=[float(x) for x in S["t"]], H_m=[float(x) for x in S["H"]]),
                ensemble=dict(floor_mean_mps=float(fl.mean()), floor_sd_mps=float(fl.std()), floor_min_mps=float(fl.min()), floor_max_mps=float(fl.max()), R_mean_ln=float(Rs.mean()),
                              regimes={k: int(v) for k, v in pd.Series([e["regime"] for e in ens]).value_counts().items()}, members=ens))


def cmd_predict(cfgv, mf_dir, out):
    ctx, info = make_ctx(cfgv); F, src = load_field(cfgv, mf_dir); rows = []
    for g in cfgv["stageg3v"]["geometry"]:
        r = predict_geometry(ctx, g, field_of(F, g, "full")); r["meanfield_source"] = src; r["constant_H_m"] = ctx.cfg["decay"]["h_fallback_m"]
        c0 = predict_geometry(ctx, g, None); r["constant_H"] = dict(floor_det_mps=c0["deterministic"]["floor_mps"], floor_ens_mps=c0["ensemble"]["floor_mean_mps"], R_ln=c0["deterministic"]["R_ln"], regime=c0["deterministic"]["regime"])
        rows.append(r); print(json.dumps({k: r[k] for k in ("id", "role")} | dict(det=r["deterministic"]["floor_mps"], ens=round(r["ensemble"]["floor_mean_mps"], 4), R=round(r["deterministic"]["R_ln"], 3), regime=r["deterministic"]["regime"], const=r["constant_H"]), default=float))
    res = dict(campaign=info, geometry_table_sha256=cfgv["stageg3v"]["table_hash"], config="stageg3v_config.toml", burn_response_ratio=ctx.ratio, geometries=rows)
    if out: pathlib.Path(out).write_text(json.dumps(res, indent=1, default=float))
    return res


# ---- mean field summary (deliverable 1) -------------------------------------------------------------------------------------------------------

def profile_H(mf_dir):
    """Local scale height H(h) = -dh / d ln rho by central differences over +-2 profile steps, per geometry and pass; returns (altitudes, array [pass-rows, altitudes], lnrho array, geometry label per row)."""
    rows = []; lab = []; alts = None
    for f in sorted(glob.glob(str(pathlib.Path(mf_dir) / "pf_*.csv"))):
        d = pd.read_csv(f)
        for k, gk in d.groupby("k"):
            gk = gk.sort_values("alt_m"); alts = gk.alt_m.values; rows.append(np.log(gk.rho_kg_m3.values)); lab.append(gk.cand.iloc[0])
    L = np.array(rows); step = alts[1] - alts[0]; Hh = np.full_like(L, np.nan); Hh[:, 2:-2] = -4 * step / (L[:, 4:] - L[:, :-4]); return alts, Hh, L, np.array(lab)


def cmd_meanfield(cfgv, mf_dir, out):
    ctx, info = make_ctx(cfgv); F, src = load_field(cfgv, mf_dir); d = pathlib.Path(src); res = dict(source=src, h_ref_m=ctx.cfg["orbit"]["h_ref_m"], constant_H_m=ctx.cfg["decay"]["h_fallback_m"], alt_offsets_m=cfgv["stageg3v"]["sampling"]["alt_offsets_m"])
    st = {}
    for g in cfgv["stageg3v"]["geometry"]:
        f = F[g["field_id"]]; H = f["H"]; l0 = f["lnrho0"]
        st[g["id"]] = dict(n_passes=len(H), H_ref_median_m=float(np.median(H)), H_ref_p10_m=float(np.percentile(H, 10)), H_ref_p90_m=float(np.percentile(H, 90)), H_ref_min_m=float(H.min()), H_ref_max_m=float(H.max()), H_ref_first_m=float(H[0]), H_ref_last_m=float(H[-1]),
                           lnrho_ref_first=float(l0[0]), lnrho_ref_range=float(np.ptp(l0)), lnrho_ref_net=float(l0[-1] - l0[0]), lnrho_ref_sd=float(np.std(l0)), rho_ref_first_kg_m3=float(np.exp(l0[0])), rho_ref_median_kg_m3=float(np.exp(np.median(l0))))
    res["h_ref_sample"] = st
    alts, Hh, L, lab = profile_H(d); tab = []
    for h in np.arange(100e3, 200.1e3, 10e3):
        j = int(np.argmin(np.abs(alts - h))); x = Hh[:, j]; tab.append(dict(alt_km=float(alts[j] / 1e3), H_median_km=float(np.nanmedian(x) / 1e3), H_p10_km=float(np.nanpercentile(x, 10) / 1e3), H_p90_km=float(np.nanpercentile(x, 90) / 1e3),
                                                                           ln_rho_median=float(np.median(L[:, j])), ln_rho_sd_over_passes=float(np.std(L[:, j])), H_over_const=float(np.nanmedian(x) / ctx.cfg["decay"]["h_fallback_m"])))
    res["profile_H_of_h"] = tab
    # profile shape: departure of ln rho from a single exponential through h_ref over the fit window of [common.fit]
    jr = int(np.argmin(np.abs(alts - ctx.cfg["orbit"]["h_ref_m"]))); below, above = cfgv["common"]["fit"]["window_below_m"], cfgv["common"]["fit"]["window_above_m"]
    w = (alts >= alts[jr] - below) & (alts <= alts[jr] + above); dev = []
    for i in range(L.shape[0]):
        p = np.polyfit(alts[w] - alts[jr], L[i, w], 1); dev.append(np.abs(L[i, w] - np.polyval(p, alts[w] - alts[jr])).max())
    res["single_exponential_fit_window"] = dict(window_km=[float((alts[w][0]) / 1e3), float(alts[w][-1] / 1e3)], max_abs_dev_lnrho_median=float(np.median(dev)), max_abs_dev_lnrho_p90=float(np.percentile(dev, 90)),
                                                H_fit_median_km=float(np.median([-1 / np.polyfit(alts[w] - alts[jr], L[i, w], 1)[0] for i in range(L.shape[0])]) / 1e3))
    # per-pass H at h_ref vs pass index (trend) and the H-versus-altitude slope in the 110 to 170 km range
    sel = (alts >= 110e3) & (alts <= 170e3); sl = [np.polyfit(alts[sel] / 1e3, Hh[i, sel] / 1e3, 1)[0] for i in range(Hh.shape[0]) if np.isfinite(Hh[i, sel]).all()]
    res["dH_dh_110_170_km_per_km_median"] = float(np.median(sl))
    # floors: constant H, sampled H(k), sampled H(k) and ln rho
    fl = []
    for g in cfgv["stageg3v"]["geometry"]:
        row = dict(id=g["id"], role=g["role"])
        for mode in ("none", "H", "full"):
            p = predict_geometry(ctx, g, field_of(F, g, mode)); row[mode] = dict(floor_det_mps=p["deterministic"]["floor_mps"], floor_ens_mps=p["ensemble"]["floor_mean_mps"], floor_ens_sd_mps=p["ensemble"]["floor_sd_mps"], R_ln=p["deterministic"]["R_ln"], N_ln=p["deterministic"]["N_ln"],
                                                                                        regime=p["deterministic"]["regime"], q_req_w_cm2=p["deterministic"]["q_req_w_cm2"], W_eff=p["deterministic"]["W_eff"], n_passes=p["deterministic"]["n_passes"], days=p["deterministic"]["days"])
        fl.append(row); print(g["id"], {m: (round(row[m]["floor_det_mps"], 3), round(row[m]["floor_ens_mps"], 3), round(row[m]["R_ln"], 3), row[m]["regime"]) for m in ("none", "H", "full")}, flush=True)
    res["floors"] = fl; txt = json.dumps(res, indent=1, default=float)
    if out: pathlib.Path(out).write_text(txt)
    print(json.dumps({k: res[k] for k in ("h_ref_sample", "profile_H_of_h", "single_exponential_fit_window", "dH_dh_110_170_km_per_km_median")}, indent=1, default=float)); return res


# ---- gate scoring -----------------------------------------------------------------------------------------------------------------------------

def sd(root): return pathlib.Path(root) / "stageg3v"
def adir(root): d = sd(root) / "analysis"; d.mkdir(parents=True, exist_ok=True); return d
def dump(root, name, res): (adir(root) / name).write_text(json.dumps(res, indent=1, default=float)); print(json.dumps(res, indent=1, default=float))


def summaries(root, phase="gate"):
    fs = sorted(glob.glob(str(sd(root) / f"summary_{phase}_venus_*.csv"))); return pd.concat([pd.read_csv(f) for f in fs], ignore_index=True) if fs else pd.DataFrame()


def outcome(vg1, vgref):
    """R8.3 outcome rule: VG1 and VGref pass -> the Venus planner claim is flown-validated at those geometries; VG1 fails -> the claim stays model-only and the failure is reported; VG2 is diagnostic (never enters).
    R8.3 does not cover VG1 passing while VGref fails: reported as its own code (ASSUMED reading: validated at VG1 only; the coordinator decides)."""
    if vg1 is None or vgref is None: return "INCOMPLETE"
    if not vg1: return "VENUS_CLAIM_MODEL_ONLY (VG1 fails; failure reported)"
    return "VENUS_CLAIM_FLOWN_VALIDATED_AT_VG1_AND_VGREF" if vgref else "VG1_PASSES_VGREF_FAILS (not covered by R8.3; ASSUMED: validated at VG1 only, coordinator decides)"


def cmd_gate(cfgv, root, mf_dir):
    ctx, _ = make_ctx(cfgv); S = summaries(root, "gate"); G3 = cfgv["stageg3v"]; ok = S[S.status == "ok"] if len(S) else S; per = {}; recs = []; F, src = load_field(cfgv, mf_dir)
    for g in G3["geometry"]:
        field = field_of(F, g, "full"); rows = ok[ok.geometry == g["id"]] if len(ok) else ok; per[g["id"]] = []
        for r in rows.itertuples():
            rec = A3.score_flight(ctx, cfgv, g, field, sd(root) / "runs" / r.run_id, float(G3["burn_response_ratio"]))
            rec.update(geometry=g["id"], role=g["role"], keeper=r.keeper, tag=r.tag, stop_reason=r.stop_reason, meanfield=src, field_mode="meanfield" if field is not None else "NO_FIELD")
            recs.append(rec); per[g["id"]].append(None if rec.get("insufficient") else rec)
    ver = {g["id"]: A3.verdict(per[g["id"]], G3["criteria"]["flights_rule"]) for g in G3["geometry"]}
    res = dict(flights=recs, verdict={gid: dict(passed=v, role=next(g["role"] for g in G3["geometry"] if g["id"] == gid)) for gid, v in ver.items()}, outcome=outcome(ver["VG1"], ver["VGref"]),
               vg2=dict(passed=ver["VG2"], note="diagnostic only; cannot block"), flights_rule=G3["criteria"]["flights_rule"], criteria=G3["criteria"], burn_response_ratio=G3["burn_response_ratio"],
               heat_cap_source=G3["heat_cap"]["source"])
    dump(root, "gate_verdict.json", res); return res


def qs_inertial(rho_peak, rp, ra, mu):
    a = 0.5 * (np.asarray(rp, float) + np.asarray(ra, float)); vp = np.sqrt(mu * (2 / np.asarray(rp, float) - 1 / a)); return 0.5 * np.asarray(rho_peak, float) * vp ** 3 / 1e4


def cmd_integrity(cfgv, root):
    S = summaries(root, "gate"); G3 = cfgv["stageg3v"]; out = dict(n_rows=len(S), expected=G3["plan_counts"]["gate"], problems=[], burn_response_ratio=G3["burn_response_ratio"], heat_cap_source=G3["heat_cap"]["source"]); P = out["problems"]
    if G3["burn_response_ratio"] != 1.0: P.append("burn_response_ratio is not 1.0")
    if not len(S): P.append("no gate summary rows"); dump(root, "integrity.json", out); return out
    out["rows_ok"] = int((S.status == "ok").sum()); (S.status != "ok").any() and P.append("rows not ok: " + ", ".join(S[S.status != "ok"].run_id))
    len(S) != G3["plan_counts"]["gate"] and P.append(f"row count {len(S)} != {G3['plan_counts']['gate']}")
    hosts = []
    for f in sorted(glob.glob(str(sd(root) / "mc_*_venus_*.toml"))):
        t = tomllib.loads(pathlib.Path(f).read_text()); t.get("mc_route") != "process" and P.append(f"route {t.get('mc_route')} in {pathlib.Path(f).name}, not process"); t.get("mc_failed", 0) and P.append(f"mc_failed={t['mc_failed']}")
        hosts.append(dict(file=pathlib.Path(f).name, host=t.get("host"), route=t.get("mc_route"), workers=t.get("mc_workers")))
    out["host_stamps"] = hosts; st = [h["host"] for h in hosts]; len(set(st)) != len(st) and P.append("host stamps are not distinct across the campaign files: " + ", ".join(map(str, st)))
    out["summary_files"] = [pathlib.Path(f).name for f in sorted(glob.glob(str(sd(root) / "summary_*_venus_*.csv")))]
    td = {g["id"]: g for g in G3["geometry"]}; tab = []; mu = G3["rotation"]["mu_m3s2"]
    for r in S[S.status == "ok"].itertuples():
        g = td[r.geometry]; ep = pd.Timestamp(g["epoch_utc"]); core = f"{ep.year}-{ep.month}-{ep.day} {ep.hour}:{ep.minute}:{ep.second}"
        not str(r.gram_core_epoch).startswith(core) and P.append(f"{r.run_id}: GRAM core epoch {r.gram_core_epoch} != {core}")
        r.keeper == "zero_burn" and r.burns != 0 and P.append(f"{r.run_id}: zero-burn flight has {r.burns} burns")
        r.keeper == "fixed" and r.stop_reason not in ("end_radius", "horizon") and P.append(f"{r.run_id}: keeper stop_reason {r.stop_reason}")
        r.keeper == "fixed" and r.stop_reason == "horizon" and P.append(f"{r.run_id}: keeper reached the horizon without the end radius")
        r.keeper == "zero_burn" and r.stop_reason in ("controller_stopped",) and P.append(f"{r.run_id}: zero-burn stop_reason {r.stop_reason}")
        D = pd.read_csv(sd(root) / "runs" / r.run_id / "passes.csv"); d = D.dropna(subset=["ra_m"])
        d[["rho_peak_kg_m3", "q_w_cm2", "vrel_mps", "rp_m"]].isna().any().any() and P.append(f"{r.run_id}: NaN in rho_peak / q / vrel / rp")
        qr = A3.qs_rel(d.rho_peak_kg_m3, d.vrel_mps); qi = qs_inertial(d.rho_peak_kg_m3, d.rp_m, d.ra_m, mu)
        tab.append(dict(run_id=r.run_id, geometry=r.geometry, role=r.role, keeper=r.keeper, passes=len(D), days=r.t_end_days, dv_mps=r.total_dv_mps, burns=r.burns, stop_reason=r.stop_reason, safety_reason=r.safety_reason,
                        qs_rel_max=float(qr.max()), qs_inertial_max=float(qi.max()), qs_rel_over_032=float((qr > 0.32).mean()), qs_rel_over_054=float((qr > 0.54).mean()), ratio_rel_over_inertial_median=float(np.median(qr / qi)), any_rel_over_065=bool((qr > 0.65).any())))
    out["scoring"] = tab; dump(root, "integrity.json", out); return out


# ---- selftest ---------------------------------------------------------------------------------------------------------------------------------

def selftest():
    cfgv = loadv(); G3 = cfgv["stageg3v"]; cr = G3["criteria"]; gs = G3["geometry"]
    assert [g["id"] for g in gs] == ["VG1", "VG2", "VGref"] and [g["role"] for g in gs] == ["VG1", "VG2", "VGref"] and [g["inc_deg"] for g in gs] == [90.0, 63.4, 89.876]
    assert all(g["epoch_utc"] == "2014-05-19T04:00:00" and g["epoch_off_d"] == 0.0 for g in gs) and G3["base_epoch_utc"] == "2014-05-19T04:00:00"
    assert G3["burn_response_ratio"] == 1.0 and G3["decision_speed"] == "relative" and G3["seeds"] == [101, 102] and G3["zero_burn_roles"] == ["VG1"]
    assert (cr["min_correlation"], cr["range_factor"], cr["range_abs_ln"], cr["floor_abs_mps"], cr["floor_rel"], cr["smooth_passes"]) == (0.85, 1.5, 0.15, 0.3, 0.30, 15)
    assert G3["safety"]["q_stop_w_cm2"] == 0.65 and G3["safety"]["n_hloc"] == 3.0 and G3["heat_cap"]["source"].startswith("NONE")
    n_keeper = len(gs) * G3["keeper_seeds_per_geometry"]; n_zb = sum(g["role"] in G3["zero_burn_roles"] for g in gs) * G3["zero_burn_seeds"]; pc = G3["plan_counts"]; assert (n_keeper, n_zb, pc["gate"]) == (6, 2, 8) == (pc["keeper"], pc["zero_burn"], pc["gate"])
    with tempfile.TemporaryDirectory() as d:
        p = pathlib.Path(d) / "bad.toml"; p.write_text(CFGV_PATH.read_text().replace(G3["table"]["sha256"], "0" * 64))
        try: loadv(p); raise SystemExit("hash mismatch was not rejected")
        except ValueError as e: assert "mismatch" in str(e)
    ctx, info = make_ctx(cfgv); cfgr = ctx.cfg; assert cfgr["floor"]["burn_response_ratio"] == G3["burn_response_ratio"] and cfgr["band"]["fast_window_passes"] == cr["smooth_passes"]
    # table floors reproduce: the constant-H deterministic floors of the three geometries equal the table
    T = json.loads((V.notes_dir() / G3["table"]["path"]).read_text())
    for g, tg in zip(gs, T["geometries"]):
        r = GR.evaluate(ctx, A3.cand_of(g), None, None)[0]; assert abs(r["floor_mps"] - tg["floor_det_mps"]) < 6e-5 and r["n_passes"] == tg["n_passes_model"], (g["id"], r["floor_mps"], tg["floor_det_mps"])
    # U15 rotation identity at Venus (independent numpy implementation inside stageg3_analysis)
    mu, w = G3["rotation"]["mu_m3s2"], G3["rotation"]["w_rot_rad_s"]; worst = worst_v = 0.0; cnt = 0; rp = ctx.Req + ctx.cfg["orbit"]["h_ref_m"]
    for inc in (39.2, 63.4, 89.876, 90.0, 116.6, 150.0):
        for raan in (0.0, 90.0, 270.0):
            for argp in (0.0, 45.0, 90.0, 200.0):
                al, wr, vr, cl = A3.rotation_identity(inc, raan, argp, rp, G3["start"]["ra0_m"], mu, w); worst = max(worst, abs(al - wr)); worst_v = max(worst_v, abs(vr - cl)); cnt += 1
    assert worst < 1e-6 and worst_v < 1e-6, (worst, worst_v)
    assert abs(float(A3.qs_rel(2 * 0.65e4 / 9500 ** 3, 9500.0)) - 0.65) < 1e-12 and abs(float(qs_inertial(1e-8, 6.2e6, 7.2e7, mu)) - 0.5e-8 * (mu * (2 / 6.2e6 - 2 / 7.82e7)) ** 1.5 / 1e4) < 1e-12
    # mean-field loader: H and ln rho from the Venus samples equal an independent secant on the raw file; the field has one row per pass of the track; the profile H(h) is finite
    F, src = load_field(cfgv); d = pathlib.Path(src); off = G3["sampling"]["alt_offsets_m"]; assert off == [-2000.0, 0.0, 2000.0]
    for g in gs:
        raw = pd.read_csv(d / f"mf_{g['id']}.csv"); k = 7; rr = raw[raw.k == k].sort_values("alt_idx").rho_kg_m3.values; f = F[g["id"]]
        assert abs(f["H"][k] - (-(off[2] - off[0]) / (np.log(rr[2]) - np.log(rr[0])))) < 1e-9 and abs(f["lnrho0"][k] - np.log(rr[1])) < 1e-12 and len(f["H"]) == T["geometries"][0]["n_passes_model"] and np.isfinite(f["H"]).all() and (f["H"] > 0).all()
    alts, Hh, L, lab = profile_H(d); assert np.isfinite(Hh[:, 2:-2]).all() and (Hh[:, 2:-2] > 0).all() and set(lab) == {"VG1", "VG2", "VGref"}
    # criteria and outcomes (the thresholds are the Mars gate's; stageg3_analysis.selftest covers them in full, repeated here on the Venus config)
    t = np.linspace(0, 1, 200); a = 0.3 * np.sin(2 * np.pi * 1.5 * t); good = A3.criteria(a, 1.1 * a, 0.5, 0.6, 0.1, 0.46, cr, 0.05); assert good["pass"] and not A3.criteria(a, -a, 0.5, 0.5, 0.1, 0.46, cr, 0.05)["a_pass"]
    assert not A3.criteria(a, 3 * a, 0.5, 0.5, 0.1, 0.46, cr, 0.05)["b_pass"] and A3.criteria(0.05 * a, 0.08 * a, 0.5, 0.5, 0.1, 0.46, cr, 0.05)["b_pass"]
    assert outcome(True, True).startswith("VENUS_CLAIM_FLOWN_VALIDATED") and outcome(False, True).startswith("VENUS_CLAIM_MODEL_ONLY") and outcome(False, False).startswith("VENUS_CLAIM_MODEL_ONLY") and outcome(True, False).startswith("VG1_PASSES") and outcome(None, True) == "INCOMPLETE"
    # model feeding and burn removal on VG1 with the sampled field: feeding the model's own epochs reproduces its path; burn removal inverts the lever weights; synthetic flights score as expected
    g1 = gs[0]; fld = field_of(F, g1, "full"); cand = A3.cand_of(g1); full = A3.model_path(ctx, cand, fld); n0 = min(full["n"], 50)
    own = A3.model_path(ctx, cand, fld, full["S"]["t"][:n0]); assert own["n"] == n0 and np.allclose(own["path"], full["path"][:n0], atol=1e-12)
    skew = A3.model_path(ctx, cand, fld, 1.03 * full["S"]["t"][:n0]); assert skew["n"] == n0
    with tempfile.TemporaryDirectory() as td:
        Sf = A3._trunc(full["S"], n0); noise = np.concatenate([[0.0], np.cumsum(np.diff(full["path"][:n0]))]); burns = {5: 0.07, 6: -0.04, 20: 0.1}; w_ = A3._fake_run(pathlib.Path(td) / "k", ctx, cand, Sf, burns, noise)
        P_, nat, wf = A3.flown_natural(ctx, pathlib.Path(td) / "k", G3["burn_response_ratio"]); assert np.allclose(nat - nat[0], noise - noise[0], atol=1e-9) and np.allclose(wf, w_)
        P2, nat2, w2 = GR.dev_natural(ctx, str(pathlib.Path(td) / "k"), G3["burn_response_ratio"]); assert np.allclose(nat2, nat) and np.allclose(w2, wf)
        A3._fake_run(pathlib.Path(td) / "z", ctx, cand, Sf, {}, noise); small = dict(cfgv); small["stageg3"] = dict(G3, criteria=dict(cr, min_passes_scored=30))
        rec = A3.score_flight(ctx, small, g1, fld, pathlib.Path(td) / "z", 1.0); assert not rec["insufficient"] and rec["a_pass"] and rec["b_pass"] and rec["corr"] > 0.999, rec
        A3._fake_run(pathlib.Path(td) / "bad", ctx, cand, Sf, {}, -noise); assert not A3.score_flight(ctx, small, g1, fld, pathlib.Path(td) / "bad", 1.0)["a_pass"]
        A3._fake_run(pathlib.Path(td) / "short", ctx, cand, A3._trunc(full["S"], 12), {}, noise[:12]); assert A3.score_flight(ctx, small, g1, fld, pathlib.Path(td) / "short", 1.0)["insufficient"]
    print(f"stageg3v selftest ok: config and table (8 flights = 6 keeper + 2 zero-burn, burn_response_ratio 1.0, heat cap flagged ASSUMED), table floors reproduced, U15 rotation identity at Venus (max error {worst:.1e} m/s over {cnt} orbits), "
          f"V_rel closed form ({worst_v:.1e}), mean-field loader on the Venus samples, criteria, R8.3 outcomes, burn removal = geometry_rank.dev_natural, epoch feeding identity, synthetic flights")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["selftest", "meanfield", "predict", "gate", "integrity"]); ap.add_argument("--results-root"); ap.add_argument("--mf"); ap.add_argument("--out"); ap.add_argument("--config", default=str(CFGV_PATH))
    a = ap.parse_args()
    if a.cmd == "selftest": return selftest()
    cfgv = loadv(a.config)
    if a.cmd == "meanfield": return cmd_meanfield(cfgv, a.mf, a.out)
    if a.cmd == "predict": return cmd_predict(cfgv, a.mf, a.out)
    a.results_root or ap.error("--results-root required")
    return cmd_gate(cfgv, a.results_root, a.mf) if a.cmd == "gate" else cmd_integrity(cfgv, a.results_root)


if __name__ == "__main__": main()
