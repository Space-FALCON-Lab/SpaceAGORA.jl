#!/usr/bin/env python3 -I
"""Stage O, WP-O1: the flown Mars Odyssey reference for the periapsis-corridor comparison.

Commands (run with `python3 -I`):
  build --dumps DIR   copy the two JAIS tables, merge the Arrow dumps from
                      stageo_extract_apsides.jl, write kernel_apsides.csv and flown_reference.json
  selftest            recompute everything from the files in stageo_data and assert the
                      plan's DERIVED numbers exactly; the Table III cross-check and the
                      nav-state cross-check are printed and stored, never asserted

Labels: SOURCED (file), DERIVED (computed here from sourced values), ASSUMED (named choice).
Data dir: $STAGEO_DATA_DIR, else <main checkout>/tasks/2026-10-15_periapsis-region-border-targeting/stageo_data.
JAIS source dir (build only): $JAIS_ODY_DIR (data/raw/odyssey_closed_loop_20261008 of the JAIS repo).
"""
import csv, hashlib, json, math, os, shutil, statistics, subprocess, sys, tomllib
from datetime import datetime, timedelta, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
TASK = "tasks/2026-10-15_periapsis-region-border-targeting"
R_KM, MU = 3396.19, 42828.375            # SOURCED plan R1.2
MP1, MP2 = (19, 126), (127, 171)         # SOURCED Smith and Bell 2005 Table 1 (plan 2.1)
PHASES = {"MP1": (19, 126), "MP2": (127, 171), "EG": (172, 248), "WO": (249, 336)}
# SOURCED Falcone and Putnam 2023 Table III via apoapsis_pi_preregistration.md (plan 2.1)
T3_RADIUS_KM = {"MP1": 28533.0, "MP2": 10038.0, "EG": 7408.0, "WO": 4904.0, "END": 3900.0}  # boundary radii: start of MP1, MP2, EG, WO, end
T3_PASSES = {"MP1": 108, "MP2": 45, "EG": 77, "WO": 88}
T3_DAYS = {"MP1": 41.8, "MP2": 7.2, "EG": 8.7, "WO": 7.4}
MP_END_R_KM = 7408.0
BINS_MM = [(20.0, 1e9), (15.0, 20.0), (10.0, 15.0), (7.4, 10.0)]   # apoapsis radius bins (Mm), SOURCED plan R1.4
ODY_TIERS = {"over_corridor_0.32": ("gt", 0.32), "over_allowable_0.54": ("gt", 0.54),
             "over_qualification_0.65": ("gt", 0.65), "low_below_0.14": ("lt", 0.14)}
PAPER_TIERS = {"over_soft_0.25": ("gt", 0.25), "over_medium_0.30": ("gt", 0.30),
               "over_hard_0.45": ("gt", 0.45), "low_below_0.05": ("lt", 0.05)}
TOL = 0.10                               # DERIVED proxy tolerance, plan 2.3 item 2


def data_dir():
    if os.environ.get("STAGEO_DATA_DIR"):
        return Path(os.environ["STAGEO_DATA_DIR"])
    common = subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "--path-format=absolute",
                                      "--git-common-dir"], text=True).strip()
    return Path(common).parent / TASK / "stageo_data"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rows(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------- loading
def load(D):
    q = {int(r["P"]): (float(r["rho_max"]), float(r["q_flown_w_cm2"])) for r in rows(D / "q_per_pass_flown.csv")}
    burns = [dict(P=int(r["P"]), dv=float(r["dv_mps"]), phase=r["phase"],
                  epoch=float(r["epoch_s"]) if r["epoch_s"] else None,
                  pre=r["pre_epoch"] == "True", term=r["terminal_raise"] == "True")
             for r in rows(D / "burn_table_flown.csv")]
    ap = {int(r["orbit"]): r for r in rows(D / "kernel_apsides.csv")}
    return q, burns, ap


def derive_epochs(ap):
    """DERIVED: Keplerian period from the kernel apsides; epochs are seconds after the
    P19 periapsis (the feather files carry no time column). Apoapsis k follows periapsis k."""
    t, out = 0.0, {}
    for k in sorted(ap):
        rp, ra = R_KM + float(ap[k]["periapsis_alt_km"]), R_KM + float(ap[k]["apoapsis_alt_km"])
        T = 2 * math.pi * math.sqrt(((rp + ra) / 2) ** 3 / MU)
        out[k] = dict(rp=rp, ra=ra, period=T, t_peri=t, t_apo=t + T / 2)
        t += T
    return out


def frac(vals, op, thr):
    n = sum(1 for v in vals if (v > thr if op == "gt" else v < thr))
    return n, n / len(vals)


def tiers(vals, spec):
    out = {}
    for name, (op, thr) in spec.items():
        n, f = frac(vals, op, thr)
        # +-10% heat-rate interval: recompute at T/1.1 and T/0.9 (plan R1.3)
        a, b = frac(vals, op, thr / 1.1)[1], frac(vals, op, thr / 0.9)[1]
        out[name] = dict(count=n, fraction=f, interval=[min(a, b), max(a, b)])
        if op == "lt":
            # The plan's R1.3 row "< 0.14" ([0.078, 0.157]) uses T*1.1 for the upper end, not the
            # T/0.9 its own rule gives (0.170). Both kept; the discrepancy is reported by the selftest.
            out[name]["interval_plan_table"] = [frac(vals, op, thr / 1.1)[1], frac(vals, op, thr * 1.1)[1]]
    return out


def lnstd(q, ks, burn_orbits):
    """std of ln(Q_k/Q_k-1) over consecutive pairs with no burn at the apoapsis between them."""
    d = [math.log(q[k][1] / q[k - 1][1]) for k in ks if k in q and k - 1 in ks and k - 1 in q and (k - 1) not in burn_orbits]
    return statistics.stdev(d), len(d)


def metrics(D):
    q, burns, ap = load(D)
    E = derive_epochs(ap)
    mp = list(range(19, 172))
    Q = [q[k][1] for k in mp]
    mpb = [b for b in burns if not b["pre"] and not b["term"] and 19 <= b["P"] <= 171]
    bo = {b["P"] for b in burns if not b["pre"]}
    s_mp, n_mp = lnstd(q, set(mp), bo)
    s_eg, n_eg = lnstd(q, set(range(172, 249)), bo)
    # proxy recomputation from kernel apsides (R1.2): checks the CSV against its own definition
    dev = 0.0
    for k in q:
        rp, ra = E[k]["rp"], E[k]["ra"]
        V = math.sqrt(MU * (2 / rp - 2 / (rp + ra))) * 1e3
        dev = max(dev, abs(0.5 * q[k][0] * V ** 3 / 1e4 / q[k][1] - 1))
    # T1-T3 from the kernel
    pairs = [(k, E[k]["ra"] - E[k + 1]["ra"], E[k + 1]["t_apo"] - E[k]["t_apo"]) for k in range(19, 171)]
    bins = {}
    for lo, hi in BINS_MM:
        sel = [(d, dt) for k, d, dt in pairs if lo * 1e3 < E[k]["ra"] <= hi * 1e3]
        key = f"{lo:g}-{hi:g}Mm" if hi < 1e8 else f">{lo:g}Mm"
        bins[key] = dict(n=len(sel), T2_median_km_per_pass=statistics.median(d for d, _ in sel),
                         T3_median_km_per_day=statistics.median(d / dt * 86400 for d, dt in sel))
    t_first = E[19]["t_peri"]
    days = lambda k: (E[k]["t_apo"] - t_first) / 86400.0
    kernel_days_to_171 = days(171)
    # Table III cross-check, per boundary (reported, not asserted)
    cross, cum_p, cum_d = [], 0, 0.0
    for name, (a, b) in PHASES.items():
        cum_p += T3_PASSES[name]; cum_d += T3_DAYS[name]
        r_tab = T3_RADIUS_KM[{"MP1": "MP2", "MP2": "EG", "EG": "WO", "WO": "END"}[name]]
        k_last = b
        first_le = next((k for k in sorted(E) if k >= 19 and E[k]["ra"] <= r_tab), None)
        cross.append(dict(phase=name, table3_end_radius_km=r_tab, table3_cum_passes=cum_p, table3_cum_days=cum_d,
                          sb_last_orbit=k_last, kernel_apo_radius_at_last_orbit_km=E[k_last]["ra"],
                          radius_dev_pct=100 * (E[k_last]["ra"] / r_tab - 1),
                          first_apo_le_table3_radius_orbit=first_le, passes_by_first_le_rule=None if first_le is None else first_le - 18,
                          passes_sb_range=k_last - 18, passes_match_exactly=(first_le is not None and first_le - 18 == cum_p),
                          radius_within_1pct=abs(E[k_last]["ra"] / r_tab - 1) <= 0.01,
                          kernel_days_peri19_to_apo_last=days(k_last)))
    # burn epoch anchors: the manifest scenario epoch is P19 periapsis + 16 min (plan 2.1)
    anch = [dict(P=b["P"], flown_epoch_s=b["epoch"], derived_epoch_s=E[b["P"]]["t_apo"] - 960.0,
                 diff_s=E[b["P"]]["t_apo"] - 960.0 - b["epoch"])
            for b in burns if b["epoch"] is not None and b["P"] in E]
    # main-phase anchors only (P <= 171): the later ones jump by ~one day (see report); the drift
    # of the Keplerian sum is ~linear in orbit number, so extend the last two anchors to P171
    am = [a for a in anch if a["P"] <= 171]
    s2 = (am[-1]["diff_s"] - am[-2]["diff_s"]) / (am[-1]["P"] - am[-2]["P"])
    d171 = am[-1]["diff_s"] + s2 * (171 - am[-1]["P"])
    return dict(
        passes=len(mp), mean=statistics.mean(Q), median=statistics.median(Q), max=max(Q),
        n_gt_025=frac(Q, "gt", .25)[0], n_gt_032=frac(Q, "gt", .32)[0], n_gt_045=frac(Q, "gt", .45)[0],
        n_gt_054=frac(Q, "gt", .54)[0], n_lt_014=frac(Q, "lt", .14)[0],
        p20_q=q[20][1], p19_q=q[19][1],
        main_burns=len(mpb), main_dv_mps=sum(abs(b["dv"]) for b in mpb),
        main_burn_orbits=[b["P"] for b in mpb],
        lnq_std_mp=s_mp, lnq_pairs_mp=n_mp, lnq_white_equiv_mp=s_mp / math.sqrt(2),
        lnq_std_eg=s_eg, lnq_pairs_eg=n_eg,
        proxy_recompute_max_rel_dev=dev,
        tiers_odyssey=tiers(Q, ODY_TIERS), tiers_paper=tiers(Q, PAPER_TIERS),
        gt_025_to_065_table={str(t): tiers(Q, {"x": ("gt", t)})["x"] for t in (.25, .30, .32, .45, .54, .65)},
        lt_table={str(t): tiers(Q, {"x": ("lt", t)})["x"] for t in (.14, .05)},
        iso_risk_eps=frac(Q, "gt", .32)[1],
        kernel_mp_end_orbit=171, kernel_mp_end_radius_km=E[171]["ra"],
        kernel_first_apo_le_7408_orbit=next(k for k in sorted(E) if E[k]["ra"] <= MP_END_R_KM),
        kernel_p19_apo_radius_km=E[19]["ra"], kernel_p19_peri_alt_km=E[19]["rp"] - R_KM,
        T1_kernel=dict(passes=len(mp), days_derived_peri19_to_apo171=kernel_days_to_171,
                       days_derived_from_scenario_epoch=(E[171]["t_apo"] - E[19]["t_peri"] - 960.0) / 86400,
                       days_anchor_corrected_from_scenario_epoch=(E[171]["t_apo"] - E[19]["t_peri"] - 960.0 - d171) / 86400,
                       anchor_drift_s_per_orbit=s2, anchor_extrapolated_diff_at_171_s=d171,
                       table3_days=T3_DAYS["MP1"] + T3_DAYS["MP2"], table3_passes=153),
        T2_T3_bins=bins,
        T2_T3_campaign_mean=dict(T2_km_per_pass=statistics.mean(d for _, d, _ in pairs),
                                 T3_km_per_day=statistics.mean(d / dt * 86400 for _, d, dt in pairs), pairs=len(pairs)),
        table3_crosscheck=cross, burn_epoch_anchor_check=anch,
        burn_epoch_anchor_max_abs_diff_s_main=max(abs(a["diff_s"]) for a in am),
        burn_epoch_anchor_max_abs_diff_s_all=max(abs(a["diff_s"]) for a in anch),
    )


# ---------------------------------------------------------------- nav-state cross-check (report only)
def nav_elements(manifest):
    sc = tomllib.loads(Path(manifest).read_text())["scenarios"][0]
    s = sc["initial_state_j2000_m"]
    r = [x / 1e3 for x in s[:3]]; v = [x / 1e3 for x in s[3:]]
    t = sc["initial_time"]
    utc = datetime(t["year"], t["month"], t["day"], t["hour"], t["minute"], int(t["second"]), tzinfo=timezone.utc)
    # TDB ~ UTC + 32 leap seconds (2001) + 32.184 s (ASSUMED: periodic TDB-TT term < 2 ms ignored)
    tdb = utc + timedelta(seconds=64.184)
    T = ((tdb - datetime(2000, 1, 1, 12, tzinfo=timezone.utc)).total_seconds() / 86400) / 36525
    # SOURCED pck00011.tpc BODY499_POLE_RA/DEC (nutation/precession terms ignored: declared)
    ra, de = map(math.radians, (317.269202 - 0.10927547 * T, 54.432516 - 0.05827105 * T))
    z = [math.cos(de) * math.cos(ra), math.cos(de) * math.sin(ra), math.sin(de)]
    nx = [-z[1], z[0], 0.0]; n = math.hypot(*nx[:2]); x = [c / n for c in nx]
    y = [z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0]]
    dot = lambda a, b: sum(p * q for p, q in zip(a, b))
    cross = lambda a, b: [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    rb = [dot(r, x), dot(r, y), dot(r, z)]; vb = [dot(v, x), dot(v, y), dot(v, z)]
    rn, vn = math.sqrt(dot(rb, rb)), math.sqrt(dot(vb, vb))
    h = cross(rb, vb); hn = math.sqrt(dot(h, h))
    ev = [(vn ** 2 - MU / rn) * rb[i] / MU - dot(rb, vb) * vb[i] / MU for i in range(3)]
    e = math.sqrt(dot(ev, ev)); a = 1 / (2 / rn - vn ** 2 / MU)
    inc = math.degrees(math.acos(h[2] / hn))
    node = [-h[1], h[0], 0.0]; nn = math.hypot(node[0], node[1])
    raan = math.degrees(math.atan2(node[1], node[0])) % 360
    aop = math.degrees(math.acos(max(-1, min(1, dot(node, ev) / (nn * e)))))
    if ev[2] < 0: aop = 360 - aop
    nu = math.degrees(math.acos(max(-1, min(1, dot(ev, rb) / (e * rn)))))
    if dot(rb, vb) < 0: nu = 360 - nu
    ref = dict(ra_km=sc["ra_m"] / 1e3, rp_alt_km=sc["rp_altitude_m"] / 1e3, i_deg=sc["i_deg"],
               raan_deg=sc["raan_deg"], aop_deg=sc["aop_deg"], ta_deg=sc["ta_deg"])
    got = dict(ra_km=a * (1 + e), rp_alt_km=a * (1 - e) - R_KM, i_deg=inc, raan_deg=raan, aop_deg=aop, ta_deg=nu)
    return dict(label="DERIVED; REPORT ONLY (plan R1.1)", frame="Mars body-equator inertial (pole from pck00011, T=%.9f cy)" % T,
                reference=ref, from_nav_state=got, difference={k: got[k] - ref[k] for k in ref},
                note="reference elements are 'reference only' in the manifest; the nav state is not used as a start state")


# ---------------------------------------------------------------- build
SRC_NAMES = ["q_per_pass_flown.csv", "burn_table_flown.csv"]
FEATHERS = ["True_Odyssey_periapsis_alts_kernel.feather", "True_Odyssey_apoapsis_alts_kernel.feather"]


def build(dumps):
    D = data_dir(); D.mkdir(parents=True, exist_ok=True)
    src = Path(os.environ["JAIS_ODY_DIR"])
    prov = {}
    for n in SRC_NAMES:
        shutil.copyfile(src / n, D / n)
        prov[n] = dict(sha256=sha(D / n), label="SOURCED (JAIS repo data/raw/odyssey_closed_loop_20261008, copied unchanged)")
    p = {int(float(r["orbit"])): float(r["altitude"]) for r in rows(Path(dumps) / "True_Odyssey_periapsis_alts_kernel.csv")}
    a = {int(float(r["orbit"])): float(r["altitude"]) for r in rows(Path(dumps) / "True_Odyssey_apoapsis_alts_kernel.csv")}
    assert sorted(p) == sorted(a)
    for n in FEATHERS:
        prov[n] = dict(sha256=sha(REPO / "data/telemetry/Odyssey" / n), label="SOURCED (read via Arrow.jl; values copied to kernel_apsides.csv)")
    E = derive_epochs({k: dict(periapsis_alt_km=repr(p[k]), apoapsis_alt_km=repr(a[k])) for k in p})
    with open(D / "kernel_apsides.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["orbit", "periapsis_alt_km", "apoapsis_alt_km", "t_peri_s_derived", "t_apo_s_derived"])
        for k in sorted(p):
            w.writerow([k, repr(p[k]), repr(a[k]), f"{E[k]['t_peri']:.3f}", f"{E[k]['t_apo']:.3f}"])
    m = metrics(D)
    nav = nav_elements(REPO / "benchmarks/studies/telemetry_validation/manifests/odyssey_tolson.toml")
    ref = dict(plan="odyssey_comparison_plan.md sha256 f38bfde788ab702360489f6c0306e336e4b2d8d8679181c4a5ab8ae34ae9a68e",
               labels="SOURCED/DERIVED/ASSUMED per plan section 2",
               inputs=prov,
               kernel_apsides_csv_sha256=sha(D / "kernel_apsides.csv"),
               epoch_note="DERIVED: the feather files carry no time column; epochs are Keplerian-period sums of the kernel apsides, anchored to the 14 flown burn epochs (burn_epoch_anchor_check)",
               scoring="Q_s = 0.5 rho_peak V_p^3 / 1e4, mu 42828.375, R 3396.19 (plan R1.2)",
               metrics=m, nav_state_cross_check=nav)
    (D / "flown_reference.json").write_text(json.dumps(ref, indent=1) + "\n")
    print("wrote", D)


# ---------------------------------------------------------------- selftest
def selftest():
    D = data_dir()
    ref = json.loads((D / "flown_reference.json").read_text())
    for n, v in ref["inputs"].items():                       # hashes recorded and still match
        f = D / n if n in SRC_NAMES else REPO / "data/telemetry/Odyssey" / n
        assert sha(f) == v["sha256"], n
    assert sha(D / "kernel_apsides.csv") == ref["kernel_apsides_csv_sha256"]
    m = metrics(D)
    assert json.loads(json.dumps(m)) == ref["metrics"], "flown_reference.json is stale; rerun build"
    ok = lambda c, s: print(("PASS " if c else "FAIL ") + s) or c
    chk = [
        m["passes"] == 153, round(m["mean"], 4) == 0.2238, round(m["median"], 4) == 0.2104, round(m["max"], 3) == 0.514,
        m["n_gt_025"] == 50, m["n_gt_032"] == 18, m["n_gt_045"] == 2, m["n_gt_054"] == 0, m["n_lt_014"] == 19,
        m["main_burns"] == 13, round(m["main_dv_mps"], 3) == 2.924,
        round(m["lnq_std_mp"], 3) == 0.424 and m["lnq_pairs_mp"] == 139,
        round(m["lnq_white_equiv_mp"], 3) == 0.300, round(m["lnq_std_eg"], 3) == 0.268,
        round(m["iso_risk_eps"], 3) == 0.118, round(m["p19_q"], 3) == 0.185, round(m["p20_q"], 3) == 0.193,
        round(m["T1_kernel"]["table3_days"], 1) == 49.0, m["kernel_mp_end_orbit"] == 171,
        m["kernel_first_apo_le_7408_orbit"] == 171, m["proxy_recompute_max_rel_dev"] < 1e-3,
    ]
    # R1.3 table (nominal fraction, count, +-10% interval), three decimals
    tab = {">0.25": (50, .327, .235, .399), ">0.30": (24, .157, .078, .255), ">0.32": (18, .118, .052, .183),
           ">0.45": (2, .013, .007, .033), ">0.54": (0, 0, 0, .007), ">0.65": (0, 0, 0, 0)}
    for t, (n, f, lo, hi) in tab.items():
        r = m["gt_025_to_065_table"][str(float(t[1:]))]
        chk.append((r["count"], round(r["fraction"], 3), round(r["interval"][0], 3), round(r["interval"][1], 3)) == (n, f, lo, hi))
        print("  R1.3 %s got %s" % (t, (r["count"], round(r["fraction"], 3), [round(x, 3) for x in r["interval"]])))
    r = m["lt_table"]["0.14"]
    chk.append((r["count"], round(r["fraction"], 3), [round(x, 3) for x in r["interval_plan_table"]]) == (19, .124, [.078, .157]))
    print("  NOTE plan R1.3 '< 0.14' row [0.078, 0.157] reproduces only with the upper threshold 0.14*1.1; the plan's own rule (T/0.9) gives %s"
          % [round(x, 3) for x in r["interval"]])
    r = m["lt_table"]["0.05"]
    chk.append((r["count"], r["interval"]) == (0, [0, 0]))
    bad = [i for i, c in enumerate(chk) if not c]
    ok(not bad, "plan DERIVED numbers reproduced (%d checks%s)" % (len(chk), "" if not bad else ", failed idx %s" % bad))
    print("  stats: mean %.4f median %.4f max %.3f lnq std %.3f (%d pairs) eg %.3f; dV %.3f over %d burns %s"
          % (m["mean"], m["median"], m["max"], m["lnq_std_mp"], m["lnq_pairs_mp"], m["lnq_std_eg"], m["main_dv_mps"], m["main_burns"], m["main_burn_orbits"]))
    print("  proxy recompute from kernel apsides: max rel dev %.2e" % m["proxy_recompute_max_rel_dev"])
    print("Table III cross-check (REPORTED, not asserted):")
    mism = 0
    for c in m["table3_crosscheck"]:
        bad_ = not (c["passes_match_exactly"] and c["radius_within_1pct"]); mism += bad_
        print("  %-3s end r Table III %7.0f km vs kernel apo@P%d %8.1f km (%+.2f%%); passes Table III %3d, first-apo<=radius rule %s (P%s), S&B range %3d; days Table III %.1f vs derived %.2f  %s"
              % (c["phase"], c["table3_end_radius_km"], c["sb_last_orbit"], c["kernel_apo_radius_at_last_orbit_km"], c["radius_dev_pct"],
                 c["table3_cum_passes"], c["passes_by_first_le_rule"] or "none", c["first_apo_le_table3_radius_orbit"] or "-", c["passes_sb_range"],
                 c["table3_cum_days"], c["kernel_days_peri19_to_apo_last"], "MISMATCH" if bad_ else "ok"))
    print("  boundaries with a mismatch (>1%% radius or pass count): %d" % mism)
    t1 = m["T1_kernel"]
    print("T1 kernel: %d passes, %.2f d (peri P19 to apo P171) / %.2f d (from scenario epoch) vs Table III 153, 49.0 d" % (t1["passes"], t1["days_derived_peri19_to_apo171"], t1["days_derived_from_scenario_epoch"]))
    print("  anchor-corrected days (scenario epoch to apo P171): %.2f d; drift %.1f s/orbit" % (t1["days_anchor_corrected_from_scenario_epoch"], t1["anchor_drift_s_per_orbit"]))
    print("  burn-epoch anchors: max |derived - flown| = %.0f s over main-phase burns (P<=171), %.0f s over all 25 (later ones jump ~1 day: derived epochs are not trusted past P171)"
          % (m["burn_epoch_anchor_max_abs_diff_s_main"], m["burn_epoch_anchor_max_abs_diff_s_all"]))
    for k, b in m["T2_T3_bins"].items():
        print("  bin %-10s n=%3d T2 median %7.1f km/pass  T3 median %8.1f km/day" % (k, b["n"], b["T2_median_km_per_pass"], b["T3_median_km_per_day"]))
    print("  campaign mean T2 %.1f km/pass, T3 %.1f km/day" % (m["T2_T3_campaign_mean"]["T2_km_per_pass"], m["T2_T3_campaign_mean"]["T3_km_per_day"]))
    nv = ref["nav_state_cross_check"]
    print("nav-state elements (report only):")
    for k in nv["reference"]:
        print("  %-10s ref %12.4f  nav %12.4f  diff %+10.4f" % (k, nv["reference"][k], nv["from_nav_state"][k], nv["difference"][k]))
    print("kernel P19 apo radius %.1f km, MP end apo radius (P171) %.1f km" % (m["kernel_p19_apo_radius_km"], m["kernel_mp_end_radius_km"]))
    assert not bad, "selftest FAILED"
    print("SELFTEST OK")


if __name__ == "__main__":
    c = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    if c == "selftest": selftest()
    elif c == "build": build(sys.argv[sys.argv.index("--dumps") + 1])
    else: sys.exit(__doc__)
