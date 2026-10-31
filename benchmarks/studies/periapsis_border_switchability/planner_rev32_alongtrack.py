"""Revision 32 (analysis only, Python, read-only on flights, no Julia, no commit): where does the multi-degree along-track position error at G0975M, G0063M and G0890 come from?
Per flight the Revision 21 telescoping split is kept: tot = A (timing) + C (history: flown (a, e, i) against the nominal plan's) + B (residual, everything J2 does not explain), vectors in the tangent plane at the target
(planner_rev21_position). The extensions: (i) drag-induced apse rotation: regression of the per-pass apse change not explained by secular J2 on the apoapsis drop of the pass (pooled per geometry, bootstrap over flights),
total = slope x campaign drop; (ii) tesseral kicks along the FLOWN tracks (tesseral_lib.gauss_kick, 12 planet phases, coherent accumulation, rms and worst-case sum); (iii) the history term C; (iv) timing A and the
closing-apoapsis definition (T from the osculating a against the post-pass orbit); (v) the target direction if the nominal plan used the drop calibration kappa fitted on these flights (V = log map of the kappa-calibrated target
at t* against the plan's). Shares of the squared error = projection of each term on the total error vector / |total|^2 (they add up to 1 for A + C + B); intervals by bootstrap over flights.
"""
import json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_rev21_position as PS
import planner_lib as PLB
import planner_rev18_offline as O

OUT = PS.TD / "rev32_data"
TOKS = ("G0975M", "G0063M", "G0890")
BLOCKS = ("rev31pos", "rev23geo", "rev28t", "rev21v")
GROUP = lambda arm: "v4+L1" if arm.startswith("OLRT3-W4") else ("v3 (default-t*)" if arm.endswith("W3T") else "v3")


def nominal_k(tok, kappa=None):
    """PS.nominal with the plan's levels kept and the drop calibration kappa replaced (None = the plan's own)."""
    RP = PS._rp(); ic = PS.headline_ic(tok); x0 = math.log(RP.QFIRST_CAL.get(tok, O.Q_FIRST)) + RP.ic_shift(ic); tab = PS.get_tab(tok); sw, P = RP.arm_def("OLR")
    inp = PLB.ol_inputs(sw, tab, P, x0); pl = PLB.dp_plan(inp["level"], inp["lo"], inp["hi"], inp["w"], P["dp_tie_eps"], pen=P["dp_burn_penalty_mps"]); x = np.concatenate(([x0], np.asarray(inp["level"]) + pl["nu"]))
    h = PS.arrival_history(tab, x, x0, PS.kappa_for(tok) if kappa is None else kappa, PS.target_radius(tok)); i0, om0, Om0 = PS.geometry_elements(tok)
    return dict(track=PS.J2Track(h["t"], h["a"], h["e"], np.full(len(h["t"]), i0), om0, Om0, i0), days=h["days"], hist=h, x=x, N=h["N"])


def flights():
    out = []
    for b in BLOCKS:
        for d in PS.block_dirs(PS.ROOT, b):
            k = PS.parse_id(d.name)
            if k["tok"] in TOKS and PS.finished(d) and GROUP(k["arm"]) != "v3 (default-t*)":
                out.append(dict(dir=d, tok=k["tok"], arm=k["arm"], block=b, seed=k["seed"], group=GROUP(k["arm"])))
    return out


def read(d, tok):
    P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); ra = P.ra_m.values.astype(float); hit = np.where(np.isfinite(ra) & (ra <= PS.target_radius(tok)))[0]
    if not len(hit):
        return None
    n = int(hit[0] + 1); Pn = P.iloc[:n]
    return dict(n=n, t=Pn.time_s.values.astype(float), a=Pn.a_m.values.astype(float), e=Pn.e.values.astype(float), inc=Pn.i_rad.values.astype(float), om=Pn.omega_rad.values.astype(float), Om=Pn.raan_rad.values.astype(float),
                ra=Pn.ra_m.values.astype(float), rp=Pn.rp_m.values.astype(float), q=Pn.q_w_cm2.values.astype(float))


def vectors(fl, nom, nomk):
    tok = fl["tok"]; D = read(fl["dir"], tok)
    if D is None:
        return None
    i0, om0, Om0 = PS.geometry_elements(tok); tr = PS.J2Track(D["t"], D["a"], D["e"], D["inc"], om0, Om0, i0); t_arr = float(D["t"][-1] + 0.5 * tr.T[-1]); nt = nom["track"]; tstar = nom["days"] * 86400.0
    dom, dOm = PS.j2_orbit(D["a"][-1], D["e"][-1], D["inc"][-1]); u_fl = PS.apo_hat(D["om"][-1] + 0.5 * dom, D["Om"][-1] + 0.5 * dOm, D["inc"][-1])
    u_tgt = nt.apo(tstar); u_nA = nt.apo(t_arr); u_fp = tr.apo(t_arr)
    tot = PS.logmap(u_tgt, u_fl); A = PS.logmap(u_tgt, u_nA); C = PS.tproj(u_tgt, PS.logmap(u_nA, u_fp)); B = PS.tproj(u_tgt, PS.logmap(u_fp, u_fl)); V = PS.logmap(u_tgt, nomk["track"].apo(tstar))
    s = PS.tproj(u_tgt, nt.rate_vec(tstar)); s = s / np.linalg.norm(s)
    # (iv): definition offset: the closing time from the osculating a at the periapsis event against the post-pass orbit (ra + rp)/2
    T_osc = 2 * math.pi * math.sqrt(D["a"][-1] ** 3 / PS.MU); T_post = 2 * math.pi * math.sqrt((0.5 * (D["ra"][-1] + D["rp"][-1])) ** 3 / PS.MU); rate_deg_s = PS.deg(float(np.linalg.norm(nt.rate_vec(tstar))))
    return dict(D=D, tr=tr, t_arr=t_arr, tstar=tstar, u_fl=u_fl, u_tgt=u_tgt, u_fp=u_fp, tot=tot, A=A, B=B, C=C, V=V, s_hat=s, def_deg=rate_deg_s * 0.5 * abs(T_osc - T_post), dt_d=(t_arr - tstar) / 86400.0, rate_deg_d=rate_deg_s * 86400.0)


def drag_fit(vs, lag=0, boot=None, rng=None):
    """Pooled OLS of the per-pass apse residual (deg) on the apoapsis drop (km): r_p = dOmega_obs - dOmega_J2(a_p, e_p, i_p), x_p = drop of pass p+1-lag. Returns slope (deg/km), intercept (deg/pass)."""
    X, Y = [], []
    for v in vs:
        D = v["D"]; n = D["n"]
        for p in range(1, n - 1):
            dobs = (D["om"][p] - D["om"][p - 1] + math.pi) % (2 * math.pi) - math.pi; dj2 = PS.j2_orbit(D["a"][p - 1], D["e"][p - 1], D["inc"][p - 1])[0]
            drop = (D["ra"][p - 1 + (0 if lag == 0 else -1)] - D["ra"][p + (0 if lag == 0 else -1)]) / 1e3 if p + (0 if lag == 0 else -1) < n else np.nan
            if np.isfinite(drop) and abs(dobs) < 1.0:
                X.append(drop); Y.append(PS.deg(dobs - dj2))
    X = np.array(X); Y = np.array(Y); b = np.polyfit(X, Y, 1)
    return float(b[0]), float(b[1]), len(X)


def tesseral_flown(v, n_phase=12, nf=64):
    import tesseral_lib as tl
    cs = tl.coefficients("Mars", 20); Cc, Ss = tl.select(cs["C"], cs["S"], lambda n, m: m >= 1); D = v["D"]; n = D["n"]; disp = np.zeros((n_phase, n, 3)); per = np.zeros((n_phase, n))
    for ip in range(n_phase):
        th0 = 2 * math.pi * ip / n_phase
        for k in range(n):
            kk = tl.gauss_kick(float(D["a"][k]), float(D["e"][k]), float(D["inc"][k]), float(D["Om"][k]), float(D["om"][k]), th0 + PS.W_ROT * float(D["t"][k]), PS.W_ROT, cs, Cc, Ss, nf=nf)
            u0 = PS.apo_hat(D["om"][k], D["Om"][k], D["inc"][k]); u1 = PS.apo_hat(D["om"][k] + kk["omega"], D["Om"][k] + kk["raan"], D["inc"][k] + kk["i"]); disp[ip, k] = PS.logmap(u0, u1); per[ip, k] = PS.ang(u0, u1)
    return dict(coherent_deg=[PS.deg(float(np.linalg.norm(disp[ip].sum(0)))) for ip in range(n_phase)], rms_deg=[PS.deg(float(math.sqrt((per[ip] ** 2).sum()))) for ip in range(n_phase)], sum_abs_deg=[PS.deg(float(per[ip].sum())) for ip in range(n_phase)])


def deg(x): return PS.deg(x)


def analyse(n_boot=2000, seed=20261031):
    from scipy.stats import spearmanr
    rng = np.random.default_rng(seed); fls = flights(); res = {"cells": {}, "flights": []}; noms = {}; vs_all = []
    kap = {}
    for tok in TOKS:
        ds = [f["dir"] for f in fls if f["tok"] == tok]; kap[tok] = PS.calibrate_kappa(tok, ds); noms[tok] = (PS.nominal(tok), nominal_k(tok, kap[tok]["kappa_med"]), nominal_k(tok))
    for f in fls:
        nom, nomk, nom2 = noms[f["tok"]]; v = vectors(f, nom, nomk)
        if v is not None:
            v.update(f); vs_all.append(v)
    for tok in TOKS:
        vs = [v for v in vs_all if v["tok"] == tok]; nom = noms[tok][0]
        # (i) drag regression (pooled per geometry), lag 0 and 1, bootstrap over flights
        fits = {}
        for lag in (0, 1):
            b0 = drag_fit(vs, lag); bs = []
            for _ in range(300):
                idx = rng.integers(0, len(vs), len(vs)); bs.append(drag_fit([vs[i] for i in idx], lag)[:2])
            bs = np.array(bs); fits[lag] = dict(slope_deg_per_km=b0[0], intercept_deg_per_pass=b0[1], n=b0[2], slope_ci=[float(np.percentile(bs[:, 0], 2.5)), float(np.percentile(bs[:, 0], 97.5))], intercept_ci=[float(np.percentile(bs[:, 1], 2.5)), float(np.percentile(bs[:, 1], 97.5))])
        lag = 0; beta = fits[lag]["slope_deg_per_km"]
        for v in vs:
            D = v["D"]; dw = math.radians(beta * (D["ra"][0] - D["ra"][-1]) / 1e3)               # drag-induced apse rotation over the campaign (slope x total drop)
            u2 = PS.apo_hat(D["om"][-1] + dw, D["Om"][-1], D["inc"][-1]); v["Dvec"] = PS.tproj(v["u_tgt"], PS.logmap(v["u_fl"], v["u_fl"] if True else u2) + (PS.logmap(v["u_fl"], u2))) if False else PS.tproj(v["u_tgt"], PS.logmap(v["u_fl"], u2))
            v["drag_deg"] = deg(float(np.linalg.norm(v["Dvec"])))
        # (ii) tesseral along the flown tracks (one flight per arm group and block is enough for the bound; here all)
        for v in vs:
            v["tess"] = tesseral_flown(v)
        def shares(sel):
            tt = sum(float(np.dot(v["tot"], v["tot"])) for v in sel); d = lambda key: sum(float(np.dot(v[key], v["tot"])) for v in sel) / tt
            return dict(A=d("A"), C=d("C"), B=d("B"), V=d("V"), drag=d("Dvec"), tess_upper=sum(max(v["tess"]["coherent_deg"]) ** 2 for v in sel) / (tt * (180 / math.pi) ** 2) if False else sum(math.radians(max(v["tess"]["coherent_deg"])) ** 2 for v in sel) / tt, tess_rms_upper=sum(math.radians(max(v["tess"]["rms_deg"])) ** 2 for v in sel) / tt)
        groups = {"all": vs, "v3": [v for v in vs if v["group"] == "v3"], "v4+L1": [v for v in vs if v["group"] == "v4+L1"]}
        cell = dict(kappa=kap[tok], drag_fit=fits, nominal_days=nom["days"], groups={})
        for gname, sel in groups.items():
            if not sel:
                continue
            s0 = shares(sel); bs = [shares([sel[i] for i in rng.integers(0, len(sel), len(sel))]) for _ in range(n_boot)]
            ci = {k: [float(np.percentile([b[k] for b in bs], 2.5)), float(np.percentile([b[k] for b in bs], 97.5))] for k in s0}
            med = lambda f: float(np.median([f(v) for v in sel]))
            cell["groups"][gname] = dict(n=len(sel), shares=s0, ci=ci, med_err=med(lambda v: deg(float(np.linalg.norm(v["tot"])))), med_A=med(lambda v: deg(float(np.linalg.norm(v["A"])))), med_C=med(lambda v: deg(float(np.linalg.norm(v["C"])))), med_B=med(lambda v: deg(float(np.linalg.norm(v["B"])))),
                                         med_V=med(lambda v: deg(float(np.linalg.norm(v["V"])))), med_drag=med(lambda v: v["drag_deg"]), med_def=med(lambda v: v["def_deg"]), med_dt=med(lambda v: v["dt_d"]),
                                         tess_coh_med=med(lambda v: float(np.median(v["tess"]["coherent_deg"]))), tess_coh_max=float(max(max(v["tess"]["coherent_deg"]) for v in sel)), tess_rms_max=float(max(max(v["tess"]["rms_deg"]) for v in sel)), tess_sum_max=float(max(max(v["tess"]["sum_abs_deg"]) for v in sel)))
        res["cells"][tok] = cell
    # P32-2: along-track error against the campaign-pace mismatch, pooled
    rows = []
    for v in vs_all:
        nom = noms[v["tok"]][0]; hist = nom["hist"]; x = None
        D = v["D"]; N = D["n"]; Nn = nom["N"]
        # integrated level difference: ln q flown against the nominal plan level at the same fraction of the campaign duration
        tf = (D["t"] - D["t"][0]) / (D["t"][-1] - D["t"][0]); tn = (hist["t"] - hist["t"][0]) / (hist["t"][-1] - hist["t"][0])
        xn = noms[v["tok"]][2]["x"]; xn = np.asarray(xn)[:len(tn)]; Ld = float(np.mean(np.log(D["q"]) - np.interp(tf, tn, xn)))
        al = deg(float(np.dot(v["tot"], v["s_hat"])))
        rows.append(dict(tok=v["tok"], run=v["dir"].name, along_deg=al, abs_along=abs(al), dt_d=v["dt_d"], dN=N - Nn, level_diff=Ld, arr_minus_nominal_d=v["t_arr"] / 86400.0 - nom["days"]))
    res["flights"] = rows
    sp = {}
    for key in ("dt_d", "dN", "level_diff"):
        for yk in ("along_deg", "abs_along"):
            r, p = spearmanr([r_[key] for r_ in rows], [r_[yk] for r_ in rows]); sp[f"{key}|{yk}"] = (float(r), float(p))
    res["spearman"] = sp; res["n_flights"] = len(rows)
    return res


if __name__ == "__main__" and len(sys.argv) == 1:
    OUT.mkdir(parents=True, exist_ok=True)
    R_ = analyse(); json.dump(R_, open(str(OUT / "alongtrack.json"), "w"), indent=1, default=float)
    for tok, c in R_["cells"].items():
        print(tok, "kappa", round(c["kappa"]["kappa_med"], 3), "iqr_rel %.3f" % c["kappa"]["iqr_rel"], "drag fit lag0", {k: (round(v, 5) if not isinstance(v, list) else [round(x, 5) for x in v]) for k, v in c["drag_fit"][0].items()})
        for g, v in c["groups"].items():
            print("  ", g, "n", v["n"], "err %.2f A %.2f C %.2f B %.2f V %.2f drag %.3f def %.4f dt %.2f | shares" % (v["med_err"], v["med_A"], v["med_C"], v["med_B"], v["med_V"], v["med_drag"], v["med_def"], v["med_dt"]), {k: (round(x, 3), [round(y, 3) for y in v["ci"][k]]) for k, x in v["shares"].items()}, "tess coh max %.2f rms max %.2f sum max %.1f" % (v["tess_coh_max"], v["tess_rms_max"], v["tess_sum_max"]))
    print({k: (round(a, 2), round(b, 3)) for k, (a, b) in R_["spearman"].items()}, R_["n_flights"])


def pass_count_model(res, n_boot=2000, seed=7):
    """Along-track error against the pass-count difference dN = N_flown - N_nominal: OLS slope (deg per pass) per geometry with a bootstrap interval, against the J2 apse rotation per orbit of the final orbits (the
    mechanism: J2 rotates the apse per ORBIT, so a campaign that reaches the same arrival time with a different number of orbits accumulates a different rotation)."""
    rng = np.random.default_rng(seed); out = {}
    for tok in TOKS:
        rows = [r for r in res["flights"] if r["tok"] == tok]; dN = np.array([r["dN"] for r in rows], float); al = np.array([r["along_deg"] for r in rows]); fl = [f for f in flights() if f["tok"] == tok]
        D = read(fl[0]["dir"], tok); i0 = PS.geometry_elements(tok)[0]; last = [PS.deg(PS.j2_orbit(D["a"][k], D["e"][k], D["inc"][k])[0]) for k in range(D["n"] - 50, D["n"])]
        b = np.polyfit(dN, al, 1); bs = [np.polyfit(dN[i], al[i], 1)[0] for i in (rng.integers(0, len(dN), len(dN)) for _ in range(n_boot)) if np.ptp(dN[i]) > 0]
        out[tok] = dict(n=len(rows), dN_med=float(np.median(dN)), dN_range=[float(dN.min()), float(dN.max())], slope_deg_per_pass=float(b[0]), slope_ci=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))], intercept=float(b[1]),
                        dom_per_orbit_final_deg=float(np.mean(last)), passes_per_deg=float(1.0 / abs(np.mean(last))), along_med=float(np.median(al)), r2=float(np.corrcoef(dN, al)[0, 1] ** 2))
    return out


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "passmodel":
    res = json.load(open(str(OUT / "alongtrack.json"))); pm = pass_count_model(res); json.dump(pm, open(str(OUT / "passcount_model.json"), "w"), indent=1)
    for k, v in pm.items():
        print(k, {a: (round(b, 3) if isinstance(b, float) else ([round(x, 3) for x in b] if isinstance(b, list) else b)) for a, b in v.items()})


def calibrated_nominal_check():
    """IN-SAMPLE check (kappa fitted on the same flights): angle between the flown apoapsis direction and the nominal's at the flown arrival time (no timing term), for the plan's own kappa and for kappa_cal; pass counts N_nom, N_cal against N_flown."""
    out = {}
    for tok in TOKS:
        fl = [f for f in flights() if f["tok"] == tok]; kc = PS.calibrate_kappa(tok, [f["dir"] for f in fl])["kappa_med"]; nom = PS.nominal(tok); nk = nominal_k(tok, kc); e0, e1, Nf, ar0, ar1, arf = [], [], [], nom["days"], nk["days"], []
        for f in fl:
            v = vectors(f, nom, nk)
            if v is None: continue
            e0.append(PS.deg(PS.ang(v["u_fl"], nom["track"].apo(v["t_arr"])))); e1.append(PS.deg(PS.ang(v["u_fl"], nk["track"].apo(v["t_arr"])))); Nf.append(v["D"]["n"]); arf.append(v["t_arr"] / 86400.0)
        out[tok] = dict(kappa_cal=kc, N_nominal=nom["N"], N_cal=nk["N"], N_flown_med=float(np.median(Nf)), t_nominal_d=ar0, t_cal_d=ar1, t_flown_med_d=float(np.median(arf)), err_vs_nominal_med=float(np.median(e0)), err_vs_calibrated_med=float(np.median(e1)), within1_cal=float(np.mean(np.array(e1) < 1)), n=len(e0))
    return out


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "calcheck":
    o = calibrated_nominal_check(); json.dump(o, open(str(OUT / "calibrated_nominal_check.json"), "w"), indent=1)
    for k, v in o.items():
        print(k, {a: round(b, 3) if isinstance(b, float) else b for a, b in v.items()})


def repeatability():
    """Flight-to-flight repeatability of the history (the upper bound of what a nominal fitted on dev flights can give): for each flight i the angle between its apoapsis direction at its arrival and the J2 track of every OTHER flight j of the same
    geometry and arm group (the apse direction j would have at i's arrival time, no timing term other than the arrival difference); summary per geometry and per group."""
    out = {}
    for tok in TOKS:
        fl = [f for f in flights() if f["tok"] == tok]; nom = PS.nominal(tok); V = []
        for f in fl:
            v = vectors(f, nom, nom)
            if v is not None: v.update(f); V.append(v)
        for gname in ("v3", "v4+L1", "all"):
            sel = [v for v in V if gname == "all" or v["group"] == gname]; errs = []
            for i, vi in enumerate(sel):
                row = [PS.deg(PS.ang(vi["u_fl"], vj["tr"].apo(vi["t_arr"]))) for j, vj in enumerate(sel) if j != i]
                errs.append(float(np.median(row)))
            if errs: out[f"{tok}|{gname}"] = dict(n=len(sel), median_of_medians=float(np.median(errs)), max=float(max(errs)), within1=float(np.mean(np.array(errs) < 1)))
    return out


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "repeat":
    o = repeatability(); json.dump(o, open(str(OUT / "repeatability.json"), "w"), indent=1)
    for k, v in o.items(): print(k, {a: round(b, 3) for a, b in v.items()})
