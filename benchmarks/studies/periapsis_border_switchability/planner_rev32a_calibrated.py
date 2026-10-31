"""Revision 32a (analysis only, Python, read-only on flights, no Julia, no commit): a nominal calibrated on DEV flights (rev32dev: OLRT3-W3E5, seeds 1-4, G0975M, G0063M, G0890) that predicts the arrival apse direction.
Method (minimal, fitted on dev only). For each dev flight j: its J2 track tr_j (planner_rev21_position.J2Track on the flown (t_k, a_k, e_k, i_k) history from the geometry's initial elements) and its residual b_j at arrival
(actual apoapsis direction minus the track's, tangent vector). The calibrated nominal is
    u_cal(t) = normalize( mean_j tr_j.apo(t) ), moved by the mean residual b (variant 'resid'), 
t* -> prediction u_cal(t*); median dev pass count N_cal and median dev arrival t_cal are reported. Variant choice (track mean alone, 'track', or with the mean residual, 'resid') is made by leave-one-out on the dev flights.
Held check: the prediction against the held flights' actual arrival apse direction, at the nominal t* (operational) and at the flight's own arrival time (history only), against the table-based nominal. Along-track = component
on the apse rate direction at t*, cross-track the rest."""
import json, math, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_rev21_position as PS
import planner_rev32_alongtrack as A32

OUT = PS.TD / "rev32_data"
TOKS = A32.TOKS


def flight_state(d, tok):
    D = A32.read(d, tok)
    if D is None:
        return None
    i0, om0, Om0 = PS.geometry_elements(tok); tr = PS.J2Track(D["t"], D["a"], D["e"], D["inc"], om0, Om0, i0); t_arr = float(D["t"][-1] + 0.5 * tr.T[-1])
    dom, dOm = PS.j2_orbit(D["a"][-1], D["e"][-1], D["inc"][-1]); u_fl = PS.apo_hat(D["om"][-1] + 0.5 * dom, D["Om"][-1] + 0.5 * dOm, D["inc"][-1]); u_fp = tr.apo(t_arr)
    return dict(tr=tr, t_arr=t_arr, u_fl=u_fl, u_fp=u_fp, b=PS.logmap(u_fp, u_fl), N=D["n"], dir=d)


def expmap(u, v):
    v = PS.tproj(u, v); n = float(np.linalg.norm(v))
    return u if n < 1e-15 else math.cos(n) * u + math.sin(n) * v / n


def predict(dev, t, variant):
    m = np.mean([f["tr"].apo(t) for f in dev], axis=0); m = m / np.linalg.norm(m)
    if variant == "resid":
        m = expmap(m, np.mean([PS.tproj(m, f["b"]) for f in dev], axis=0))
    return m


def dev_fit(tok):
    ds = sorted((PS.ROOT / "rev32dev").glob(f"{tok}_PL_OLRT3-W3E5_pert_ODY_s*_L0")); dev = [s for s in (flight_state(d, tok) for d in ds if PS.finished(d)) if s]
    return dev


def components(u_pred, u_act, u_tgt, s_hat):
    e = PS.logmap(u_tgt, u_act) - PS.logmap(u_tgt, u_pred)       # actual minus predicted, in the tangent plane at the nominal target direction
    al = PS.deg(float(np.dot(e, s_hat))); cr = PS.deg(float(np.linalg.norm(e - np.dot(e, s_hat) * s_hat))); return PS.deg(PS.ang(u_pred, u_act)), al, cr


def evaluate(flts, nom, fn_pred, tok):
    tstar = nom["days"] * 86400.0; u_tgt = nom["track"].apo(tstar); s = PS.tproj(u_tgt, nom["track"].rate_vec(tstar)); s = s / np.linalg.norm(s); rows = []
    for f in flts:
        r_t = components(fn_pred(tstar), f["u_fl"], u_tgt, s); r_a = components(fn_pred(f["t_arr"]), f["u_fl"], u_tgt, s)
        rows.append(dict(run=f["dir"].name, ang_tstar=r_t[0], along_tstar=r_t[1], cross_tstar=r_t[2], ang_arr=r_a[0], along_arr=r_a[1], cross_arr=r_a[2], dt_d=(f["t_arr"] - tstar) / 86400.0))
    return rows


def loo(tok, dev, nom):
    out = {}
    for variant in ("track", "resid"):
        rows = []
        for i in range(len(dev)):
            rest = dev[:i] + dev[i + 1:]; rows += evaluate([dev[i]], nom, lambda t, rest=rest: predict(rest, t, variant), tok)
        out[variant] = dict(median_arr=float(np.median([r["ang_arr"] for r in rows])), median_tstar=float(np.median([r["ang_tstar"] for r in rows])), within1_arr=float(np.mean([r["ang_arr"] < 1 for r in rows])), within1_tstar=float(np.mean([r["ang_tstar"] < 1 for r in rows])))
    return out


def fit_all():
    fit = {}
    for tok in TOKS:
        dev = dev_fit(tok); nom = PS.nominal(tok); l = loo(tok, dev, nom)
        variant = min(l, key=lambda v: l[v]["median_tstar"])         # chosen by dev leave-one-out only
        fit[tok] = dict(variant=variant, loo=l, n_dev=len(dev), N_cal=float(np.median([f["N"] for f in dev])), t_cal_d=float(np.median([f["t_arr"] for f in dev]) / 86400.0), t_nominal_d=nom["days"], dev_runs=[f["dir"].name for f in dev])
    return fit


if __name__ == "__main__" and sys.argv[1] == "fit":
    OUT.mkdir(parents=True, exist_ok=True); fit = fit_all(); json.dump(fit, open(str(OUT / "rev32a_fit.json"), "w"), indent=1)
    for k, v in fit.items():
        print(k, "variant", v["variant"], "N_cal %.0f t_cal %.1f d (table nominal %.1f)" % (v["N_cal"], v["t_cal_d"], v["t_nominal_d"]), "LOO", {a: {b: round(c, 3) for b, c in x.items()} for a, x in v["loo"].items()})


def held_flights(tok):
    out = []
    for b in ("rev31pos", "rev23geo", "rev28t"):
        for d in PS.block_dirs(PS.ROOT, b):
            k = PS.parse_id(d.name)
            if k["tok"] == tok and PS.finished(d) and A32.GROUP(k["arm"]) != "v3 (default-t*)":
                s = flight_state(d, tok)
                if s: s.update(group=A32.GROUP(k["arm"]), block=b); out.append(s)
    return out


if __name__ == "__main__" and sys.argv[1] == "held":
    fit = json.load(open(str(OUT / "rev32a_fit.json"))); res = {}
    for tok in TOKS:
        dev = dev_fit(tok); nom = PS.nominal(tok); held = held_flights(tok); var = fit[tok]["variant"]
        cal = evaluate(held, nom, lambda t: predict(dev, t, var), tok); tab = evaluate(held, nom, lambda t: nom["track"].apo(t), tok)
        for r, h in zip(cal + tab, held + held): r["group"] = h["group"]
        res[tok] = dict(variant=var, n=len(held), cal=cal, table=tab)
        for lab, rows in (("calibrated", cal), ("table nominal", tab)):
            a = np.array([r["ang_tstar"] for r in rows]); b = np.array([r["ang_arr"] for r in rows])
            print(tok, lab, "n", len(rows), "at t*: median %.2f within1 %.2f | along med |.| %.2f cross med %.2f | at own arrival: median %.2f within1 %.2f" % (np.median(a), (a < 1).mean(), np.median([abs(r["along_tstar"]) for r in rows]), np.median([r["cross_tstar"] for r in rows]), np.median(b), (b < 1).mean()))
    json.dump(res, open(str(OUT / "rev32a_held.json"), "w"), indent=1, default=float)
