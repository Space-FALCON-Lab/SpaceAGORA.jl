"""Revision 28B (analysis only, Python replay, no Julia, no flights): the terminal-aim margin term_aim_below_m (12.5 km ASSUMED in the config) with guard v3 on DEV paths, drop noise of the last passes drawn from a pool of the
dev flights' own drop-noise values (eps = flown drop / model drop), common draws for all margins. Metrics per path and draw: final apoapsis error (km), low-side miss (< -25 km), passes, window, campaign dV.
`local_check`: the flown misses re-simulated locally from the first terminal-aim pass with the flown per-pass drop ratio and prediction error (any target radius)."""
import json, math, pathlib, sys, time
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import planner_rev20_replay as R
import planner_rev21_guard as G

R.DEV_EXT = True
OUT = R.TD / "rev28_data"
MARGINS = (12500.0, 10000.0, 7500.0, 5000.0)


def eps_pool(last=40):
    v = []
    for g in ("ODY", "G1"):
        for rid in R.dev_paths_ext(g):
            F = R.get_flight(rid); n = F["n"]; v += [float(x) for x in F["eps"][max(1, n - last):n]]
    return np.array(v)


MARGIN_FROZEN_M = 5000.0      # Revision 28B choice, frozen after the dev run: the smallest candidate margin with 0 low-side misses on dev (0 of 108 draws) and no window or dV cost; it also avoids all four flown misses in the local check
HELD_MARGINS = (12500.0, 7500.0, 5000.0)


def dev_run(margins=MARGINS, draws=4, geoms=("ODY", "G1"), name="dev_margin", paths_fn=None):
    pool = eps_pool(); OUT.mkdir(parents=True, exist_ok=True)
    print("eps pool n %d, mean %.3f, sd of ln %.3f, max %.2f" % (len(pool), pool.mean(), np.std(np.log(pool)), pool.max()), flush=True)
    rows = []
    for g in geoms:
        ic = R.get_ic(g, "S-fixed"); ts = R.tstar_for(g, ic)
        for rid in (paths_fn or R.dev_paths_ext)(g):
            F = R.get_flight(rid); sd = R.split_id(rid)["seed"] * 100
            for k in range(draws):
                for mg in margins:
                    r = R.replay(F, g, "OLRT3", ic=ic, guard=G.GUARD_V3, tstar_d=ts, seed=sd, term_below_m=mg, eps_pool=pool, eps_seed=k + 1)
                    if not r["reached"]:
                        continue
                    rows.append(dict(geom=g, rid=rid, draw=k, margin=mg, err=r["final_err_km"], n=r["n"], win=bool(abs(r["t_arr_d"] - ts) <= r["w_max_d"]), dv=r["dv_campaign"], hard=bool(r["hard_violation"]), dt=r["t_arr_d"] - ts))
        json.dump(rows, open(str(OUT / f"{name}.json"), "w"), indent=1, default=float)
    return rows


def summarize(rows):
    out = {}
    for mg in sorted(set(r["margin"] for r in rows), reverse=True):
        rr = [r for r in rows if r["margin"] == mg]; base = {(r["rid"], r["draw"]): r for r in rows if r["margin"] == 12500.0}
        e = np.array([r["err"] for r in rr]); extra = [r["n"] - base[(r["rid"], r["draw"])]["n"] for r in rr]; ddv = [r["dv"] - base[(r["rid"], r["draw"])]["dv"] for r in rr]
        winlost = sum(1 for r in rr if base[(r["rid"], r["draw"])]["win"] and not r["win"]); wingain = sum(1 for r in rr if (not base[(r["rid"], r["draw"])]["win"]) and r["win"])
        out[mg] = dict(n=len(rr), low=int((e < -25).sum()), low_frac=float((e < -25).mean()), med=float(np.median(e)), p10=float(np.percentile(e, 10)), p90=float(np.percentile(e, 90)), minerr=float(e.min()), extra_pass_frac=float(np.mean([x > 0 for x in extra])),
                       fewer_pass_frac=float(np.mean([x < 0 for x in extra])), mean_extra_pass=float(np.mean(extra)), ddv_med=float(np.median(ddv)), ddv_mean=float(np.mean(ddv)), win=int(sum(r["win"] for r in rr)), win_lost=winlost, win_gain=wingain, hard=int(sum(r["hard"] for r in rr)), dv_med=float(np.median([r["dv"] for r in rr])))
    return out


if __name__ == "__main__" and sys.argv[1] == "dev":
    rows = dev_run(draws=int(sys.argv[2]) if len(sys.argv) > 2 else 4)
    for mg, v in summarize(rows).items():
        print("margin %.1f km:" % (mg / 1e3), {k: (round(x, 3) if isinstance(x, float) else x) for k, x in v.items()})


# ------------------------------------------------------------------------------------------------ the flown misses, re-simulated locally from the first terminal-aim pass
import pandas as pd
import planner_lib as PLB
import planner_rev21_position as PS


def local_sim(dirpath, tok, margin_m, term_step=None, steps=None):
    """Flown flight dir -> final apoapsis error (km) and pass count when the terminal aim uses `margin_m`. Before the aim is active the flown level is used; once active (planner_terminal_aim on the simulated state) the level is
    ln q = aim + the flown one-pass prediction error e_p; the realized drop per pass = the flown drop-per-q ratio of that pass x the simulated q (so the drop noise of the flight is kept)."""
    d = pathlib.Path(dirpath); P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True); Rp = pd.read_csv(d / "replan.csv")
    tgt = PS.target_radius(tok); ra = P.ra_m.values.astype(float); n = int(np.where(ra <= tgt)[0][0] + 1); q = P.q_w_cm2.values.astype(float)[:n]; x = np.log(q); ra = ra[:n]
    ra0 = float(ra[0] + (P.ra_m.values[0] - ra[0]))                                          # first pass: ra before = initial apoapsis of the table is not needed (pass >= 2 only)
    ratio = np.zeros(n); ratio[1:] = (ra[:-1] - ra[1:]) / q[1:]                              # drop per unit q of passes 2..n (m per W/cm2)
    row = {int(p): i for i, p in enumerate(Rp["pass"].values)}; e = np.zeros(n); term = np.zeros(n, bool)
    for p in range(2, n + 1):
        i = row.get(p)
        if i is not None:
            dv = float(Rp.executed_dv_mps.values[i]); w = float(Rp.w0.values[i]); z = float(Rp.z0_pred_lnq.values[i]); e[p - 1] = x[p - 1] + dv / w - z; term[p - 1] = (Rp.reason.values[i] == "terminal_aim")
    cfg = PLB.load_cfg(); Pp = PLB.flat_P(cfg, "OLRT3"); Pp = dict(Pp, term_aim_below_m=float(margin_m))
    if term_step: Pp["term_step_m"] = term_step
    p0 = int(np.argmax(term)) + 1 if term.any() else n + 1
    pst = max(3, p0 - 3); ras = list(ra[:pst - 1]); c_hist = [ratio[k] for k in range(1, pst - 1) if ratio[k] > 0]          # c entries for passes 2 .. pst-2 (the decision at pst uses passes up to pst-1)
    xs = list(x[:pst - 1]); ra_prev = ras[-1]; p = pst; fin = None
    while p <= n + 12:
        i = min(p, n) - 1
        if p > pst:
            c_hist.append(ratio[min(p - 1, n) - 1] if False else (ras[-2] - ras[-1]) / math.exp(xs[-1]) if len(ras) >= 2 else ratio[1])
        ta = PLB.terminal_aim(ra_prev - tgt, c_hist, Pp)
        xp = (ta["lnq"] + e[i]) if ta is not None else x[i]
        drop = ratio[i] * math.exp(xp)
        ra_new = ra_prev - drop; ras.append(ra_new); xs.append(xp); ra_prev = ra_new
        if ra_new <= tgt:
            fin = (ra_new - tgt) / 1e3; return dict(err_km=fin, n=p, flown_err_km=(ra[-1] - tgt) / 1e3, flown_n=n)
        p += 1
    return dict(err_km=float("nan"), n=p, flown_err_km=(ra[-1] - tgt) / 1e3, flown_n=n)


FLOWN_MISSES = (("rev21v", "ODY_PL_OLRT3-W3_pert_ODY_s18_L0", "ODY"), ("rev26val2", "ODY_PL_OLRT3-W4L_pert_ODY_s13_L0", "ODY"), ("rev26val2", "ODY_PL_OLRT3-W4L_pert_ODY_s13_L2", "ODY"), ("rev27low", "G1038_PL_OLRT3-W3E_pert_ODY_s16_L0", "G1038"))


def check_misses(margins=MARGINS):
    out = {}
    for blk, name, tok in FLOWN_MISSES:
        d = R.ROOT / blk / name
        out[name] = {mg / 1e3: local_sim(d, tok, mg) for mg in margins}
    return out


if __name__ == "__main__" and sys.argv[1] == "misses":
    for name, v in check_misses().items():
        print(name, {k: (round(x["err_km"], 2), x["n"]) for k, x in v.items()}, "flown", round(list(v.values())[0]["flown_err_km"], 2), list(v.values())[0]["flown_n"])


if __name__ == "__main__" and sys.argv[1] == "held":
    geoms = tuple(sys.argv[2].split(",")) if len(sys.argv) > 2 else ("ODY", "G1")
    rows = dev_run(margins=HELD_MARGINS, draws=4, geoms=geoms, name="held_margin_" + "_".join(geoms), paths_fn=R.held_paths)
    for mg, v in summarize(rows).items():
        print("held", "/".join(geoms), "margin", mg / 1e3, {k: (round(x, 3) if isinstance(x, float) else x) for k, x in v.items()})
