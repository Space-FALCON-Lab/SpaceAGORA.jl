"""Stage N UNFROZEN supplement (written AFTER the frozen results; changes no verdict). Oracle persistent offset plus a CAUSAL drift tracker:
the Stage N Kalman filter (dev-fit r, q) started at the cell's true mean ysys with prior variance r*1e-6 (offset known), so only drift (and white noise) remain.
python3 -B stagen_unfrozen.py --stagei ROOT --f20 ROOT [--out f.json]"""
import argparse, json, pathlib, sys
import numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import stagej_analysis as sj, stagen_analysis as sn

def run(root, sub, pl, rng):
    P = sj.load(root, sub, pl, "fixed_pert"); P = P[P.inw].reset_index(drop=True); D = P[P.split == "dev"].reset_index(drop=True); H = P[P.split == "heldout"].reset_index(drop=True)
    kp = sn.fit_kalman(D); E, T = [], []
    for t, G in H.groupby("tag", sort=False):
        b = G.ysys.mean(); _, _, e, tg = sn.kalman(G, kp["r"], kp["q"], b, kp["r"] * 1e-6); E += list(e); T += list(tg)   # q from the dev-fit rho: drift allowed
    E, T = np.array(E), np.array(T); hl = float(H.h_local_m.median()); sp = sn.spread(E)
    return dict(spread_ln=sp, ci=sn.boot(T, dict(e=E), lambda d: sn.spread(d["e"]), rng), altitude_shift_m=sp * hl, corridor_ln_left_paper=sn.LN5 - sp, q_hi={str(p): float(np.quantile(E, 1 - p)) for p in sn.EPS})

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--f20"); ap.add_argument("--stagei"); ap.add_argument("--out"); a = ap.parse_args(); rng = np.random.default_rng(sn.SEED)
    R = {n: {pl: run(r, s, pl, rng) for pl in sj.PLANETS} for n, r, s in (("I", a.stagei, "stagei"), ("F20", a.f20, "stagef20"))}
    if a.out: json.dump(R, open(a.out, "w"), indent=1, default=float)
    print(json.dumps(R, indent=1))
