#!/usr/bin/env python3
"""Cross-tool comparison for the SpaceAGORA / Basilisk MuJoCo cross-validation (see README.md).

  python3 compare.py <case.toml> <sagora_dir> <bsk_dir> <out_prefix> [--k 10]

Reads sagora_dt*.csv and bsk_dt*.csv (same columns, same sample times), and for each compared quantity writes
  D        max over samples of the norm of (finest SpaceAGORA run - finest Basilisk run)
  e_S, e_B max over samples of the norm of (finest run - second-finest run) of the same tool: the
           self-convergence error of each side
  bound    k * max(e_S, e_B)
  pass     D <= bound
Nothing is floored or loosened: a zero self-convergence error gives a zero bound.
"""
import argparse, csv, glob, math, os, re, sys
import numpy as np

try:
    import tomllib
except ImportError:  # pragma: no cover
    import tomli as tomllib


def load_ladder(d, tool):
    runs = []
    for f in glob.glob(os.path.join(d, f"{tool}_*.csv")):
        m = re.search(r"_(?:dt|tol)([0-9.eE+-]+)\.csv$", f)
        h = float(m.group(1))
        data = np.genfromtxt(f, delimiter=",", names=True)
        runs.append((h, data))
    runs.sort(key=lambda r: -r[0])  # coarse -> fine
    if len(runs) < 3:
        sys.exit(f"{tool}: need at least 3 levels in {d}, found {len(runs)}")
    t0 = runs[0][1]["t"]
    for h, dat in runs:
        if dat["t"].shape != t0.shape or np.max(np.abs(dat["t"] - t0)) > 1e-9:
            sys.exit(f"{tool}: sample times differ at level {h}")
    return runs


def vec(dat, prefix, comps):
    return np.stack([dat[prefix + c] for c in comps], axis=1)


def quantities(dat, cfg):
    """name -> (T, n) array of the compared quantities, for one run."""
    out = {}
    free = [b["name"] for b in cfg["body"]]
    ref = free[0]
    names = [c[:-3] for c in dat.dtype.names if c.endswith("_rx")]
    for b in names:
        if b == ref:
            continue
        out[("rel_pos", b)] = vec(dat, b + "_r", "xyz") - vec(dat, ref + "_r", "xyz")
        out[("rel_vel", b)] = vec(dat, b + "_v", "xyz") - vec(dat, ref + "_v", "xyz")
    for b in free:
        out[("attitude_q", b)] = vec(dat, b + "_q", "wxyz")
        out[("body_rate_w", b)] = vec(dat, b + "_w", "xyz")
    for c in dat.dtype.names:
        if c.startswith("q_"):
            out[("joint_angle", c[2:])] = dat[c][:, None]
        if c.startswith("qd_"):
            out[("joint_rate", c[3:])] = dat[c][:, None]
    return out


def align_sign(q, qref):
    s = np.sign(np.sum(q * qref, axis=1, keepdims=True))
    s[s == 0] = 1.0
    return q * s


def diff(kind, a, b):
    if kind == "attitude_q":
        a = align_sign(a, b)
    d = np.linalg.norm(a - b, axis=1)
    return float(d.max()), int(d.argmax())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case"); ap.add_argument("sagora"); ap.add_argument("bsk"); ap.add_argument("out")
    ap.add_argument("--k", type=float, default=10.0)
    a = ap.parse_args()
    cfg = tomllib.load(open(a.case, "rb"))
    S, Bk = load_ladder(a.sagora, "sagora"), load_ladder(a.bsk, "bsk")
    qs = [quantities(d, cfg) for _, d in S]
    qb = [quantities(d, cfg) for _, d in Bk]
    t = S[0][1]["t"]
    rows = []
    for key in qs[0]:
        kind, ent = key
        D, iD = diff(kind, qs[-1][key], qb[-1][key])
        eS, _ = diff(kind, qs[-1][key], qs[-2][key]); eS_prev, _ = diff(kind, qs[-2][key], qs[-3][key])
        eB, _ = diff(kind, qb[-1][key], qb[-2][key]); eB_prev, _ = diff(kind, qb[-2][key], qb[-3][key])
        order = lambda p, e: math.log2(p / e) if e > 0 and p > 0 else float("nan")
        # Sharper, informational checks (not used for pass/fail):
        #  D_by_level: D of every SpaceAGORA level against the finest Basilisk run; order_D from its last two entries
        #  D_extrap_S: SpaceAGORA first-order Richardson extrapolation of the two finest levels vs Basilisk
        #  D_extrap2_S: second-order (error c1 h + c2 h^2) extrapolation of the three finest levels vs Basilisk
        xs = [q[key] for q in qs]
        if kind == "attitude_q":
            xs = [align_sign(x, xs[-1]) for x in xs]
        Dl = [diff(kind, x, qb[-1][key])[0] for x in xs]
        Dx, _ = diff(kind, 2 * xs[-1] - xs[-2], qb[-1][key])
        Dx2, _ = diff(kind, (8 * xs[-1] - 6 * xs[-2] + xs[-3]) / 3, qb[-1][key])
        orderD = math.log2(Dl[-2] / Dl[-1]) if Dl[-1] > 0 and Dl[-2] > 0 else float("nan")
        bound = a.k * max(eS, eB)
        rows.append(dict(case=cfg["name"], quantity=kind, entity=ent, D=D, D_extrap_S=Dx, D_extrap2_S=Dx2, order_D=orderD, D_by_level=";".join(f"{d:.3e}" for d in Dl), t_of_D=float(t[iD]), e_S=eS, e_B=eB,
                         k=a.k, bound=bound, ratio=(D / bound if bound > 0 else float("inf")), passed=D <= bound,
                         order_S=order(eS_prev, eS), order_B=order(eB_prev, eB),
                         sagora_levels=";".join(f"{h:g}" for h, _ in S), bsk_levels=";".join(f"{h:g}" for h, _ in Bk)))
    cols = list(rows[0])
    with open(a.out + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, cols); w.writeheader(); w.writerows(rows)
    with open(a.out + ".md", "w") as f:
        f.write(f"# {cfg['name']}: SpaceAGORA ProximityScene vs Basilisk MJScene\n\n")
        f.write(f"SpaceAGORA dt ladder [s]: {rows[0]['sagora_levels']}; Basilisk step ladder [s]: {rows[0]['bsk_levels']}; "
                f"k = {a.k:g}. D = max over samples of the cross-tool difference between the finest runs; e_S, e_B = "
                f"difference between the two finest runs of each tool; pass if D <= k max(e_S, e_B). D_extrap_S is D with the finest SpaceAGORA run replaced by its first-order Richardson extrapolation (informational only, not used for pass/fail); D_extrap2_S uses the three finest levels at second order; order D is the observed order of D over the two finest SpaceAGORA levels (D proportional to dt^order).\n\n")
        f.write("| quantity | entity | D | e_S | e_B | bound | D/bound | order S | order B | order D | D_extrap_S | D_extrap2_S | pass |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            f.write(f"| {r['quantity']} | {r['entity']} | {r['D']:.3e} | {r['e_S']:.3e} | {r['e_B']:.3e} | {r['bound']:.3e} | "
                    f"{r['ratio']:.3g} | {r['order_S']:.2f} | {r['order_B']:.2f} | {r['order_D']:.2f} | {r['D_extrap_S']:.2e} | {r['D_extrap2_S']:.2e} | {'PASS' if r['passed'] else 'FAIL'} |\n")
        f.write(f"\n{sum(r['passed'] for r in rows)} of {len(rows)} quantities pass.\n")
        if cfg.get("verdict"):
            f.write(f"\nVerdict: {cfg['verdict']}\n")
    print(open(a.out + ".md").read())
    return 0 if all(r["passed"] for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
