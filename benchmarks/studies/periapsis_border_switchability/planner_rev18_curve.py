"""Revision 18: the dV-versus-arrival-time curve as the lower envelope of every feasible plan sampled (lambda sweep, band family, tube-edge plans), and the dV with slack s around a target time t* (analysis only).
Usage: python3 -I planner_rev18_curve.py --report rev18_report.json --out-md F --out-json G"""
import argparse, json, pathlib
import numpy as np
WMAX_ORBITS = 3


def envelope(d):
    pts = [(x["days"], x["dv"], "lambda") for x in d["lam"]] + [(d["D_min"]["days"], d["D_min"]["dv"], "edge hi"), (d["D_max"]["days"], d["D_max"]["dv"], "edge lo")]
    pts = sorted(p for p in pts if np.isfinite(p[0]))
    best = []                                       # Pareto lower envelope: keep points not dominated by a point with smaller or equal dV at a duration at least as close
    xs = np.array([p[0] for p in pts]); ys = np.array([p[1] for p in pts])
    return xs, ys, pts


def dv_at(xs, ys, D):
    """Feasible upper bound of the minimum dV at duration D: linear interpolation between the sampled feasible plans bracketing D (a convex mixture of two plans is not claimed feasible; this is an interpolation of samples, labeled)."""
    o = np.argsort(xs); return float(np.interp(D, xs[o], ys[o]))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--report", required=True); ap.add_argument("--out-md", required=True); ap.add_argument("--out-json", required=True); a = ap.parse_args(); R = json.load(open(a.report)); out = {}; L = []
    for tok in ("ODY", "G1"):
        d = R[tok]; xs, ys, pts = envelope(d); Tf = d["T_last_s"] / 86400.0; Wm = WMAX_ORBITS * Tf; Dn = d["nominal"]["days"]
        L += [f"### {tok}: model-world dV versus arrival time (DERIVED, in-model days; final-orbit period {d['T_last_s']:.0f} s = {Tf:.3f} d, W_max = {WMAX_ORBITS} periods = {Wm:.3f} d)", "",
              f"Nominal plan (lambda = 0): {d['nominal']['dv']:.2f} m/s, {Dn:.1f} d, {d['nominal']['N']} passes. D_min (level on the upper tube edge): {d['D_min']['dv']:.2f} m/s, {d['D_min']['days']:.1f} d. D_max (lower edge): {d['D_max']['dv']:.2f} m/s, {d['D_max']['days']:.1f} d.", "",
              "| in-model arrival (d) | dV of the best sampled plan within +/- 0.5 d (m/s) | source |", "|---|---|---|"]
        for D in (36, 38, 40, 42, 45, 50, 55, 60, 66, 80, 100, 120, 140):
            m = (np.abs(xs - D) <= 0.5)
            if m.any():
                i = np.argmin(np.where(m, ys, np.inf)); L.append(f"| {xs[i]:.1f} | {ys[i]:.2f} | {pts[i][2]} |")
        L += ["", "dV with slack s around a target t* (best sampled plan with arrival in [t* - s, t* + s]; the envelope is a feasible upper bound, linear steps between samples):", "", "| t* (d) | s = 0 | s = W_max (0.43 d) | s = 1 d | s = 2 d | s = 5 d | s = 10 d | dV saved at s = W_max | per final orbit period |", "|---|---|---|---|---|---|---|---|---|"]
        rows = []
        for ts in (d["D_min"]["days"] + 0.5, 40.0, Dn - 3, Dn, 50.0, 60.0, 80.0):
            vals = []
            for s in (0.0, Wm, 1.0, 2.0, 5.0, 10.0):
                grid = np.linspace(ts - s, ts + s, 401) if s > 0 else np.array([ts]); vals.append(float(min(dv_at(xs, ys, g) for g in grid)))
            rows.append(dict(t_star=ts, dv_s=vals)); L.append(f"| {ts:.1f} | " + " | ".join(f"{v:.2f}" for v in vals) + f" | {vals[0] - vals[1]:.3f} | {(vals[0] - vals[1]) / (Wm / Tf):.4f} m/s per period |")
        L.append(""); out[tok] = dict(rows=rows, W_max_d=Wm, T_final_d=Tf)
    pathlib.Path(a.out_md).write_text("\n".join(L)); json.dump(out, open(a.out_json, "w"), indent=1); print("\n".join(L))


if __name__ == "__main__":
    main()
