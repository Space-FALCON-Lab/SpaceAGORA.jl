"""Revision 27 sensitivity re-read (analysis only; Python, read-only on flight data): every earlier claim that rests on a per-flight aerobraking dV difference, re-stated with n, a bootstrap 95% CI (10,000 resamples, fixed seed 20261030; pairs resampled
jointly for paired claims, each group resampled for unpaired ones), a sign test (exact two-sided binomial on the nonzero paired differences; unpaired claims get an exact permutation test on the difference of medians), and a verdict by a fixed rule.
Verdict rule (fixed here, ASSUMED): still holds = the CI excludes 0 and (sign p < 0.05, or, with n = 4 where p cannot go below 0.125, all pairs have the claimed sign and the CI lower bound exceeds the one-pair sensitivity 1.6 m/s; unpaired 4 + 4: complete separation of the two groups instead of unanimity);
holds only in sign = the point estimate has the claimed sign and at least 2/3 of the pairs agree, but the rule above is not met; not resolved = otherwise. dV = aerobraking dV to the last counted pass, ve ln(m0/m_last) (planner_rev20_score.flight).
Usage: python3 -I planner_rev27_sens.py --out-dir DIR"""
import argparse, itertools, json, math, os, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev21_val as V
import planner_rev20_score as SC
ROOT = V.ROOT; NB = 10000; SEED = 20261030; S_PAIR = 1.6


def sign_p(d):
    d = [x for x in d if x != 0]; n = len(d); k = sum(1 for x in d if x > 0)
    if n == 0: return 1.0
    kk = min(k, n - k); return min(1.0, 2 * sum(math.comb(n, i) for i in range(kk + 1)) / 2 ** n)


def paired(d, claim_sign):
    d = np.asarray(d, float); rng = np.random.default_rng(SEED); bs = np.median(rng.choice(d, size=(NB, len(d)), replace=True), axis=1); lo, hi = np.percentile(bs, [2.5, 97.5])
    agree = float(np.mean(np.sign(d) == claim_sign)); p = sign_p(d); excl = lo > 0 or hi < 0; ok_sign = np.sign(np.median(d)) == claim_sign
    strong = excl and ok_sign and (p < 0.05 or (len(d) <= 4 and agree == 1.0 and (lo if claim_sign > 0 else -hi) > S_PAIR))
    verdict = "still holds" if strong else ("holds only in sign" if ok_sign and agree >= 2 / 3 else "not resolved")
    return dict(n=len(d), est=float(np.median(d)), mean=float(np.mean(d)), ci=[float(lo), float(hi)], sign_p=p, agree=agree, pos=int((d > 0).sum()), neg=int((d < 0).sum()), verdict=verdict, kind="paired", d=[float(x) for x in d])


def unpaired(a, b, claim_sign):
    a = np.asarray(a, float); b = np.asarray(b, float); rng = np.random.default_rng(SEED)
    bs = np.median(rng.choice(a, size=(NB, len(a)), replace=True), axis=1) - np.median(rng.choice(b, size=(NB, len(b)), replace=True), axis=1); lo, hi = np.percentile(bs, [2.5, 97.5])
    z = np.concatenate((a, b)); n = len(z); obs = np.median(a) - np.median(b); cnt = tot = 0
    for idx in itertools.combinations(range(n), len(a)):
        m = np.zeros(n, bool); m[list(idx)] = True; dd = np.median(z[m]) - np.median(z[~m]); tot += 1; cnt += abs(dd) >= abs(obs) - 1e-12
    p = cnt / tot; excl = lo > 0 or hi < 0; ok_sign = np.sign(obs) == claim_sign
    sep = (claim_sign > 0 and a.min() > b.max()) or (claim_sign < 0 and a.max() < b.min()); strong = excl and ok_sign and (p < 0.05 or (min(len(a), len(b)) <= 4 and sep and (lo if claim_sign > 0 else -hi) > S_PAIR))
    verdict = "still holds" if strong else ("holds only in sign" if ok_sign and (np.sign(np.median(a) - np.median(b)) == claim_sign) and ((claim_sign > 0 and a.min() > b.max()) or (claim_sign < 0 and a.max() < b.min())) else "not resolved")
    return dict(n=(len(a), len(b)), est=float(obs), ci=[float(lo), float(hi)], sign_p=p, agree=float("nan"), verdict=verdict, kind="unpaired (permutation p)", a=[float(x) for x in a], b=[float(x) for x in b])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); args = ap.parse_args(); out = pathlib.Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    V3 = V.load("rev21v", "summary_rev21v_*.csv"); C = V.load("rev21c", "summary_rev21c_*.csv"); V2 = V.load("rev20val", "summary_rev20val*.csv"); R17 = V.load("rev17", "summary_rev17_*.csv"); GE = V.load("rev23geo", "summary_rev23geo_*.csv")
    EI = V.load("rev27ei", "summary_rev27ei_*.csv"); LO = V.load("rev25low", "summary_rev25low_*.csv"); G1 = V.load("gen1", "summary_gen1_*.csv")
    def pair(A, B, filt, armA, armB, seeds=(13, 14, 15, 16)):
        d = []; keys = []
        for k, fa in A.items():
            if k[0] != armA or not filt(k): continue
            kb = (armB,) + k[1:]
            if kb in B: d.append(fa["dv_sum"] - B[kb]["dv_sum"]); keys.append(k[1:])
        return d, keys
    R = {}; S = lambda x: 1 if x > 0 else -1
    # (a) v3 - v2 (claim: +, v3 higher)
    R["a ODY v3-v2"] = (paired(pair(V3, V2, lambda k: k[1] == "ODY" and k[3] <= 16, "OLRT3W3", "OLRT3W")[0], 1), "+0.37 m/s (20 pairs)")
    R["a G1 v3-v2"] = (paired(pair(V3, V2, lambda k: k[1] == "G1" and k[2] == "L0", "OLRT3W3", "OLRT3W")[0], 1), "+0.24 m/s (4 pairs)")
    dall = pair(V3, V2, lambda k: k[1] in ("ODY", "G1") and k[3] <= 16 and k[2] in ("L0", "L1", "L2", "L3", "L4"), "OLRT3W3", "OLRT3W")[0]; R["a all v3-v2"] = (paired(dall, 1), "+0.37 m/s (24 pairs)")
    # (b) timing price: guarded - unguarded
    def extra(Fg, armg, Fc, geom, lvl_ok=lambda l: True, seeds=(13, 14, 15, 16)):
        d = []
        for k, f in Fg.items():
            if k[0] != armg or k[1] != geom or not lvl_ok(k[2]) or k[3] not in seeds: continue
            c = Fc.get(("OLRT3",) + k[1:])
            if c: d.append(f["dv_sum"] - c["dv_sum"])
        return d
    R["b ODY v2 - OLRT3 (Rev 17)"] = (paired(extra(V2, "OLRT3W", R17, "ODY"), 1), "+2.13 m/s (20 pairs)")
    R["b ODY v3 - OLRT3 (Rev 17)"] = (paired(extra(V3, "OLRT3W3", R17, "ODY"), 1), "+3.10 m/s (20 pairs)")
    R["b ODY v3 - OLRT3 (rev21c, s17-20, L0/L2)"] = (paired(extra(V3, "OLRT3W3", C, "ODY", seeds=(17, 18, 19, 20)), 1), "(not in the earlier note as a median; 8 pairs)")
    R["b G1 v2 - OLRT3"] = (paired(extra(V2, "OLRT3W", V2, "G1"), 1), "+1.58 m/s (4 pairs)")
    R["b G1 v3 - OLRT3"] = (paired(extra(V3, "OLRT3W3", V2, "G1"), 1), "+1.95 m/s (4 pairs)")
    R["b G0063 v3 - OLRT3"] = (paired(extra(V3, "OLRT3W3", C, "G0063"), 1), "+0.45 m/s (4 pairs)")
    R["b G0975 L0 v3 - OLRT3"] = (paired(extra(V3, "OLRT3W3", C, "G0975", lambda l: l == "L0"), 1), "+11.41 m/s (4 pairs)")
    R["b G0975 L2p v3 - OLRT3"] = (paired(extra(V3, "OLRT3W3", C, "G0975", lambda l: l == "L2p"), 1), "+6.71 m/s (4 pairs)")
    # (c) geometry spread (unpaired) and G0063M vs preview
    g = lambda F, geom, arm="OLRT3W3": [f["dv_sum"] for k, f in F.items() if k[0] == arm and k[1] == geom and k[2] == "L0" and k[3] in (13, 14, 15, 16)]
    R["c G0890 - G0528 (unpaired)"] = (unpaired(g(GE, "G0890"), g(GE, "G0528"), 1), "8.15 vs 2.39 m/s (difference 5.76; P23-6 needs >= 2.0)")
    R["c G0890 - G0528 (paired by seed)"] = (paired([a - b for a, b in zip(g(GE, "G0890"), g(GE, "G0528"))], 1), "same, seeds paired (same perturbation stream)")
    R["c G0063M - G0063 preview (unpaired)"] = (unpaired(g(GE, "G0063M"), g(V3, "G0063"), -1), "2.67 vs 3.94 m/s (mean-field table lower)")
    R["c G0063M - G0063 preview (paired by seed)"] = (paired([a - b for a, b in zip(g(GE, "G0063M"), g(V3, "G0063"))], -1), "same, seeds paired")
    # (d) S-free - S-fixed
    R["d ODY WF3-W3"] = (paired(pair(V3, V3, lambda k: k[1] == "ODY" and k[2] == "L0", "OLRT3WF3", "OLRT3W3")[0], 1), "+0.50 m/s (4 pairs); model-world saving would be negative")
    R["d G1 WF3-W3"] = (paired(pair(V3, V3, lambda k: k[1] == "G1" and k[2] == "L0", "OLRT3WF3", "OLRT3W3")[0], -1), "-0.34 m/s (4 pairs)")
    # (e) EI 150 - EI 2000
    dei = [EI[("OLRT3W3E", "ODY", "L0", s)]["dv_sum"] - V3[("OLRT3W3", "ODY", "L0", s)]["dv_sum"] for s in (13, 14, 15, 16)]
    R["e rev27ei EI150-EI2000 (perturbed)"] = (paired(dei, -1), "median -0.12 m/s; |median| 1.04 > 0.3 (smallness rule)")
    Sb = pd.read_csv(next(ROOT.glob("summary_rev27b_*.csv"))).set_index("run_id"); d27b = float(Sb.loc["ODY_PL_OLRT3-W3E_pert_ODY_s13_mean"].total_dv_mps - Sb.loc["ODY_PL_OLRT3-W3_pert_ODY_s13_mean"].total_dv_mps)
    R["e rev27b EI150-EI2000 (unperturbed, 1 pair)"] = (dict(n=1, est=d27b, ci=[float("nan")] * 2, sign_p=float("nan"), agree=float("nan"), verdict="not resolved (n = 1)", kind="single pair"), "-1.59 m/s (30 vs 49 burns)")
    R["e pooled 5 pairs (rev27ei x4 + rev27b)"] = (paired(dei + [d27b], -1), "pooled EI 150 - EI 2000 (sign of the unperturbed pair: less dV at EI 150)")
    # (f) other stated dV findings
    def stop_dv(geom):
        o = []
        for k, f in sorted(LO.items()):
            if k[0] == "OLRT3W3" and k[1] == geom and k[2] == "L0":
                m = pd.read_csv(pathlib.Path(f["dir"]) / "passes.csv").sort_values("pass").mass_kg.values.astype(float); o.append(SC.VE * math.log(m[0] / m[-1]))      # dV to the stop (the scorer's flight() stops at the 7,408 km pass)
        return o
    R["f rev25low G1038 - G0486 (unpaired)"] = (unpaired(stop_dv("G1038"), stop_dv("G0486"), 1), "G1038 9.08 vs G0486 4.17 m/s (P25-4: median <= 4)")
    R["f rev25low G1038 - G1230 (unpaired)"] = (unpaired(stop_dv("G1038"), stop_dv("G1230"), 1), "G1038 9.08 vs G1230 2.19 m/s")
    for arm, lab in (("BODYT", "B-ODY-T"), ("N2IT", "N2-I-T")):
        R[f"f G0975 L2p {lab} - OLRT3 (rev21c, paired by seed)"] = (paired([a - b for a, b in zip([f["dv_sum"] for k, f in sorted(C.items(), key=lambda t: t[0][3]) if k[0] == arm and k[1] == "G0975" and k[2] == "L2p"], [f["dv_sum"] for k, f in sorted(C.items(), key=lambda t: t[0][3]) if k[0] == "OLRT3" and k[1] == "G0975" and k[2] == "L2p"])], 1), f"baseline vs OLRT3 dV (Revision 21 note: baselines 'lower dV' claim for {lab} was 3.7/5.1 vs OLRT3 3.1)")
        R[f"f G0063 L0 {lab} - OLRT3 (rev21c, paired by seed)"] = (paired([a - b for a, b in zip([f["dv_sum"] for k, f in sorted(C.items(), key=lambda t: t[0][3]) if k[0] == arm and k[1] == "G0063"], [f["dv_sum"] for k, f in sorted(C.items(), key=lambda t: t[0][3]) if k[0] == "OLRT3" and k[1] == "G0063"])], -1), f"baseline lower than OLRT3 at G0063 (1.67/3.50 vs 3.72)")
    # sensitivity estimate
    dei_a = np.array(dei); sd_ei = float(np.std(dei_a, ddof=1)); sens = dict(rev27b_pair=d27b, rev27ei_diffs=dei, rev27ei_sd=sd_ei, rev27ei_rms=float(np.sqrt(np.mean(dei_a ** 2))), sigma_flight_from_ei=sd_ei / math.sqrt(2), sigma_flight_from_ei_rms=float(np.sqrt(np.mean(dei_a ** 2)) / math.sqrt(2)), sigma_flight_from_27b=abs(d27b) / math.sqrt(2),
                                                                      within_cell_sd=float(np.sqrt(np.mean([np.var([f['dv_sum'] for k2, f in V3.items() if k2[0] == 'OLRT3W3' and k2[1] == 'ODY' and k2[2] == lv and k2[3] <= 16], ddof=1) for lv in ('L0', 'L1', 'L2', 'L3', 'L4')]))), pooled_rms_5=float(np.sqrt(np.mean(np.append(dei_a, d27b) ** 2))), v3v2_sd=float(np.std(dall, ddof=1)), v3v2_n=len(dall))
    res = {k: dict(v[0], claim=v[1]) for k, v in R.items()}; (out / "rev27_sensitivity.json").write_text(json.dumps(dict(claims=res, sensitivity=sens, nb=NB, seed=SEED), indent=1, default=float))
    for k, v in res.items(): print(f"{k:55s} n={v['n']} est={v['est']:+.2f} ci=[{v['ci'][0]:+.2f},{v['ci'][1]:+.2f}] p={v['sign_p']:.3f} {v['verdict']}")
    print(sens)


if __name__ == "__main__":
    main()
