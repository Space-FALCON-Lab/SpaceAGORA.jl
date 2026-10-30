"""Half-lever attribution test (open issue I1): is the 0.5 periapsis response a pass-attribution / selection artifact?

Usage: python3 -I halflever_attribution.py [--selftest] [--dev DIR]   (DIR = tasks/.../stageo_data/dev)
Tests (per arm, 48 Stage O dev flights): (1) re-pair each burn by epoch (first periapsis event after the request time);
(2) lag scan of the per-pass response; (3) distributed-lag regression with drift/drag covariates, matched non-burn
passes, a threshold-style instrument; (4) arms separately.  Lever = 4 a sqrt((1-e)/(1+e))/sqrt(mu/a) from the previous
pass's (ra, rp) (the pre-burn orbit), dv from decisions.csv, apsis rows from passes.csv.
"""
import sys, math, pathlib, csv, glob
import numpy as np
MU = 4.282837e13  # Mars GM, only used to scale the lever; ASSUMED standard value (same as stagek_lib to ~1e-6)

def lever(ra, rp):
    a = 0.5 * (ra + rp); e = (ra - rp) / (ra + rp)
    return 4 * a * math.sqrt((1 - e) / (1 + e)) / math.sqrt(MU / a)

def load_run(run):
    P = list(csv.DictReader(open(run / "passes.csv"))); D = list(csv.DictReader(open(run / "decisions.csv")))
    f = lambda k: np.array([float(p[k]) for p in P])
    run_d = dict(t=f("time_s"), rp=f("rp_m"), ra=f("ra_m"), w=f("omega_rad"), alt=f("alt_m"), q=f("q_w_cm2"),
                 rho=f("rho"), rhof=f("rho_fit"), n=len(P))
    dec = [(float(d["request_time_s"]), int(d["target_pass"]), float(d["delta_v_mps"]),
            0.5 * (float(d["lower_m"]) + float(d["upper_m"]))) for d in D]
    return run_d, dec

def design(R, dec):
    """Per pass index i>=1: dv_i (0 if none), epoch-paired index, lever*dv, band center. Pass number k <-> index k-1."""
    n = R["n"]; dv = np.zeros(n); ctr = np.full(n, np.nan); repair = np.full(n, -1); req = np.full(n, np.nan)
    for t, k, d, c in dec:
        i = k - 1
        if 1 <= i < n:
            ctr[i] = c
            if d != 0:
                dv[i] = d; req[i] = t
                repair[i] = int(np.searchsorted(R["t"], t, side="right"))  # first periapsis event after the request
    x = np.zeros(n)
    for i in range(1, n):
        if dv[i] != 0: x[i] = lever(R["ra"][i - 1], R["rp"][i - 1]) * dv[i]
    return dv, x, ctr, repair, req

def wrap(a): return (a + math.pi) % (2 * math.pi) - math.pi

def lagstats(name, X, Y):
    m = np.abs(X) > 0; x, y = X[m], Y[m]
    beta = (x @ y) / (x @ x); c = np.corrcoef(x, y)[0, 1] if len(x) > 2 else float("nan")
    print(f"  {name:34s} n={m.sum():4d} median ratio {np.median(y / x):7.3f} slope {beta:7.3f} corr {c:6.3f}")
    return (int(m.sum()), float(np.median(y / x)), float(beta), float(c))

def test12(runs, label, win=6):
    """Re-pairing by epoch, and the lag scan (baseline = median non-burn pair change within +-win, excluding pairs near any burn)."""
    mism = tot = 0; out = {}; fr = []
    rows = {l: ([], [], [], []) for l in (-1, 0, 1, 2)}; cum = {l: ([], []) for l in (0, 1, 2)}
    for R, dec in runs:
        dv, x, ctr, repair, req = design(R, dec); n = R["n"]; dr = np.diff(R["rp"], prepend=np.nan); dw = np.diff(R["w"], prepend=np.nan)
        for i in np.nonzero(dv)[0]:
            tot += 1; mism += int(repair[i] != i); fr.append((req[i] - R["t"][i - 1]) / (R["t"][i] - R["t"][i - 1]))
        free = [j for j in range(1, n) if all(dv[jj] == 0 for jj in range(max(1, j - 2), min(n, j + 3)))]  # no burn within +-2 passes
        fa = np.array(free)
        for i in np.nonzero(dv)[0]:
            near = fa[np.abs(fa - i) <= win]
            if len(near) < 3: continue
            b_r, b_w = np.median(dr[near]), np.median(wrap_arr(dw[near]))
            for l in rows:
                j = i + l
                if 1 <= j < n:
                    rows[l][0].append(x[i]); rows[l][1].append(dr[j] - b_r); rows[l][2].append(dv[i]); rows[l][3].append(wrap(dw[j] - b_w))
            for l in cum:
                j = i + l
                if j < n: cum[l][0].append(x[i]); cum[l][1].append(R["rp"][j] - R["rp"][i - 1] - (l + 1) * b_r)
    print(f"{label}: epoch re-pairing: burns {tot}, burn target pass differs from first periapsis after request in {mism}")
    out["repair_mismatch"] = (tot, mism)
    print(f"  request time as a fraction of the periapsis-to-periapsis interval (apoapsis is near 0.5): min {min(fr):.3f} median {np.median(fr):.3f} max {max(fr):.3f}")
    for l in rows:
        X, Y, V, W = (np.array(a) for a in rows[l])
        out[f"lag{l:+d}"] = lagstats(f"lag {l:+d} single pair", X, Y)
        out[f"lag{l:+d}_omega"] = float((V @ W) / (V @ V)) if len(V) else float("nan")
        print(f"  {'':34s} omega shift per m/s {out[f'lag{l:+d}_omega']:+.3e} rad")
    for l in cum:
        X, Y = (np.array(a) for a in cum[l]); out[f"cum{l}"] = lagstats(f"cumulative thru lag {l} (baseline x{l+1})", X, Y)
    return out

def wrap_arr(a): return (np.asarray(a) + math.pi) % (2 * math.pi) - math.pi

def build_reg(runs, covariates=True, lags=(-2, -1, 0, 1, 2)):
    """Stacked per-pass rows; run-demeaned. y = d rp_i; regressors x_{i-l} (l=0 current-pass burn), plus drift/drag covariates."""
    Xs, ys, ws, groups, zs, burn = [], [], [], [], [], []
    for g, (R, dec) in enumerate(runs):
        dv, x, ctr, repair, req = design(R, dec); n = R["n"]; dr = np.diff(R["rp"], prepend=np.nan); dw = wrap_arr(np.diff(R["w"], prepend=np.nan))
        for i in range(3, n - 3):
            row = [x[i - l] for l in lags]
            if covariates:
                row += [dr[i - 1], math.log(max(R["q"][i - 1], 1e-30)),
                        math.log(R["rhof"][i - 1] / R["rho"][i - 1]) if R["rho"][i - 1] > 0 and R["rhof"][i - 1] > 0 else 0.0,
                        (R["alt"][i - 1] - ctr[i]) / 1000 if not np.isnan(ctr[i]) else 0.0]
            Xs.append(row); ys.append(dr[i]); ws.append(dw[i]); groups.append(g); zs.append(R["alt"][i - 1] - ctr[i] if not np.isnan(ctr[i]) else np.nan); burn.append(dv[i] != 0)
    return np.array(Xs), np.array(ys), np.array(ws), np.array(groups), np.array(zs), np.array(burn)

def ols_cluster(X, y, g):
    Xd, yd = X.copy(), y.copy()
    for u in np.unique(g):  # run fixed effects
        m = g == u; Xd[m] -= Xd[m].mean(0); yd[m] -= yd[m].mean()
    b = np.linalg.lstsq(Xd, yd, rcond=None)[0]; r = yd - Xd @ b
    XtXi = np.linalg.pinv(Xd.T @ Xd); meat = np.zeros((X.shape[1],) * 2)
    for u in np.unique(g):
        m = g == u; s = Xd[m].T @ r[m]; meat += np.outer(s, s)
    V = XtXi @ meat @ XtXi
    return b, np.sqrt(np.diag(V))

def test3(runs, label, lags=(-2, -1, 0, 1, 2)):
    out = {}
    for cov in (False, True):
        X, y, w, g, z, burn = build_reg(runs, cov, lags)
        if (np.abs(X[:, lags.index(0)]) > 0).sum() < 5: print(f"{label}: too few burns"); return out
        for tgt, name in ((y, "d rp"), (w, "d omega")):
            b, se = ols_cluster(X, tgt, g)
            s = "  ".join(f"x[{l:+d}] {b[i]:+.3f}({se[i]:.3f})" for i, l in enumerate(lags))
            print(f"{label} {'with covariates' if cov else 'lags only      '} {name:8s}: {s}")
            out[(cov, name)] = (b[lags.index(0)], se[lags.index(0)])
    # matched non-burn passes: same arm, nearest on (z, previous change); effect = (burn pair - matched pair)/lever*dv
    X, y, w, g, z, burn = build_reg(runs, True, lags)
    x0 = X[:, lags.index(0)]; prev = X[:, len(lags)]; clean = np.all(np.abs(X[:, :len(lags)]) == 0, axis=1) & ~np.isnan(z)
    bi = np.nonzero((np.abs(x0) > 0) & ~np.isnan(z))[0]; ci = np.nonzero(clean)[0]
    if len(bi) >= 5 and len(ci) > 20:
        sz, sp = np.nanstd(z[ci]), np.std(prev[ci]); ratios = []
        for i in bi:
            d = ((z[ci] - z[i]) / sz) ** 2 + ((prev[ci] - prev[i]) / sp) ** 2
            kk = ci[np.argsort(d)[:5]]; ratios.append((y[i] - np.mean(y[kk])) / x0[i])
        ratios = np.array(ratios); xs = x0[bi]; ex = np.array([(y[i] - np.mean(y[ci[np.argsort(((z[ci] - z[i]) / sz) ** 2 + ((prev[ci] - prev[i]) / sp) ** 2)[:5]]])) for i in bi])
        beta = (xs @ ex) / (xs @ xs)
        print(f"{label} matched (5 nearest clean non-burn passes on z, prev change): n={len(bi)} median ratio {np.median(ratios):.3f} slope {beta:.3f}")
        out["matched"] = (len(bi), float(np.median(ratios)), float(beta))
    # NOTE: log Q or rho at the burned pass itself is post-treatment (it encodes the new periapsis altitude), so only pass i-1 values are used.
    # control model fit on clean (no burn within +-2 passes) passes only, predicted on burn passes: avoids the x~z collinearity
    # that makes the all-pass regression with z unidentified (dv is almost a function of z)
    for names, cols in (("z+prev", [len(lags) + 0, len(lags) + 3]), ("z+prev+Q+rho-err (pre-burn only)", list(range(len(lags), len(lags) + 4))),
                        ("prev+Q+rho-err (no z)", list(range(len(lags), len(lags) + 3)))):
        if len(bi) < 5 or len(ci) < 50: break
        C = np.c_[np.ones(len(ci)), X[ci][:, cols]]; coef = np.linalg.lstsq(C, y[ci], rcond=None)[0]
        pred = np.c_[np.ones(len(bi)), X[bi][:, cols]] @ coef; ex = y[bi] - pred; xs = x0[bi]
        cr = (len(ci), float(np.sqrt(np.mean((y[ci] - C @ coef) ** 2))))
        print(f"{label} control model [{names}] fit on {cr[0]} clean passes (rmse {cr[1]:.0f} m): burn n={len(bi)} median ratio {np.median(ex / xs):.3f} slope {(xs @ ex) / (xs @ xs):.3f} (unadjusted slope {(xs @ y[bi]) / (xs @ xs):.3f})")
        out["ctrl_" + names] = (len(bi), float(np.median(ex / xs)), float((xs @ ex) / (xs @ xs)))
    zc = X[ci][:, len(lags) + 3]; cc = np.linalg.lstsq(np.c_[np.ones(len(ci)), zc], y[ci], rcond=None)[0] if len(ci) > 20 else [np.nan, np.nan]
    print(f"{label} natural mean reversion on clean passes: d rp = {cc[0]:+.0f} {cc[1]:+.0f}*z[km] m (z = alt before minus band center); burn-pass mean z {np.nanmean(z[bi]) / 1000:+.2f} km, |z| {np.nanmean(np.abs(z[bi])) / 1000:.2f} km vs clean |z| {np.nanmean(np.abs(z[ci])) / 1000:.2f} km")
    # threshold-style instrument: sign-aligned |z| beyond the keeper margin (z known pre-burn); IV slope of dRp on dv-lever
    ok = ~np.isnan(z)
    zin = (np.abs(z[ok]) > 750).astype(float); xx = x0[ok]; yy = y[ok]
    if zin.sum() > 5 and (xx != 0).sum() > 5:
        # first stage: x ~ z (signed, dv corrects toward center) ; instrument = -z (clipped by margin); 2SLS with run FE not applied (pooled)
        Z = np.clip(-z[ok], -1500, 1500) * zin
        fs = (Z @ xx) / (Z @ Z); xh = fs * Z; beta = (xh @ yy) / (xh @ xx) if abs(xh @ xx) > 0 else float("nan")
        print(f"{label} IV (instrument = clipped -z beyond 750 m of band center; z = alt before burn minus center): first-stage slope {fs:+.3e}, IV ratio {beta:.3f}  [INVALID as an instrument if natural mean reversion depends on z, which the clean-pass fit tests; shown for completeness]")
        out["iv"] = float(beta)
    return out


def level_shift(runs, label):
    """Selection-free level estimate: the keeper decides from pass i-1 (its z is selected), so skip it. pre = mean rp over passes i-4..i-2,
    post = mean over i..i+2 (or i+1..i+3), minus the same statistic on clean passes (no burn in i-4..i+3) of the same runs."""
    out = {}
    for name, (a, b) in (("post i..i+2", (0, 3)), ("post i+1..i+3", (1, 4))):
        bx, by, cs, bv, bw, cw = [], [], [], [], [], []
        for R, dec in runs:
            dv, x, *_ = design(R, dec); n = R["n"]
            for i in range(5, n - 4):
                if np.count_nonzero(dv[i - 4:i + 4]) > (1 if dv[i] else 0) or (dv[i] == 0 and np.count_nonzero(dv[i - 4:i + 4])): continue
                st = R["rp"][i + a:i + b].mean() - R["rp"][i - 4:i - 1].mean()
                sw = R["w"][i + a:i + b].mean() - R["w"][i - 4:i - 1].mean()
                if dv[i]: bx.append(x[i]); by.append(st); bv.append(dv[i]); bw.append(sw)
                else: cs.append(st); cw.append(sw)
        if len(bx) < 5: print(f"{label} level shift [{name}]: only {len(bx)} isolated burns"); continue
        base = float(np.median(cs)); X = np.array(bx); Y = np.array(by) - base
        beta = (X @ Y) / (X @ X)
        print(f"{label} level shift [{name}]: isolated burns n={len(bx)}, clean n={len(cs)} (baseline median {base:+.0f} m): median ratio {np.median(Y / X):.3f} slope {beta:.3f} corr {np.corrcoef(X, Y)[0, 1]:.3f}")
        V = np.array(bv); W = np.array(bw) - float(np.median(cw)); ow = (V @ W) / (V @ V)
        print(f"{label} level shift [{name}]: omega level shift {ow:+.3e} rad per m/s (DRIFT-CONFOUNDED, not used; corr {np.corrcoef(V, W)[0, 1]:.2f}; clean omega drift median {np.median(cw):+.2e} rad)")
        out[name] = (len(bx), float(np.median(Y / X)), float(beta), float(ow))
    return out

def selftest():
    rng = np.random.default_rng(1); ok = True
    for beta_true, lagshift in ((0.5, 0), (1.0, 0), (1.0, 1)):
        runs = []
        for s in range(6):
            n = 200; t = np.cumsum(np.full(n, 6e4)) + 3e4; ra = np.full(n, 2.8e7); rp = np.full(n, 3.48e6); w = np.zeros(n)
            rpv = np.zeros(n); dec = []; dvs = {}
            for k in rng.choice(np.arange(10, 190), 15, replace=False):
                dvs[k] = rng.uniform(-0.3, 0.4)
            drift = np.cumsum(rng.normal(0, 30, n)); rp = rp + drift
            for k, dv in dvs.items():
                L = lever(ra[k - 2], rp[k - 2]) * dv; rp[k - 1 + lagshift:] += beta_true * L  # level shift (burn pass is index k-1)
                dec.append((t[k - 2] + 2e4, k, dv, 0.0))
            for k in range(2, n):
                if k not in dvs: dec.append((t[k - 2] + 2e4, k, 0.0, 0.0))
            dec.sort(key=lambda d: d[1])
            R = dict(t=t, rp=rp, ra=ra, w=w, alt=rp - 3.39e6, q=np.full(n, 0.1), rho=np.full(n, 1e-8), rhof=np.full(n, 1e-8), n=n)
            runs.append((R, dec))
        _, y, _, _, _, _ = build_reg(runs, True)
        X, y, w, g, z, burn = build_reg(runs, True)
        b, se = ols_cluster(X, y, g); i0 = 2
        want = (beta_true if lagshift == 0 else 0.0); got = b[i0]
        print(f"selftest beta {beta_true} lag {lagshift}: lag0 coef {got:.3f} (expect {want}), lag-1 coef {b[1]:.3f}")
        ok &= abs(got - want) < 0.15
        if lagshift == 1: ok &= abs(b[3] - beta_true) < 0.15  # effect shows at lag +1 = burn at previous pass index
        mism = sum(int(design(R, dec)[3][i] != i) for R, dec in runs for i in np.nonzero(design(R, dec)[0])[0])
        ok &= mism == 0
    print("SELFTEST", "PASS" if ok else "FAIL"); return ok

def main():
    if "--selftest" in sys.argv: sys.exit(0 if selftest() else 1)
    dev = pathlib.Path(sys.argv[sys.argv.index("--dev") + 1]); arms = {}
    for run in sorted(glob.glob(str(dev / "runs" / "*"))):
        run = pathlib.Path(run)
        try: arms.setdefault(run.name.split("_")[0], []).append(load_run(run))
        except Exception as ex: print("skip", run.name, ex)
    for arm, runs in arms.items():
        print(f"\n=== {arm}: {len(runs)} runs, burns {sum(int((design(R, d)[0] != 0).sum()) for R, d in runs)} ===")
        test12(runs, arm); level_shift(runs, arm); test3(runs, arm)
    # pooled N2
    n2 = [r for a, rs in arms.items() if a != "FO" for r in rs]
    print("\n=== pooled N2 ==="); level_shift(n2, "N2"); test3(n2, "N2")

if __name__ == "__main__": main()
