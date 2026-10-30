"""Probe GC-2 analysis: realism of the cached GRAM and 10x10 gravity (plan: tasks/2026-10-15_periapsis-region-border-targeting/probe_gc2_plan.md).

Usage: python3 -I probe_gc2_analysis.py --results-root DIR [--out FILE]      reads DIR/probe_gc2/
       python3 -I probe_gc2_analysis.py --selftest                           synthetic data, no flights
Thresholds are the frozen ones of the plan; this script only evaluates them.
"""
import argparse, csv, glob, gzip, io, math, os, sys
import numpy as np

BANDS_KM = [0, 150, 250, 500, 1000, 2000]; GATED_BANDS_UP_TO_KM = 500; MIN_SAMPLES = 200
R1_BIAS, R1_RMS = 0.1, 0.25; CORE_FRAC = 0.05; WINDOW_ALT_M = 250e3
R2_FRAC = 0.8; BOOT = 2000; SPLITS = 500
R4_BIAS, R4_RMS = 0.25, 0.5
CFGS = [("direct", 20), ("cache", 20), ("direct", 10), ("cache", 10)]


def read_csv(path):
    with open(path) as f: return list(csv.DictReader(f))


def col(rows, k): return np.array([float(r[k]) for r in rows])


def lag1(x):
    x = np.asarray(x, float); x = x - x.mean(); d = float(np.sum(x * x))
    return float(np.sum(x[:-1] * x[1:]) / d) if len(x) > 2 and d > 0 else float("nan")


def load_knots(path):
    if not os.path.isfile(path): return None
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        rows = list(csv.reader(f))[1:]
    return np.array([[float(x) for x in r] for r in rows]) if rows else None   # pass, t, alt, r


def load_flight(d):
    rows = read_csv(os.path.join(d, "passes.csv"))
    kn = None
    for n in ("knots.csv.gz", "knots.csv"):
        kn = load_knots(os.path.join(d, n))
        if kn is not None: break
    return dict(q=col(rows, "q_eff_w_cm2"), alt=col(rows, "alt_m"), apo=col(rows, "apo_alt_m"), r=col(rows, "r_pert"), knots=kn)


# ---------------------------------------------------------------- R1
def r1_stats(samples, sigma_of, label):
    """samples: dict of arrays pass, alt_m, rho_cache, rho_direct, vrel_mps. sigma_of(lo_km, hi_km) -> sigma of ln r in that altitude range."""
    e = np.log(samples["rho_cache"] / samples["rho_direct"]); alt = samples["alt_m"] / 1e3; q = samples["rho_direct"] * samples["vrel_mps"] ** 3
    core = np.zeros(len(e), bool)
    for p in np.unique(samples["pass"]):
        m = samples["pass"] == p; core[m] = q[m] > CORE_FRAC * q[m].max()
    out = {}; ok = True
    def block(name, m, sig, gated):
        nonlocal ok
        if m.sum() == 0: return
        b = dict(n=int(m.sum()), bias=float(e[m].mean()), rms=float(np.sqrt(np.mean(e[m] ** 2))), max=float(np.max(np.abs(e[m]))), sigma=sig)
        if sig and sig > 0: b["bias_over_sigma"] = b["bias"] / sig; b["rms_over_sigma"] = b["rms"] / sig
        if gated and m.sum() >= MIN_SAMPLES and sig and sig > 0:
            b["pass"] = bool(abs(b["bias"]) <= R1_BIAS * sig and b["rms"] <= R1_RMS * sig); ok &= b["pass"]
        out[name] = b
    core_hi = float(np.percentile(alt[core], 95)) if core.any() else 0.0
    block("core", core, sigma_of(0, core_hi), True)
    for lo, hi in zip(BANDS_KM[:-1], BANDS_KM[1:]):
        block(f"{lo}-{hi} km", (alt >= lo) & (alt < hi), sigma_of(lo, hi), hi <= GATED_BANDS_UP_TO_KM)
    return dict(label=label, bands=out, passed=bool(ok), n_total=int(len(e)))


def sigma_function(direct_flights):
    al = np.concatenate([f["knots"][:, 2] for f in direct_flights if f["knots"] is not None]) / 1e3
    lr = np.log(np.concatenate([f["knots"][:, 3] for f in direct_flights if f["knots"] is not None]))
    def sig(lo, hi):
        m = (al >= lo) & (al < hi) if hi > lo else np.zeros(len(al), bool)
        return float(np.std(lr[m], ddof=1)) if m.sum() > 50 else float("nan")
    return sig


def seed_diff_rms(flights, passes=(1, 2, 3)):
    """rms of ln r differences between different seeds at the same pass and the same time since the first knot at or below 300 km."""
    tabs = []
    for f in flights:
        kn = f["knots"]
        if kn is None: continue
        d = {}
        for p in passes:
            m = (kn[:, 0] == p) & (kn[:, 2] <= 300e3)
            if m.sum() == 0: continue
            t = kn[m, 1]; d[p] = {int(round(a)): math.log(b) for a, b in zip(t - t.min(), kn[m, 3])}
        tabs.append(d)
    diffs = []
    for i in range(len(tabs) - 1):
        for p in passes:
            a, b = tabs[i].get(p), tabs[i + 1].get(p)
            if a and b: diffs += [a[k] - b[k] for k in a.keys() & b.keys()]
    return float(np.sqrt(np.mean(np.square(diffs)))) if diffs else float("nan")


# ---------------------------------------------------------------- R2 / R3
def per_flight(f, K):
    q, alt, apo, r = f["q"][:K], f["alt"][:K], f["apo"][:K], f["r"][:K]
    lr = np.log(r); within = dict(var=float("nan"), lag=float("nan"))
    kn = f["knots"]
    if kn is not None:
        ss = dof = num = den = 0.0
        for p in np.unique(kn[:, 0]):
            m = (kn[:, 0] == p) & (kn[:, 2] <= WINDOW_ALT_M); x = np.log(kn[m, 3])
            if len(x) < 3: continue
            x = x - x.mean(); ss += float(np.sum(x * x)); dof += len(x) - 1; num += float(np.sum(x[:-1] * x[1:])); den += float(np.sum(x * x))
        if dof: within = dict(var=ss / dof, lag=num / den)
    return dict(q=q, dalt=np.diff(alt), decay=-np.diff(apo), lnr=lr, var_peri=float(np.var(lr, ddof=1)), lag_peri=lag1(lr), var_within=within["var"], lag_within=within["lag"])


def stats_of(pf):
    """pf: list of per-flight dicts. Returns {stat name: value} by family."""
    s = {}
    for fam in ("q", "dalt", "decay"):
        v = np.concatenate([p[fam] for p in pf])
        s[f"{fam}:mean"] = float(v.mean()); s[f"{fam}:std"] = float(v.std(ddof=1))
        for qq in (5, 50, 95): s[f"{fam}:q{qq:02d}"] = float(np.percentile(v, qq))
    for k in ("var_peri", "lag_peri", "var_within", "lag_within"): s[f"lnr:{k}"] = float(np.nanmean([p[k] for p in pf]))
    return s


def w1(a, b, n=400):
    u = (np.arange(n) + 0.5) / n; return float(np.mean(np.abs(np.quantile(a, u) - np.quantile(b, u))))


def family_pool(pf, fam):
    return np.concatenate([p["lnr"] if fam == "lnr" else p[fam] for p in pf])


def compare(D, C, rng_seed=20261022):
    """D, C: lists of per-flight dicts (same seeds, same order). Containment (bootstrap over D) and Wasserstein tests per family."""
    rng = np.random.default_rng(rng_seed); nD = len(D); sD = stats_of(D); sC = stats_of(C)
    boot = {k: [] for k in sD}
    for _ in range(BOOT):
        idx = rng.integers(0, nD, nD); s = stats_of([D[i] for i in idx])
        for k, v in s.items(): boot[k].append(v)
    lo = {k: float(np.nanpercentile(v, 5)) for k, v in boot.items()}; hi = {k: float(np.nanpercentile(v, 95)) for k, v in boot.items()}
    inside = {k: bool(lo[k] <= sC[k] <= hi[k]) for k in sD}
    fams = {"q": "q", "dalt": "dalt", "decay": "decay", "lnr": "lnr"}; out = {}
    half = nD // 2; splits = [rng.permutation(nD) for _ in range(SPLITS)]
    for fam in fams:
        keys = [k for k in sD if k.startswith(fam + ":")]; frac = float(np.mean([inside[k] for k in keys]))
        cc = []; dd = []
        for pm in splits:
            A, B = pm[:half], pm[half:2 * half]
            cc.append(w1(family_pool([C[i] for i in A], fam), family_pool([D[i] for i in B], fam))); dd.append(w1(family_pool([D[i] for i in A], fam), family_pool([D[i] for i in B], fam)))
        med_c, p95_d = float(np.median(cc)), float(np.percentile(dd, 95))
        out[fam] = dict(contained_fraction=frac, w1_cmp_median=med_c, w1_direct_halves_p95=p95_d, w1_direct_halves_median=float(np.median(dd)),
                        passed=bool(frac >= R2_FRAC and med_c <= p95_d), stats={k: dict(cmp=sC[k], direct=sD[k], lo=lo[k], hi=hi[k], inside=inside[k]) for k in keys})
    return dict(families=out, passed=bool(all(v["passed"] for v in out.values())))


# ---------------------------------------------------------------- R4
def cfg_name(r):
    np_ = int(float(r["np"])) if r.get("np") not in (None, "") else 20
    return f"{r['cfg']}{r['degree']}" + ("" if np_ == 20 else f"n{np_}")


def r4_analysis(rows, ref="direct20"):
    ok = [r for r in rows if r["status"] == "ok"]
    val = {}
    for r in ok: val[(cfg_name(r), int(r["restart"]), r["arm"], int(r["seed"]))] = (math.log(float(r["q_eff"])), float(r["alt_m"]), float(r["apo_change_m"]))
    names = ("ln q_eff", "periapsis altitude (m)", "apoapsis change (m)")
    res = {}
    byr = {}
    for k, v in val.items():
        if k[0] == ref and k[2] == "pert": byr.setdefault(k[1], []).append(v)
    sig = [math.sqrt(np.mean([np.var(np.array(v)[:, j], ddof=1) for v in byr.values() if len(v) > 1])) for j in range(3)]
    res["sigma_gram"] = dict(zip(names, sig))
    for cn in sorted({k[0] for k in val} - {ref}):
        blk = {}; allok = True
        for arm in ("pert", "mean"):
            ks = sorted({k[1:] for k in val if k[0] == cn and k[2] == arm and (ref,) + k[1:] in val})
            if not ks: continue
            diff = np.array([np.array(val[(cn,) + k]) - np.array(val[(ref,) + k]) for k in ks])
            for j, nm in enumerate(names):
                b = dict(n=len(ks), bias=float(diff[:, j].mean()), rms=float(np.sqrt(np.mean(diff[:, j] ** 2))), sigma=sig[j])
                b["bias_over_sigma"] = b["bias"] / sig[j]; b["rms_over_sigma"] = b["rms"] / sig[j]
                if arm == "pert":
                    b["pass"] = bool(abs(b["bias"]) <= R4_BIAS * sig[j] and b["rms"] <= R4_RMS * sig[j]); allok &= b["pass"]
                blk[f"{arm}:{nm}"] = b
        res[cn] = dict(items=blk, passed=bool(allok))
    return res


# ---------------------------------------------------------------- driver
def fmt(o, ind=0):
    return "\n".join(" " * ind + l for l in repr(o).split("\n"))


def run(root):
    pg = os.path.join(root, "probe_gc2"); L = []; rep = {}
    runs = os.path.join(pg, "runs")
    ens = {c: [load_flight(os.path.join(runs, f"e{k:02d}_{c[0]}{c[1]}_pert")) for k in range(16) if os.path.isfile(os.path.join(runs, f"e{k:02d}_{c[0]}{c[1]}_pert", "passes.csv"))] for c in CFGS}
    for c, v in ens.items(): L.append(f"ensemble {c[0]}{c[1]}: {len(v)} flights, passes {sorted({len(f['q']) for f in v})}")
    nflights = min(len(v) for v in ens.values()); K = min(len(f["q"]) for v in ens.values() for f in v)
    pf = {c: [per_flight(f, K) for f in v[:nflights]] for c, v in ens.items()}
    D = pf[("direct", 20)]; rep["K_passes"] = K; rep["n_seeds"] = nflights
    L.append(f"K = {K} passes, {nflights} seeds per configuration")
    # R1
    sig = sigma_function(ens[("direct", 20)]); sd = seed_diff_rms(ens[("direct", 20)][:8]); rep["R1"] = {}; rep["R1_seed_diff_rms_ln_r_passes1to3"] = sd
    L.append(f"R1 seed-to-seed rms(ln r difference) at the same pass and knot, passes 1-3, seeds 0-7 consecutive pairs: {sd:.4f}")
    for cn in ("cache20", "cache10", "cache20n40", "cache20n80"):
        p = os.path.join(runs, f"dens_{cn}_mean", "r1_samples.csv")
        if not os.path.isfile(p): continue
        rows = read_csv(p); s = {k: col(rows, k) for k in ("pass", "alt_m", "rho_cache", "rho_direct", "vrel_mps")}
        good = (s["rho_direct"] > 0) & (s["rho_cache"] > 0); rep["R1_dropped_nonpositive"] = rep.get("R1_dropped_nonpositive", {}); rep["R1_dropped_nonpositive"][cn] = int((~good).sum()); s = {k: v[good] for k, v in s.items()}
        rep["R1"][cn] = r1_stats(s, sig, cn); L.append(f"R1 {cn}: pass={rep['R1'][cn]['passed']} n={rep['R1'][cn]['n_total']}")
        for b, v in rep["R1"][cn]["bands"].items(): L.append(f"   {b}: {v}")
    # R2 / R3a
    rep["R2"] = {}
    for c in CFGS[1:]:
        r = compare(D, pf[c]); rep["R2"][f"{c[0]}{c[1]}"] = r
        L.append(f"R2/R3a {c[0]}{c[1]} vs direct20: passed={r['passed']}")
        for fam, v in r["families"].items(): L.append(f"   {fam}: contained {v['contained_fraction']:.2f}, W1 cmp median {v['w1_cmp_median']:.4g} vs direct-halves p95 {v['w1_direct_halves_p95']:.4g} (median {v['w1_direct_halves_median']:.4g}) -> {v['passed']}")
    # R3b
    p20, p10 = os.path.join(runs, "m_direct20_mean"), os.path.join(runs, "m_direct10_mean")
    if os.path.isfile(os.path.join(p20, "passes.csv")) and os.path.isfile(os.path.join(p10, "passes.csv")):
        a, b = load_flight(p20), load_flight(p10); n = 10
        d20, d10 = np.diff(a["alt"][:n]), np.diff(b["alt"][:n]); frac = float(np.sqrt(np.mean((d10 - d20) ** 2)) / np.sqrt(np.mean(d20 ** 2)))
        altrms = float(np.sqrt(np.mean((b["alt"][:n] - a["alt"][:n]) ** 2)))
        sp = float(np.sqrt(np.mean([np.var([f["alt"][i] for f in ens[("direct", 20)]], ddof=1) for i in range(n)])))
        rep["R3b"] = dict(drift_change_rms_ratio=frac, alt_rms_m=altrms, seed_spread_std_m=sp, alt_rms_over_seed_spread=altrms / sp)
        L.append(f"R3b (first {n} passes, mean arm): rms(change10 - change20)/rms(change20) = {frac:.3f}; periapsis altitude rms diff {altrms:.1f} m vs seed spread std {sp:.1f} m (ratio {altrms/sp:.3f})")
    # R4 (frozen), exploratory R4v (variants, one orbit) and E1 R4b (two orbits)
    allf = sorted(glob.glob(os.path.join(pg, "restarts_*.csv")))
    base = [f for f in allf if "restarts_r4b_" not in f and "restarts_r4v_" not in f]; v1 = [f for f in allf if "restarts_r4v_" in f]; v2 = [f for f in allf if "restarts_r4b_" in f]
    def show(tag, rows):
        r4 = r4_analysis(rows); L.append(f"{tag} sigma_GRAM: {r4['sigma_gram']}")
        for cn, v in r4.items():
            if cn == "sigma_gram": continue
            L.append(f"{tag} {cn}: passed={v['passed']}")
            for k, b in v["items"].items(): L.append(f"   {k}: n={b['n']} bias {b['bias']:.4g} ({b['bias_over_sigma']:.2f} sigma) rms {b['rms']:.4g} ({b['rms_over_sigma']:.2f} sigma){' pass=' + str(b['pass']) if 'pass' in b else ''}")
        return r4
    if base:
        rows = read_csv(base[-1]); L.append(f"R4 file {os.path.basename(base[-1])}: {len(rows)} rows, not ok {sum(r['status'] != 'ok' for r in rows)}")
        rep["R4"] = show("R4", rows)
        if v1: rep["R4v"] = show("R4v (exploratory variants, one orbit)", rows + [r for f in v1 for r in read_csv(f)])
    if v2: rep["R4b"] = show("R4b (exploratory E1, two orbits)", [r for f in v2 for r in read_csv(f)])
    # verdicts
    def v(cn, r1k):
        return dict(R1=rep["R1"].get(r1k, {}).get("passed"), R2=rep["R2"][cn]["passed"], R4=rep.get("R4", {}).get(cn, {}).get("passed"))
    ver = {"cache20": v("cache20", "cache20"), "cache10": v("cache10", "cache10"), "direct10": dict(R3a=rep["R2"]["direct10"]["passed"], R4=rep.get("R4", {}).get("direct10", {}).get("passed"))}
    for k in ("cache20", "cache10"): ver[k]["realistic"] = all(x is True for x in ver[k].values())
    ver["direct10"]["realistic"] = all(x is True for x in ver["direct10"].values())
    rep["verdicts"] = ver; L.append(f"VERDICTS: {ver}")
    # timing
    tp = glob.glob(os.path.join(pg, "summary_timed_*.csv"))
    if tp:
        t = [r for r in read_csv(tp[0]) if r["run_id"].startswith("timed_")]
        for r in t: L.append(f"TIMED {r['run_id']}: wall {float(r['run_wall_s']):.1f} s, passes {r['passes']}, per pass {float(r['run_wall_s'])/int(r['passes']):.3f} s, load {r['load1_before']} -> {r['load1_after']}")
    tv = glob.glob(os.path.join(pg, "summary_timedv_*.csv"))
    if tv:
        for r in [r for r in read_csv(tv[0]) if r["run_id"].startswith("timedv_")]: L.append(f"TIMEDV {r['run_id']}: wall {float(r['run_wall_s']):.1f} s, passes {r['passes']}, per pass {float(r['run_wall_s'])/int(r['passes']):.3f} s, load {r['load1_before']} -> {r['load1_after']}")
    return rep, "\n".join(L)


def selftest():
    rng = np.random.default_rng(3)
    # R1 synthetic: error 0 -> pass; bias 0.5 -> fail
    n = 3000; s = dict(**{"pass": np.repeat([1, 2, 3], n // 3), "alt_m": np.tile(np.linspace(100e3, 1500e3, n // 3), 3), "vrel_mps": np.full(n, 4000.0)})
    s["rho_direct"] = np.exp(-(s["alt_m"] - 100e3) / 8e3) * 1e-8
    sg = lambda lo, hi: 0.3
    s["rho_cache"] = s["rho_direct"] * np.exp(rng.normal(0, 0.02, n)); ok = r1_stats(s, sg, "t"); assert ok["passed"], ok
    s["rho_cache"] = s["rho_direct"] * np.exp(0.5); assert not r1_stats(s, sg, "t")["passed"]
    # R2 synthetic: C = D distribution passes, C shifted fails
    def flight(shift=0.0, scale=1.0):
        K = 30; lnr = rng.normal(0, 0.35, K); alt = 1.1e5 + np.cumsum(rng.normal(0, 300, K)); apo = 2.4e7 - 1e4 * np.arange(K) + rng.normal(0, 3e3, K) + 0 * shift
        return dict(q=0.15 * np.exp(rng.normal(shift, 0.3 * scale, K)), alt=alt, apo=apo, r=np.exp(lnr), knots=None)
    D = [per_flight(flight(), 30) for _ in range(16)]; C = [per_flight(flight(), 30) for _ in range(16)]; C2 = [per_flight(flight(shift=1.0), 30) for _ in range(16)]
    r_same, r_shift = compare(D, C), compare(D, C2)
    assert r_same["families"]["q"]["passed"] or r_same["families"]["q"]["contained_fraction"] >= 0.6, r_same["families"]["q"]
    assert not r_shift["families"]["q"]["passed"], r_shift["families"]["q"]
    # R4 synthetic
    rows = []
    for rs in range(1, 6):
        for sd in (1, 2, 3):
            base = (0.15 * math.exp(rng.normal(0, 0.3)), 1.1e5 + rng.normal(0, 500), -3e4 + rng.normal(0, 2e3))
            for c, d in CFGS:
                e = 0.0 if (c, d) == ("direct", 20) else 0.01
                rows.append(dict(restart=rs, cfg=c, degree=d, arm="pert", seed=sd, status="ok", q_eff=base[0] * (1 + e), alt_m=base[1] + 10 * e, apo_change_m=base[2] + 20 * e))
    r4 = r4_analysis(rows); assert r4["cache20"]["passed"], r4["cache20"]
    for r in rows:
        if r["cfg"] == "cache" and r["degree"] == 20: r["q_eff"] *= 3
    assert not r4_analysis(rows)["cache20"]["passed"]
    assert abs(lag1(np.sin(0.2 * np.arange(300))) - math.cos(0.2)) < 0.03
    assert abs(w1(np.arange(100.0), np.arange(100.0) + 2) - 2) < 0.05
    print("SELFTEST PASS")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--results-root"); ap.add_argument("--out"); ap.add_argument("--selftest", action="store_true"); a = ap.parse_args()
    if a.selftest: selftest(); sys.exit(0)
    rep, text = run(os.path.abspath(a.results_root)); print(text)
    if a.out: open(a.out, "w").write(text + "\n")
