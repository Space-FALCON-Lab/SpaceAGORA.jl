"""Stage M cell lists (written and run BEFORE any production run). Python only; no flights.
Per planet it reproduces the Stage F20/I cell list (same draw as stagef20_cells.py: numpy default_rng(cell_seed), distinct (RAAN, omega) per (epoch, inclination) group, first n_dev_per_group
of each group development, rest held out) and asserts it is byte-identical to the Stage I list (pinned sha256). From it, per scope:
  se (Mars)      dev = first dev cell of each (epoch, inclination) group (12); held out = 24 cells x 4 seed streams (seed_tag = tag + _s1.._s4)
  broad (Mars)   dev 12, held out 24, ra0 10,038 km
  full (Mars)    optional: the first held-out cell of each (epoch, inclination) group (12), ra0 = the config orbit
  g_earth, g_venus  dev 12, held out 24, Pi-matched ra0
and the seeded paper Table IV draw T4 (48 cells: 24 dev, 24 held out; plan A.2 and Addendum B).
Usage: python3 -I stagem_cells.py --results-root DIR [--set=base|se|broad|full|g|t4|all]  -> DIR/stagem/<planet>/cells_*.csv
       python3 -I stagem_cells.py --selftest"""
import argparse, datetime as dt, hashlib, pathlib, subprocess, sys, tempfile, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent
C = tomllib.load(open(HERE / "stagem_config.toml", "rb"))
S, M = C["stagef"], C["stagem"]
# Pinned sha256 of the Stage I (= Stage F20) cells_base.csv per planet (project log, Stage I entry), and of the T4 list drawn with the seed in the config.
BASE_SHA = dict(mars="254bc49f73f9b34f34a626f37b38b145e25fc5be960ce6e9c263a4ad73f0be83", earth="096b26d002d6a35f2d10520f7fc2af13df645235a99c17aebcd4820876901106",
                venus="1d7c0821db29532e199f01fcb64ec86218768dc80ec1a88c8c7e73b2c075d6db")
T4_SHA = "96dd8d6816432783249e0236301ba946195238dc91f7252f25a3e9988f5c443d"
SCOPES = dict(se="mars", broad="mars", full="mars", g_earth="earth", g_venus="venus")

def tag(i, r, w, e): return f"I{round(10*i):04d}_R{round(r):03d}_W{round(w):03d}_E{e}"
def sha(path): return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()

def base_cells():
    """{planet: DataFrame}: the Stage F20/I draw, same rng call order as stagef20_cells.py."""
    rng = np.random.default_rng(S["cell_seed"]); out = {}
    for pl in S["planets"]:
        g = C["planets"][pl]["grid"]; combos = [(r, w) for r in g["raan_deg"] for w in g["omega_deg"]]; rows = []
        for e in range(1, S["n_epochs"] + 1):
            for inc in g["inclination_deg"]:
                pick = rng.choice(len(combos), size=S["n_per_group"], replace=False)
                for k, j in enumerate(pick):
                    r, w = combos[j]
                    rows.append(dict(tag=tag(inc, r, w, e), inc=inc, raan=r, omega=w, epoch_idx=e, split="dev" if k < S["n_dev_per_group"] else "heldout"))
        out[pl] = pd.DataFrame(rows)
    return out

def dev_half(D):
    """First dev cell of each (epoch, inclination) group in file order (the Stage L rule)."""
    d = D[D.split == "dev"]; return d[~d.duplicated(["epoch_idx", "inc"], keep="first")]

def held_first(D):
    d = D[D.split == "heldout"]; return d[~d.duplicated(["epoch_idx", "inc"], keep="first")]

def scope_table(scope, D):
    sc = M["scope"][scope]; ra0 = sc["ra0_m"]; hz = round(M["horizon_factor"] * sc["predicted_days"], 6)
    dev, held = dev_half(D), D[D.split == "heldout"]
    if scope == "full": dev, held = D.iloc[0:0], held_first(D)
    parts = [dev.assign(seed_tag=dev.tag)]
    if scope == "se":
        parts += [held.assign(seed_tag=held.tag + s) for s in M["se_seed_suffixes"]]
    else:
        parts.append(held.assign(seed_tag=held.tag))
    T = pd.concat(parts, ignore_index=True); T["ra0_m"] = ra0; T["horizon_days"] = hz
    return T

def t4_cells():
    """Paper Table IV draw (plan A.2). One default_rng(t4_cell_seed); one uniform vector of length n per field, in the order below."""
    T, n = M["t4"], M["n_t4"]; rng = np.random.default_rng(M["t4_cell_seed"]); u = lambda lo, hi: rng.uniform(lo, hi, size=n)
    ra = u(T["ra_center_km"] - T["ra_half_width_km"], T["ra_center_km"] + T["ra_half_width_km"])
    alt = u(T["vacuum_alt_center_km"] - T["vacuum_alt_half_width_km"], T["vacuum_alt_center_km"] + T["vacuum_alt_half_width_km"])
    inc, raan, aop = u(*T["inclination_deg"]), u(*T["raan_deg"]), u(*T["aop_deg"])
    ta = u(T["true_anomaly_center_deg"] - T["true_anomaly_half_width_deg"], T["true_anomaly_center_deg"] + T["true_anomaly_half_width_deg"])
    t0, t1 = (dt.datetime.fromisoformat(T[k]) for k in ("epoch_start_utc", "epoch_end_utc")); span = (t1 - t0).total_seconds()
    sec = np.floor(u(0.0, span)).astype(int); cdf = u(*T["cd_factor"])
    rows = [dict(tag=f"T4_{i + 1:02d}", split="dev" if i < M["n_t4_dev"] else "heldout", ra_m=round(ra[i] * 1e3, 3), rp_m=round(T["r_eq_m"] + alt[i] * 1e3, 3),
                 vacuum_alt_km=round(alt[i], 6), inc=round(inc[i], 6), raan=round(raan[i], 6), omega=round(aop[i], 6), ta_deg=round(ta[i], 6),
                 epoch_utc=(t0 + dt.timedelta(seconds=int(sec[i]))).strftime("%Y-%m-%dT%H:%M:%S"), cd_factor=round(cdf[i], 6)) for i in range(n)]
    D = pd.DataFrame(rows); D["ra0_m"] = D.ra_m
    D["horizon_days"] = round(M["horizon_factor"] * M["scope"]["t4"]["predicted_days"], 6); D["seed_tag"] = D.tag
    return D

def write_all(root, which="all"):
    """Writes the lists; returns {relative path: sha256}."""
    root = pathlib.Path(root) / S["results_subdir"]; B = base_cells(); out = {}
    def put(pl, name, D):
        p = root / pl / name; p.parent.mkdir(parents=True, exist_ok=True); D.to_csv(p, index=False); out[f"{pl}/{name}"] = sha(p)
    for pl, D in B.items():
        if which in ("all", "base", "se", "broad", "full", "g"): put(pl, "cells_base.csv", D)
    for sc, pl in SCOPES.items():
        if which in ("all", sc, "g") and not (which == "g" and pl == "mars") and (which != "se" or sc == "se") and (which != "broad" or sc == "broad") and (which != "full" or sc == "full"):
            put(pl, f"cells_{sc}.csv", scope_table(sc, B[pl]))
    if which in ("all", "t4"): put("mars", "cells_t4.csv", t4_cells())
    return out

def selftest():
    ok = []
    def chk(name, cond): ok.append(bool(cond)); print(("PASS " if cond else "FAIL ") + name)
    with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2, tempfile.TemporaryDirectory() as d3:
        h1, h2 = write_all(d1), write_all(d2)
        chk("determinism: two generations give identical sha256 for every file", h1 == h2)
        for pl in S["planets"]: chk(f"base list {pl} sha256 equals the pinned Stage I hash", h1[f"{pl}/cells_base.csv"] == BASE_SHA[pl])
        subprocess.run([sys.executable, "-I", str(HERE / "stagef20_cells.py"), "--results-root", d3], check=True, capture_output=True)
        for pl in S["planets"]:
            chk(f"base list {pl} byte-identical to stagef20_cells.py output", (pathlib.Path(d3) / "stagef20" / pl / "cells_base.csv").read_bytes() == (pathlib.Path(d1) / S["results_subdir"] / pl / "cells_base.csv").read_bytes())
        si = tomllib.load(open(HERE / "stagei_config.toml", "rb")); chk("[stagef] and [planets] equal stagei_config.toml (except results_subdir)",
            {k: v for k, v in S.items() if k != "results_subdir"} == {k: v for k, v in si["stagef"].items() if k != "results_subdir"} and C["planets"] == si["planets"] and C["common"] == si["common"] and C["analysis"] == si["analysis"])
        B = base_cells()
        for pl in S["planets"]:
            D = B[pl]; chk(f"{pl}: 48 cells, 24 dev, 24 heldout, tags unique", len(D) == 48 and (D.split == "dev").sum() == 24 and (D.split == "heldout").sum() == 24 and D.tag.is_unique)
            chk(f"{pl}: dev half = 12, one per (epoch, inclination)", len(dev_half(D)) == 12 and not dev_half(D).duplicated(["epoch_idx", "inc"]).any())
        for sc, pl in SCOPES.items():
            T = scope_table(sc, B[pl]); nd, nh = (T.split == "dev").sum(), (T.split == "heldout").sum()
            exp = dict(se=(12, 96), broad=(12, 24), full=(0, 12), g_earth=(12, 24), g_venus=(12, 24))[sc]
            chk(f"scope {sc}: dev/held rows {exp}", (nd, nh) == exp and T.seed_tag.is_unique)
            chk(f"scope {sc}: held-out tags are Stage I held-out cells, dev tags Stage I dev cells", set(T[T.split == "heldout"].tag) <= set(B[pl][B[pl].split == "heldout"].tag) and set(T[T.split == "dev"].tag) <= set(B[pl][B[pl].split == "dev"].tag))
        Ts = scope_table("se", B["mars"]); chk("se held seeds: 24 cells x _s1.._s4", Ts[Ts.split == "heldout"].seed_tag.str.extract(r"(_s\d)$")[0].value_counts().to_dict() == {s: 24 for s in M["se_seed_suffixes"]})
        T = t4_cells(); T2 = t4_cells(); R = M["t4"]
        chk("T4: 48 cells, 24 dev, 24 heldout, dev first", len(T) == 48 and (T.split == "dev").sum() == 24 and (T.split.iloc[:24] == "dev").all() and T.tag.is_unique)
        chk("T4: same seed gives the same table; another seed differs", T.equals(T2) and not np.allclose(T.inc, np.random.default_rng(1).uniform(*R["inclination_deg"], size=48)))
        rng_ok = lambda x, lo, hi: bool(((x >= lo) & (x <= hi)).all())
        chk("T4 ranges: ra, rp altitude, i, RAAN, AOP, true anomaly", rng_ok(T.ra_m, 1e3 * (R["ra_center_km"] - R["ra_half_width_km"]), 1e3 * (R["ra_center_km"] + R["ra_half_width_km"])) and
            rng_ok(T.rp_m - R["r_eq_m"], 1e3 * (R["vacuum_alt_center_km"] - R["vacuum_alt_half_width_km"]), 1e3 * (R["vacuum_alt_center_km"] + R["vacuum_alt_half_width_km"])) and
            rng_ok(T.inc, *R["inclination_deg"]) and rng_ok(T.raan, *R["raan_deg"]) and rng_ok(T.omega, *R["aop_deg"]) and rng_ok(T.ta_deg, 180 - R["true_anomaly_half_width_deg"], 180 + R["true_anomaly_half_width_deg"]))
        chk("T4 ranges: epoch in 1-31 December 2001, C_D factor in [0.9, 1.1]", T.epoch_utc.between(R["epoch_start_utc"], "2001-12-31T23:59:59").all() and rng_ok(T.cd_factor, *R["cd_factor"]))
        chk("T4: vacuum_alt_km column equals (rp - R_eq)/1e3 (R_eq 3396.2 km)", np.allclose(T.vacuum_alt_km, (T.rp_m - R["r_eq_m"]) / 1e3, atol=1e-3))
        chk("T4: periapsis below apoapsis", (T.rp_m < T.ra_m).all())
        chk("T4 file sha256 equals the pinned value", h1["mars/cells_t4.csv"] == T4_SHA)
    print("SELFTEST", "PASSED" if all(ok) else "FAILED", f"({sum(ok)}/{len(ok)})"); return all(ok)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--results-root"); ap.add_argument("--set", default="all", choices=["all", "base", "se", "broad", "full", "g", "t4"]); ap.add_argument("--selftest", action="store_true"); a = ap.parse_args()
    if a.selftest: sys.exit(0 if selftest() else 1)
    if not a.results_root: ap.error("--results-root is required")
    H = write_all(a.results_root, a.set)
    for k, v in sorted(H.items()): print(v, k)
    for pl in S["planets"]:
        k = f"{pl}/cells_base.csv"
        if k in H: assert H[k] == BASE_SHA[pl], f"{k} differs from the Stage I list"
    if "mars/cells_t4.csv" in H: assert H["mars/cells_t4.csv"] == T4_SHA, "T4 list differs from the pinned hash"

if __name__ == "__main__": main()
