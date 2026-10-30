"""Stage E, part 1: does picking the geometry with estimator v2 save propellant versus a baseline geometry?

No truth simulation. Inputs (under the results root):
  oracle/oracle_cell_labels.csv           per-cell open-loop dV for every policy (superposition, perfect forecast)
  estimator_v2/selection_cell_table.csv   estimator-v2 predictions per cell (columns prefixed v1_ are the v2 values in this directory)
  blind/blind_cell_table.csv              45 blind cells (top-5 picks per epoch are group a_top5)
Output: <results_root>/compare_baselines/{set_policy_stats.csv, set_membership.csv, drift_floor.csv, ratios.csv, fig_*.png, assumptions.txt}

Usage: python3 -I compare_baselines.py [--results-root R] [--boot 10000] [--seed 20261016]
Completed and decayed cells are NEVER pooled: a decayed cell only has the cost to termination, which understates a 60-day cost.
"""
import argparse, math, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from results_root import add_results_root_argument, results_root
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

POLICIES = {"lazy_s0": "lazy_s0", "near_passive": "near_passive", "hold_nearest": "hold_nearest",
            "reset_opposite": "reset_opposite", "dp": "dp"}   # policy -> column prefix (<p>_dv, <p>_burns)
QUIET = [(63.4, 180.0), (116.6, 0.0)]      # SOURCED: selection_skill_config.toml baselines.quiet_set
TOPFRACS = (0.01, 0.05, 0.10)              # ASSUMED reading of "k = 1, 5, 10%"
ASSUMPTIONS = """ASSUMPTIONS (all numbers in this folder inherit them)
1. Superposition: a burn shifts the next periapsis by lever x dV and leaves the later passive drift unchanged (open-loop; the superposition check showed it fails for open-loop execution).
2. Perfect forecast: every policy sees the whole passive history; the estimator-v2 error is NOT in the dV, only in which cells are selected.
3. Unbounded burns, exact execution, one burn opportunity per pass, band shrunk 1 km per side, 60-day horizon.
4. dV of a decayed cell is the cost to termination (<= 60 d). It understates a full 60-day cost and is reported in its own rows, never pooled with completed cells.
5. The 336 cells are a designed sweep (4 inclinations x RAAN x omega x 3 epochs, plus 48 at 50 deg), not a natural prior on mission geometry. Baseline sets are defined by grid membership, which is an assumption.
6. 'Selected' picks are drawn from the same 336 cells, ranked by estimator-v2 predicted lazy_s0 dV (predicted decays last; ties by predicted burns, then tag). Its true cost is then read from the oracle table.
7. B-median = completed cells whose DP dV lies between the 25th and 75th percentile of completed cells (reading of 'middle quartile': the interquartile middle half).
8. Drift floor: lower bound dV >= L_min * max(0, passive_alt_max - passive_alt_min - w), L_min = min(Keplerian lever first, last), w = shrunk_width_km_min. The width varies per pass, so this is approximate, not rigorous.
9. Bootstrap: independent resampling of both sets with replacement, ratio of medians, fixed seed; the 95% interval is the 2.5/97.5 percentile. Selected sets are nested subsets of the same grid; intervals for small k are wide.
"""


def load(root):
    L = pd.read_csv(os.path.join(root, "oracle", "oracle_cell_labels.csv"))
    T = pd.read_csv(os.path.join(root, "estimator_v2", "selection_cell_table.csv"))
    B = pd.read_csv(os.path.join(root, "blind", "blind_cell_table.csv"))
    return L, T, B


def quiet_flag(L):
    return np.array([any(abs(i - qi) < 0.05 and abs(r - qr) < 0.5 for qi, qr in QUIET) for i, r in zip(L.inc_deg, L.raan_deg)])


def build_sets(L, T):
    """Return {set name: tag array}. Selected: rank all 336 by estimator-v2 predicted dV."""
    S = {}
    comp = L.status == "completed"
    S["B-ref50"] = L.tag[np.isclose(L.inc_deg, 50.0)].values
    S["B-grid(336)"] = L.tag.values
    S["B-grid(completed)"] = L.tag[comp].values
    q1, q3 = L.dp_dv[comp].quantile([0.25, 0.75])
    S["B-median"] = L.tag[comp & (L.dp_dv >= q1) & (L.dp_dv <= q3)].values
    S["B-quiet"] = L.tag[quiet_flag(L)].values
    T = T.set_index("tag").loc[L.tag]
    score = np.where(T.v1_decayed.values, np.inf, T.v1_lazy_dv.values.astype(float))
    order = pd.DataFrame(dict(tag=L.tag.values, s=score, b=T.v1_lazy_burns.values)).sort_values(["s", "b", "tag"], kind="stable")
    for f in TOPFRACS:
        k = math.ceil(f * len(L))
        S[f"Selected-top{int(round(f*100))}%(k={k})"] = order.tag.values[:k]
    return S, order


def stats(v):
    v = np.asarray(v, float)
    if len(v) == 0:
        return dict(n=0, median=np.nan, q1=np.nan, q3=np.nan, min=np.nan, max=np.nan)
    return dict(n=len(v), median=np.median(v), q1=np.quantile(v, .25), q3=np.quantile(v, .75), min=v.min(), max=v.max())


def boot_ratio(a, b, n, rng):
    """median(a)/median(b) 95% interval by independent bootstrap."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return np.nan, np.nan
    ia = rng.integers(0, len(a), (n, len(a))); ib = rng.integers(0, len(b), (n, len(b)))
    ma = np.median(a[ia], axis=1); mb = np.median(b[ib], axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = ma / mb
    r = r[np.isfinite(r)]
    return (np.quantile(r, .025), np.quantile(r, .975)) if len(r) else (np.nan, np.nan)


def table(L, S, nboot, seed):
    Li = L.set_index("tag")
    rows, ratios = [], []
    rng = np.random.default_rng(seed)
    for name, tags in S.items():
        sub = Li.loc[tags]
        frac_dec = float((sub.status == "decayed").mean())
        for stat_name in ("completed", "decayed"):
            g = sub[sub.status == stat_name]
            for p in POLICIES:
                d = stats(g[f"{p}_dv"]); bn = stats(g[f"{p}_burns"])
                rows.append(dict(set=name, stratum=stat_name, policy=p, **d, burns_median=bn["median"], burns_min=bn["min"], burns_max=bn["max"],
                                 frac_decayed_in_set=frac_dec, note=("cost_to_termination: understates a 60 d cost" if stat_name == "decayed" else "")))
    R = pd.DataFrame(rows)
    sel = [n for n in S if n.startswith("Selected")]
    bas = [n for n in S if n.startswith("B-")]
    for sn in sel:
        for bn_ in bas:
            for p in POLICIES:
                a = Li.loc[S[sn]]; a = a[a.status == "completed"][f"{p}_dv"].values
                b = Li.loc[S[bn_]]; b = b[b.status == "completed"][f"{p}_dv"].values
                lo, hi = boot_ratio(a, b, nboot, rng)
                ma, mb = (np.median(a) if len(a) else np.nan), (np.median(b) if len(b) else np.nan)
                ratios.append(dict(selected=sn, baseline=bn_, policy=p, n_sel_completed=len(a), n_base_completed=len(b), median_sel=ma, median_base=mb,
                                   ratio_sel_over_base=ma / mb if mb else np.nan, ci95_lo=lo, ci95_hi=hi, saving_mps=mb - ma,
                                   ci_excludes_1=bool(np.isfinite(lo) and (hi < 1 or lo > 1))))
    return R, pd.DataFrame(ratios)


def blind_rows(L, B, nboot, seed):
    """Blind picks (top-5 per epoch, group a_top5): true costs from the A1 pipeline; only lazy_s0, near and dp exist."""
    rng = np.random.default_rng(seed + 1)
    pick = B[B.group == "a_top5"]
    cols = {"lazy_s0": "true_lazy_s0_mps", "near_passive": "true_near_mps", "dp": "true_dp_mps"}
    comp = L[L.status == "completed"]
    base_sets = {"B-grid(completed)": comp, "B-ref50(completed)": comp[np.isclose(comp.inc_deg, 50)]}
    rows = []
    for p, c in cols.items():
        a = pick[pick.status == "completed"][c].dropna().values
        for bn_, bd in base_sets.items():
            b = bd[f"{p}_dv"].values
            lo, hi = boot_ratio(a, b, nboot, rng)
            rows.append(dict(selected="Blind top-5/epoch (15 cells)", baseline=bn_, policy=p, n_sel_completed=len(a), n_base_completed=len(b),
                             n_sel_decayed=int((pick.status == "decayed").sum()), median_sel=np.median(a) if len(a) else np.nan, median_base=np.median(b),
                             ratio_sel_over_base=np.median(a) / np.median(b), ci95_lo=lo, ci95_hi=hi,
                             saving_mps=np.median(b) - np.median(a),
                             note="hold_nearest, reset_opposite not computed for blind cells"))
    return pd.DataFrame(rows)


def drift_floor(L, S):
    Li = L.set_index("tag").copy()
    Li["excursion_km"] = Li.passive_alt_max_km - Li.passive_alt_min_km
    Li["lever_mean"] = 0.5 * (Li.lever_kepler_first_mps_per_km + Li.lever_kepler_last_mps_per_km)
    Li["lever_min"] = Li[["lever_kepler_first_mps_per_km", "lever_kepler_last_mps_per_km"]].min(axis=1)
    Li["drift_dv_mean_lever"] = Li.excursion_km * Li.lever_mean              # cost to cancel the full excursion at the mean lever
    Li["floor_dv"] = Li.lever_min * np.maximum(0.0, Li.excursion_km - Li.shrunk_width_km_min)   # approximate lower bound
    Li["dp_over_floor"] = Li.dp_dv / Li.floor_dv.replace(0, np.nan)
    Li["dp_over_drift"] = Li.dp_dv / Li.drift_dv_mean_lever
    Li["dp_below_floor"] = Li.dp_dv < Li.floor_dv
    rows = []
    for name, tags in S.items():
        sub = Li.loc[tags]
        for st in ("completed", "decayed"):
            g = sub[sub.status == st]
            if len(g) == 0:
                continue
            rows.append(dict(set=name, stratum=st, n=len(g), excursion_km_median=g.excursion_km.median(), dp_dv_median=g.dp_dv.median(),
                             drift_dv_mean_lever_median=g.drift_dv_mean_lever.median(), floor_dv_median=g.floor_dv.median(),
                             dp_over_drift_median=g.dp_over_drift.median(), dp_over_floor_median=g.dp_over_floor.median(),
                             n_dp_below_approx_floor=int(g.dp_below_floor.sum()),
                             spearman_dp_vs_excursion=g.dp_dv.corr(g.excursion_km, method="spearman") if len(g) > 2 else np.nan))
    return pd.DataFrame(rows), Li


def figures(L, S, R, RT, out):
    comp = L[L.status == "completed"].set_index("tag")
    names = [n for n in S]
    fig, ax = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    for j, p in enumerate(("lazy_s0", "near_passive")):
        data = [comp.loc[[t for t in S[n] if t in comp.index], f"{p}_dv"].values for n in names]
        ax[j].boxplot(data, vert=False, labels=[f"{n} (n={len(d)})" for n, d in zip(names, data)], showfliers=True)
        ax[j].set_xlabel(f"{p} dV, completed cells only (m/s, open loop)"); ax[j].set_title(p)
    fig.savefig(os.path.join(out, "fig_dv_by_set.png"), dpi=130); plt.close(fig)
    r = RT[(RT.policy == "near_passive") & (RT.baseline.isin(["B-ref50", "B-grid(completed)", "B-median", "B-quiet"]))]
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    sels = sorted(r.selected.unique()); bases = list(r.baseline.unique()); w = 0.8 / len(bases)
    for k, b in enumerate(bases):
        g = r[r.baseline == b].set_index("selected").loc[sels]
        x = np.arange(len(sels)) + k * w
        ax.errorbar(x, g.ratio_sel_over_base, yerr=[g.ratio_sel_over_base - g.ci95_lo, g.ci95_hi - g.ratio_sel_over_base], fmt="o", capsize=3, label=f"vs {b}")
    ax.axhline(1, color="k", lw=0.8); ax.set_xticks(np.arange(len(sels)) + 0.4 - w / 2); ax.set_xticklabels(sels, rotation=15, fontsize=8)
    ax.set_ylabel("median dV ratio Selected / baseline (near_passive)"); ax.legend(fontsize=8)
    fig.savefig(os.path.join(out, "fig_ratio_near_passive.png"), dpi=130); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(); add_results_root_argument(ap)
    ap.add_argument("--boot", type=int, default=10000); ap.add_argument("--seed", type=int, default=20261016)
    a = ap.parse_args()
    root = results_root(a.results_root); out = os.path.join(root, "compare_baselines"); os.makedirs(out, exist_ok=True)
    L, T, B = load(root)
    S, order = build_sets(L, T)
    R, RT = table(L, S, a.boot, a.seed)
    RB = blind_rows(L, B, a.boot, a.seed)
    D, Li = drift_floor(L, S)
    mem = pd.concat([pd.DataFrame(dict(set=n, tag=t)) for n, t in S.items()])
    mem.to_csv(os.path.join(out, "set_membership.csv"), index=False)
    R.to_csv(os.path.join(out, "set_policy_stats.csv"), index=False)
    pd.concat([RT, RB], ignore_index=True).to_csv(os.path.join(out, "ratios.csv"), index=False)
    D.to_csv(os.path.join(out, "drift_floor.csv"), index=False)
    Li[["excursion_km", "lever_mean", "lever_min", "drift_dv_mean_lever", "floor_dv", "dp_dv", "dp_over_floor", "dp_over_drift", "status"]].to_csv(os.path.join(out, "drift_floor_cells.csv"))
    open(os.path.join(out, "assumptions.txt"), "w").write(ASSUMPTIONS)
    figures(L, S, R, RT, out)
    pd.set_option("display.width", 250, "display.max_columns", 30, "display.max_rows", 500)
    print(ASSUMPTIONS)
    print("== completed-cell stats (decayed rows are in set_policy_stats.csv; cost_to_termination only)")
    print(R[(R.stratum == "completed") & R.policy.isin(["lazy_s0", "near_passive", "dp"])].round(3).drop(columns=["note"]).to_string())
    print("== decayed strata (cost to termination, NOT comparable)")
    print(R[(R.stratum == "decayed") & (R.n > 0) & R.policy.isin(["near_passive"])].round(3).drop(columns=["note"]).to_string())
    print("== ratios Selected/baseline")
    print(pd.concat([RT, RB]).round(3).drop(columns=["note"], errors="ignore").to_string())
    print("== drift / floor"); print(D.round(3).to_string())
    print("Selected top-10 picks:", list(order.tag.values[:10]))


if __name__ == "__main__":
    main()
