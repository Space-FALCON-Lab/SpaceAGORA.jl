"""dV correction (analysis only, no flights). Every run stops at the first closing apoapsis at or below the target, but the keeper has already commanded a burn for the NEXT pass there and the run ends part way through it, so the summary
total_dv_mps (ve ln(m0/m_end)) and the replan.csv sum |dv| both include burn that comes after the target is reached. Corrected dV = ve ln(m0 / m_n) on passes.csv mass_kg at the last pass n that counts (the first closing apoapsis <= target;
else the last finite pass), ve = 220 s x 9.80665 (Isp SOURCED mars_generality_config.toml). Cross-check: it equals the sum |delta_v| of the decisions with target_pass <= n.
Builds corrected COPIES of the summary files (run directories symlinked, originals untouched) so the unchanged analysis scripts can be re-run on them.
Usage: python3 -I planner_dv_correction.py build --planner-root R/planner --stageo-root S/stageo --out-dir O [--json F]"""
import argparse, glob, json, os, pathlib, sys
import numpy as np, pandas as pd
TARGET = 7.408e6; VE = 220.0 * 9.80665


def corrected(d):
    """(dv_corrected, n, dv_decisions_through_n, dv_logged_all) for a run directory, or None."""
    P = pd.read_csv(d / "passes.csv").sort_values("pass").reset_index(drop=True)
    if "mass_kg" not in P:
        return None
    P = P[np.isfinite(P.rho_peak_kg_m3) & np.isfinite(P.vrel_mps)].reset_index(drop=True); ra = P.ra_m.values.astype(float); hit = np.where(np.isfinite(ra) & (ra <= TARGET))[0]
    n = int(hit[0] + 1) if len(hit) else len(P); m = P.mass_kg.values.astype(float); dv = VE * np.log(m[0] / m[n - 1])
    dd = None; all_ = None; f = d / "decisions.csv"
    if f.exists() and f.stat().st_size > 2:
        D = pd.read_csv(f); dd = float(D[D.target_pass <= n].delta_v_mps.abs().sum()); all_ = float(D.delta_v_mps.abs().sum())
    return dv, n, dd, all_


def find_dir(root, phase, rid):
    for c in (root / phase / rid, root / "runs" / rid, root / "mars" / "runs" / rid, root / "dev" / rid, root / "sens" / rid, root / "smoke" / rid, root / "gen1" / rid):
        if (c / "passes.csv").exists():
            return c
    return None


def build(roots, out_dir, rows_out):
    for root, sub in roots:
        dst = out_dir / sub; dst.mkdir(parents=True, exist_ok=True)
        for p in sorted(root.iterdir()):
            if p.is_dir() and not (dst / p.name).exists():
                (dst / p.name).symlink_to(p.resolve())
            elif p.is_file() and p.suffix != ".csv":
                pass
        for f in sorted(glob.glob(str(root / "summary_*.csv"))):
            S = pd.read_csv(f); phase = pathlib.Path(f).name.split("_")[1]; new = S.total_dv_mps.copy() if "total_dv_mps" in S else None
            for i, r in S.iterrows():
                if "run_id" not in S or "total_dv_mps" not in S or str(r.get("status", "ok")) != "ok":
                    continue
                d = find_dir(root, phase, r.run_id); c = corrected(d) if d is not None else None
                if c is None:
                    rows_out.append(dict(root=str(root), summary=pathlib.Path(f).name, run_id=r.run_id, arm=r.get("arm"), status="not_corrected")); continue
                dv, n, dd, all_ = c
                rows_out.append(dict(root=str(root), summary=pathlib.Path(f).name, phase=phase, run_id=r.run_id, arm=str(r.get("arm")), geometry=str(r.get("geometry", "ODY")), dv_summary=float(r.total_dv_mps), dv_corrected=dv,
                                     dv_decisions_to_n=dd, dv_logged_all=all_, post_target_flown=float(r.total_dv_mps) - dv, n_counted=n, n_passes=int(r.passes) if "passes" in r else None))
                S.loc[i, "total_dv_mps"] = dv
            S.to_csv(dst / pathlib.Path(f).name, index=False)


def screen_main(argv):
    """Re-run planner_screen_score.py with a flight() that counts only the decisions up to the last counted pass (sum |dv|, net dv, burns, reversals); every other quantity is the original."""
    HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE)); import planner_screen_score as SS
    orig = SS.flight

    def flight(d):
        r = orig(d); c = corrected(d)
        if c is None:
            return r
        n = c[1]; R = pd.read_csv(d / "replan.csv"); ex = R[R["pass"] <= n].executed_dv_mps.values; b = np.abs(ex) >= 1e-3; sg = np.sign(ex[b])
        r.update(sumdv=float(np.abs(ex).sum()), netdv=float(abs(ex.sum())), burns=int(b.sum()), rev=int((np.diff(sg) != 0).sum()) if b.sum() > 1 else 0)
        assert abs(r["sumdv"] - c[0]) < 1e-6, (d, r["sumdv"], c[0])
        return r
    SS.flight = flight; sys.argv = ["planner_screen_score.py"] + argv; SS.main()


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "screen":
        return screen_main(sys.argv[2:])
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["build"]); ap.add_argument("--planner-root", required=True); ap.add_argument("--stageo-root"); ap.add_argument("--out-dir", required=True); ap.add_argument("--json"); a = ap.parse_args()
    out = pathlib.Path(a.out_dir); rows = []; roots = [(pathlib.Path(a.planner_root), "planner")] + ([(pathlib.Path(a.stageo_root), "stageo")] if a.stageo_root else [])
    build(roots, out, rows)
    df = pd.DataFrame(rows); ok = df[df.get("status", pd.Series(dtype=str)).ne("not_corrected")] if "status" in df else df
    x = ok.dropna(subset=["dv_decisions_to_n"]); print(f"flights corrected {len(ok)}, not corrected {len(df) - len(ok)}; with decision log {len(x)}; max |corrected - decisions through n| = {np.abs(x.dv_corrected - x.dv_decisions_to_n).max():.2e} m/s; "
          f"flights with post-target flown dV > 1e-6: {(ok.post_target_flown > 1e-6).sum()} of {len(ok)}; max {ok.post_target_flown.max():.3f} m/s")
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(rows, indent=1, default=float))


if __name__ == "__main__":
    main()
