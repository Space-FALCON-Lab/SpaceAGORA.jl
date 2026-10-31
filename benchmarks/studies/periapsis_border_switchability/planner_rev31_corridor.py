"""Revision 31.4: the status table with the final-orbit corridor counted and not counted, and every flight in every block whose last-pass q lies outside [0.05, 0.25] W/cm2. Analysis only; Python, read-only on flight data.
q = passes.csv q_w_cm2 (mean-field pass heat rate) of the closing pass (the pass that brings the apoapsis to the target), as planner_rev20_score.flight (q_last). Usage: python3 -I planner_rev31_corridor.py --out-dir DIR"""
import argparse, json, pathlib, sys
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_rev26_val2 as V
import planner_rev30_seed as S30
import planner_rev20_score as SC
import planner_rev21_position as PS


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args(); out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    # (1) status-table rows with and without the corridor
    spec = [("Odyssey 93.5 deg, 7,408 km (v4 + L1, 5 km)", "rev30seed", "OLRT3W4L5", "ODY"), ("Odyssey 93.5 deg (v3, 5 km)", "rev29aim", "OLRT3W3E5", "ODY"), ("G0528 63.4 deg, 7,408 km", "rev30seed", "OLRT3W3L5", "G0528"), ("G1230 116.6 deg, 1,000 km, 5 km aim", "rev30seed", "OLRT3W3E5", "G1230"),
            ("G1038 93.5 deg, 500 km", "rev29aim", "OLRT3W3E5", "G1038"), ("G0486 63.4 deg, 2,000 km", "rev27low", "OLRT3W3E", "G0486"), ("G1 116.6 deg, 7,408 km", "rev26val2", "OLRT3W4L", "G1"), ("G0975M (nominal t*)", "rev28t", "OLRT3W3E", "G0975M"), ("G0975M (default t*)", "rev28t", "OLRT3W3T", "G0975M"),
            ("G0063M 25.2 deg", "rev23geo", "OLRT3W3", "G0063M"), ("G0890 93.5 deg", "rev23geo", "OLRT3W3", "G0890"), ("G1230, earlier 12.5 km aim (rev27low s13-16 + rev28t W3E s17-20)", None, None, "G1230")]
    rows = []
    for lab, blk, arm, geom in spec:
        if blk is None:
            fs = [f for b in ("rev27low", "rev28t") for k, f in V.load(b, None).items() if k[0] == "OLRT3W3E" and k[1] == geom and k[2] == "L0"]
        else: fs = [f for k, f in V.load(blk, None).items() if k[0] == arm and k[1] == geom and k[2] == "L0"]
        n = len(fs); l1 = lambda f: bool(f["pos"] and f["pos"]["l1"] < 1.0); nc = sum(bool(f["ok"] and not f["hard"] and f["window"] and l1(f)) for f in fs); cc = sum(bool(f["ok"] and not f["hard"] and f["window"] and l1(f) and f["corr"]) for f in fs)
        ql = [float(f["q_last"]) for f in fs]; rows.append(dict(label=lab, n=n, all_not_counted=nc, all_counted=cc, corridor=int(sum(f["corr"] for f in fs)), below=int(sum(q < 0.05 for q in ql)), above=int(sum(q > 0.25 for q in ql)), wil_nc=S30.wilson(nc, n), wil_c=S30.wilson(cc, n), q_min=min(ql), q_max=max(ql)))
    # (2) every flight in every block outside the corridor
    root = V.ROOT; flagged = []; nflights = 0; skipped = []
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        try: F = V.load(d.name, None)
        except Exception as ex: skipped.append((d.name, type(ex).__name__)); continue
        for k, f in F.items():
            nflights += 1
            if f["reached"] and (f["q_last"] < 0.05 or f["q_last"] > 0.25): flagged.append(dict(block=d.name, arm=k[0], geom=k[1], level=k[2], seed=k[3], q_last=float(f["q_last"]), qs_last=float(f["qs_last"]), err=float(f["err"]), ei=float(f["summ"].get("ei_km", 2000.0)) if f["summ"] else 2000.0))
    (out / "rev31_corridor.json").write_text(json.dumps(SC.nan_none(dict(rows=rows, flagged=flagged, n_scanned=nflights, skipped=skipped)), indent=1))
    for r in rows: print(f"{r['label']:70s} n={r['n']} ALL not counted {r['all_not_counted']} counted {r['all_counted']} corridor {r['corridor']} below {r['below']} above {r['above']} q [{r['q_min']:.3f},{r['q_max']:.3f}]")
    print("scanned", nflights, "skipped", skipped, "flagged", len(flagged))
    df = pd.DataFrame(flagged); print(df.groupby(["block", "arm", "geom"]).agg(n=("q_last", "size"), qmin=("q_last", "min"), qmax=("q_last", "max")).round(3).to_string() if len(df) else "none")


if __name__ == "__main__":
    main()
