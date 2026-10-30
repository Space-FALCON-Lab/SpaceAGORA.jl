"""Revision 25 (WP25-P) bit-identity digest of key Python outputs with every new switch off: sha256 of a canonical JSON of (a) the target radius constants, (b) planner_rev20_ic.select_ic plans for ODY, G1, G0975
(ref, S-fixed, S-free x horizon arrival/full, no terminal term), (c) the ODY and G1 ground tables rebuilt from the mean field (planner_rev23_geometry), (d) planner_rev19_offline.model_plan, (e) planner_rev18_offline.arrival
on the stored ODY table. Usage: python3 -I planner_p25_digest.py --write F | --check F | --print.  Env SPACEAGORA_RA_TARGET_M must be unset for the default digest."""
import hashlib, json, math, pathlib, sys
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_analysis as A
import planner_rev18_offline as O
import planner_rev19_offline as R19
import planner_rev20_ic as R


def _c(x):
    if isinstance(x, dict): return {str(k): _c(v) for k, v in sorted(x.items(), key=lambda kv: str(kv[0]))}
    if isinstance(x, (list, tuple)): return [_c(v) for v in x]
    if isinstance(x, np.ndarray): return [_c(v) for v in x.tolist()]
    if isinstance(x, (np.floating, float)): return float(x).hex() if math.isfinite(float(x)) else str(x)
    if isinstance(x, (np.integer,)): return int(x)
    return x


def compute():
    d = dict(ra_target=A.RA_TARGET, ra_t_ic=R.RA_T, ra_t_18=O.RA_T)
    d["select_ic"] = {}
    for tok in ("ODY", "G1", "G0975"):
        for sc in ("ref", "S-fixed", "S-free"):
            for hz in ("arrival", "full"):
                ic = R.select_ic(tok, sc, horizon=hz); d["select_ic"]["%s|%s|%s" % (tok, sc, hz)] = {k: v for k, v in ic.items() if k != "labels"}
    import planner_rev23_geometry as X
    G = X.init(); geoms = G["PP"].geometries(G["cfg"], G["cfg3"])
    d["tables"] = {}
    for tok in ("ODY", "G1"):
        g = geoms[tok]; field, _ = G["G3A"].meanfield_for(G["cfg3"], G["ctx"], None, g["g"]); t = X.table_from(X.build(g["cand"], field))
        d["tables"][tok] = {k: hashlib.sha256(np.ascontiguousarray(np.asarray(t[k], float)).tobytes()).hexdigest() for k in ("path", "slow", "w", "H", "a", "e", "t", "rp", "ra")}
    d["model_plan"] = {}
    for tok in ("ODY", "G1"):
        tab = O.load_tab(tok)
        for term in (False, True): d["model_plan"]["%s|%s" % (tok, term)] = R19.model_plan(tab, term)
    d["arrival_ody"] = O.arrival(O.load_tab("ODY"), np.full(300, math.log(0.2)))
    return d


def digest(d=None):
    return hashlib.sha256(json.dumps(_c(compute() if d is None else d), sort_keys=True).encode()).hexdigest()


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[:1] == ["--write"]:
        h = digest(); pathlib.Path(a[1]).write_text(h + "\n"); print(h)
    elif a[:1] == ["--check"]:
        h = digest(); ref = pathlib.Path(a[1]).read_text().strip(); print(h, "IDENTICAL" if h == ref else "DIFFERENT"); sys.exit(0 if h == ref else 1)
    else:
        print(digest())
