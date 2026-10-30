"""Venus gate V-G(V): predicted periapsis tracks of the gate geometries for the mean-field sampler stageg3v_meanfield.jl (the Venus counterpart of geometry_rank.write_tracks, which is not edited).
Per geometry of geometry_gate_table_venus.json (sha256 verified against stageg3v_config.toml): the constant-H model schedule (Venus ranking context), body-fixed longitude and geodetic latitude of every periapsis.
Usage: python3 -I stageg3v_tracks.py OUT_CSV"""
import hashlib, json, pathlib, sys, tomllib
import numpy as np, pandas as pd
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import geometry_rank as G
import geometry_rank_venus as V


def main(out):
    cfg3 = tomllib.loads((HERE / "stageg3v_config.toml").read_text()); t = cfg3["stageg3v"]["table"]; raw = (V.notes_dir() / t["path"]).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == t["sha256"], "geometry table sha256 mismatch"
    T = json.loads(raw); ctx, info = V.make_ctx(); rows = []
    for g in T["geometries"]:
        c = dict(id=g["id"], inc=g["inclination_deg"], epoch_off=g["epoch_offset_days"], omega=g["omega0_deg"], raan=g["raan0_deg"])
        q = G.solve_q_req(ctx, c, ctx.cfg["band"]["d_ref_days"]); S = G.schedule(ctx, c, q); lon, lat, _ = G.lon_lat(ctx, c, S)
        for k in range(S["n"]): rows.append((g["id"], g["epoch_utc"], g["epoch_offset_days"], g["omega0_deg"], g["raan0_deg"], g["inclination_deg"], k, S["t"][k], lat[k], lon[k], ctx.cfg["orbit"]["h_ref_m"]))
    df = pd.DataFrame(rows, columns="cand epoch_utc epoch_off_d omega0_deg raan0_deg inc_deg k t_s lat_gd_rad lon_rad h_center_m".split())
    pathlib.Path(out).parent.mkdir(parents=True, exist_ok=True); df.to_csv(out, index=False); print(df.groupby("cand").agg(n=("k", "size"), t_end_d=("t_s", lambda x: x.max() / 86400), lat_min=("lat_gd_rad", "min"), lat_max=("lat_gd_rad", "max")).round(3).to_string())


if __name__ == "__main__":
    main(sys.argv[1])
