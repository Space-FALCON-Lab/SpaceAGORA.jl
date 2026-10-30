"""Ground tables of the Revision 16 generic orbits (GA to GD) for the planner arms: planner_ground_tables_generic.json (+ .sha256) in the task dir, pinned in planner_config.toml [planner.ground_table_generic].
Same construction as planner_predictions.ground_table (geometry_rank model at the calibrated context: q_cap 0.2503, planning D_ref 70.4 d) but WITHOUT the GRAM mean-field term (no sampling job was run for these geometries: the
preview model; a stated handicap of the planner arms relative to Odyssey and G1). Usage: python3 -I planner_generic_tables.py"""
import json, pathlib, sys
HERE = pathlib.Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import planner_lib as PLB, planner_predictions as PP

cfg = PLB.load_cfg(); GR, G3A = PP._imports(); cfg3 = G3A.load3(); ctx = PP.make_ctx_cal(cfg3, cfg); out = {}
for gid, o in cfg["planner"]["geometry"].items():
    if gid == "ODY":
        continue
    cand = dict(id=gid, inc=float(o["inc_deg"]), epoch_off=float(o["epoch_off_d"]), omega=float(o["omega_deg"]), raan=float(o["raan_deg"]))
    tab, m = PP.ground_table(ctx, cand, None, cfg); tab.update(token=o["role"], field_mode="preview_no_field", inc_deg=cand["inc"], omega_deg=cand["omega"], raan_deg=cand["raan"], epoch_off_d=cand["epoch_off"]); out[gid] = tab
    print(gid, o["role"], "n", tab["n"], "q_req %.4f" % tab["q_req"], "W %.4f" % tab["W"])
doc = dict(schema="planner_ground_tables/1", note="Revision 16 generic orbits; no mean-field term", geometries=out)
p = PLB.task_dir() / "planner_ground_tables_generic.json"; p.write_text(json.dumps(doc, indent=0, default=float)); h = PLB.sha256_file(p)
(PLB.task_dir() / "planner_ground_tables_generic.json.sha256").write_text(h + "  planner_ground_tables_generic.json\n"); print("sha256", h)
