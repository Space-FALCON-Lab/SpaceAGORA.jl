# Estimator v3 Earth regression (the c0 pattern): the Earth cells of c0_earth_regression.jl (12 attribution cells plus 12 further, 22 distinct) predicted by
#   v2            predict_history2 / calibrate_seed2                                  -> <out>/v2/<tag>_{truthseed,selfcal}.csv
#   v3 empty set  predict_history3 with the carried set forced empty (force = [])      -> <out>/v3_empty/...   (must be byte-identical to v2)
#   v3 criterion  predict_history3 with the budget criterion (default 1 km)           -> <out>/v3_crit/...    (identical to v2 where the carried set is empty)
# and a table <out>/carried.csv with the carried degrees and the criterion values per cell. Compare the trees with cmp.
#   julia --project=. c0_v3_regression.jl --out=DIR --results-root=R
include(joinpath(@__DIR__, "estimator_v3.jl"))
using Statistics
const CLI0 = Dict(m.captures[1] => m.captures[2] for a in ARGS for m in [match(r"^--([^=]+)=(.*)$", a)] if m !== nothing)
const OUT0 = abspath(CLI0["out"]); const ROOT0 = abspath(CLI0["results-root"])
const STUDY = load_study(); const CFG2 = TOML.parsefile(joinpath(STUDY_DIR, "estimator_v2_config.toml"))
const SETT = settings_from(CFG2)
const OPTS = (v = CFG2["v2"]; V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"])))
const TAB = SPTable(OPTS.n_sp)
cell_dir(tag, inc) = inc == 50.0 ? joinpath(ROOT0, "border", "grid", tag[7:end]) : joinpath(ROOT0, "border2", "grid", tag)
pdf(p::Prediction) = DataFrame(k=1:length(p.t_s), t_s=p.t_s, alt_m=p.alt_m, radius_m=p.radius_m, lat_rad=p.lat_rad, vrel_mps=p.vrel_mps, rho=p.rho, q_w_cm2=p.q_w_cm2,
    a_m=p.a_m, e=p.e, i_rad=p.i_rad, raan_rad=p.raan_rad, omega_rad=p.omega_rad, h_qhigh_m=p.h_qhigh_m, h_qlow_m=p.h_qlow_m)

function main0()
    nbody_planet(); for d in ("v2", "v3_empty", "v3_crit"); mkpath(joinpath(OUT0, d)); end
    rec = let r = earth_record(STUDY), cfg = load_corridor_config(CORRIDOR_CONFIG_PATH)
        PlanetRecord(; (f => getfield(r, f) for f in fieldnames(PlanetRecord) if f != :harm)...,
            harm=GravitationalHarmonicsModel(20, 20, joinpath(REPO_ROOT, String(cfg.values["gravity"]["coefficient_file"])), r.planet))
    end
    ctx = make_context(rec, SETT); zf = zonal_field(rec.harm)
    zo_empty = ZonalOpts(force=Int[]); zo_crit = ZonalOpts()
    labels = CSV.read(joinpath(ROOT0, "oracle", "oracle_cell_labels.csv"), DataFrame)
    att = TOML.parsefile(joinpath(STUDY_DIR, "attribution_config.toml"))["cells"]["list"]
    atags = ["I$(lpad(round(Int, 10 * c[1]), 4, '0'))_R$(lpad(round(Int, c[2]), 3, '0'))_W$(lpad(round(Int, c[3]), 3, '0'))_E$(Int(c[4]))" for c in att]
    tags = unique(vcat(atags, labels.tag[1:28:end])); rows = NamedTuple[]
    for tag in tags
        c = labels[labels.tag .== tag, :][1, :]
        inc, raan, om = deg2rad(c.inc_deg), deg2rad(c.raan_deg), deg2rad(c.omega_deg); et0 = epoch_et(STUDY, Int(c.epoch))
        seed = CSV.read(joinpath(cell_dir(tag, c.inc_deg), "calibration.csv"), DataFrame).seed_alt_m[end]
        pa = predict_history2(ctx, SETT, OPTS, TAB, seed, inc, raan, om, et0); sa = calibrate_seed2(ctx, SETT, OPTS, TAB, inc, raan, om, et0; tol_m=1.0, max_iter=8)
        pas = predict_history2(ctx, SETT, OPTS, TAB, sa, inc, raan, om, et0)
        CSV.write(joinpath(OUT0, "v2", tag * "_truthseed.csv"), pdf(pa)); CSV.write(joinpath(OUT0, "v2", tag * "_selfcal.csv"), pdf(pas))
        for (nm, zo) in (("v3_empty", zo_empty), ("v3_crit", zo_crit))
            pb, info = predict_history3(ctx, SETT, OPTS, TAB, zf, zo, seed, inc, raan, om, et0); sb = calibrate_seed3(ctx, SETT, OPTS, TAB, zf, zo, inc, raan, om, et0; tol_m=1.0, max_iter=8)
            pbs, _ = predict_history3(ctx, SETT, OPTS, TAB, zf, zo, sb, inc, raan, om, et0)
            CSV.write(joinpath(OUT0, nm, tag * "_truthseed.csv"), pdf(pb)); CSV.write(joinpath(OUT0, nm, tag * "_selfcal.csv"), pdf(pbs))
            if nm == "v3_crit"
                n = min(length(pa.alt_m), length(pb.alt_m))
                push!(rows, (tag=tag, carried=join(info[1], ";"), D3_km=info[2][1], D4_km=info[2][2], D5_km=info[2][3], D6_km=info[2][4],
                    max_alt_change_m=n == 0 ? NaN : maximum(abs.(pb.alt_m[1:n] .- pa.alt_m[1:n])), n_v2=length(pa.alt_m), n_v3=length(pb.alt_m)))
            end
        end
    end
    CSV.write(joinpath(OUT0, "carried.csv"), DataFrame(rows)); println("wrote $(length(tags)) cells")
end
main0()
