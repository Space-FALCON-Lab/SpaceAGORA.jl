# Stage C, C0: Earth regression test of the planet-record parametrization.
# Runs the Earth estimator (v1 truth seed, v2 truth seed, v2 self-calibrated seed, B1 hand model) on the 12 attribution cells plus 12 further cells
# of the 336 (every 28th label row) with the sources found in --src (a directory holding common.jl, estimator_v1.jl, estimator_v2.jl,
# shortperiod_analytic.jl, study_config.toml, estimator_v2_config.toml, attribution_config.toml) and writes one CSV per cell and predictor to --out.
# Run it once with the frozen pre-parametrization copy and once with this directory, then compare the two output trees byte for byte (cmp).
#   julia --project=. c0_earth_regression.jl --src=DIR --out=DIR --results-root=R
const CLI0 = Dict(m.captures[1] => m.captures[2] for a in ARGS for m in [match(r"^--([^=]+)=(.*)$", a)] if m !== nothing)
const SRC = abspath(CLI0["src"]); const OUT0 = abspath(CLI0["out"]); const ROOT0 = abspath(CLI0["results-root"])
include(joinpath(SRC, "shortperiod_analytic.jl"))
using Statistics

const STUDY = load_study(joinpath(SRC, "study_config.toml"))
const CFG2 = TOML.parsefile(joinpath(SRC, "estimator_v2_config.toml"))
const SETT = settings_from(CFG2)
const OPTS = (v = CFG2["v2"]; V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"])))
const TAB = SPTable(OPTS.n_sp)
cell_dir(tag, inc) = inc == 50.0 ? joinpath(ROOT0, "border", "grid", tag[7:end]) : joinpath(ROOT0, "border2", "grid", tag)
epoch_et0(idx) = SimulationModel.ephemerides_time_seconds(parse_epoch(STUDY["grid"]["epochs_utc"][idx]), SpiceEphemeridesModel())
pdf(p::Prediction) = DataFrame(k=1:length(p.t_s), t_s=p.t_s, alt_m=p.alt_m, radius_m=p.radius_m, lat_rad=p.lat_rad, vrel_mps=p.vrel_mps, rho=p.rho, q_w_cm2=p.q_w_cm2,
    a_m=p.a_m, e=p.e, i_rad=p.i_rad, raan_rad=p.raan_rad, omega_rad=p.omega_rad, h_qhigh_m=p.h_qhigh_m, h_qlow_m=p.h_qlow_m)

function main0()
    nbody_planet(); mkpath(OUT0)
    ctx = make_context(STUDY, SETT)
    labels = CSV.read(joinpath(ROOT0, "oracle", "oracle_cell_labels.csv"), DataFrame)
    att = TOML.parsefile(joinpath(SRC, "attribution_config.toml"))["cells"]["list"]
    atags = ["I$(lpad(round(Int, 10 * c[1]), 4, '0'))_R$(lpad(round(Int, c[2]), 3, '0'))_W$(lpad(round(Int, c[3]), 3, '0'))_E$(Int(c[4]))" for c in att]
    tags = unique(vcat(atags, labels.tag[1:28:end]))
    for tag in tags
        c = labels[labels.tag .== tag, :][1, :]
        inc, raan, om = deg2rad(c.inc_deg), deg2rad(c.raan_deg), deg2rad(c.omega_deg); et0 = epoch_et0(Int(c.epoch))
        seed = CSV.read(joinpath(cell_dir(tag, c.inc_deg), "calibration.csv"), DataFrame).seed_alt_m[end]
        p1 = predict_history(ctx, SETT, seed, inc, raan, om, et0)
        p2 = predict_history2(ctx, SETT, OPTS, TAB, seed, inc, raan, om, et0)
        seed_s = calibrate_seed2(ctx, SETT, OPTS, TAB, inc, raan, om, et0; tol_m=1.0, max_iter=8)
        p3 = predict_history2(ctx, SETT, OPTS, TAB, seed_s, inc, raan, om, et0)
        h = hand_history(ctx, SETT, seed, inc, raan, om, et0, p1.alt_m[1], p1.vrel_mps[1])
        CSV.write(joinpath(OUT0, tag * "_v1.csv"), pdf(p1)); CSV.write(joinpath(OUT0, tag * "_v2.csv"), pdf(p2))
        CSV.write(joinpath(OUT0, tag * "_v2_selfcal.csv"), pdf(p3)); CSV.write(joinpath(OUT0, tag * "_hand.csv"), DataFrame(h))
        # closed-form series pieces that use the primary and the frame (bodystate): first-pass mean state and 5 windows
        m0, = (m = SV5(p2.a_m[1], p2.e[1] * cos(p2.omega_rad[1]), p2.e[1] * sin(p2.omega_rad[1]), p2.i_rad[1], p2.raan_rad[1]); (m,))
        tp = p2.t_s[1:min(end, 40)]
        an = analytic_series(ctx, m0, et0, tp, ["Moon", "Sun"])
        op = oscillation_params(ctx, m0, et0, tp[1], tp[end], "Moon")
        CSV.write(joinpath(OUT0, tag * "_analytic.csv"), DataFrame(t_s=tp, rbar=an.rbar, sp=an.sp))
        open(joinpath(OUT0, tag * "_osc.txt"), "w") do io
            for k in sort(collect(keys(op))); println(io, k, " ", op[k]); end
        end
    end
    println("wrote $(length(tags)) cells x 6 files to $OUT0")
end
main0()
