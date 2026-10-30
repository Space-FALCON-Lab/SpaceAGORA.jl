# Estimator v5 driver (Julia half): identity self-tests, Mars and Earth predictions, timing. Serial, single process (SPICE is not thread-safe). One Julia job at a time, memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. run_estimator_v5.jl --phase=selftest|mars|earth|timing
#        [--results-root=DIR] [--config=MARS_TOML] [--name=NAME] [--cells=FILE] [--window=integrated|kepler] [--nwindow=N] [--nback=N] [--cluster=X] [--postdrag=true|false]
#        [--kickpd=true|false] [--zsp=true|false]
#   selftest   identities: (1) window=:kepler and postdrag=false returns exactly v3 where the carry rule is false and exactly v4 where it is true (Earth cells, Mars cells);
#              (2) the corrected stepper with a zero tesseral field equals the v4 zero-field stepper where only the drag window differs (not run: see diagnostics)
#   mars       v5 predictions for the cells of <root>/mars (truth seeds from truth/<tag>/calibration.csv) -> pred_<name>/, pred_<name>_selfcal/, pred_cells_<name>.csv (default name v5)
#   earth      v5 predictions for the Earth cells of <root>/oracle/oracle_cell_labels.csv -> <root>/estimator_<name>/{pred,pred_selfcal,pred_cells.csv}
#   timing     warm cost per 60-day Mars prediction (v3, v4, v5), single job alone
using Statistics, OrdinaryDiffEq
include(joinpath(@__DIR__, "mars_common.jl"))
include(joinpath(@__DIR__, "estimator_v5.jl"))

function parse_cli5()
    kv = Dict{String, String}()
    for a in ARGS
        m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2])
    end
    return kv
end
const CLI5 = parse_cli5()
const STUDY5 = load_study(abspath(get(CLI5, "study", STUDY_CONFIG_PATH)))
const ROOT5 = haskey(CLI5, "results-root") ? abspath(CLI5["results-root"]) : results_root(STUDY5)
const CFG25 = TOML.parsefile(joinpath(STUDY_DIR, "estimator_v2_config.toml"))
const OPTS5 = (v = CFG25["v2"]; V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"])))
const TAB5 = SPTable(OPTS5.n_sp)
const SETT5 = settings_from(CFG25)
const ZO5 = ZonalOpts()
const NAME5 = get(CLI5, "name", "v5")
const TT5 = TesseralOpts(force=haskey(CLI5, "force") ? parse(Bool, CLI5["force"]) : nothing)
const TO5 = TimingOpts(window=Symbol(get(CLI5, "window", "integrated")), n_window=parse(Int, get(CLI5, "nwindow", "32")), n_back=parse(Int, get(CLI5, "nback", "8")),
    cluster=parse(Float64, get(CLI5, "cluster", "1.0")), post_drag_period=parse(Bool, get(CLI5, "postdrag", "true")),
    kick_post_drag=parse(Bool, get(CLI5, "kickpd", "true")), zonal_sp=parse(Bool, get(CLI5, "zsp", "false")))

pred_df5(p::Prediction) = DataFrame(k=1:length(p.t_s), t_s=p.t_s, alt_m=p.alt_m, radius_m=p.radius_m, lat_rad=p.lat_rad, vrel_mps=p.vrel_mps,
    rho=p.rho, q_w_cm2=p.q_w_cm2, a_m=p.a_m, e=p.e, i_rad=p.i_rad, raan_rad=p.raan_rad, omega_rad=p.omega_rad, h_qhigh_m=p.h_qhigh_m, h_qlow_m=p.h_qlow_m)
same_pred(p, q) = all(isequal(getfield(p, f), getfield(q, f)) for f in fieldnames(Prediction))

function earth_record_deg(study, degree::Int, order::Int)
    rec = earth_record(study); cfg = load_corridor_config(CORRIDOR_CONFIG_PATH)
    harm = GravitationalHarmonicsModel(degree, order, joinpath(REPO_ROOT, String(cfg.values["gravity"]["coefficient_file"])), rec.planet)
    return PlanetRecord(; (f => getfield(rec, f) for f in fieldnames(PlanetRecord) if f != :harm)..., harm=harm)
end
function mars_numerics(out)
    nodes_f = joinpath(out, "estimator_numerics.toml")
    return isfile(nodes_f) ? with_settings(SETT5; n_drag_half=Int(TOML.parsefile(nodes_f)["n_drag_half"])) : SETT5
end

function cell_row5(tag, seed_t, seed_s, p1, p2, info1)
    ns, D = info1; d(n) = n - 2 <= length(D) ? D[n - 2] : NaN
    return (tag=tag, seed_truth_m=seed_t, seed_selfcal_m=seed_s, n_pred=length(p1.t_s), decayed_pred=p1.decayed, t_decay_pred_days=p1.t_decay_s / SECONDS_PER_DAY, ended=String(p1.ended),
        n_pred_selfcal=length(p2.t_s), decayed_pred_selfcal=p2.decayed, t_decay_pred_selfcal_days=p2.t_decay_s / SECONDS_PER_DAY,
        n_pred_v1=0, decayed_pred_v1=false, t_decay_pred_v1_days=NaN, wall_s_first_call=NaN,
        carried=join(ns, ";"), D3_km=d(3), D4_km=d(4), D5_km=d(5), D6_km=d(6), D8_km=d(8), D12_km=d(12), D20_km=d(20), DT_km=NaN, tesseral=false, DT_km_selfcal=NaN)
end

function phase_mars()
    mcp = abspath(get(CLI5, "config", MARS_CONFIG_PATH)); mc = load_mars(mcp); ms = mars_setup(mc); out = mars_results_dir(ROOT5, mc)
    ep(i) = parse_epoch(mc["grid"]["epochs_utc"][i]); etm(i) = SimulationModel.ephemerides_time_seconds(ep(i), SpiceEphemeridesModel())
    cells = CSV.read(joinpath(out, get(CLI5, "cells", "cells.csv")), DataFrame); g = Float64(mc["orbit"]["initial_periapsis_altitude_guess_m"])
    S = mars_numerics(out); zf = zonal_field(ms.harm); tf = tesseral_field(ms.harm); sfx = NAME5
    for d in ("pred_$sfx", "pred_$(sfx)_selfcal"); mkpath(joinpath(out, d)); end
    println("estimator v5 at Mars: window=$(TO5.window) nwindow=$(TO5.n_window) nback=$(TO5.n_back) cluster=$(TO5.cluster) postdrag=$(TO5.post_drag_period) force=$(TT5.force)")
    rows = NamedTuple[]; ctxs = Dict{Int, EstContext}(); t0 = time()
    for c in eachrow(cells)
        ei = Int(c.epoch_idx); ctx = get!(ctxs, ei) do; make_context(mars_record(ms, ep(ei)), S); end
        calf = joinpath(out, "truth", c.tag, "calibration.csv"); isfile(calf) || continue
        seed_t = CSV.read(calf, DataFrame).seed_alt_m[end]; inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = etm(ei)
        tth = haskey(CLI5, "oracle-time") ? CSV.read(joinpath(out, "truth", c.tag, "passes_3b.csv"), DataFrame).time_s : nothing     # UNFROZEN diagnostic (name it v5x_*)
        p1, info1 = predict_history5(ctx, S, OPTS5, TAB5, zf, ZO5, tf, TT5, TO5, seed_t, inc, raan, om, et0; theta_times=tth)
        seed_s = calibrate_seed5(ctx, S, OPTS5, TAB5, zf, ZO5, tf, TT5, TO5, inc, raan, om, et0; seed0=g, tol_m=Float64(CFG25["estimator"]["calibration_tol_m"]), max_iter=Int(CFG25["estimator"]["calibration_max_iter"]))
        p2, info2 = predict_history5(ctx, S, OPTS5, TAB5, zf, ZO5, tf, TT5, TO5, seed_s, inc, raan, om, et0; theta_times=tth)
        CSV.write(joinpath(out, "pred_$sfx", c.tag * ".csv"), pred_df5(p1)); CSV.write(joinpath(out, "pred_$(sfx)_selfcal", c.tag * ".csv"), pred_df5(p2))
        r = cell_row5(c.tag, seed_t, seed_s, p1, p2, info1); push!(rows, merge(r, (DT_km=info1[3], tesseral=info1[4])))
    end
    CSV.write(joinpath(out, "pred_cells_$sfx.csv"), DataFrame(rows)); println("mars v5 done: $(length(rows)) cells in $(round(time() - t0; digits=1)) s (includes compilation and SPICE; not a timing result)")
end

function phase_earth()
    nbody_planet(); degree = parse(Int, get(CLI5, "degree", "20")); order = parse(Int, get(CLI5, "order", "20"))
    rec = earth_record_deg(STUDY5, degree, order); zf = zonal_field(rec.harm); tf = tesseral_field(rec.harm)
    out = haskey(CLI5, "out") ? abspath(CLI5["out"]) : joinpath(ROOT5, "estimator_" * NAME5)
    for d in ("pred", "pred_selfcal"); mkpath(joinpath(out, d)); end
    ctx = make_context(rec, SETT5)
    labels = CSV.read(joinpath(ROOT5, "oracle", "oracle_cell_labels.csv"), DataFrame)
    cell_dir5(tag, inc) = inc == 50.0 ? joinpath(ROOT5, "border", "grid", tag[7:end]) : joinpath(ROOT5, "border2", "grid", tag)
    println("estimator v5 at Earth: window=$(TO5.window) nwindow=$(TO5.n_window) nback=$(TO5.n_back) cluster=$(TO5.cluster) postdrag=$(TO5.post_drag_period) force=$(TT5.force)")
    rows = NamedTuple[]; t0 = time()
    for c in eachrow(labels)
        inc, raan, om = deg2rad(c.inc_deg), deg2rad(c.raan_deg), deg2rad(c.omega_deg); et0 = epoch_et(STUDY5, Int(c.epoch))
        seed_t = CSV.read(joinpath(cell_dir5(c.tag, c.inc_deg), "calibration.csv"), DataFrame).seed_alt_m[end]
        p1, info1 = predict_history5(ctx, SETT5, OPTS5, TAB5, zf, ZO5, tf, TT5, TO5, seed_t, inc, raan, om, et0)
        seed_s = calibrate_seed5(ctx, SETT5, OPTS5, TAB5, zf, ZO5, tf, TT5, TO5, inc, raan, om, et0; tol_m=Float64(CFG25["estimator"]["calibration_tol_m"]), max_iter=Int(CFG25["estimator"]["calibration_max_iter"]))
        p2, info2 = predict_history5(ctx, SETT5, OPTS5, TAB5, zf, ZO5, tf, TT5, TO5, seed_s, inc, raan, om, et0)
        CSV.write(joinpath(out, "pred", c.tag * ".csv"), pred_df5(p1)); CSV.write(joinpath(out, "pred_selfcal", c.tag * ".csv"), pred_df5(p2))
        r = cell_row5(c.tag, seed_t, seed_s, p1, p2, info1); push!(rows, merge(r, (DT_km=info1[3], tesseral=info1[4])))
    end
    CSV.write(joinpath(out, "pred_cells.csv"), DataFrame(rows)); println("earth v5 done: $(length(rows)) cells in $(round(time() - t0; digits=1)) s -> $out")
end

"""Identities of the options (regression): with `window = :kepler` and `post_drag_period = false`, v5 is v3 where the carry rule is false and v4 where it is true."""
function phase_selftest()
    n_ok = 0; n_bad = 0; legacy = TimingOpts(window=:kepler, post_drag_period=false, kick_post_drag=false, zonal_sp=false)
    # Mars (carried) and Earth (uncarried) cells
    mc = load_mars(abspath(get(CLI5, "config", MARS_CONFIG_PATH))); ms = mars_setup(mc); out = mars_results_dir(ROOT5, mc); S = mars_numerics(out)
    ep(i) = parse_epoch(mc["grid"]["epochs_utc"][i]); etm(i) = SimulationModel.ephemerides_time_seconds(ep(i), SpiceEphemeridesModel())
    zf = zonal_field(ms.harm); tf = tesseral_field(ms.harm); cells = CSV.read(joinpath(out, "cells.csv"), DataFrame)
    for c in eachrow(cells[1:8:end, :])
        ei = Int(c.epoch_idx); ctx = make_context(mars_record(ms, ep(ei)), S); inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = etm(ei)
        seed = CSV.read(joinpath(out, "truth", c.tag, "calibration.csv"), DataFrame).seed_alt_m[end]
        p4, i4 = predict_history4(ctx, S, OPTS5, TAB5, zf, ZO5, tf, TesseralOpts(), seed, inc, raan, om, et0)
        p5, i5 = predict_history5(ctx, S, OPTS5, TAB5, zf, ZO5, tf, TesseralOpts(), legacy, seed, inc, raan, om, et0)
        s4 = calibrate_seed4(ctx, S, OPTS5, TAB5, zf, ZO5, tf, TesseralOpts(), inc, raan, om, et0; seed0=1.09e5); s5 = calibrate_seed5(ctx, S, OPTS5, TAB5, zf, ZO5, tf, TesseralOpts(), legacy, inc, raan, om, et0; seed0=1.09e5)
        ok = same_pred(p4, p5) && i4[4] == i5[4] && s4 == s5; ok ? (n_ok += 1) : (n_bad += 1); ok || println("MISMATCH Mars ", c.tag, " carried=", i4[4], i5[4], " seeds ", s4, " ", s5)
    end
    println("selftest (Mars, legacy options == v4): $n_ok identical, $n_bad different")
    nbody_planet(); rec = earth_record_deg(STUDY5, 20, 20); zfe = zonal_field(rec.harm); tfe = tesseral_field(rec.harm); ctxe = make_context(rec, SETT5)
    labels = CSV.read(joinpath(ROOT5, "oracle", "oracle_cell_labels.csv"), DataFrame); n_ok = 0; n_bad = 0
    cell_dir5(tag, inc) = inc == 50.0 ? joinpath(ROOT5, "border", "grid", tag[7:end]) : joinpath(ROOT5, "border2", "grid", tag)
    for c in eachrow(labels[1:28:end, :])
        inc, raan, om = deg2rad(c.inc_deg), deg2rad(c.raan_deg), deg2rad(c.omega_deg); et0 = epoch_et(STUDY5, Int(c.epoch))
        seed = CSV.read(joinpath(cell_dir5(c.tag, c.inc_deg), "calibration.csv"), DataFrame).seed_alt_m[end]
        p3, _ = predict_history3(ctxe, SETT5, OPTS5, TAB5, zfe, ZO5, seed, inc, raan, om, et0)
        p5, i5 = predict_history5(ctxe, SETT5, OPTS5, TAB5, zfe, ZO5, tfe, TesseralOpts(), legacy, seed, inc, raan, om, et0)
        s3 = calibrate_seed3(ctxe, SETT5, OPTS5, TAB5, zfe, ZO5, inc, raan, om, et0); s5 = calibrate_seed5(ctxe, SETT5, OPTS5, TAB5, zfe, ZO5, tfe, TesseralOpts(), legacy, inc, raan, om, et0)
        ok = same_pred(p3, p5) && !i5[4] && s3 == s5; ok ? (n_ok += 1) : (n_bad += 1); ok || println("MISMATCH Earth ", c.tag, " carried=", i5[4], " seeds ", s3, " ", s5)
    end
    println("selftest (Earth, legacy options == v3 where the rule is false): $n_ok identical, $n_bad different")
end

"""Warm per-prediction cost (ms) of v3, v4 and v5 on every 5th Mars cell, 5 repeats, single thread, truth seed, 60-day horizon. Run ONLY with the machine otherwise idle."""
function phase_timing()
    mc = load_mars(abspath(get(CLI5, "config", MARS_CONFIG_PATH))); ms = mars_setup(mc); out = mars_results_dir(ROOT5, mc); S = mars_numerics(out)
    ep(i) = parse_epoch(mc["grid"]["epochs_utc"][i]); etm(i) = SimulationModel.ephemerides_time_seconds(ep(i), SpiceEphemeridesModel())
    zf = zonal_field(ms.harm); tf = tesseral_field(ms.harm); cells = CSV.read(joinpath(out, "cells.csv"), DataFrame)
    ctxs = Dict(i => make_context(mars_record(ms, ep(i)), S) for i in 1:length(mc["grid"]["epochs_utc"]))
    jobs = [(ctxs[Int(c.epoch_idx)], deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega), etm(Int(c.epoch_idx)),
             CSV.read(joinpath(out, "truth", c.tag, "calibration.csv"), DataFrame).seed_alt_m[end]) for c in eachrow(cells[1:parse(Int, get(CLI5, "stride", "5")):end, :])]
    rows = NamedTuple[]
    variants = Any[("v3", j -> predict_history3(j[1], S, OPTS5, TAB5, zf, ZO5, j[6], j[2], j[3], j[4], j[5])),
        ("v4", j -> predict_history4(j[1], S, OPTS5, TAB5, zf, ZO5, tf, TesseralOpts(), j[6], j[2], j[3], j[4], j[5])),
        ("v5_kepler_window", j -> predict_history5(j[1], S, OPTS5, TAB5, zf, ZO5, tf, TesseralOpts(), TimingOpts(window=:kepler), j[6], j[2], j[3], j[4], j[5])),
        ("v5", j -> predict_history5(j[1], S, OPTS5, TAB5, zf, ZO5, tf, TesseralOpts(), TO5, j[6], j[2], j[3], j[4], j[5]))]
    for (name, f) in variants
        for j in jobs; f(j); end        # warm-up
        ts = Float64[]
        for rep in 1:5, j in jobs; push!(ts, 1e3 * @elapsed f(j)); end
        push!(rows, (variant=name, n_predictions=length(ts), median_ms=median(ts), q1_ms=quantile(ts, 0.25), q3_ms=quantile(ts, 0.75), mean_ms=mean(ts)))
        println(rows[end])
    end
    mkpath(joinpath(ROOT5, "v5")); CSV.write(joinpath(ROOT5, "v5", get(CLI5, "timing-file", "estimator_v5_timing.csv")), DataFrame(rows))
end

function main5()
    phase = get(CLI5, "phase", "selftest")
    phase == "selftest" && phase_selftest()
    phase == "mars" && phase_mars()
    phase == "earth" && phase_earth()
    phase == "timing" && phase_timing()
end
main5()
