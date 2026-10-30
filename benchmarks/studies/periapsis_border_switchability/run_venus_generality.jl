# Venus generality driver (Julia half): truth sweep, resolution rule, frozen estimators v3 and v5. One Julia job at a time, memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. run_venus_generality.jl --phase=cells|pitable|convergence|pilot|truth|predict [--results-root=DIR] [--workers=N] [--only=tag,tag] [--config=TOML]
# Plan and pre-registered expectations: tasks/2026-10-15_periapsis-region-border-targeting/venus_generality_plan.md. Writes under <results root>/venus.
#   cells        cells.csv (written BEFORE any truth run)
#   pitable      planet_constants.csv, body_distances.csv, pole_drift.toml (no simulation)
#   convergence  numerical-resolution rule of plan section 5 (no truth): drag nodes 32/64/128 against a 4x reference; v5 window (32, 8) against (192, 16) -> estimator_numerics.toml
#   pilot        one cell, serial (truth; first cell or --only)
#   truth        the 32-cell truth sweep on the process route (run_monte_carlo via run_process_campaign; route asserted :process)
#   predict      frozen v3 and v5 (truth seed and self-calibrated seed) for every cell -> pred_v3, pred_v3_selfcal, pred_v5, pred_v5_selfcal, pred_cells_v3.csv, pred_cells_v5.csv
using Statistics, OrdinaryDiffEq
include(joinpath(@__DIR__, "venus_common.jl"))
include(joinpath(@__DIR__, "estimator_v3.jl"))
include(joinpath(@__DIR__, "estimator_v5.jl"))

function parse_cli()
    kv = Dict{String, String}()
    for a in ARGS
        m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2])
    end
    return kv
end
const CLI = parse_cli()
const VC = load_venus(abspath(get(CLI, "config", VENUS_CONFIG_PATH)))
const ROOT = haskey(CLI, "results-root") ? abspath(CLI["results-root"]) : abspath(get(ENV, "SPACEAGORA_BORDER_RESULTS", joinpath(REPO_ROOT, "output", "periapsis_border_switchability")))
const OUT = mars_results_dir(ROOT, VC)
const CFG2 = TOML.parsefile(joinpath(STUDY_DIR, "estimator_v2_config.toml"))
const OPTS = (v = CFG2["v2"]; V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"])))
const TAB = SPTable(OPTS.n_sp)
const SETT0 = settings_from(CFG2)
const ZO = ZonalOpts()
const TT = TesseralOpts()
const NUMERICS_PATH = joinpath(OUT, "estimator_numerics.toml")

epoch_of(idx) = parse_epoch(VC["grid"]["epochs_utc"][idx])
et_of(idx) = SimulationModel.ephemerides_time_seconds(epoch_of(idx), SpiceEphemeridesModel())
guess_alt() = Float64(VC["orbit"]["initial_periapsis_altitude_guess_m"])
cal_kw() = (tol_m=Float64(CFG2["estimator"]["calibration_tol_m"]), max_iter=Int(CFG2["estimator"]["calibration_max_iter"]))

pred_df(p::Prediction) = DataFrame(k=1:length(p.t_s), t_s=p.t_s, alt_m=p.alt_m, radius_m=p.radius_m, lat_rad=p.lat_rad, vrel_mps=p.vrel_mps,
    rho=p.rho, q_w_cm2=p.q_w_cm2, a_m=p.a_m, e=p.e, i_rad=p.i_rad, raan_rad=p.raan_rad, omega_rad=p.omega_rad, h_qhigh_m=p.h_qhigh_m, h_qlow_m=p.h_qlow_m)

"""Estimator numerics chosen by `convergence` (defaults if the file is absent): (drag nodes per half, v5 n_window, v5 n_back)."""
function numerics()
    nd = Int(SETT0.n_drag_half); nw, nb = Int.(VC["estimator_numerics"]["v5_window_default"])
    if isfile(NUMERICS_PATH)
        t = TOML.parsefile(NUMERICS_PATH); nd = Int(t["n_drag_half"]); nw = Int(t["n_window"]); nb = Int(t["n_back"])
    end
    return nd, nw, nb
end
timing_opts(nw, nb) = TimingOpts(window=:integrated, n_window=nw, n_back=nb, cluster=1.0, post_drag_period=true, kick_post_drag=true, zonal_sp=false)

# ---------------------------------------------------------------------------------------------------------------

function phase_cells()
    mkpath(OUT); cells = mars_cells(VC)
    CSV.write(joinpath(OUT, "cells.csv"), DataFrame(cells)); println("wrote $(length(cells)) cells to $(joinpath(OUT, "cells.csv"))")
end

function phase_pitable(ms::MarsSetup)
    rec = mars_record(ms, epoch_of(1)); h = rec.harm
    CSV.write(joinpath(OUT, "planet_constants.csv"), DataFrame([(planet="Venus", mu_m3s2=h.gm_m3s2, Re_m=h.reference_radius_m, J2=-sqrt(5.0) * h.C[3, 1], J3=-sqrt(7.0) * h.C[4, 1],
        Rp_e_m=rec.planet.Rp_e, Rp_p_m=rec.planet.Rp_p, spin_rad_s=rec.planet.ω[3], beta_kg_m2=rec.beta, scale_height_m=Float64(rec.fit.scale_heights_m[1]),
        rho_ref_kg_m3=rec.planet.ρ_ref, h_ref_m=rec.planet.h_ref, apo_alt_m=rec.apo_alt_m, area_m2=ms.area_m2, wet_mass_kg=ms.wet_mass_kg)]))
    drows = NamedTuple[]
    for ei in eachindex(VC["grid"]["epochs_utc"]), b in rec.body_names
        et = et_of(ei); st = spkezr(uppercase(b), et, "J2000", "NONE", rec.center)[1]; r = norm(st[1:3]) * 1e3; v = norm(st[4:6]) * 1e3
        push!(drows, (planet="Venus", epoch=ei, body=b, r_b_m=r, mu_b_m3s2=MU_BODY[b], n_b_rad_s=v / r))
    end
    CSV.write(joinpath(OUT, "body_distances.csv"), DataFrame(drows))
    l0 = SimulationModel.planet_frame_lpi(ms.planet, et_of(1), SpiceEphemeridesModel()); l1 = SimulationModel.planet_frame_lpi(ms.planet, et_of(1) + 60 * SECONDS_PER_DAY, SpiceEphemeridesModel())
    open(joinpath(OUT, "pole_drift.toml"), "w") do io; println(io, "pole_drift_deg_over_60d = $(rad2deg(acos(clamp(dot(l0[3, :], l1[3, :]), -1.0, 1.0))))"); end
    println("pitable done")
end

# ---------------------------------------------------------------------------------------------------------------

"""Plan section 5 resolution rule; no truth used; self-calibrated seeds."""
function phase_convergence(ms::MarsSetup)
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame); rec(ei) = mars_record(ms, epoch_of(ei)); ctxs = Dict{Tuple{Int, Int}, EstContext}()
    getctx(ei, nd) = get!(ctxs, (ei, nd)) do; make_context(rec(ei), with_settings(SETT0; n_drag_half=nd)); end
    zf = zonal_field(ms.harm); tf = tesseral_field(ms.harm); g = guess_alt(); nw0, nb0 = Int.(VC["estimator_numerics"]["v5_window_default"]); nw1, nb1 = Int.(VC["estimator_numerics"]["v5_window_refined"])
    rows = NamedTuple[]; vrows = NamedTuple[]
    for c in eachrow(cells)
        ei = Int(c.epoch_idx); inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = et_of(ei)
        for nd in Int.(VC["estimator_numerics"]["drag_nodes_candidates"])
            s_c = with_settings(SETT0; n_drag_half=nd); s_r = with_settings(SETT0; n_drag_half=4nd); cc = getctx(ei, nd); cr = getctx(ei, 4nd)
            seed = calibrate_seed3(cr, s_r, OPTS, TAB, zf, ZO, inc, raan, om, et0; seed0=g, cal_kw()...)
            pc, _ = predict_history3(cc, s_c, OPTS, TAB, zf, ZO, seed, inc, raan, om, et0); pr, _ = predict_history3(cr, s_r, OPTS, TAB, zf, ZO, seed, inc, raan, om, et0)
            n = min(length(pc.alt_m), length(pr.alt_m))
            push!(rows, (tag=c.tag, n_drag=nd, max_alt_diff_m=n == 0 ? NaN : maximum(abs.(pc.alt_m[1:n] .- pr.alt_m[1:n])), n_c=length(pc.alt_m), n_r=length(pr.alt_m)))
        end
        S = SETT0; ctx = getctx(ei, Int(SETT0.n_drag_half))   # the v5 window rule is tested at the default drag nodes
        seed = calibrate_seed5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, timing_opts(nw0, nb0), inc, raan, om, et0; seed0=g, cal_kw()...)
        p0, _ = predict_history5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, timing_opts(nw0, nb0), seed, inc, raan, om, et0)
        p1, _ = predict_history5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, timing_opts(nw1, nb1), seed, inc, raan, om, et0)
        n = min(length(p0.alt_m), length(p1.alt_m))
        push!(vrows, (tag=c.tag, max_alt_diff_m=n == 0 ? NaN : maximum(abs.(p0.alt_m[1:n] .- p1.alt_m[1:n])), max_dt_diff_s=n == 0 ? NaN : maximum(abs.(p0.t_s[1:n] .- p1.t_s[1:n])), n_0=length(p0.alt_m), n_1=length(p1.alt_m)))
    end
    D = DataFrame(rows); CSV.write(joinpath(OUT, "convergence_drag_nodes.csv"), D); V = DataFrame(vrows); CSV.write(joinpath(OUT, "convergence_v5_window.csv"), V)
    tol = Float64(VC["estimator_numerics"]["convergence_tolerance_m"]); chosen = 0
    summ = combine(groupby(D, :n_drag), :max_alt_diff_m => (x -> maximum(filter(isfinite, x); init=0.0)) => :worst_m, :max_alt_diff_m => (x -> median(filter(isfinite, x))) => :median_m,
        [:n_c, :n_r] => ((a, b) -> count(a .!= b)) => :n_cells_with_different_pass_count)
    println(summ)
    for r in eachrow(summ); r.worst_m < tol && chosen == 0 && (chosen = r.n_drag); end
    chosen == 0 && (chosen = maximum(Int.(VC["estimator_numerics"]["drag_nodes_candidates"])); println("WARNING: no drag-node candidate passed; using the largest"))
    wa = maximum(filter(isfinite, V.max_alt_diff_m); init=0.0); wt = maximum(filter(isfinite, V.max_dt_diff_s); init=0.0)
    ok = wa < Float64(VC["estimator_numerics"]["v5_altitude_tolerance_m"]) && wt < Float64(VC["estimator_numerics"]["v5_passtime_tolerance_s"])
    nw, nb = ok ? (nw0, nb0) : (nw1, nb1)
    println("v5 window: worst altitude change $wa m, worst pass-time change $wt s -> ", ok ? "default ($nw0, $nb0) accepted" : "refined ($nw1, $nb1) used")
    open(NUMERICS_PATH, "w") do io
        println(io, "n_drag_half = $chosen\nn_window = $nw\nn_back = $nb\nworst_drag_node_diff_m = $(maximum(filter(isfinite, D[D.n_drag .== Int(SETT0.n_drag_half), :max_alt_diff_m]); init=0.0))")
        println(io, "v5_window_worst_alt_diff_m = $wa\nv5_window_worst_dt_diff_s = $wt\nrule = \"plan section 5; no truth used\"")
    end
end

"""Addendum A1 (no truth): is the refined window (192, 16) itself converged? Compares it with (384, 32) on all cells; reports only, writes convergence_v5_window_refined.csv."""
function phase_convergence_refined(ms::MarsSetup)
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame); zf = zonal_field(ms.harm); tf = tesseral_field(ms.harm); g = guess_alt(); S = SETT0; ctxs = Dict{Int, EstContext}()
    nwa, nba = Int.(VC["estimator_numerics"]["v5_window_refined"]); rows = NamedTuple[]
    for c in eachrow(cells)
        ei = Int(c.epoch_idx); ctx = get!(ctxs, ei) do; make_context(mars_record(ms, epoch_of(ei)), S); end
        inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = et_of(ei)
        seed = calibrate_seed5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, timing_opts(nwa, nba), inc, raan, om, et0; seed0=g, cal_kw()...)
        p0, _ = predict_history5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, timing_opts(nwa, nba), seed, inc, raan, om, et0)
        p1, _ = predict_history5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, timing_opts(2nwa, 2nba), seed, inc, raan, om, et0)
        n = min(length(p0.alt_m), length(p1.alt_m))
        push!(rows, (tag=c.tag, max_alt_diff_m=n == 0 ? NaN : maximum(abs.(p0.alt_m[1:n] .- p1.alt_m[1:n])), max_dt_diff_s=n == 0 ? NaN : maximum(abs.(p0.t_s[1:n] .- p1.t_s[1:n])), n_0=length(p0.alt_m), n_1=length(p1.alt_m)))
    end
    V = DataFrame(rows); CSV.write(joinpath(OUT, "convergence_v5_window_refined.csv"), V)
    println("refined (192,16) vs (384,32): worst altitude change $(maximum(V.max_alt_diff_m)) m, worst pass-time change $(maximum(V.max_dt_diff_s)) s; cells with different pass count: $(count(V.n_0 .!= V.n_1))")
end

# ---------------------------------------------------------------------------------------------------------------

mk_seeds(cells, out, hd) = [(; config_path=abspath(get(CLI, "config", VENUS_CONFIG_PATH)), out, tag=r.tag, inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch_idx), horizon_days=hd) for r in eachrow(cells)]

function report(res)
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1))")
    foreach(f -> println("FAILED $(f.seed): $(first(sprint(showerror, f.error), 400))"), res.failed)
end

function phase_pilot()
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame)
    only = haskey(CLI, "only") ? split(CLI["only"], ",") : [cells.tag[1]]; cells = cells[in.(cells.tag, Ref(only)), :]
    for r in mk_seeds(cells, joinpath(OUT, "pilot"), Float64(VC["grid"]["horizon_days"]))
        t0 = time(); row = venus_cell_task(r); println("PILOT ", row, "\nwall total ", round(time() - t0; digits=1), " s (not a timing result)")
    end
end

function phase_truth()
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame)
    haskey(CLI, "only") && (cells = cells[in.(cells.tag, Ref(split(CLI["only"], ","))), :])
    workers = haskey(CLI, "workers") ? parse(Int, CLI["workers"]) : Int(VC["grid"]["workers"]); hd = Float64(VC["grid"]["horizon_days"])
    out = joinpath(OUT, "truth"); mkpath(out)
    res = run_process_campaign(venus_cell_task, mk_seeds(cells, out, hd), workers; mission_time_s=hd * SECONDS_PER_DAY, worker_files=[joinpath(@__DIR__, "venus_common.jl")])
    report(res); write_rows_m(joinpath(OUT, "truth_run_summary.csv"), [s.value for s in res.successful])
    open(joinpath(OUT, "truth_mc.toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_elapsed_s = $(res.elapsed_s)")
    end
end

# ---------------------------------------------------------------------------------------------------------------

function phase_predict(ms::MarsSetup)
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame); nd, nw, nb = numerics(); S = with_settings(SETT0; n_drag_half=nd); TO = timing_opts(nw, nb)
    zf = zonal_field(ms.harm); tf = tesseral_field(ms.harm); g = guess_alt(); ctxs = Dict{Int, EstContext}()
    for d in ("pred_v3", "pred_v3_selfcal", "pred_v5", "pred_v5_selfcal"); mkpath(joinpath(OUT, d)); end
    println("estimator numerics: drag nodes per half $nd, v5 window ($nw, $nb); v3 zonal budget $(ZO.budget_km) km; v5 timing $(TO)")
    rows3 = NamedTuple[]; rows5 = NamedTuple[]; t0 = time()
    for c in eachrow(cells)
        ei = Int(c.epoch_idx); ctx = get!(ctxs, ei) do; make_context(mars_record(ms, epoch_of(ei)), S); end
        calf = joinpath(OUT, "truth", c.tag, "calibration.csv"); isfile(calf) || continue
        seed_t = CSV.read(calf, DataFrame).seed_alt_m[end]; inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = et_of(ei)
        p1, i1 = predict_history3(ctx, S, OPTS, TAB, zf, ZO, seed_t, inc, raan, om, et0)
        s3 = calibrate_seed3(ctx, S, OPTS, TAB, zf, ZO, inc, raan, om, et0; seed0=g, cal_kw()...); p2, _ = predict_history3(ctx, S, OPTS, TAB, zf, ZO, s3, inc, raan, om, et0)
        CSV.write(joinpath(OUT, "pred_v3", c.tag * ".csv"), pred_df(p1)); CSV.write(joinpath(OUT, "pred_v3_selfcal", c.tag * ".csv"), pred_df(p2))
        push!(rows3, (tag=c.tag, seed_truth_m=seed_t, seed_selfcal_m=s3, n_pred=length(p1.t_s), decayed_pred=p1.decayed, t_decay_pred_days=p1.t_decay_s / SECONDS_PER_DAY, ended=String(p1.ended),
            n_pred_selfcal=length(p2.t_s), decayed_pred_selfcal=p2.decayed, t_decay_pred_selfcal_days=p2.t_decay_s / SECONDS_PER_DAY, carried=join(i1[1], ";")))
        q1, j1 = predict_history5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, TO, seed_t, inc, raan, om, et0)
        s5 = calibrate_seed5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, TO, inc, raan, om, et0; seed0=g, cal_kw()...); q2, _ = predict_history5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, TO, s5, inc, raan, om, et0)
        CSV.write(joinpath(OUT, "pred_v5", c.tag * ".csv"), pred_df(q1)); CSV.write(joinpath(OUT, "pred_v5_selfcal", c.tag * ".csv"), pred_df(q2))
        push!(rows5, (tag=c.tag, seed_truth_m=seed_t, seed_selfcal_m=s5, n_pred=length(q1.t_s), decayed_pred=q1.decayed, t_decay_pred_days=q1.t_decay_s / SECONDS_PER_DAY, ended=String(q1.ended),
            n_pred_selfcal=length(q2.t_s), decayed_pred_selfcal=q2.decayed, t_decay_pred_selfcal_days=q2.t_decay_s / SECONDS_PER_DAY, carried=join(j1[1], ";"), DT_km=j1[3], tesseral=j1[4]))
    end
    CSV.write(joinpath(OUT, "pred_cells_v3.csv"), DataFrame(rows3)); CSV.write(joinpath(OUT, "pred_cells_v5.csv"), DataFrame(rows5))
    println("predict done: $(length(rows3)) cells in $(round(time() - t0; digits=1)) s (includes compilation and SPICE; not a timing result)")
end

# ---------------------------------------------------------------------------------------------------------------

function main()
    phase = get(CLI, "phase", "cells"); mkpath(OUT); ms = venus_setup(VC)
    println("phase=$phase out=$OUT mu=$(ms.harm.gm_m3s2) Re=$(ms.harm.reference_radius_m) beta=$(ms.beta) area=$(ms.area_m2) wet=$(ms.wet_mass_kg) H=$(ms.planet.H) rho_ref=$(ms.planet.ρ_ref) h_ref=$(ms.planet.h_ref)")
    phase == "cells" && phase_cells()
    phase == "pitable" && phase_pitable(ms)
    phase == "convergence" && phase_convergence(ms)
    phase == "convergence_refined" && phase_convergence_refined(ms)
    phase == "pilot" && phase_pilot()
    phase == "truth" && phase_truth()
    phase == "predict" && phase_predict(ms)
end
main()
