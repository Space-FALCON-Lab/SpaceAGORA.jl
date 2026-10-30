# Venus generality with the native Venus-GRAM atmosphere as truth. One Julia job at a time, memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. venus_gram_generality.jl --phase=guess|convergence|convergence_refined|pilot|truth|predict --results-root=DIR [--workers=N] [--only=tag,tag]
# Plan (frozen before any GRAM truth run): tasks/2026-10-15_periapsis-region-border-targeting/venus_gram_plan.md. Writes under <results root>/venus.
#   guess        cells.csv (grid + Pi-rule altitude + per-cell exponential fit to GRAM mean density + fit residuals/H_local), planet_constants.csv, body_distances.csv (GRAM density only; no truth)
#   convergence, convergence_refined   the exponential-Venus resolution rule (no truth) on the GRAM-fit records -> estimator_numerics.toml
#   pilot        one cell, serial (truth)
#   truth        32-cell truth sweep via run_monte_carlo on the process route (route asserted :process)
#   predict      frozen v3 and v5 with the per-cell exponential fit (truth seed and self-calibrated seed) -> pred_* folders and pred_cells_v3/v5.csv
using Statistics, OrdinaryDiffEq
include(joinpath(@__DIR__, "venus_gram_common.jl"))
include(joinpath(@__DIR__, "estimator_v3.jl"))
include(joinpath(@__DIR__, "estimator_v5.jl"))

function parse_cli()
    kv = Dict{String, String}()
    for a in ARGS; m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2]); end
    return kv
end
const CLI = parse_cli()
const VGC = load_venus_gram(abspath(get(CLI, "config", VENUS_GRAM_CONFIG_PATH)))
const ROOT = haskey(CLI, "results-root") ? abspath(CLI["results-root"]) : abspath(get(ENV, "SPACEAGORA_BORDER_RESULTS", joinpath(REPO_ROOT, "output", "periapsis_border_switchability_gram")))
const OUT = mars_results_dir(ROOT, VGC)
const CFG2 = TOML.parsefile(joinpath(STUDY_DIR, "estimator_v2_config.toml"))
const OPTS = (v = CFG2["v2"]; V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"])))
const TAB = SPTable(OPTS.n_sp)
const SETT0 = settings_from(CFG2)
const ZO = ZonalOpts()
const TT = TesseralOpts()
const NUMERICS_PATH = joinpath(OUT, "estimator_numerics.toml")

epoch_of(idx) = parse_epoch(VGC["grid"]["epochs_utc"][idx])
et_of(idx) = SimulationModel.ephemerides_time_seconds(epoch_of(idx), SpiceEphemeridesModel())
cal_kw() = (tol_m=Float64(CFG2["estimator"]["calibration_tol_m"]), max_iter=Int(CFG2["estimator"]["calibration_max_iter"]))
pred_df(p::Prediction) = DataFrame(k=1:length(p.t_s), t_s=p.t_s, alt_m=p.alt_m, radius_m=p.radius_m, lat_rad=p.lat_rad, vrel_mps=p.vrel_mps,
    rho=p.rho, q_w_cm2=p.q_w_cm2, a_m=p.a_m, e=p.e, i_rad=p.i_rad, raan_rad=p.raan_rad, omega_rad=p.omega_rad, h_qhigh_m=p.h_qhigh_m, h_qlow_m=p.h_qlow_m)

function numerics()
    nd = Int(SETT0.n_drag_half); nw, nb = Int.(VGC["estimator_numerics"]["v5_window_default"])
    if isfile(NUMERICS_PATH); t = TOML.parsefile(NUMERICS_PATH); nd = Int(t["n_drag_half"]); nw = Int(t["n_window"]); nb = Int(t["n_back"]); end
    return nd, nw, nb
end
timing_opts(nw, nb) = TimingOpts(window=:integrated, n_window=nw, n_back=nb, cluster=1.0, post_drag_period=true, kick_post_drag=true, zonal_sp=false)

const FITKIND = get(CLI, "fit", "piecewise")   # piecewise (primary) | single (secondary)
FITKIND in ("piecewise", "single") || error("--fit must be piecewise or single")
const NODES = isfile(joinpath(OUT, "fit_nodes.csv")) ? CSV.read(joinpath(OUT, "fit_nodes.csv"), DataFrame) : nothing
"""Per-cell setup: the GRAM-fit estimator record (the fit uses data available at the cell's start only)."""
function cell_ms(ms0::MarsSetup, c)
    temp = Float64(VGC["atmosphere"]["temperature_k"])
    FITKIND == "single" && return with_fit(ms0, exp_fit(Float64(c.fit_rho_g), Float64(c.h_guess_m), Float64(c.fit_H_m), Float64(VGC["atmosphere"]["fit_top_altitude_m"]), temp, "cell fit single"))
    n = NODES[NODES.tag .== c.tag, :]; return with_fit(ms0, pw_fit(n.h_m, n.rho_kg_m3, temp, "cell fit piecewise"))
end

# ---------------------------------------------------------------------------------------------------------------

function phase_guess(ms::MarsSetup)
    mkpath(OUT); cells = mars_cells(VGC); epoch = epoch_of(1); et0 = et_of(1); rot = mars_rot(ms.planet, epoch)
    model = gram_model(VGC, epoch); pargs = gram_pargs(ms.planet, Float64(VGC["atmosphere_gram"]["entry_interface_km"])); rows = NamedTuple[]; nrows = NamedTuple[]
    println("GRAM model core epoch: ", model.core.initial_time)
    for c in cells
        g = pi_rule_altitude(VGC, ms, rot, et0, model, pargs, c.inc, c.raan, c.omega); f = gram_exp_fit(VGC, model, pargs, g.h_g, g.site)
        push!(rows, merge(c, (h_guess_m=g.h_g, site_lat_rad=g.site.lat, site_lon_rad=g.site.lon, site_t_s=g.site.t_p, site_vrel_mps=g.site.vrel, q_target_w_m2=g.q_target_w_m2,
            fit_H_m=f.H, fit_rho_g=f.rho_g, fit_rms_lnrho=f.rms_lnrho, fit_max_abs_lnrho=f.max_abs_lnrho, fit_res_at_hg=f.res_at_hg, hloc_min_m=f.hloc_min, hloc_max_m=f.hloc_max, hloc_at_hg_m=f.hloc_at_hg, fit_n=f.n_window, pw_max_abs_lnrho=f.pw_max_abs_lnrho, pw_rms_lnrho=f.pw_rms_lnrho, pw_n_nodes=length(f.pw_h))))
        for (k, (hh, rr)) in enumerate(zip(f.pw_h, f.pw_rho)); push!(nrows, (tag=c.tag, k=k, h_m=hh, rho_kg_m3=rr)); end
    end
    D = DataFrame(rows); CSV.write(joinpath(OUT, "cells.csv"), D); CSV.write(joinpath(OUT, "fit_nodes.csv"), DataFrame(nrows)); println("wrote $(nrow(D)) cells to $(joinpath(OUT, "cells.csv"))")
    show(stdout, "text/plain", describe(D[:, [:h_guess_m, :fit_H_m, :fit_rms_lnrho, :fit_max_abs_lnrho, :pw_rms_lnrho, :pw_max_abs_lnrho, :hloc_min_m, :hloc_max_m]], :min, :median, :max)); println()
    # planet constants / bodies (as venus pitable; atmosphere columns = median of the per-cell fits)
    rec = mars_record(ms, epoch); h = rec.harm
    CSV.write(joinpath(OUT, "planet_constants.csv"), DataFrame([(planet="Venus", mu_m3s2=h.gm_m3s2, Re_m=h.reference_radius_m, J2=-sqrt(5.0) * h.C[3, 1], J3=-sqrt(7.0) * h.C[4, 1],
        Rp_e_m=ms.planet.Rp_e, Rp_p_m=ms.planet.Rp_p, spin_rad_s=ms.planet.ω[3], beta_kg_m2=ms.beta, scale_height_m=median(D.fit_H_m), rho_ref_kg_m3=median(D.fit_rho_g), h_ref_m=median(D.h_guess_m),
        apo_alt_m=ms.apo_alt_m, area_m2=ms.area_m2, wet_mass_kg=ms.wet_mass_kg)]))
    drows = NamedTuple[]
    for b in rec.body_names
        st = spkezr(uppercase(b), et0, "J2000", "NONE", rec.center)[1]; r = norm(st[1:3]) * 1e3; v = norm(st[4:6]) * 1e3
        push!(drows, (planet="Venus", epoch=1, body=b, r_b_m=r, mu_b_m3s2=MU_BODY[b], n_b_rad_s=v / r))
    end
    CSV.write(joinpath(OUT, "body_distances.csv"), DataFrame(drows)); println("guess done")
end

# ---------------------------------------------------------------------------------------------------------------

function phase_convergence(ms0::MarsSetup)
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame); zf = zonal_field(ms0.harm); tf = tesseral_field(ms0.harm)
    nw0, nb0 = Int.(VGC["estimator_numerics"]["v5_window_default"]); nw1, nb1 = Int.(VGC["estimator_numerics"]["v5_window_refined"]); rows = NamedTuple[]; vrows = NamedTuple[]
    for c in eachrow(cells)
        ei = Int(c.epoch_idx); inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = et_of(ei); ms = cell_ms(ms0, c); rec = mars_record(ms, epoch_of(ei)); g = Float64(c.h_guess_m)
        for nd in Int.(VGC["estimator_numerics"]["drag_nodes_candidates"])
            s_c = with_settings(SETT0; n_drag_half=nd); s_r = with_settings(SETT0; n_drag_half=4nd); cc = make_context(rec, s_c); cr = make_context(rec, s_r)
            seed = calibrate_seed3(cr, s_r, OPTS, TAB, zf, ZO, inc, raan, om, et0; seed0=g, cal_kw()...)
            pc, _ = predict_history3(cc, s_c, OPTS, TAB, zf, ZO, seed, inc, raan, om, et0); pr, _ = predict_history3(cr, s_r, OPTS, TAB, zf, ZO, seed, inc, raan, om, et0)
            n = min(length(pc.alt_m), length(pr.alt_m))
            push!(rows, (tag=c.tag, n_drag=nd, max_alt_diff_m=n == 0 ? NaN : maximum(abs.(pc.alt_m[1:n] .- pr.alt_m[1:n])), n_c=length(pc.alt_m), n_r=length(pr.alt_m)))
        end
        S = SETT0; ctx = make_context(rec, S)
        seed = calibrate_seed5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, timing_opts(nw0, nb0), inc, raan, om, et0; seed0=g, cal_kw()...)
        p0, _ = predict_history5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, timing_opts(nw0, nb0), seed, inc, raan, om, et0)
        p1, _ = predict_history5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, timing_opts(nw1, nb1), seed, inc, raan, om, et0)
        n = min(length(p0.alt_m), length(p1.alt_m))
        push!(vrows, (tag=c.tag, max_alt_diff_m=n == 0 ? NaN : maximum(abs.(p0.alt_m[1:n] .- p1.alt_m[1:n])), max_dt_diff_s=n == 0 ? NaN : maximum(abs.(p0.t_s[1:n] .- p1.t_s[1:n])), n_0=length(p0.alt_m), n_1=length(p1.alt_m)))
    end
    D = DataFrame(rows); CSV.write(joinpath(OUT, "convergence_drag_nodes.csv"), D); V = DataFrame(vrows); CSV.write(joinpath(OUT, "convergence_v5_window.csv"), V)
    tol = Float64(VGC["estimator_numerics"]["convergence_tolerance_m"]); chosen = 0
    summ = combine(groupby(D, :n_drag), :max_alt_diff_m => (x -> maximum(filter(isfinite, x); init=0.0)) => :worst_m, :max_alt_diff_m => (x -> median(filter(isfinite, x))) => :median_m,
        [:n_c, :n_r] => ((a, b) -> count(a .!= b)) => :n_cells_with_different_pass_count)
    println(summ)
    for r in eachrow(summ); r.worst_m < tol && chosen == 0 && (chosen = r.n_drag); end
    chosen == 0 && (chosen = maximum(Int.(VGC["estimator_numerics"]["drag_nodes_candidates"])); println("WARNING: no drag-node candidate passed; using the largest"))
    wa = maximum(filter(isfinite, V.max_alt_diff_m); init=0.0); wt = maximum(filter(isfinite, V.max_dt_diff_s); init=0.0)
    ok = wa < Float64(VGC["estimator_numerics"]["v5_altitude_tolerance_m"]) && wt < Float64(VGC["estimator_numerics"]["v5_passtime_tolerance_s"])
    nw, nb = ok ? (nw0, nb0) : (nw1, nb1)
    println("v5 window: worst altitude change $wa m, worst pass-time change $wt s -> ", ok ? "default ($nw0, $nb0) accepted" : "refined ($nw1, $nb1) used")
    open(NUMERICS_PATH, "w") do io
        println(io, "n_drag_half = $chosen\nn_window = $nw\nn_back = $nb\nworst_drag_node_diff_m = $(maximum(filter(isfinite, D[D.n_drag .== Int(SETT0.n_drag_half), :max_alt_diff_m]); init=0.0))")
        println(io, "v5_window_worst_alt_diff_m = $wa\nv5_window_worst_dt_diff_s = $wt\nrule = \"plan section 5; no truth used\"")
    end
end

function phase_convergence_refined(ms0::MarsSetup)
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame); zf = zonal_field(ms0.harm); tf = tesseral_field(ms0.harm); S = SETT0
    nwa, nba = Int.(VGC["estimator_numerics"]["v5_window_refined"]); rows = NamedTuple[]
    for c in eachrow(cells)
        ei = Int(c.epoch_idx); ctx = make_context(mars_record(cell_ms(ms0, c), epoch_of(ei)), S)
        inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = et_of(ei)
        seed = calibrate_seed5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, timing_opts(nwa, nba), inc, raan, om, et0; seed0=Float64(c.h_guess_m), cal_kw()...)
        p0, _ = predict_history5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, timing_opts(nwa, nba), seed, inc, raan, om, et0)
        p1, _ = predict_history5(ctx, S, OPTS, TAB, zf, ZO, tf, TT, timing_opts(2nwa, 2nba), seed, inc, raan, om, et0)
        n = min(length(p0.alt_m), length(p1.alt_m))
        push!(rows, (tag=c.tag, max_alt_diff_m=n == 0 ? NaN : maximum(abs.(p0.alt_m[1:n] .- p1.alt_m[1:n])), max_dt_diff_s=n == 0 ? NaN : maximum(abs.(p0.t_s[1:n] .- p1.t_s[1:n])), n_0=length(p0.alt_m), n_1=length(p1.alt_m)))
    end
    V = DataFrame(rows); CSV.write(joinpath(OUT, "convergence_v5_window_refined.csv"), V)
    println("refined vs 2x: worst altitude change $(maximum(V.max_alt_diff_m)) m, worst pass-time change $(maximum(V.max_dt_diff_s)) s; cells with different pass count: $(count(V.n_0 .!= V.n_1))")
end

# ---------------------------------------------------------------------------------------------------------------

mk_seeds(cells, out, hd) = [(; config_path=abspath(get(CLI, "config", VENUS_GRAM_CONFIG_PATH)), out, tag=r.tag, inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch_idx), horizon_days=hd,
    h_guess_m=r.h_guess_m, fit_rho_g=r.fit_rho_g, fit_H_m=r.fit_H_m, pw_h=NODES[NODES.tag .== r.tag, :h_m], pw_rho=NODES[NODES.tag .== r.tag, :rho_kg_m3]) for r in eachrow(cells)]

function report(res)
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1))")
    foreach(f -> println("FAILED $(f.seed.tag): $(first(sprint(showerror, f.error), 600))"), res.failed)
end

function phase_pilot()
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame)
    only = haskey(CLI, "only") ? split(CLI["only"], ",") : [cells.tag[1]]; cells = cells[in.(cells.tag, Ref(only)), :]
    for r in mk_seeds(cells, joinpath(OUT, "pilot"), Float64(VGC["grid"]["horizon_days"]))
        t0 = time(); row = venus_gram_cell_task(r); println("PILOT ", row, "\nwall total ", round(time() - t0; digits=1), " s (not a timing result)")
    end
end

function phase_truth()
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame)
    haskey(CLI, "only") && (cells = cells[in.(cells.tag, Ref(split(CLI["only"], ","))), :])
    workers = haskey(CLI, "workers") ? parse(Int, CLI["workers"]) : Int(VGC["grid"]["workers"]); hd = Float64(VGC["grid"]["horizon_days"])
    out = joinpath(OUT, "truth"); mkpath(out)
    res = run_gram_campaign(venus_gram_cell_task, mk_seeds(cells, out, hd), workers; mission_time_s=hd * SECONDS_PER_DAY, worker_files=[joinpath(@__DIR__, "venus_gram_common.jl")])
    report(res); write_rows_m(joinpath(OUT, "truth_run_summary.csv"), [s.value for s in res.successful])
    open(joinpath(OUT, "truth_mc.toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_elapsed_s = $(res.elapsed_s)\nmc_density_family = \"gram\"")
    end
end

# ---------------------------------------------------------------------------------------------------------------

function phase_predict(ms0::MarsSetup)
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame); nd, nw, nb = numerics(); S = with_settings(SETT0; n_drag_half=nd); TO = timing_opts(nw, nb)
    zf = zonal_field(ms0.harm); tf = tesseral_field(ms0.harm)
    for d in ("pred_v3", "pred_v3_selfcal", "pred_v5", "pred_v5_selfcal"); mkpath(joinpath(OUT, d)); end
    println("estimator numerics: drag nodes per half $nd, v5 window ($nw, $nb); v3 zonal budget $(ZO.budget_km) km; v5 timing $(TO)")
    rows3 = NamedTuple[]; rows5 = NamedTuple[]; t0 = time()
    for c in eachrow(cells)
        ei = Int(c.epoch_idx); ctx = make_context(mars_record(cell_ms(ms0, c), epoch_of(ei)), S); g = Float64(c.h_guess_m)
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
    phase = get(CLI, "phase", "guess"); mkpath(OUT); ms = venus_gram_setup(VGC)
    println("phase=$phase out=$OUT mu=$(ms.harm.gm_m3s2) Re=$(ms.harm.reference_radius_m) beta=$(ms.beta) area=$(ms.area_m2) wet=$(ms.wet_mass_kg) density=GRAM mean EI_km=$(VGC["atmosphere_gram"]["entry_interface_km"])")
    phase == "guess" && phase_guess(ms)
    phase == "convergence" && phase_convergence(ms)
    phase == "convergence_refined" && phase_convergence_refined(ms)
    phase == "pilot" && phase_pilot()
    phase == "truth" && phase_truth()
    phase == "predict" && phase_predict(ms)
end
main()
