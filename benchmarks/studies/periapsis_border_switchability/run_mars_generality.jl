# Stage C driver (Julia half): Mars generality of the frozen periapsis estimator.
# One Julia job at a time, memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. run_mars_generality.jl --phase=cells|pitable|pi|convergence|pilot|truth|predict|diag [--results-root=DIR] [--workers=N] [--only=tag,tag] [--variant=NAME]
# Writes under <results root>/mars. The Python half (analyze_mars_generality.py) scores the predictions.
#   cells        cells.csv (the cell list; written BEFORE any truth run)
#   pitable      planet_constants.csv, body_distances.csv (Earth and Mars constants and perturber distances at the epochs; no simulation)
#   convergence  numerical resolution of the estimator's drag quadrature at Mars against its own refined reference (uses no truth)
#   pilot        one cell, serial (truth, with timing)
#   truth        the 128-cell truth sweep on the process route (run_monte_carlo; route asserted :process)
#   predict      frozen estimator v1 / v2 / v2 self-calibrated seed / B1 hand model for every completed-calibration cell
#   diag         truth variants (nosun, nodrag, j2only) for the cells in diag_cells.csv (attribution subset)
#   predict_diag frozen-estimator counterparts of those variants (tides off / drag off / unchanged) for the same cells -> pred_diag/<variant>/<tag>.csv
#   predict_j3   UNFROZEN DIAGNOSTIC: estimator v2 plus the averaged J3 term, all cells -> pred_j3, pred_j3_selfcal, pred_cells_j3.csv
include(joinpath(@__DIR__, "mars_common.jl"))
include(joinpath(@__DIR__, "estimator_v2_j3_diag.jl"))   # UNFROZEN diagnostic only (phases predict_j3, predict_diag); the frozen path never calls it
using Statistics

function parse_cli()
    kv = Dict{String, String}(); flags = Set{String}()
    for a in ARGS
        m = match(r"^--([^=]+)=(.*)$", a); m === nothing ? (startswith(a, "--") && push!(flags, a[3:end])) : (kv[m.captures[1]] = m.captures[2])
    end
    return kv, flags
end
const CLI, FLAGS = parse_cli()
const MC = load_mars(get(CLI, "config", MARS_CONFIG_PATH))
const STUDY = load_study()
const ROOT = haskey(CLI, "results-root") ? abspath(CLI["results-root"]) : results_root(STUDY)
const OUT = mars_results_dir(ROOT, MC)
const CFG2 = TOML.parsefile(joinpath(STUDY_DIR, "estimator_v2_config.toml"))
const OPTS = (v = CFG2["v2"]; V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"])))
const TAB = SPTable(OPTS.n_sp)
const SETT0 = settings_from(CFG2)

pred_df(p::Prediction) = DataFrame(k=1:length(p.t_s), t_s=p.t_s, alt_m=p.alt_m, radius_m=p.radius_m, lat_rad=p.lat_rad, vrel_mps=p.vrel_mps,
    rho=p.rho, q_w_cm2=p.q_w_cm2, a_m=p.a_m, e=p.e, i_rad=p.i_rad, raan_rad=p.raan_rad, omega_rad=p.omega_rad, h_qhigh_m=p.h_qhigh_m, h_qlow_m=p.h_qlow_m)

epoch_of(idx) = parse_epoch(MC["grid"]["epochs_utc"][idx])
et_of(idx) = SimulationModel.ephemerides_time_seconds(epoch_of(idx), SpiceEphemeridesModel())
const SELECTED_NODES_PATH = joinpath(OUT, "estimator_numerics.toml")
function settings_final()
    nodes = isfile(SELECTED_NODES_PATH) ? Int(TOML.parsefile(SELECTED_NODES_PATH)["n_drag_half"]) : SETT0.n_drag_half
    return with_settings(SETT0; n_drag_half=nodes)
end

# ---------------------------------------------------------------------------------------------------------------

function phase_cells()
    mkpath(OUT)
    cells = mars_cells(MC)
    CSV.write(joinpath(OUT, "cells.csv"), DataFrame(cells))
    println("wrote $(length(cells)) cells to $(joinpath(OUT, "cells.csv"))")
end

# ---------------------------------------------------------------------------------------------------------------

function phase_pitable(ms::MarsSetup)
    ph = MC["phobos_deimos"]
    # constants (Earth from its record, Mars from the Mars record)
    srec = load_study(); erec = earth_record(srec); mrec = mars_record(ms, epoch_of(1))
    rows = NamedTuple[]
    for (rec, nm) in ((erec, "Earth"), (mrec, "Mars"))
        h = rec.harm; push!(rows, (planet=nm, mu_m3s2=h.gm_m3s2, Re_m=h.reference_radius_m, J2=-sqrt(5.0) * h.C[3, 1], J3=-sqrt(7.0) * h.C[4, 1],
            Rp_e_m=rec.planet.Rp_e, Rp_p_m=rec.planet.Rp_p, flattening=(rec.planet.Rp_e - rec.planet.Rp_p) / rec.planet.Rp_e, spin_rad_s=rec.planet.ω[3],
            beta_kg_m2=rec.beta, scale_height_m=nm == "Mars" ? Float64(rec.fit.scale_heights_m[1]) : NaN, apo_alt_m=rec.apo_alt_m))
    end
    CSV.write(joinpath(OUT, "planet_constants.csv"), DataFrame(rows))
    # perturber distances and mean motions at the epochs (SPICE), in the J2000 frame
    drows = NamedTuple[]
    for (rec, nm, nep, epfun) in ((erec, "Earth", length(srec["grid"]["epochs_utc"]), i -> epoch_et(srec, i)), (mrec, "Mars", length(MC["grid"]["epochs_utc"]), et_of))
        for ei in 1:nep, b in rec.body_names
            et = epfun(ei); st = spkezr(uppercase(b), et, "J2000", "NONE", rec.center)[1]; r = norm(st[1:3]) * 1e3; v = norm(st[4:6]) * 1e3
            push!(drows, (planet=nm, epoch=ei, body=b, r_b_m=r, mu_b_m3s2=MU_BODY[b], n_b_rad_s=v / r))
        end
    end
    CSV.write(joinpath(OUT, "body_distances.csv"), DataFrame(drows))
    # pole drift over the horizon
    l0 = SimulationModel.planet_frame_lpi(ms.planet, et_of(1), SpiceEphemeridesModel()); l1 = SimulationModel.planet_frame_lpi(ms.planet, et_of(1) + 60 * SECONDS_PER_DAY, SpiceEphemeridesModel())
    drift = rad2deg(acos(clamp(dot(l0[3, :], l1[3, :]), -1.0, 1.0)))
    # Phobos / Deimos validity guard (outer-perturber expansion needs a / r_b << 1) at the Mars initial orbit
    mu = ms.harm.gm_m3s2
    rp0 = ms.planet.Rp_e + Float64(MC["orbit"]["initial_periapsis_altitude_guess_m"]); ra0 = ms.planet.Rp_e + ms.apo_alt_m; a0 = 0.5 * (rp0 + ra0)
    open(joinpath(OUT, "pole_and_moons.toml"), "w") do io
        println(io, "pole_drift_deg_over_60d = $drift\na_initial_m = $a0\nrp_initial_m = $rp0\nra_initial_m = $ra0")
        for (nm, gm, ab) in (("phobos", ph["phobos_gm_m3s2"], ph["phobos_a_m"]), ("deimos", ph["deimos_gm_m3s2"], ph["deimos_a_m"]))
            v_enc = sqrt(mu * (2 / ab - 1 / a0))   # speed where the spacecraft orbit crosses the moon's orbital radius
            println(io, "$(nm)_a_over_rb = $(a0 / ab)\n$(nm)_rp_below_moon_orbit = $(rp0 < ab)\n$(nm)_ra_above_moon_orbit = $(ra0 > ab)\n$(nm)_GM_over_mu = $(gm / mu)\n$(nm)_crossing_speed_mps = $v_enc")
            # acceleration of the moon at distance d relative to the Sun's tide at the spacecraft semi-major axis (order of magnitude)
            sun_tide = MU_BODY["Sun"] * a0 / (norm(spkpos("SUN", et_of(1), "J2000", "NONE", "MARS")[1]) * 1e3)^3
            println(io, "$(nm)_accel_at_1000km_mps2 = $(gm / (1.0e6)^2)\n$(nm)_accel_at_3000km_mps2 = $(gm / (3.0e6)^2)\nsun_tide_accel_at_a_mps2 = $sun_tide")
            # impulse per fly-by at impact parameter b: dv = 2 GM / (b v); effect on the periapsis radius for a tangential-at-apoapsis equivalent lever (dr_p = 4 a^2 v_a dv / mu)
            va = sqrt(mu * (2 / ra0 - 1 / a0)); lev = 4 * a0^2 * va / mu
            for b in (1.0e6, 3.0e6)
                dv = 2 * gm / (b * v_enc); println(io, "$(nm)_flyby_b$(Int(b / 1000))km_dv_mps = $dv\n$(nm)_flyby_b$(Int(b / 1000))km_equiv_rp_shift_m = $(dv * lev)")
            end
        end
        println(io, "kepler_lever_m_per_mps = $(4 * a0^2 * sqrt(mu * (2 / ra0 - 1 / a0)) / mu)")
    end
    println("pitable done; pole drift over 60 d = $drift deg")
end

# ---------------------------------------------------------------------------------------------------------------

function pi_row(rec::PlanetRecord, ctx::EstContext, planet_name, tag, epoch_idx, et0, alt, a, e, rho, rp; bodies_r)
    p = a * (1 - e^2); n = sqrt(ctx.mu / a^3)
    pij2 = ctx.J2 * (ctx.Re / p)^2
    p3 = [ctx.body_mu[k] / ctx.mu * (a / bodies_r[k])^3 for k in eachindex(bodies_r)]
    H = corridor_scale_height(ctx.fit, alt)
    return (planet=planet_name, tag=tag, epoch=epoch_idx, a_m=a, e=e, rp_m=rp, Pi_J2=pij2, Pi_3B=sum(p3), Gamma_3B_J2=sum(p3) / pij2,
        Pi_drag=0.5 * rho * rp * (1 + e) / ctx.beta, Pi_H=H / rp, a_over_rb_max=maximum(a ./ bodies_r), n_rad_s=n)
end

function phase_pi(ms::MarsSetup)
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame); rows = NamedTuple[]
    srec = load_study(); erec = earth_record(srec); ectx = make_context(erec, SETT0)
    # Earth cells: first-pass truth values; epoch distances from SPICE
    labels = CSV.read(joinpath(ROOT, "oracle", "oracle_cell_labels.csv"), DataFrame)
    for c in eachrow(labels)
        dir = c.inc_deg == 50.0 ? joinpath(ROOT, "border", "grid", c.tag[7:end]) : joinpath(ROOT, "border2", "grid", c.tag)
        t1 = CSV.read(joinpath(dir, "passes_3b.csv"), DataFrame)[1, :]
        et0 = epoch_et(srec, Int(c.epoch)); br = [norm(spkpos(uppercase(b), et0, "J2000", "NONE", "EARTH")[1]) * 1e3 for b in erec.body_names]
        push!(rows, pi_row(erec, ectx, "Earth", c.tag, Int(c.epoch), et0, t1.alt_m, t1.a_m, t1.e, t1.rho, t1.radius_m; bodies_r=br))
    end
    for c in eachrow(cells)
        f = joinpath(OUT, "truth", c.tag, "passes_3b.csv"); isfile(f) || continue
        t1 = CSV.read(f, DataFrame)[1, :]; mrec = mars_record(ms, epoch_of(Int(c.epoch_idx))); mctx = make_context(mrec, SETT0)
        et0 = et_of(Int(c.epoch_idx)); br = [norm(spkpos(uppercase(b), et0, "J2000", "NONE", "MARS")[1]) * 1e3 for b in mrec.body_names]
        push!(rows, pi_row(mrec, mctx, "Mars", c.tag, Int(c.epoch_idx), et0, t1.alt_m, t1.a_m, t1.e, t1.rho, t1.radius_m; bodies_r=br))
    end
    CSV.write(joinpath(OUT, "pi_cells.csv"), DataFrame(rows)); println("pi_cells.csv: $(length(rows)) rows")
end

# ---------------------------------------------------------------------------------------------------------------

"""Numerical resolution of the drag quadrature at Mars, estimator against its own 4x refined reference (no truth). Self-calibrated seeds."""
function phase_convergence(ms::MarsSetup)
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame); rows = NamedTuple[]
    ctxs = Dict{Tuple{Int, Int}, EstContext}()
    getctx(ei, nd) = get!(ctxs, (ei, nd)) do; make_context(mars_record(ms, epoch_of(ei)), with_settings(SETT0; n_drag_half=nd)); end
    for nd in Int.(MC["estimator_numerics"]["drag_nodes_candidates"]), c in eachrow(cells)
        ei = Int(c.epoch_idx); s_c = with_settings(SETT0; n_drag_half=nd); s_r = with_settings(SETT0; n_drag_half=4nd)
        cc = getctx(ei, nd); cr = getctx(ei, 4nd); inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = et_of(ei)
        g = Float64(MC["orbit"]["initial_periapsis_altitude_guess_m"])
        seed = calibrate_seed2(cr, s_r, OPTS, TAB, inc, raan, om, et0; seed0=g)     # one seed for the pair (the reference's)
        pc = predict_history2(cc, s_c, OPTS, TAB, seed, inc, raan, om, et0); pr = predict_history2(cr, s_r, OPTS, TAB, seed, inc, raan, om, et0)
        n = min(length(pc.alt_m), length(pr.alt_m))
        push!(rows, (tag=c.tag, n_drag=nd, max_alt_diff_m=n == 0 ? NaN : maximum(abs.(pc.alt_m[1:n] .- pr.alt_m[1:n])), n_c=length(pc.alt_m), n_r=length(pr.alt_m),
            decayed_c=pc.decayed, decayed_r=pr.decayed))
    end
    D = DataFrame(rows); CSV.write(joinpath(OUT, "convergence_drag_nodes.csv"), D)
    tol = Float64(MC["estimator_numerics"]["convergence_tolerance_m"]); chosen = 0
    summ = combine(groupby(D, :n_drag), :max_alt_diff_m => (x -> maximum(filter(isfinite, x); init=0.0)) => :worst_m, :max_alt_diff_m => (x -> median(filter(isfinite, x))) => :median_m,
        [:n_c, :n_r] => ((a, b) -> count(a .!= b)) => :n_cells_with_different_pass_count)
    println(summ)
    for r in eachrow(summ)
        if r.worst_m < tol && chosen == 0; chosen = r.n_drag; end
    end
    chosen == 0 && (chosen = maximum(Int.(MC["estimator_numerics"]["drag_nodes_candidates"])); println("WARNING: no candidate passed; using the largest"))
    open(SELECTED_NODES_PATH, "w") do io; println(io, "n_drag_half = $chosen\nrule = \"smallest candidate with worst altitude difference to a 4x refined reference < $(tol) m (no truth used)\""); end
    println("selected n_drag_half = $chosen")
end

# ---------------------------------------------------------------------------------------------------------------

function mk_seeds(cells, variant, out, hd)
    return [(; config_path=abspath(get(CLI, "config", MARS_CONFIG_PATH)), out, tag=r.tag, inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch_idx), horizon_days=hd, variant) for r in eachrow(cells)]
end

report_m(res) = (println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1))");
    foreach(f -> println("FAILED $(f.seed): $(first(sprint(showerror, f.error), 400))"), res.failed))

function save_mc_m(res, out, name)
    open(joinpath(out, "$(name)_mc.toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_elapsed_s = $(res.elapsed_s)")
    end
end

function phase_pilot()
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame)
    only = haskey(CLI, "only") ? split(CLI["only"], ",") : [cells.tag[1]]
    cells = cells[in.(cells.tag, Ref(only)), :]
    for r in mk_seeds(cells, get(CLI, "variant", "full"), joinpath(OUT, "pilot"), Float64(MC["grid"]["horizon_days"]))
        t0 = time(); row = mars_cell_task(r)
        println("PILOT ", row, "\nwall total ", round(time() - t0; digits=1), " s")
    end
end

function phase_truth()
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame)
    haskey(CLI, "only") && (cells = cells[in.(cells.tag, Ref(split(CLI["only"], ","))), :])
    workers = haskey(CLI, "workers") ? parse(Int, CLI["workers"]) : Int(MC["grid"]["workers"]); hd = Float64(MC["grid"]["horizon_days"])
    out = joinpath(OUT, "truth"); mkpath(out)
    seeds = mk_seeds(cells, "full", out, hd)
    res = run_process_campaign(mars_cell_task, seeds, workers; mission_time_s=hd * SECONDS_PER_DAY, worker_files=[joinpath(@__DIR__, "mars_common.jl")])
    report_m(res); write_rows_m(joinpath(OUT, "truth_run_summary.csv"), [s.value for s in res.successful]); save_mc_m(res, OUT, "truth")
end

function phase_diag()
    cells = CSV.read(joinpath(OUT, "diag_cells.csv"), DataFrame)
    variants = split(get(CLI, "variants", "nosun,nodrag,j2only"), ",")
    workers = haskey(CLI, "workers") ? parse(Int, CLI["workers"]) : Int(MC["grid"]["workers"]); hd = Float64(MC["grid"]["horizon_days"])
    out = joinpath(OUT, "truth")
    seeds = vcat([mk_seeds(cells, String(v), out, hd) for v in variants]...)
    res = run_process_campaign(mars_cell_task, seeds, workers; mission_time_s=hd * SECONDS_PER_DAY, worker_files=[joinpath(@__DIR__, "mars_common.jl")])
    report_m(res); write_rows_m(joinpath(OUT, "diag_run_summary.csv"), [s.value for s in res.successful]); save_mc_m(res, OUT, "diag")
end

# ---------------------------------------------------------------------------------------------------------------

function phase_predict(ms::MarsSetup)
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame); S = settings_final()
    for d in ("pred_v1", "pred", "pred_selfcal", "hand"); mkpath(joinpath(OUT, d)); end
    println("estimator drag nodes per half: $(S.n_drag_half) (default $(SETT0.n_drag_half))")
    rows = NamedTuple[]; ctxs = Dict{Int, EstContext}(); t0 = time()
    for c in eachrow(cells)
        ei = Int(c.epoch_idx); ctx = get!(ctxs, ei) do; make_context(mars_record(ms, epoch_of(ei)), S); end
        calf = joinpath(OUT, "truth", c.tag, "calibration.csv"); isfile(calf) || continue
        seed_t = CSV.read(calf, DataFrame).seed_alt_m[end]
        inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = et_of(ei)
        g = Float64(MC["orbit"]["initial_periapsis_altitude_guess_m"])
        p0 = predict_history(ctx, S, seed_t, inc, raan, om, et0)                                                 # v1
        w1 = @elapsed p1 = predict_history2(ctx, S, OPTS, TAB, seed_t, inc, raan, om, et0)                       # v2, truth seed
        seed_s = calibrate_seed2(ctx, S, OPTS, TAB, inc, raan, om, et0; seed0=g, tol_m=Float64(CFG2["estimator"]["calibration_tol_m"]), max_iter=Int(CFG2["estimator"]["calibration_max_iter"]))
        p2 = predict_history2(ctx, S, OPTS, TAB, seed_s, inc, raan, om, et0)
        CSV.write(joinpath(OUT, "pred_v1", c.tag * ".csv"), pred_df(p0)); CSV.write(joinpath(OUT, "pred", c.tag * ".csv"), pred_df(p1))
        CSV.write(joinpath(OUT, "pred_selfcal", c.tag * ".csv"), pred_df(p2))
        if !isempty(p1.alt_m)
            h = hand_history(ctx, S, seed_t, inc, raan, om, et0, p1.alt_m[1], p1.vrel_mps[1]); CSV.write(joinpath(OUT, "hand", c.tag * ".csv"), DataFrame(h))
        end
        push!(rows, (tag=c.tag, seed_truth_m=seed_t, seed_selfcal_m=seed_s, n_pred=length(p1.t_s), decayed_pred=p1.decayed, t_decay_pred_days=p1.t_decay_s / SECONDS_PER_DAY, ended=String(p1.ended),
            n_pred_selfcal=length(p2.t_s), decayed_pred_selfcal=p2.decayed, t_decay_pred_selfcal_days=p2.t_decay_s / SECONDS_PER_DAY,
            n_pred_v1=length(p0.t_s), decayed_pred_v1=p0.decayed, t_decay_pred_v1_days=p0.t_decay_s / SECONDS_PER_DAY, wall_s_first_call=w1))
    end
    CSV.write(joinpath(OUT, "pred_cells.csv"), DataFrame(rows)); println("predict done: $(length(rows)) cells in $(round(time() - t0; digits=1)) s")
end

function phase_predict_diag(ms::MarsSetup)
    cells = CSV.read(joinpath(OUT, "diag_cells.csv"), DataFrame); S = settings_final(); ctxs = Dict{Int, EstContext}()
    for (variant, sv) in (("nosun", with_settings(S; tides=false)), ("nodrag", with_settings(S; drag=false)), ("j2only", S), ("full", S))
        mkpath(joinpath(OUT, "pred_diag", variant))
        for c in eachrow(cells)
            ei = Int(c.epoch_idx); ctx = get!(ctxs, ei) do; make_context(mars_record(ms, epoch_of(ei)), S); end
            seed_t = CSV.read(joinpath(OUT, "truth", c.tag, "calibration.csv"), DataFrame).seed_alt_m[end]
            p = predict_history2(ctx, sv, OPTS, TAB, seed_t, deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega), et_of(ei))
            CSV.write(joinpath(OUT, "pred_diag", variant, c.tag * ".csv"), pred_df(p))
        end
    end
    println("predict_diag done for $(nrow(cells)) cells")
end

function phase_predict_j3(ms::MarsSetup)
    cells = CSV.read(joinpath(OUT, "cells.csv"), DataFrame); S = settings_final(); J3 = -sqrt(7.0) * ms.harm.C[4, 1]
    for d in ("pred_j3", "pred_j3_selfcal"); mkpath(joinpath(OUT, d)); end
    println("UNFROZEN DIAGNOSTIC: J3 = $J3 (from the truth gravity file)")
    rows = NamedTuple[]; ctxs = Dict{Int, EstContext}(); g = Float64(MC["orbit"]["initial_periapsis_altitude_guess_m"])
    function calibrate_j3(ctx, inc, raan, om, et0)   # the same self-calibration rule as calibrate_seed2, with the J3 predictor
        q_target = sqrt(ctx.q_low * ctx.q_high); seed = g
        for _ in 1:8
            p = predict_history2_j3(ctx, S, OPTS, TAB, J3, seed, inc, raan, om, et0; max_passes=1); isempty(p.t_s) && return NaN
            h_target = altitude_of_density_fast(ctx.fit, 2 * q_target * ctx.w2m / p.vrel_mps[1]^3); abs(h_target - p.alt_m[1]) < 1.0 && return seed; seed += h_target - p.alt_m[1]
        end
        return seed
    end
    for c in eachrow(cells)
        ei = Int(c.epoch_idx); ctx = get!(ctxs, ei) do; make_context(mars_record(ms, epoch_of(ei)), S); end
        calf = joinpath(OUT, "truth", c.tag, "calibration.csv"); isfile(calf) || continue
        seed_t = CSV.read(calf, DataFrame).seed_alt_m[end]; inc, raan, om = deg2rad(c.inc), deg2rad(c.raan), deg2rad(c.omega); et0 = et_of(ei)
        p1 = predict_history2_j3(ctx, S, OPTS, TAB, J3, seed_t, inc, raan, om, et0)
        seed_s = calibrate_j3(ctx, inc, raan, om, et0); p2 = predict_history2_j3(ctx, S, OPTS, TAB, J3, seed_s, inc, raan, om, et0)
        CSV.write(joinpath(OUT, "pred_j3", c.tag * ".csv"), pred_df(p1)); CSV.write(joinpath(OUT, "pred_j3_selfcal", c.tag * ".csv"), pred_df(p2))
        push!(rows, (tag=c.tag, seed_truth_m=seed_t, seed_selfcal_m=seed_s, n_pred=length(p1.t_s), decayed_pred=p1.decayed, t_decay_pred_days=p1.t_decay_s / SECONDS_PER_DAY, ended=String(p1.ended),
            n_pred_selfcal=length(p2.t_s), decayed_pred_selfcal=p2.decayed, t_decay_pred_selfcal_days=p2.t_decay_s / SECONDS_PER_DAY,
            n_pred_v1=0, decayed_pred_v1=false, t_decay_pred_v1_days=NaN, wall_s_first_call=NaN))
    end
    CSV.write(joinpath(OUT, "pred_cells_j3.csv"), DataFrame(rows)); println("predict_j3 done: $(length(rows)) cells")
end

# ---------------------------------------------------------------------------------------------------------------

function main()
    phase = get(CLI, "phase", "cells"); mkpath(OUT)
    ms = mars_setup(MC)
    println("phase=$phase root=$ROOT out=$OUT mu=$(ms.harm.gm_m3s2) Re=$(ms.harm.reference_radius_m) beta=$(ms.beta) area=$(ms.area_m2) wet=$(ms.wet_mass_kg) apo_alt=$(ms.apo_alt_m)")
    phase == "cells" && phase_cells()
    phase == "pitable" && phase_pitable(ms)
    phase == "pi" && phase_pi(ms)
    phase == "convergence" && phase_convergence(ms)
    phase == "pilot" && phase_pilot()
    phase == "truth" && phase_truth()
    phase == "diag" && phase_diag()
    phase == "predict" && phase_predict(ms)
    phase == "predict_diag" && phase_predict_diag(ms)
    phase == "predict_j3" && phase_predict_j3(ms)
end
main()
