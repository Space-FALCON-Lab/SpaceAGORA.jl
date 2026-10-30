# Option 3 at Mars: a NUMERICAL predictor N(field) in place of the closed-form estimator. Plan and thresholds:
# tasks/2026-10-15_periapsis-region-border-targeting/numerical_predictor_plan.md (written before the runs); results: numerical_predictor_results.md.
# N(field) is the Mars truth engine itself (mars_common.jl: SPICE Mars, Sun, exponential atmosphere + constant-beta drag, truth tolerances), passive, restarted from the
# initial or the measured osculating state, with Mars50c harmonics of a chosen degree x order. No engine code is re-implemented; no src/ or examples/ change.
# One Julia job at a time, memory capped:
#   systemd-run --user --scope -p MemoryMax=16G julia --project=. run_numerical_predictor.jl --phase=probe|ladder|timing|verify|closedloop
#        --results-root=DIR --config=MARS_TOML [--fields=2x0,3x0,...] [--tol=T1|T2] [--workers=N] [--only=tag,tag] [--band=fixed|perpass] [--field=20x20]
#   probe       one cell, serial: wall time of N for a few fields (scale of the campaign; not a table number)
#   ladder      history predictions of N for the fields on every cell of <root>/mars/cells.csv -> pred_num<DxO>[L]/, pred_num<DxO>[L]_selfcal/, pred_cells_num<DxO>[L].csv
#               (process route via run_monte_carlo, route asserted :process; wall times contended)
#   timing      serial warm wall time per 60-day candidate for each field and tolerance on a fixed cell subset (run only with the machine otherwise idle) -> num_timing.csv
#   verify      restart consistency: N(field) from the MEASURED apoapsis state, passive, predicted next periapsis vs truth next periapsis -> num_verify_<tag>.csv
#   closedloop  Arm B controller (near_passive) with N(field) as the finite-thrust predictor on <root>/mars/<cells-file> -> <root>/mars/closedloop_num<DxO>[_perpass]/
include(joinpath(@__DIR__, "run_mars_closedloop.jl"))   # mars_common.jl, PredictorPolicy, mars_corridor_config, mars_closedloop_args (its cl_main does not run when included)
using Statistics

# ASSUMED looser integration setting for the second tolerance of the ladder (plan section 2): not a truth setting.
const LOOSE_TOL = (reltol=1.0e-6, abstol=1.0e-8, dt_max_orbit_s=120.0, dt_max_atmosphere_s=5.0)

parse_field(s::AbstractString) = (m = match(r"^(\d+)x(\d+)$", s); m === nothing ? throw(ArgumentError("field must be DxO, got $s")) : (parse(Int, m.captures[1]), parse(Int, m.captures[2])))
num_name(deg, ord, tol) = "num$(deg)x$(ord)" * (tol == "T2" ? "L" : tol == "T0" ? "T" : "")
# ASSUMED tighter setting T0 (reltol 1e-11, abstol 1e-13, truth step caps): the reference for the integrator error of the truth settings T1 (plan C1c).
const TIGHT_TOL = (reltol=1.0e-11, abstol=1.0e-13)

"""MarsSetup of the numerical predictor: the truth setup with its harmonics swapped for degree x order, and (tol == "T2") the loose integration setting."""
function num_setup(ms::MarsSetup, deg::Int, ord::Int, tol::AbstractString)
    mc = deepcopy(ms.mc)
    if tol == "T2"
        ig = mc["integration"]; ig["reltol"] = LOOSE_TOL.reltol; ig["abstol"] = LOOSE_TOL.abstol
        ig["dt_max_orbit_s"] = LOOSE_TOL.dt_max_orbit_s; ig["dt_max_atmosphere_s"] = LOOSE_TOL.dt_max_atmosphere_s
    elseif tol == "T0"
        ig = mc["integration"]; ig["reltol"] = TIGHT_TOL.reltol; ig["abstol"] = TIGHT_TOL.abstol
    elseif tol != "T1"
        throw(ArgumentError("tol must be T0, T1 or T2"))
    end
    harm = GravitationalHarmonicsModel(deg, ord, joinpath(REPO_ROOT, String(mc["gravity"]["coefficient_file"])), ms.planet)
    return MarsSetup(mc, ms.planet, harm, ms.fit, ms.beta, ms.area_m2, ms.wet_mass_kg, ms.dry_mass_kg, ms.apo_alt_m)
end

pred_frame(rows) = DataFrame(k=1:length(rows), t_s=[r.time_s for r in rows], alt_m=[r.alt_m for r in rows], radius_m=[r.radius_m for r in rows], lat_rad=[r.lat_rad for r in rows],
    vrel_mps=[r.vrel_mps for r in rows], rho=[r.rho for r in rows], q_w_cm2=[r.q_w_cm2 for r in rows], a_m=[r.a_m for r in rows], e=[r.e for r in rows], i_rad=[r.i_rad for r in rows],
    raan_rad=[r.raan_rad for r in rows], omega_rad=[r.omega_rad for r in rows], h_qhigh_m=[r.h_qhigh_m for r in rows], h_qlow_m=[r.h_qlow_m for r in rows])

"""One ladder task: N(field) passive history from the truth-calibrated seed and from the self-calibrated seed (the predictor's own peak-Q calibration). Never throws."""
function num_cell_task(seed)
    t0 = time()
    try
        mc = load_mars(seed.config_path); ms = mars_setup(mc); msN = num_setup(ms, seed.deg, seed.ord, seed.tol)
        epoch = parse_epoch(mc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms.planet, epoch)
        dir = joinpath(seed.out, "scratch", seed.name, seed.tag); mkpath(dir)
        seed_t = CSV.read(joinpath(seed.calib_dir, "calibration.csv"), DataFrame).seed_alt_m[end]
        horizon = Float64(seed.horizon_days) * SECONDS_PER_DAY
        r1 = mars_run(msN, rot, epoch, horizon, dir; seed_alt_m=seed_t, inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega)
        t_cal = time()
        seed_s, _ = mars_calibrate(msN, rot, epoch, dir; inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega); cal_wall = time() - t_cal
        r2 = mars_run(msN, rot, epoch, horizon, dir; seed_alt_m=seed_s, inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega)
        mkpath(joinpath(seed.out, "pred_$(seed.name)")); mkpath(joinpath(seed.out, "pred_$(seed.name)_selfcal"))
        CSV.write(joinpath(seed.out, "pred_$(seed.name)", seed.tag * ".csv"), pred_frame(r1.rows)); CSV.write(joinpath(seed.out, "pred_$(seed.name)_selfcal", seed.tag * ".csv"), pred_frame(r2.rows))
        td(r, res) = isempty(r) ? NaN : r[end].time_s / SECONDS_PER_DAY
        return (; tag=seed.tag, name=seed.name, status="ok", seed_truth_m=seed_t, seed_selfcal_m=seed_s, n_pred=length(r1.rows), decayed_pred=r1.retcode != "Success",
            t_decay_pred_days=r1.retcode != "Success" ? td(r1.rows, r1) : NaN, ended=r1.retcode, n_pred_selfcal=length(r2.rows), decayed_pred_selfcal=r2.retcode != "Success",
            t_decay_pred_selfcal_days=r2.retcode != "Success" ? td(r2.rows, r2) : NaN, n_pred_v1=0, decayed_pred_v1=false, t_decay_pred_v1_days=NaN,
            wall_s_first_call=r1.wall_s, wall_s_selfcal_run=r2.wall_s, wall_s_calibration=cal_wall, wall_s_task=time() - t0, error="")
    catch err
        return (; tag=seed.tag, name=seed.name, status="run_error", seed_truth_m=NaN, seed_selfcal_m=NaN, n_pred=0, decayed_pred=false, t_decay_pred_days=NaN, ended="",
            n_pred_selfcal=0, decayed_pred_selfcal=false, t_decay_pred_selfcal_days=NaN, n_pred_v1=0, decayed_pred_v1=false, t_decay_pred_v1_days=NaN,
            wall_s_first_call=NaN, wall_s_selfcal_run=NaN, wall_s_calibration=NaN, wall_s_task=time() - t0, error=replace(first(sprint(showerror, err), 300), "\n" => " "))
    end
end

# ---- finite-thrust next-periapsis predictor (Arm A pattern) with the numerical field ------------------------------------------------------

"""Next geodetic periapsis altitude and relative speed after a tangential `dv` at the next apoapsis, from the measured J2000 state (r0, v0), mass and elapsed time:
the truth engine with the field of `msN`, the example's fixed-burn guidance and a thruster, terminated at the next periapsis (the Arm A predictor of
run_baseline_closedloop.jl with the Mars setup)."""
function predict_numerical(msN::MarsSetup, epoch, r0, v0, mass_kg, dv, elapsed_s, dir)
    mu = msN.harm.gm_m3s2; oe = _corridor_element_summary(r0, v0, mu)
    oe.a > 0.0 && 0.0 <= oe.e < 1.0 || throw(ArgumentError("predictor state is not elliptic"))
    horizon = 1.2 * 2π * sqrt(oe.a^3 / mu)
    ep = _shift_initial_time(epoch, elapsed_s)
    fixed = CorridorFixedBurnGuidance(dv, 1, false)
    v = msN.mc["vehicle"]
    thruster = BaseThrusterModel([Float64(v["thrust_n"])], [0.0], [0.0], [-1.0], [-1.0], [Float64(v["specific_impulse_s"])])
    args = mars_args(msN, ep, horizon, dir; state=(SVector{3, Float64}(r0), SVector{3, Float64}(v0)), mass_kg=mass_kg, controller=fixed, thruster=thruster)
    scheduled = Ref(false); res = Ref((NaN, NaN)); geodetic_of = _geodetic_of(args)
    env = args.environment_model; planet = env.planet; et0 = SimulationModel.ephemerides_time_seconds(args.initial_time, env.ephemerides_model)
    function record!(integrator)
        if scheduled[]
            pos = integrator.u.sc[1].pos; vel = integrator.u.sc[1].vel
            lpi = SimulationModel.planet_frame_lpi(planet, et0 + Float64(integrator.t), env.ephemerides_model)
            res[] = (geodetic_of(Float64(integrator.t), pos)[1], norm(lpi * vel - cross(planet.ω, lpi * pos)))
            terminate!(integrator)
        end
        return nothing
    end
    peri = ContinuousCallback((u, t, integrator) -> dot(u.sc[1].pos, u.sc[1].vel), record!, nothing; save_positions=(false, false))
    sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=(_fixed_burn_schedule_callback(fixed, scheduled), peri))
    isfinite(res[][1]) && return (alt=res[][1], vrel=res[][2], impact=false)
    # The candidate burn lowered the orbit so far that the engine ended the run at the impact event before the next periapsis: report the altitude and relative speed at the
    # end of the run (at or below the impact altitude), a monotone lower bound for the controller's search, and flag it (counted as `predictor_impacts`).
    pos = sol.u[end].sc[1].pos; vel = sol.u[end].sc[1].vel; t_end = Float64(sol.t[end])
    lpi = SimulationModel.planet_frame_lpi(planet, et0 + t_end, env.ephemerides_model)
    (string(sol.retcode) == "Success" || !all(isfinite, pos)) && throw(ErrorException("numerical predictor did not reach the next periapsis (ret $(sol.retcode))"))
    return (alt=geodetic_of(t_end, pos)[1], vrel=norm(lpi * vel - cross(planet.ω, lpi * pos)), impact=true)
end

"""Closed-loop `make_predictor` for N: memoized per decision (the controller's search repeats dv values); counts engine calls. With `band == :perpass`, the corridor edges
and target are recomputed at every decision from the predicted relative speed at the next periapsis (study-local variant, plan C5)."""
function num_maker(msN::MarsSetup, epoch, dir, cfg, wrapper_ref, counter::Ref{Int}, impacts::Ref{Int}, band::Symbol)
    return (p, r, v, mass, t) -> begin
        memo = Dict{Float64, Any}()
        get_pred(dv) = get!(memo, Float64(dv)) do
            counter[] += 1; pr = predict_numerical(msN, epoch, r, v, mass, Float64(dv), t, dir); pr.impact && (impacts[] += 1); pr
        end
        if band === :perpass
            vrel = get_pred(0.0).vrel; b = msN.mc["band"]; w2m = Float64(b["w_cm2_to_w_m2"])
            hq = altitude_of_density_fast(msN.fit, 2 * Float64(b["q_high_w_cm2"]) * w2m / vrel^3); hl = altitude_of_density_fast(msN.fit, 2 * Float64(b["q_low_w_cm2"]) * w2m / vrel^3)
            cfg.values["mission"]["corridor_min_altitude_m"] = hq; cfg.values["mission"]["corridor_max_altitude_m"] = hl
            wrapper_ref[].inner.target_altitude_m = (hq + hl) / 2
        end
        return dv -> get_pred(dv).alt
    end
end

function num_closedloop_task(seed)
    t0 = time()
    try
        mc = load_mars(seed.config_path); ms = mars_setup(mc); msN = num_setup(ms, seed.deg, seed.ord, seed.tolname)
        epoch = parse_epoch(mc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms.planet, epoch)
        calib = CSV.read(joinpath(seed.calib_dir, "calibration.csv"), DataFrame); seed_alt = calib.seed_alt_m[end]
        cfg, hq, hl = mars_corridor_config(ms, calib.vrel_mps[end], Float64(mc["closedloop"]["margin_m"]))
        dir = joinpath(seed.out, seed.tag); mkpath(dir); pdir = joinpath(dir, "pred_scratch"); mkpath(pdir)
        counter = Ref(0); impacts = Ref(0); wref = Ref{Any}(nothing)
        maker = num_maker(msN, epoch, pdir, cfg, wref, counter, impacts, seed.band)
        state0 = mars_initial_state(ms, rot, seed_alt, seed.inc, seed.raan, seed.omega)
        horizon = Float64(seed.horizon_days) * SECONDS_PER_DAY
        args, wrapper = mars_closedloop_args(ms, epoch, rot, horizon, dir, cfg, state0, maker); wref[] = wrapper
        policy = wrapper.inner; obs = NamedTuple[]
        cbs = (_apsis_observation_callback(obs, args; controller=policy), _preapoapsis_guidance_callback(policy))
        tr = time(); sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=cbs); run_wall = time() - tr
        rows = mars_periapsis_table(sol, args, ms, rot)
        write_rows_m(joinpath(dir, "passes_cl.csv"), rows); write_rows_m(joinpath(dir, "decisions.csv"), policy.decisions)
        ve = Float64(mc["vehicle"]["specific_impulse_s"]) * 9.80665; m0 = ms.wet_mass_kg; mf = Float64(sol.u[end].sc[1].mass)
        plans = policy.burn_plans; bs = band_stats(rows, mc); ret = string(sol.retcode); stopped = !isempty(policy.failures)
        return (; tag=seed.tag, stratum=seed.stratum, status="ok", retcode=ret, t_end_days=sol.t[end] / SECONDS_PER_DAY, horizon_days=seed.horizon_days,
            completed=ret == "Success" && !stopped, controller_stopped=stopped, controller_failures=join(policy.failures, " | "),
            total_dv_mps=ve * log(m0 / mf), dv_commanded_mps=sum(p.requested_delta_v_mps for p in plans; init=0.0), burns=length(plans),
            first_burn_days=isempty(plans) ? NaN : plans[1].execution_start_time_s / SECONDS_PER_DAY, passes=length(rows),
            peak_q_w_cm2=bs.passes == 0 ? NaN : bs.q_max, q_min_w_cm2=bs.passes == 0 ? NaN : bs.q_min, exceeds_hard=bs.passes == 0 ? false : bs.exceeds_hard,
            minutes_in_band=bs.passes == 0 ? NaN : bs.time_in_band_days * 1440.0, passes_in_band=bs.passes == 0 ? 0 : bs.passes_in_band,
            min_alt_km=isempty(rows) ? NaN : minimum(r.alt_m for r in rows) / 1e3, corridor_min_alt_m=hq, corridor_max_alt_m=hl,
            predictor_evals=counter[], predictor_impacts=impacts[], decisions=length(policy.decisions), run_wall_s=run_wall, wall_s=time() - t0)
    catch err
        return (; tag=seed.tag, stratum=seed.stratum, status="run_error", retcode="", t_end_days=NaN, horizon_days=seed.horizon_days, completed=false,
            controller_stopped=false, controller_failures=replace(first(sprint(showerror, err), 300), "\n" => " "), total_dv_mps=NaN, dv_commanded_mps=NaN,
            burns=0, first_burn_days=NaN, passes=0, peak_q_w_cm2=NaN, q_min_w_cm2=NaN, exceeds_hard=false, minutes_in_band=NaN, passes_in_band=0,
            min_alt_km=NaN, corridor_min_alt_m=NaN, corridor_max_alt_m=NaN, predictor_evals=0, predictor_impacts=0, decisions=0, run_wall_s=NaN, wall_s=time() - t0)
    end
end

"""Restart consistency (plan C1b): passive closed-loop run of the truth engine with N as the decision predictor; at every decision N predicts the next periapsis altitude
from the measured state, compared with the next truth periapsis of that same run."""
function num_verify(seed)
    mc = load_mars(seed.config_path); ms = mars_setup(mc); msN = num_setup(ms, seed.deg, seed.ord, seed.tolname)
    epoch = parse_epoch(mc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms.planet, epoch)
    calib = CSV.read(joinpath(seed.calib_dir, "calibration.csv"), DataFrame); seed_alt = calib.seed_alt_m[end]
    cfg, hq, hl = mars_corridor_config(ms, calib.vrel_mps[end], Float64(mc["closedloop"]["margin_m"]))
    cfg.values["mission"]["corridor_min_altitude_m"] = 60_000.0; cfg.values["mission"]["corridor_max_altitude_m"] = 400_000.0      # never burns
    cfg.values["comparison"]["target_periapsis_altitude_m"] = 230_000.0; cfg.values["mission"]["reference_altitude_m"] = 230_000.0
    dir = joinpath(seed.out, seed.tag); mkpath(dir); pdir = joinpath(dir, "pred_scratch"); mkpath(pdir)
    log = NamedTuple[]
    maker = (p, r, v, mass, t) -> begin
        walls = @elapsed pr = predict_numerical(msN, epoch, r, v, mass, 0.0, t, pdir)
        push!(log, (; t_s=t, pred_num_m=pr.alt, pred_wall_s=walls)); dv -> pr.alt
    end
    state0 = mars_initial_state(ms, rot, seed_alt, seed.inc, seed.raan, seed.omega)
    args, wrapper = mars_closedloop_args(ms, epoch, rot, Float64(seed.horizon_days) * SECONDS_PER_DAY, dir, cfg, state0, maker)
    policy = wrapper.inner; obs = NamedTuple[]
    sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=(_apsis_observation_callback(obs, args; controller=policy), _preapoapsis_guidance_callback(policy)))
    rows = mars_periapsis_table(sol, args, ms, rot); out = NamedTuple[]
    for (k, l) in enumerate(log)
        k + 1 <= length(rows) || break
        push!(out, merge(l, (; decision=k, truth_alt_k_m=rows[k].alt_m, truth_alt_next_m=rows[k + 1].alt_m, err_m=l.pred_num_m - rows[k + 1].alt_m)))
    end
    write_rows_m(joinpath(dir, "verify_num.csv"), out); return out
end

# ---- driver ------------------------------------------------------------------------------------------------------------------------------

function parse_cli_n()
    kv = Dict{String, String}()
    for a in ARGS
        m = match(r"^--([^=]+)=(.*)$", a); m === nothing || (kv[m.captures[1]] = m.captures[2])
    end
    return kv
end
loadavg() = (try; split(read("/proc/loadavg", String))[1:3]; catch; ["?", "?", "?"]; end)

function num_main()
    cli = parse_cli_n(); phase = get(cli, "phase", "probe")
    cfg_path = abspath(get(cli, "config", MARS_CONFIG_PATH)); mc = load_mars(cfg_path); study = load_study()
    root = haskey(cli, "results-root") ? abspath(cli["results-root"]) : results_root(study); base = mars_results_dir(root, mc)
    cells = CSV.read(joinpath(base, "cells.csv"), DataFrame)
    haskey(cli, "only") && (cells = cells[in.(cells.tag, Ref(split(cli["only"], ","))), :])
    fields = parse_field.(split(get(cli, "fields", "2x0,3x0,8x0,20x0,4x4,8x8,12x12,20x20"), ","))
    tol = get(cli, "tol", "T1"); hd = Float64(mc["grid"]["horizon_days"]); workers = parse(Int, get(cli, "workers", "4"))
    calib_dir(tag) = joinpath(base, "truth", tag)
    if phase == "ladder"
        seeds = [(; config_path=cfg_path, out=base, tag=c.tag, inc=c.inc, raan=c.raan, omega=c.omega, epoch_idx=Int(c.epoch_idx), deg=d, ord=o, tol, name=num_name(d, o, tol),
                  calib_dir=calib_dir(c.tag), horizon_days=hd) for (d, o) in fields for c in eachrow(cells)]
        println("N ladder: $(length(seeds)) tasks ($(length(fields)) fields x $(nrow(cells)) cells), tol $tol, workers $workers; load $(loadavg())")
        res = run_process_campaign(num_cell_task, seeds, workers; mission_time_s=hd * SECONDS_PER_DAY, worker_files=[joinpath(@__DIR__, "run_numerical_predictor.jl")])
        println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1))")
        foreach(f -> println("FAILED $(f.seed.tag): $(first(sprint(showerror, f.error), 400))"), res.failed)
        df = DataFrame([s.value for s in res.successful])
        for (d, o) in fields
            nm = num_name(d, o, tol); sub = df[df.name .== nm, :]
            write_rows_m(joinpath(base, "pred_cells_$nm.csv"), select(sub, Not(:name, :status, :error)))
            println("$nm: $(nrow(sub)) cells, errors $(sum(sub.status .!= "ok")), median task wall $(round(median(sub.wall_s_task); digits=1)) s (contended)")
        end
        write_rows_m(joinpath(base, "num_ladder_summary_$(tol).csv"), df)
        open(joinpath(base, "num_ladder_$(tol)_mc.toml"), "w") do io
            println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_elapsed_s = $(res.elapsed_s)")
        end
    elseif phase == "probe"
        c = first(eachrow(cells)); println("probe cell $(c.tag), load $(loadavg())")
        for (d, o) in fields, tl in split(get(cli, "tols", "T1"), ",")
            seed = (; config_path=cfg_path, out=joinpath(root, "num_probe"), tag=c.tag, inc=c.inc, raan=c.raan, omega=c.omega, epoch_idx=Int(c.epoch_idx), deg=d, ord=o, tol=String(tl),
                    name=num_name(d, o, String(tl)), calib_dir=calib_dir(c.tag), horizon_days=hd)
            r = num_cell_task(seed); println("  $(seed.name): $(r.status) first-call run $(round(r.wall_s_first_call; digits=2)) s, calibration $(round(r.wall_s_calibration; digits=2)) s, passes $(r.n_pred), ended $(r.ended) $(r.error)")
        end
    elseif phase == "timing"
        nsub = parse(Int, get(cli, "ncells", "8")); idx = round.(Int, range(1, nrow(cells), length=min(nsub, nrow(cells)))); sub = cells[unique(idx), :]; reps = parse(Int, get(cli, "reps", "3"))
        ms = mars_setup(mc); rows = NamedTuple[]; la0 = loadavg(); println("timing: $(nrow(sub)) cells, $reps repeats; load average at start $la0")
        for tl in split(get(cli, "tols", "T1,T2"), ","), (d, o) in fields
            msN = num_setup(ms, d, o, String(tl))
            jobs = [(c, parse_epoch(mc["grid"]["epochs_utc"][Int(c.epoch_idx)]), CSV.read(joinpath(calib_dir(c.tag), "calibration.csv"), DataFrame).seed_alt_m[end]) for c in eachrow(sub)]
            dir = mktempdir()
            run1(j) = (rot = mars_rot(ms.planet, j[2]); mars_run(msN, rot, j[2], hd * SECONDS_PER_DAY, dir; seed_alt_m=j[3], inc_deg=j[1].inc, raan_deg=j[1].raan, omega_deg=j[1].omega))
            run1(jobs[1])                                      # warm-up (compilation, SPICE, harmonic setup)
            for j in jobs, rep in 1:reps
                t = @elapsed r = run1(j)
                push!(rows, (; tol=String(tl), field="$(d)x$(o)", tag=j[1].tag, rep, wall_s=t, run_wall_s=r.wall_s, passes=length(r.rows), retcode=r.retcode))
            end
            w = [x.wall_s for x in rows if x.tol == tl && x.field == "$(d)x$(o)"]; println("  $tl $(d)x$(o): median $(round(median(w); digits=3)) s over $(length(w)) runs")
        end
        la1 = loadavg(); df = DataFrame(rows); df.load_start .= la0[1]; df.load_end .= la1[1]; write_rows_m(joinpath(base, "num_timing.csv"), df); println("load average at end $la1")
    elseif phase == "verify"
        tag = get(cli, "tag", cells.tag[1]); c = first(eachrow(cells[cells.tag .== tag, :])); d, o = parse_field(get(cli, "field", "20x20"))
        out = joinpath(base, "num_verify"); mkpath(out)
        seed = (; config_path=cfg_path, out, tag=c.tag, inc=c.inc, raan=c.raan, omega=c.omega, epoch_idx=Int(c.epoch_idx), deg=d, ord=o, tolname=get(cli, "tol", "T1"),
                calib_dir=calib_dir(c.tag), horizon_days=parse(Float64, get(cli, "horizon-days", "6")))
        foreach(println, num_verify(seed))
    elseif phase == "closedloop"
        cl = CSV.read(joinpath(base, get(cli, "cells-file", "closedloop_cells.csv")), DataFrame); d, o = parse_field(get(cli, "field", "20x20")); band = Symbol(get(cli, "band", "fixed"))
        haskey(cli, "only") && (cl = cl[in.(cl.tag, Ref(split(cli["only"], ","))), :])
        chd = parse(Float64, get(cli, "horizon-days", string(Float64(mc["closedloop"]["horizon_days"]))))
        out = joinpath(base, "closedloop_$(num_name(d, o, get(cli, "tol", "T1")))" * (band === :perpass ? "_perpass" : "") * (haskey(cli, "horizon-days") ? "_h$(Int(chd))" : "")); mkpath(out)
        seeds = [(; config_path=cfg_path, out, tag=r.tag, stratum=r.stratum, inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch), deg=d, ord=o, tolname=get(cli, "tol", "T1"),
                  band, calib_dir=calib_dir(r.tag), horizon_days=chd) for r in eachrow(cl)]
        println("N closed loop: $(length(seeds)) cells, field $(d)x$(o), band $band, horizon $chd d, workers $workers; load $(loadavg())")
        res = run_process_campaign(num_closedloop_task, seeds, workers; mission_time_s=chd * SECONDS_PER_DAY, worker_files=[joinpath(@__DIR__, "run_numerical_predictor.jl")])
        println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1))")
        foreach(f -> println("FAILED $(f.seed.tag): $(first(sprint(showerror, f.error), 400))"), res.failed)
        write_rows_m(joinpath(out, "closedloop_rows.csv"), [s.value for s in res.successful])
        open(joinpath(out, "closedloop_mc.toml"), "w") do io
            println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_elapsed_s = $(res.elapsed_s)")
        end
    else
        throw(ArgumentError("unknown phase $phase"))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    num_main()
end
