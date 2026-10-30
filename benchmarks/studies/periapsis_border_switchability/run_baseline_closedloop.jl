# Stage E part 2: closed-loop baseline comparison, two controller arms (both include Sun and Moon in the predicted periapsis):
#   A  "example predictor + third body": the example's finite-thrust engine predictor, built study-locally with SpiceEphemeridesModel +
#      NBodyGravityModel(Sun, Moon) exactly as the truth run (build_args, third_bodies=true).
#   B  "estimator v2": predictor(dv) = estimator_v2.jl re-initialized from the measured osculating state at the decision (impulsive dv).
# Both use select_corridor_maneuver through the example's CorridorGuidancePolicy state, the same corridor, margin, burn model (finite-thrust
# burn plan scheduled by the example's callbacks) and the same truth engine. Cells: <root>/compare_baselines/closedloop_cells.csv.
#   C = the example's own predictor (no third bodies). --cells-file=closedloop_cells_followup.csv|closedloop_cells_all.csv selects another cell list (output dir suffix _followup/_all).
#   julia --project=. run_baseline_closedloop.jl --arm=A|B|C [--smoke] [--cells-file=NAME] [--results-root=DIR] [--workers=N] [--only=tag,tag] [--horizon-days=D]
#   julia --project=. run_baseline_closedloop.jl --verify-tide=TAG [--horizon-days=D] [--results-root=DIR]   (predictor check, serial)
# Writes <root>/compare_baselines/closedloop_<arm>[_smoke]/{<tag>/passes_cl.csv, decisions.csv, closedloop_rows.csv, closedloop_mc.toml}.
# The altitude corridor given to the controller is the heat-rate band at the first-pass relative speed (fixed for the run).
include(joinpath(@__DIR__, "run_static.jl"))
include(joinpath(@__DIR__, "estimator_v2.jl"))

const CL_CONFIG_PATH = joinpath(@__DIR__, "baseline_closedloop_config.toml")

# ---- predictor-injectable policy ---------------------------------------------------------------------------------------------------
# The example's CorridorGuidancePolicy hard-wires its predictor. This wrapper keeps the example's policy object (`inner`) as the state the
# example's callbacks read and write, and replaces only the decision hook so the predictor can be injected. The hook body is the example's
# (AGORA_Earth_Periapsis_Corridor.jl, calcGuidanceEffect!(::CorridorGuidancePolicy)) with the predictor line changed.

mutable struct PredictorPolicy <: SimulationModel.AbstractGuidanceModel
    inner::CorridorGuidancePolicy
    make_predictor::Any     # (p, r, v, mass, t) -> (dv -> next periapsis geodetic altitude, m)
end

function SimulationModel.GuidanceHooks.calcGuidanceEffect!(wrapper::PredictorPolicy, u, p, t::Float64, sat_idx::Int64)
    policy = wrapper.inner
    commands = p.shared_buffers.maneuver_commands
    (sat_idx < 1 || sat_idx > length(commands)) && return nothing
    if !policy.armed
        commands[sat_idx] = SimulationModel.CommandTypes.PropulsiveManeuverCommand(valid=true, source_orbit=policy.periapsis_count + 1)
        return nothing
    end
    r = SVector{3, Float64}(u.sc[sat_idx].pos); v = SVector{3, Float64}(u.sc[sat_idx].vel); mass = Float64(u.sc[sat_idx].mass)
    vehicle = policy.config.values["vehicle"]
    exhaust_velocity = Float64(vehicle["specific_impulse_s"]) * 9.80665
    consumed_at_limit = Float64(vehicle["thrust_n"]) * policy.max_burn_duration_s / exhaust_velocity
    duration_limited_dv = consumed_at_limit < mass ? exhaust_velocity * log(mass / (mass - consumed_at_limit)) : Inf
    wet_mass = Float64(vehicle["dry_mass_kg"]) + Float64(vehicle["propellant_mass_kg"])
    delivered_so_far = mass < wet_mass ? exhaust_velocity * log(wet_mass / mass) : 0.0
    policy.remaining_budget_mps = max(0.0, Float64(policy.config.values["comparison"]["total_delta_v_budget_mps"]) - delivered_so_far)
    max_allowed_dv = min(policy.max_dv_mps, policy.remaining_budget_mps, duration_limited_dv)
    predictor = wrapper.make_predictor(p, r, v, mass, Float64(t))
    mission = policy.config.values["mission"]
    choice = select_corridor_maneuver(predictor;
        lower_m=Float64(mission["corridor_min_altitude_m"]), upper_m=Float64(mission["corridor_max_altitude_m"]),
        margin_m=Float64(mission["safety_margin_m"]), target_m=policy.target_altitude_m, policy=policy.policy,
        max_dv_mps=max_allowed_dv, deadband_m=policy.deadband_m,
        monotone_tolerance_m=Float64(policy.config.values["comparison"]["predictor_monotone_tolerance_m"]))
    policy.last_selection = choice
    target_pass = policy.periapsis_count + 1
    push!(policy.decisions, (request_time_s=Float64(t), source_orbit=target_pass, target_pass=target_pass, policy=String(policy.policy),
        delta_v_mps=choice.delta_v_mps, predicted_periapsis_altitude_m=choice.predicted_periapsis_altitude_m,
        target_periapsis_altitude_m=choice.target_periapsis_altitude_m, trigger_reason=String(choice.reason), feasible=choice.feasible,
        predictor_evaluations=choice.evaluations))
    commands[sat_idx] = SimulationModel.CommandTypes.PropulsiveManeuverCommand(valid=true,
        delta_v_mps=choice.feasible ? abs(choice.delta_v_mps) : 0.0, direction_rad=choice.delta_v_mps >= 0.0 ? 0.0 : π, source_orbit=target_pass)
    policy.armed = false
    return nothing
end

# ---- Arm A predictor: the example's finite-thrust engine predictor + Sun and Moon ----------------------------------------------------

"""The example's `_predict_next_physical_periapsis` with the truth run's third-body setup (SPICE ephemerides, NBodyGravityModel Sun+Moon)."""
function predict_engine_3b(cfg, study, source_args, r0, v0, mass_kg, dv, elapsed_s, dir)
    planet = source_args.environment_model.planet
    oe = _corridor_element_summary(r0, v0, planet.μ)
    oe.a > 0.0 && 0.0 <= oe.e < 1.0 || throw(ArgumentError("predictor state is not elliptic"))
    horizon = 1.2 * 2π * sqrt(oe.a^3 / planet.μ)
    epoch = _shift_initial_time(source_args.initial_time, elapsed_s)   # run epoch + elapsed: Sun and Moon at the true date
    args0 = build_args(cfg, study, epoch, horizon, dir; third_bodies=true, state=(r0, v0), mass_kg=mass_kg)
    fixed = CorridorFixedBurnGuidance(dv, 1, false)
    vh = cfg.values["vehicle"]
    thruster = BaseThrusterModel([Float64(vh["thrust_n"])], [0.0], [0.0], [-1.0], [-1.0], [Float64(vh["specific_impulse_s"])])
    args = SimulationConfiguration(file_paths=args0.file_paths, simulation_settings=args0.simulation_settings,
        mission_configuration=args0.mission_configuration, environment_model=args0.environment_model, dynamics_model=args0.dynamics_model,
        guidance_model=GuidanceModel((fixed,), [1.0e12]), navigation_model=args0.navigation_model,
        control_model=ControlModel((thruster,), [1.0e12]), initial_time=args0.initial_time,
        integration_tolerances=args0.integration_tolerances, solver_config=args0.solver_config)
    scheduled = Ref(false); altitude = Ref(NaN); geodetic_of = _geodetic_of(args)
    function record_next_periapsis!(integrator)
        if scheduled[]
            altitude[] = geodetic_of(Float64(integrator.t), integrator.u.sc[1].pos)[1]
            terminate!(integrator)
        end
        return nothing
    end
    peri = ContinuousCallback((u, t, integrator) -> dot(u.sc[1].pos, u.sc[1].vel), record_next_periapsis!, nothing; save_positions=(false, false))
    run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=(_fixed_burn_schedule_callback(fixed, scheduled), peri))
    isfinite(altitude[]) || throw(ErrorException("3-body engine predictor did not reach the next periapsis"))
    return altitude[]
end

# ---- Arm B predictor: estimator v2 from the measured osculating state --------------------------------------------------------------

const EST_STATE = Ref{Any}(nothing)
function est_state(study)
    if EST_STATE[] === nothing
        c2 = TOML.parsefile(joinpath(@__DIR__, "estimator_v2_config.toml")); s = settings_from(c2); v = c2["v2"]
        o = V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"]))
        EST_STATE[] = (ctx=make_context(study, s), s=s, o=o, tab=SPTable(o.n_sp))
    end
    return EST_STATE[]
end

"""Next geodetic periapsis altitude after a tangential `dv` at the next apoapsis, from the measured osculating state (r0, v0) at ephemeris time `et_now`.
The burn is modeled as the Keplerian change (a, e) of a tangential impulse at the apoapsis of the measured osculating orbit, added to the measured
elements at the measured true anomaly; the estimator then propagates the (mean-element, Sun, Moon, J2, drag) history from that state, so dv = 0 returns the
estimator's own passive prediction from the measured state. The impulsive model ignores the finite-burn duration (Arm A models it)."""
function predict_estimator(study, r0, v0, dv, et_now)
    st = est_state(study); ctx = st.ctx
    oe = _corridor_element_summary(r0, v0, ctx.mu)
    ra = oe.a * (1 + oe.e); va = sqrt(ctx.mu * (2 / ra - 1 / oe.a))
    a2 = 1 / (2 / ra - (va + dv)^2 / ctx.mu); e2 = ra / a2 - 1
    init = (a=oe.a + (a2 - oe.a), e=oe.e + (e2 - oe.e), i=oe.i, raan=oe.raan, omega=oe.omega, nu=oe.nu)
    p = predict_history2(ctx, st.s, st.o, st.tab, 0.0, 0.0, 0.0, 0.0, et_now; max_passes=1, init_override=init)
    isempty(p.alt_m) && throw(ErrorException("estimator predicted no next periapsis (ended=$(p.ended))"))
    return p.alt_m[1]
end

"""Truth args of `cfg` with the policy wrapper and a thruster in place of the passive guidance; returns `(args, wrapper)`."""
function closedloop_args(cfg, study, epoch, horizon_s, dir, policy_sym::Symbol, make_predictor)
    args0 = build_args(cfg, study, epoch, horizon_s, dir; third_bodies=true)
    v = cfg.values["vehicle"]; c = cfg.values["comparison"]
    inner = CorridorGuidancePolicy(cfg, :favorable, policy_sym, Float64(c["max_burn_delta_v_mps"]), Float64(c["max_burn_duration_s"]),
        Float64(c["total_delta_v_budget_mps"]), Float64(c["near_passive_deadband_m"]), Float64(c["target_periapsis_altitude_m"]),
        deg2rad(Float64(c["planning_lead_deg"])), 0, false, nothing, NamedTuple[], NamedTuple[], String[])
    wrapper = PredictorPolicy(inner, make_predictor)
    thruster = BaseThrusterModel([Float64(v["thrust_n"])], [0.0], [0.0], [-1.0], [-1.0], [Float64(v["specific_impulse_s"])])
    args = SimulationConfiguration(
        file_paths=args0.file_paths, simulation_settings=args0.simulation_settings, mission_configuration=args0.mission_configuration,
        environment_model=args0.environment_model, dynamics_model=args0.dynamics_model,
        guidance_model=GuidanceModel((wrapper,), [1.0e12]), navigation_model=args0.navigation_model,
        control_model=ControlModel((thruster,), [1.0e12]), initial_time=args0.initial_time,
        integration_tolerances=args0.integration_tolerances, solver_config=args0.solver_config)
    return args, wrapper
end

"""Cell setup shared by the closed-loop and verification tasks: config with the heat-rate-band corridor, epoch, calibration."""
function cell_setup(seed, cc)
    study = load_study(seed.study_path); nbody_planet()
    calib = CSV.read(joinpath(seed.calib_dir, "calibration.csv"), DataFrame)
    cfg = cell_config(study, seed.inc, seed.raan, seed.omega; seed_alt=calib.seed_alt_m[end])   # the stored accepted seed = the open-loop initial state
    epoch = parse_epoch(study["grid"]["epochs_utc"][seed.epoch_idx])
    h_hi_q, h_lo_q = band_altitudes(fitted_atmosphere_parameters(cfg), calib.vrel_mps[end], study)   # (altitude at Q=0.32, altitude at Q=0.14)
    m = cfg.values["mission"]; c = cfg.values["comparison"]
    m["corridor_min_altitude_m"] = h_hi_q; m["corridor_max_altitude_m"] = h_lo_q
    m["safety_margin_m"] = Float64(cc["margin_m"]); m["reference_altitude_m"] = (h_hi_q + h_lo_q) / 2
    c["target_periapsis_altitude_m"] = (h_hi_q + h_lo_q) / 2
    return (; study, cfg, epoch, h_hi_q, h_lo_q)
end

function predictor_maker(arm::AbstractString, cfg, study, epoch, dir)
    et0 = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel())
    if arm == "A"
        return (p, r, v, mass, t) -> (dv -> predict_engine_3b(cfg, study, p.args, r, v, mass, Float64(dv), t, dir))
    elseif arm == "C"   # the example's OWN predictor (no Sun, no Moon), unmodified
        return (p, r, v, mass, t) -> (dv -> _predict_next_physical_periapsis(cfg, :favorable, p.args, r, v, mass, Float64(dv); elapsed_s=t))
    elseif arm == "B"
        return (p, r, v, mass, t) -> (dv -> predict_estimator(study, r, v, Float64(dv), et0 + t))
    end
    throw(ArgumentError("arm must be A, B or C"))
end

"""One closed-loop cell; plain-data row. Never throws: a failure becomes a row with status = "run_error"."""
function closedloop_task(seed)
    t0 = time()
    try
        cc = TOML.parsefile(CL_CONFIG_PATH)["controller"]
        S = cell_setup(seed, cc)
        dir = joinpath(seed.out, seed.tag); mkpath(dir)
        horizon = Float64(seed.horizon_days) * SECONDS_PER_DAY
        args, wrapper = closedloop_args(S.cfg, S.study, S.epoch, horizon, dir, Symbol(cc["policy"]), predictor_maker(seed.arm, S.cfg, S.study, S.epoch, dir))
        policy = wrapper.inner
        obs = NamedTuple[]   # the example's apsis observer also counts periapses (policy.periapsis_count), as in run_corridor_case
        cbs = (_apsis_observation_callback(obs, args; controller=policy), _preapoapsis_guidance_callback(policy))
        tr = time()
        sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=cbs)
        run_wall = time() - tr
        rows = with_band(periapsis_table(sol, args, S.cfg, S.study), S.cfg, S.study)
        write_rows(joinpath(dir, "passes_cl.csv"), rows)
        write_rows(joinpath(dir, "decisions.csv"), policy.decisions)
        v = S.cfg.values["vehicle"]; ve = Float64(v["specific_impulse_s"]) * 9.80665
        m0 = Float64(v["dry_mass_kg"]) + Float64(v["propellant_mass_kg"]); mf = Float64(sol.u[end].sc[1].mass)
        plans = policy.burn_plans; bs = band_stats(rows, S.study)
        t_end = sol.t[end]; ret = string(sol.retcode)
        stopped = !isempty(policy.failures)
        return (; arm=seed.arm, tag=seed.tag, stratum=seed.stratum, status="ok", retcode=ret, t_end_days=t_end / SECONDS_PER_DAY, horizon_days=seed.horizon_days,
            completed=ret == "Success" && !stopped, controller_stopped=stopped, controller_failures=join(policy.failures, " | "),
            total_dv_mps=ve * log(m0 / mf), dv_commanded_mps=sum(p.requested_delta_v_mps for p in plans; init=0.0), burns=length(plans),
            first_burn_days=isempty(plans) ? NaN : plans[1].execution_start_time_s / SECONDS_PER_DAY, passes=length(rows),
            peak_q_w_cm2=bs.passes == 0 ? NaN : bs.q_max, q_min_w_cm2=bs.passes == 0 ? NaN : bs.q_min, exceeds_hard=bs.passes == 0 ? false : bs.exceeds_hard,
            minutes_in_band=bs.passes == 0 ? NaN : bs.time_in_band_days * 1440.0, passes_in_band=bs.passes == 0 ? 0 : bs.passes_in_band,
            min_alt_km=isempty(rows) ? NaN : minimum(r.alt_m for r in rows) / 1e3, corridor_min_alt_m=S.h_hi_q, corridor_max_alt_m=S.h_lo_q,
            predictor_evals=sum((d.predictor_evaluations for d in policy.decisions); init=0), run_wall_s=run_wall, wall_s=time() - t0)
    catch err
        return (; arm=seed.arm, tag=seed.tag, stratum=seed.stratum, status="run_error", retcode="", t_end_days=NaN, horizon_days=seed.horizon_days, completed=false,
            controller_stopped=false, controller_failures=replace(first(sprint(showerror, err), 300), "\n" => " "), total_dv_mps=NaN, dv_commanded_mps=NaN,
            burns=0, first_burn_days=NaN, passes=0, peak_q_w_cm2=NaN, q_min_w_cm2=NaN, exceeds_hard=false, minutes_in_band=NaN, passes_in_band=0,
            min_alt_km=NaN, corridor_min_alt_m=NaN, corridor_max_alt_m=NaN, predictor_evals=0, run_wall_s=NaN, wall_s=time() - t0)
    end
end

# ---- predictor verification (serial): predicted next periapsis with and without Sun/Moon against the truth pass-to-pass change --------------------

"""Passive truth run (no burns: the corridor is made wide) in which, at every decision, three dv = 0 predictions of the next periapsis are recorded:
A (example predictor, no third bodies), A3 (Arm A, with Sun and Moon), B (estimator v2). Compared with the next truth periapsis altitude."""
function verify_task(seed)
    cc = TOML.parsefile(CL_CONFIG_PATH)["controller"]
    S = cell_setup(seed, cc)
    dir = joinpath(seed.out, seed.tag); mkpath(dir)
    S.cfg.values["mission"]["corridor_min_altitude_m"] = 60_000.0; S.cfg.values["mission"]["corridor_max_altitude_m"] = 400_000.0   # never burns
    S.cfg.values["comparison"]["target_periapsis_altitude_m"] = 230_000.0; S.cfg.values["mission"]["reference_altitude_m"] = 230_000.0
    et0 = SimulationModel.ephemerides_time_seconds(S.epoch, SpiceEphemeridesModel())
    log = NamedTuple[]
    function maker(p, r, v, mass, t)
        function pred(dv)
            if dv == 0.0 && (isempty(log) || log[end].t_s != t)
                y_a = _predict_next_physical_periapsis(S.cfg, :favorable, p.args, r, v, mass, 0.0; elapsed_s=t)
                y_a3 = predict_engine_3b(S.cfg, S.study, p.args, r, v, mass, 0.0, t, dir)
                y_b = predict_estimator(S.study, r, v, 0.0, et0 + t)
                push!(log, (; t_s=t, pred_example_no3b_m=y_a, pred_armA_3b_m=y_a3, pred_armB_est_m=y_b))
            end
            return log[end].pred_armA_3b_m
        end
        return pred
    end
    args, wrapper = closedloop_args(S.cfg, S.study, S.epoch, Float64(seed.horizon_days) * SECONDS_PER_DAY, dir, :near_passive, maker)
    policy = wrapper.inner; obs = NamedTuple[]
    sol = run_simulation(args; return_solution=true, isolate_state=false,
        extra_callbacks=(_apsis_observation_callback(obs, args; controller=policy), _preapoapsis_guidance_callback(policy)))
    rows = periapsis_table(sol, args, S.cfg, S.study)
    out = NamedTuple[]
    for (k, l) in enumerate(log)   # decision k precedes truth periapsis k+1 (decision 1 follows pass 1)
        k + 1 <= length(rows) || break
        a0 = rows[k].alt_m; a1 = rows[k + 1].alt_m
        push!(out, merge(l, (; decision=k, truth_alt_k_m=a0, truth_alt_next_m=a1, truth_shift_m=a1 - a0,
            shift_example_no3b_m=l.pred_example_no3b_m - a0, shift_armA_3b_m=l.pred_armA_3b_m - a0, shift_armB_est_m=l.pred_armB_est_m - a0)))
    end
    write_rows(joinpath(dir, "verify_tide.csv"), out)
    return out
end

function cl_main(args=ARGS)
    smoke = "--smoke" in args; root = nothing; workers = nothing; only = nothing; arm = nothing; hd_arg = nothing; verify = nothing; cells_file = "closedloop_cells.csv"
    for a in args
        startswith(a, "--results-root=") && (root = abspath(split(a, "=", limit=2)[2]))
        startswith(a, "--workers=") && (workers = parse(Int, split(a, "=", limit=2)[2]))
        startswith(a, "--only=") && (only = split(split(a, "=", limit=2)[2], ","))
        startswith(a, "--arm=") && (arm = String(split(a, "=", limit=2)[2]))
        startswith(a, "--horizon-days=") && (hd_arg = parse(Float64, split(a, "=", limit=2)[2]))
        startswith(a, "--cells-file=") && (cells_file = String(split(a, "=", limit=2)[2]))
        startswith(a, "--verify-tide=") && (verify = split(a, "=", limit=2)[2])
    end
    study = load_study(); cl = TOML.parsefile(CL_CONFIG_PATH)
    root = something(root, results_root(study)); base = joinpath(root, "compare_baselines")
    cells = CSV.read(joinpath(base, cells_file), DataFrame)
    mkseed(r, out, hd, arm) = (; study_path=STUDY_CONFIG_PATH, out, arm, tag=r.tag, stratum=r.stratum, inc=r.inc_deg, raan=r.raan_deg, omega=r.omega_deg,
        epoch_idx=Int(r.epoch), calib_dir=joinpath(root, r.calib_dir), horizon_days=hd)
    if verify !== nothing
        r = only_row = first(eachrow(cells[cells.tag .== verify, :]))
        out = joinpath(base, "verify_tide"); mkpath(out)
        res = verify_task(mkseed(r, out, something(hd_arg, 6.0), "verify"))
        foreach(println, res); return
    end
    arm in ("A", "B", "C") || throw(ArgumentError("--arm=A, B or C is required"))
    only === nothing || (cells = cells[in.(cells.tag, Ref(only)), :])
    smoke && only === nothing && (cells = cells[in.(cells.tag, Ref(String.(cl["run"]["smoke_tags"]))), :])
    hd = something(hd_arg, Float64(smoke ? cl["run"]["smoke_horizon_days"] : cl["run"]["horizon_days"]))
    workers = something(workers, Int(smoke ? cl["run"]["smoke_workers"] : cl["run"]["workers"]))
    suffix = cells_file == "closedloop_cells.csv" ? "" : "_" * replace(replace(cells_file, "closedloop_cells_" => ""), ".csv" => "")   # e.g. _followup, _all
    out = joinpath(base, "closedloop_$(arm)$(suffix)" * (smoke ? "_smoke" : "")); mkpath(out)
    seeds = [mkseed(r, out, hd, arm) for r in eachrow(cells)]
    println("closed-loop arm $(arm): $(length(seeds)) cells, horizon $(hd) d, workers $(workers), smoke=$(smoke)")
    res = run_process_campaign(closedloop_task, seeds, workers; mission_time_s=hd * SECONDS_PER_DAY, worker_files=[joinpath(@__DIR__, "run_baseline_closedloop.jl")])
    report(res); write_rows(joinpath(out, "closedloop_rows.csv"), [s.value for s in res.successful]); save_mc(res, out, "closedloop")
end

if abspath(PROGRAM_FILE) == @__FILE__
    cl_main()
end
