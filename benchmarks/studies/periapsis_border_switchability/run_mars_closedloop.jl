# Stage C3 (Mars): closed-loop transfer test, Arm B only (estimator v2 predictor with the Sun, as the Earth Arm B).
# The example's select_corridor_maneuver (policy near_passive) decides through a study-local PredictorPolicy wrapper (a copy of the Earth study's
# wrapper in run_baseline_closedloop.jl: it replaces only the decision hook of the example's CorridorGuidancePolicy so the predictor can be injected);
# the example's callbacks (apsis observer, pre-apoapsis scheduler, finite-thrust burn plans) are reused unchanged. Truth: the Mars engine run of
# mars_common.jl (Sun, J2/J3, exponential atmosphere drag). No src/ or examples/ change.
#   julia --project=. run_mars_closedloop.jl [--smoke] [--cells-file=NAME] [--results-root=DIR] [--config=MARS_TOML] [--workers=N] [--only=tag,tag] [--horizon-days=D] [--verify=TAG]
# Reads <root>/mars/<cells-file> (default closedloop_cells.csv; columns tag, stratum, inc, raan, omega, epoch), <root>/mars/truth/<tag>/calibration.csv.
# Writes <root>/mars/closedloop[_smoke]/{<tag>/passes_cl.csv, decisions.csv, closedloop_rows.csv, closedloop_mc.toml}.
include(joinpath(@__DIR__, "mars_common.jl"))
include(joinpath(@__DIR__, "estimator_v2_j3_diag.jl"))   # UNFROZEN diagnostic predictor, used only with --predictor=j3diag
include(joinpath(@__DIR__, "estimator_v3.jl"))           # estimator v3 (conditional zonal terms), used with --predictor=v3
include(joinpath(@__DIR__, "estimator_v4.jl"))           # estimator v4 (v3 plus the tesseral kick), used with --predictor=v4 (kick forced on: a one-pass history cannot evaluate the carry rule)
include(joinpath(@__DIR__, "estimator_v5.jl"))           # estimator v5 (v4 plus the corrected pass timing), used with --predictor=v5 (kick forced on; a ONE-pass prediction hardly sees the timing correction)

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

# ---- Arm B predictor: estimator v2 from the measured osculating state ------------------------------------------------------------------

const EST_MARS = Ref{Any}(nothing)
"""Frozen estimator context (default settings plus the numerically selected drag nodes) per (process, epoch)."""
function est_state_mars(ms::MarsSetup, epoch_idx::Int, epoch, root; predictor::String="frozen")
    key = (epoch_idx, predictor)
    if EST_MARS[] === nothing || get(EST_MARS[], :key, 0) != key
        c2 = TOML.parsefile(joinpath(STUDY_DIR, "estimator_v2_config.toml")); s = settings_from(c2); v = c2["v2"]
        np = joinpath(mars_results_dir(root, ms.mc), "estimator_numerics.toml")
        isfile(np) && (s = with_settings(s; n_drag_half=Int(TOML.parsefile(np)["n_drag_half"])))
        o = V2Opts(sp_tide=Symbol(v["sp_tide"]), sp_j2=Symbol(v["sp_j2"]), n_sp=Int(v["n_sp_nodes"]))
        rec = mars_record(ms, epoch)
        EST_MARS[] = (; key, ctx=make_context(rec, s), s, o, tab=SPTable(o.n_sp), rot=rec.rot, j3=predictor == "j3diag" ? -sqrt(7.0) * ms.harm.C[4, 1] : nothing,
            zf=predictor in ("v3", "v4", "v5") ? zonal_field(ms.harm) : nothing, zo=ZonalOpts(), tf=predictor in ("v4", "v5") ? tesseral_field(ms.harm) : nothing, predictor)
    end
    return EST_MARS[]
end

"""Next geodetic periapsis altitude after a tangential `dv` at the next apoapsis, from the measured osculating state (J2000 `r0`, `v0`), as the Earth Arm B
but in the Mars body-equator frame (the state is rotated with the run epoch's pole)."""
function predict_estimator_mars(st, r0, v0, dv, et_now)
    ctx = st.ctx
    oe = _corridor_element_summary(SVector{3, Float64}(st.rot' * r0), SVector{3, Float64}(st.rot' * v0), ctx.mu)
    ra = oe.a * (1 + oe.e); va = sqrt(ctx.mu * (2 / ra - 1 / oe.a))
    a2 = 1 / (2 / ra - (va + dv)^2 / ctx.mu); e2 = ra / a2 - 1
    init = (a=oe.a + (a2 - oe.a), e=oe.e + (e2 - oe.e), i=oe.i, raan=oe.raan, omega=oe.omega, nu=oe.nu)
    p = st.predictor == "v5" ? predict_history5(ctx, st.s, st.o, st.tab, st.zf, st.zo, st.tf, TesseralOpts(force=true), TimingOpts(), 0.0, 0.0, 0.0, 0.0, et_now; max_passes=1, init_override=init)[1] :
        st.tf !== nothing ? predict_history4(ctx, st.s, st.o, st.tab, st.zf, st.zo, st.tf, TesseralOpts(force=true), 0.0, 0.0, 0.0, 0.0, et_now; max_passes=1, init_override=init)[1] :
        st.zf !== nothing ? predict_history3(ctx, st.s, st.o, st.tab, st.zf, st.zo, 0.0, 0.0, 0.0, 0.0, et_now; max_passes=1, init_override=init)[1] :
        st.j3 === nothing ? predict_history2(ctx, st.s, st.o, st.tab, 0.0, 0.0, 0.0, 0.0, et_now; max_passes=1, init_override=init) :
        predict_history2_j3_init(ctx, st.s, st.o, st.tab, st.j3, init, et_now)
    isempty(p.alt_m) && throw(ErrorException("estimator predicted no next periapsis (ended=$(p.ended))"))
    return p.alt_m[1]
end

"""UNFROZEN diagnostic: next-periapsis prediction with the averaged J3 term from an osculating state `init` (a, e, i, raan, omega, nu), one pass."""
function predict_history2_j3_init(ctx, s, o, tab, J3, init, et0)
    out = Prediction()
    ivec = SV5(init.a, init.e * cos(init.omega), init.e * sin(init.omega), init.i, init.raan)
    rbs0 = s.tides ? body_positions(ctx, et0) : SV3[]
    d0 = zero(SV5); m = ivec; mm = init
    for _ in 1:(o.sp_tide === :none ? 1 : 3)
        el = vec_el(ivec - d0)
        mm = osculating_to_mean((a=el.a, e=el.e, i=el.i, raan=el.raan, omega=el.omega, nu=init.nu); mu=ctx.mu, radius=ctx.Re, j2=ctx.J2).elements
        m = mean_state(mm, ctx); d0 = tide_sp(m, ctx, o, tab, rbs0, :apo)
    end
    M0 = true_to_mean_anomaly(mm.e, mm.nu); mdot = mean_motion_rate(m, ctx); T = 2π / mdot
    dt_in = (2π - M0) / mdot; dt_out = 0.5 * T; t_k = dt_in
    rbs = body_positions(ctx, et0 + t_k)
    i1, dO1, dw1 = mean_rates_r_j3(m, ctx, s, rbs, J3)
    i2, dO2, dw2 = mean_rates_r_j3(advance(m, i1, dO1, dw1, dt_in), ctx, s, rbs, J3)
    m_pre = advance(m, i2, dO2, dw2, dt_in)
    d3 = tide_sp(m_pre, ctx, o, tab, s.tides ? rbs : SV3[], :peri)
    ov = osc_vec2(m_pre, ctx, o, tab, d3); osc = vec_el(ov); din = zero(SV5)
    if s.drag; din, _ = drag_deltas(osc, ctx, s, fit_rho(ctx)); end
    fin = rp_form(ov, din); m_peri = shift_mean2(m_pre, fin, ctx, o, tab, d3)
    osc = vec_el(osc_vec2(m_peri, ctx, o, tab, d3)); rp = osc.a * (1 - osc.e)
    rv = elements_to_cartesian((a=osc.a, e=osc.e, i=osc.i, raan=osc.raan, omega=osc.omega, nu=0.0); mu=ctx.mu)
    alt, _ = geodetic_altitude_latitude(SV3(rv[1]), ctx.planet); push!(out.alt_m, alt)
    return out
end

# ---- corridor config of one cell -----------------------------------------------------------------------------------------------------

"""CorridorConfig-compatible values for the example's controller classes: heat-rate band as altitudes at the first-pass relative speed, margin, limits."""
function mars_corridor_config(ms::MarsSetup, vrel_first, margin_m)
    mc = ms.mc; b = mc["band"]; w2m = Float64(b["w_cm2_to_w_m2"]); v = mc["vehicle"]
    hq = altitude_of_density_fast(ms.fit, 2 * Float64(b["q_high_w_cm2"]) * w2m / vrel_first^3)   # altitude at Q = 0.32 (lower edge)
    hl = altitude_of_density_fast(ms.fit, 2 * Float64(b["q_low_w_cm2"]) * w2m / vrel_first^3)    # altitude at Q = 0.14 (upper edge)
    mid = (hq + hl) / 2
    values = Dict{String, Any}(
        "mission" => Dict{String, Any}("corridor_min_altitude_m" => hq, "corridor_max_altitude_m" => hl, "reference_altitude_m" => mid, "safety_margin_m" => margin_m,
            "horizon_s" => Float64(mc["grid"]["horizon_days"]) * SECONDS_PER_DAY),
        "vehicle" => Dict{String, Any}("dry_mass_kg" => ms.dry_mass_kg, "propellant_mass_kg" => Float64(v["prop_mass_kg"]), "thrust_n" => Float64(v["thrust_n"]),
            "specific_impulse_s" => Float64(v["specific_impulse_s"])),
        "comparison" => Dict{String, Any}("target_periapsis_altitude_m" => mid, "max_burn_delta_v_mps" => Float64(v["max_burn_delta_v_mps"]),
            "max_burn_duration_s" => Float64(v["max_burn_duration_s"]), "total_delta_v_budget_mps" => Float64(v["total_delta_v_budget_mps"]),
            "planning_lead_deg" => Float64(v["planning_lead_deg"]), "predictor_monotone_tolerance_m" => Float64(v["predictor_monotone_tolerance_m"]),
            "near_passive_deadband_m" => margin_m))
    return CorridorConfig("mars_generality_config.toml", values), hq, hl
end

function mars_closedloop_args(ms::MarsSetup, epoch, rot, horizon_s, dir, cfg::CorridorConfig, state, make_predictor)
    c = cfg.values["comparison"]
    inner = CorridorGuidancePolicy(cfg, :favorable, Symbol(ms.mc["closedloop"]["policy"]), Float64(c["max_burn_delta_v_mps"]), Float64(c["max_burn_duration_s"]),
        Float64(c["total_delta_v_budget_mps"]), Float64(c["near_passive_deadband_m"]), Float64(c["target_periapsis_altitude_m"]),
        deg2rad(Float64(c["planning_lead_deg"])), 0, false, nothing, NamedTuple[], NamedTuple[], String[])
    wrapper = PredictorPolicy(inner, make_predictor)
    v = cfg.values["vehicle"]
    thruster = BaseThrusterModel([Float64(v["thrust_n"])], [0.0], [0.0], [-1.0], [-1.0], [Float64(v["specific_impulse_s"])])
    args = mars_args(ms, epoch, horizon_s, dir; state=state, controller=wrapper, thruster=thruster)
    return args, wrapper
end

"""One closed-loop Mars cell; plain-data row. Never throws: a failure becomes a row with status = "run_error"."""
function mars_closedloop_task(seed)
    t0 = time()
    try
        mc = load_mars(seed.config_path); ms = mars_setup(mc)
        epoch = parse_epoch(mc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms.planet, epoch)
        calib = CSV.read(joinpath(seed.calib_dir, "calibration.csv"), DataFrame)
        seed_alt = calib.seed_alt_m[end]
        cfg, hq, hl = mars_corridor_config(ms, calib.vrel_mps[end], Float64(mc["closedloop"]["margin_m"]))
        dir = joinpath(seed.out, seed.tag); mkpath(dir)
        st = est_state_mars(ms, seed.epoch_idx, epoch, seed.root; predictor=seed.predictor)
        et0 = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel())
        maker = (p, r, v, mass, t) -> (dv -> predict_estimator_mars(st, r, v, Float64(dv), et0 + t))
        state0 = mars_initial_state(ms, rot, seed_alt, seed.inc, seed.raan, seed.omega)
        horizon = Float64(seed.horizon_days) * SECONDS_PER_DAY
        args, wrapper = mars_closedloop_args(ms, epoch, rot, horizon, dir, cfg, state0, maker)
        policy = wrapper.inner; obs = NamedTuple[]
        cbs = (_apsis_observation_callback(obs, args; controller=policy), _preapoapsis_guidance_callback(policy))
        tr = time(); sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=cbs); run_wall = time() - tr
        rows = mars_periapsis_table(sol, args, ms, rot)
        write_rows_m(joinpath(dir, "passes_cl.csv"), rows); write_rows_m(joinpath(dir, "decisions.csv"), policy.decisions)
        ve = Float64(mc["vehicle"]["specific_impulse_s"]) * 9.80665; m0 = ms.wet_mass_kg; mf = Float64(sol.u[end].sc[1].mass)
        plans = policy.burn_plans; bs = band_stats(rows, mc); t_end = sol.t[end]; ret = string(sol.retcode); stopped = !isempty(policy.failures)
        return (; tag=seed.tag, stratum=seed.stratum, status="ok", retcode=ret, t_end_days=t_end / SECONDS_PER_DAY, horizon_days=seed.horizon_days,
            completed=ret == "Success" && !stopped, controller_stopped=stopped, controller_failures=join(policy.failures, " | "),
            total_dv_mps=ve * log(m0 / mf), dv_commanded_mps=sum(p.requested_delta_v_mps for p in plans; init=0.0), burns=length(plans),
            first_burn_days=isempty(plans) ? NaN : plans[1].execution_start_time_s / SECONDS_PER_DAY, passes=length(rows),
            peak_q_w_cm2=bs.passes == 0 ? NaN : bs.q_max, q_min_w_cm2=bs.passes == 0 ? NaN : bs.q_min, exceeds_hard=bs.passes == 0 ? false : bs.exceeds_hard,
            minutes_in_band=bs.passes == 0 ? NaN : bs.time_in_band_days * 1440.0, passes_in_band=bs.passes == 0 ? 0 : bs.passes_in_band,
            min_alt_km=isempty(rows) ? NaN : minimum(r.alt_m for r in rows) / 1e3, corridor_min_alt_m=hq, corridor_max_alt_m=hl,
            predictor_evals=sum((d.predictor_evaluations for d in policy.decisions); init=0), run_wall_s=run_wall, wall_s=time() - t0)
    catch err
        return (; tag=seed.tag, stratum=seed.stratum, status="run_error", retcode="", t_end_days=NaN, horizon_days=seed.horizon_days, completed=false,
            controller_stopped=false, controller_failures=replace(first(sprint(showerror, err), 300), "\n" => " "), total_dv_mps=NaN, dv_commanded_mps=NaN,
            burns=0, first_burn_days=NaN, passes=0, peak_q_w_cm2=NaN, q_min_w_cm2=NaN, exceeds_hard=false, minutes_in_band=NaN, passes_in_band=0,
            min_alt_km=NaN, corridor_min_alt_m=NaN, corridor_max_alt_m=NaN, predictor_evals=0, run_wall_s=NaN, wall_s=time() - t0)
    end
end

# ---- predictor check (serial): estimator next-periapsis shift against the truth pass-to-pass change, passive ---------------------------------------

function mars_verify_task(seed)
    mc = load_mars(seed.config_path); ms = mars_setup(mc)
    epoch = parse_epoch(mc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms.planet, epoch)
    calib = CSV.read(joinpath(seed.calib_dir, "calibration.csv"), DataFrame); seed_alt = calib.seed_alt_m[end]
    cfg, hq, hl = mars_corridor_config(ms, calib.vrel_mps[end], Float64(mc["closedloop"]["margin_m"]))
    cfg.values["mission"]["corridor_min_altitude_m"] = 60_000.0; cfg.values["mission"]["corridor_max_altitude_m"] = 400_000.0      # never burns
    cfg.values["comparison"]["target_periapsis_altitude_m"] = 230_000.0; cfg.values["mission"]["reference_altitude_m"] = 230_000.0
    dir = joinpath(seed.out, seed.tag); mkpath(dir)
    st = est_state_mars(ms, seed.epoch_idx, epoch, seed.root; predictor=seed.predictor); et0 = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel())
    log = NamedTuple[]
    maker = (p, r, v, mass, t) -> (dv -> (dv == 0.0 && (isempty(log) || log[end].t_s != t) && push!(log, (; t_s=t, pred_armB_est_m=predict_estimator_mars(st, r, v, 0.0, et0 + t))); log[end].pred_armB_est_m))
    state0 = mars_initial_state(ms, rot, seed_alt, seed.inc, seed.raan, seed.omega)
    args, wrapper = mars_closedloop_args(ms, epoch, rot, Float64(seed.horizon_days) * SECONDS_PER_DAY, dir, cfg, state0, maker)
    policy = wrapper.inner; obs = NamedTuple[]
    sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=(_apsis_observation_callback(obs, args; controller=policy), _preapoapsis_guidance_callback(policy)))
    rows = mars_periapsis_table(sol, args, ms, rot); out = NamedTuple[]
    for (k, l) in enumerate(log)
        k + 1 <= length(rows) || break
        a0 = rows[k].alt_m; a1 = rows[k + 1].alt_m
        push!(out, merge(l, (; decision=k, truth_alt_k_m=a0, truth_alt_next_m=a1, truth_shift_m=a1 - a0, shift_armB_est_m=l.pred_armB_est_m - a0, err_m=l.pred_armB_est_m - a1)))
    end
    write_rows_m(joinpath(dir, "verify_predictor.csv"), out); return out
end

function cl_main(args=ARGS)
    smoke = "--smoke" in args; predictor = "frozen"; root = nothing; workers = nothing; only = nothing; hd_arg = nothing; verify = nothing; cells_file = "closedloop_cells.csv"; cfg_path = MARS_CONFIG_PATH
    for a in args
        startswith(a, "--results-root=") && (root = abspath(split(a, "=", limit=2)[2]))
        startswith(a, "--workers=") && (workers = parse(Int, split(a, "=", limit=2)[2]))
        startswith(a, "--only=") && (only = split(split(a, "=", limit=2)[2], ","))
        startswith(a, "--horizon-days=") && (hd_arg = parse(Float64, split(a, "=", limit=2)[2]))
        startswith(a, "--cells-file=") && (cells_file = String(split(a, "=", limit=2)[2]))
        startswith(a, "--verify=") && (verify = split(a, "=", limit=2)[2])
        startswith(a, "--config=") && (cfg_path = abspath(split(a, "=", limit=2)[2]))
        startswith(a, "--predictor=") && (predictor = String(split(a, "=", limit=2)[2]))
    end
    mc = load_mars(cfg_path); study = load_study(); root = something(root, results_root(study)); base = mars_results_dir(root, mc)
    cells = CSV.read(joinpath(base, cells_file), DataFrame)
    mkseed(r, out, hd) = (; config_path=cfg_path, root, out, tag=r.tag, stratum=r.stratum, inc=r.inc, raan=r.raan, omega=r.omega, epoch_idx=Int(r.epoch),
        calib_dir=joinpath(base, "truth", r.tag), horizon_days=hd, predictor)
    if verify !== nothing
        r = first(eachrow(cells[cells.tag .== verify, :])); out = joinpath(base, "closedloop_verify"); mkpath(out)
        foreach(println, mars_verify_task(mkseed(r, out, something(hd_arg, 6.0)))); return
    end
    only === nothing || (cells = cells[in.(cells.tag, Ref(only)), :])
    smoke && only === nothing && (cells = first(cells, 2))
    hd = something(hd_arg, smoke ? 5.0 : Float64(mc["closedloop"]["horizon_days"])); workers = something(workers, smoke ? 2 : Int(mc["closedloop"]["workers"]))
    suffix = (cells_file == "closedloop_cells.csv" ? "" : "_" * replace(replace(cells_file, "closedloop_cells_" => ""), ".csv" => "")) * (predictor == "frozen" ? "" : "_" * predictor)
    out = joinpath(base, "closedloop" * suffix * (smoke ? "_smoke" : "")); mkpath(out)
    seeds = [mkseed(r, out, hd) for r in eachrow(cells)]
    println("Mars closed loop (Arm B): $(length(seeds)) cells, horizon $(hd) d, workers $(workers), smoke=$(smoke)")
    res = run_process_campaign(mars_closedloop_task, seeds, workers; mission_time_s=hd * SECONDS_PER_DAY, worker_files=[joinpath(@__DIR__, "run_mars_closedloop.jl")])
    println("MC route=$(res.route) workers=$(res.threads) samples=$(length(res.samples)) failed=$(length(res.failed)) elapsed_s=$(round(res.elapsed_s; digits=1))")
    foreach(f -> println("FAILED $(f.seed): $(first(sprint(showerror, f.error), 400))"), res.failed)
    write_rows_m(joinpath(out, "closedloop_rows.csv"), [s.value for s in res.successful])
    open(joinpath(out, "closedloop_mc.toml"), "w") do io
        println(io, "mc_route = \"$(res.route)\"\nmc_workers = $(res.threads)\nmc_samples = $(length(res.samples))\nmc_failed = $(length(res.failed))\nmc_elapsed_s = $(res.elapsed_s)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    cl_main()
end
