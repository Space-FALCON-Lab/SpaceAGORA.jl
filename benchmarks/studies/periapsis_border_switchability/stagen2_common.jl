# Stage N2: edge-aimed chance-constrained controller (plan: tasks/2026-10-15_periapsis-region-border-targeting/stagen2_plan.md). Reuses stagen_flight_common.jl (shim, Stage H files, online Kalman filter, included, not edited).
# Changes ONLY: (1) the edge multiplier t_mult may differ from the Stage N value (calibrated, plan section 3); (2) the aim point: target = the altitude where the predicted Q equals Q_high exp(-(c + g_a)), i.e. g_a >= 0 in ln Q below the
# chance-constrained edge (g_a = 0 aims AT the edge), clamped to [lo + 10 m, hi - 10 m] (a target exactly on the edge can land 1 cm outside it and the burn is then dropped as infeasible); the Stage N target was the band center. Filter constants, low-heat edge (Stage G rule capped to lo + w_m) and everything else as Stage N-flight.
include(joinpath(@__DIR__, "stagen_flight_common.jl"))

"""Arm code; must equal arm_code() in stagen2_analysis.py: cc<eps%>_a<100 g_a>_t<100 t_mult>_w<w>."""
stagen2_arm_code(eps::Real, t_mult::Real, aim::Real, w_m::Real) = "cc" * lpad(string(round(Int, 100eps)), 2, '0') * "_a" * lpad(string(round(Int, 100aim)), 3, '0') * "_t" * lpad(string(round(Int, 100t_mult)), 3, '0') * "_w" * lpad(string(round(Int, w_m)), 5, '0')

mutable struct StageN2Policy <: SimulationModel.AbstractGuidanceModel
    f::StageFPolicy
    kr::Float64; kq::Float64; km0::Float64; kv0::Float64
    t_mult::Float64               # edge multiplier c = t_mult s
    w_m::Float64                  # band-width cap (m); <= 0 means none
    aim_ln::Float64               # aim offset below the edge in ln Q (>= 0)
end

"""Edges as stagen_edges plus the edge-side aim target."""
function stagen2_edges(pol)
    f = pol.f; V3 = f.vrel_first^3; n = length(f.meas)
    pr = stagen_predict(Float64[m.y for m in f.meas], pol.kr, pol.kq, pol.km0, pol.kv0); s = sqrt(pr.s2); c = pol.t_mult * s
    lo = altitude_of_density_fast(f.fit, 2 * f.q_high_w_m2 / V3 * exp(-(pr.m + c)))
    hi = altitude_of_density_fast(f.fit, 2 * f.q_low_w_m2 / V3 * exp(-pr.m))
    pol.w_m > 0 && (hi = min(hi, lo + pol.w_m))
    narrow = !(hi - lo >= f.min_band_m)
    narrow && (hi = lo + f.min_band_m)
    ta = altitude_of_density_fast(f.fit, 2 * f.q_high_w_m2 / V3 * exp(-(pr.m + c + pol.aim_ln)))
    return (lo=lo, hi=hi, margin=0.0, target=clamp(ta, lo + 10.0, hi - 10.0), yhat=pr.m, collapsed=narrow, n=n, m=pr.m, s=s, c=c)
end

function SimulationModel.GuidanceHooks.calcGuidanceEffect!(wrapper::StageN2Policy, u, p, t::Float64, sat_idx::Int64)
    f = wrapper.f; policy = f.inner; commands = p.shared_buffers.maneuver_commands
    (sat_idx < 1 || sat_idx > length(commands)) && return nothing
    if !policy.armed
        commands[sat_idx] = SimulationModel.CommandTypes.PropulsiveManeuverCommand(valid=true, source_orbit=policy.periapsis_count + 1)
        return nothing
    end
    r = SVector{3, Float64}(u.sc[sat_idx].pos); v = SVector{3, Float64}(u.sc[sat_idx].vel); mass = Float64(u.sc[sat_idx].mass)
    vehicle = policy.config.values["vehicle"]; exhaust_velocity = Float64(vehicle["specific_impulse_s"]) * 9.80665
    consumed_at_limit = Float64(vehicle["thrust_n"]) * policy.max_burn_duration_s / exhaust_velocity
    duration_limited_dv = consumed_at_limit < mass ? exhaust_velocity * log(mass / (mass - consumed_at_limit)) : Inf
    wet_mass = Float64(vehicle["dry_mass_kg"]) + Float64(vehicle["propellant_mass_kg"])
    delivered_so_far = mass < wet_mass ? exhaust_velocity * log(wet_mass / mass) : 0.0
    policy.remaining_budget_mps = max(0.0, Float64(policy.config.values["comparison"]["total_delta_v_budget_mps"]) - delivered_so_far)
    max_allowed_dv = min(policy.max_dv_mps, policy.remaining_budget_mps, duration_limited_dv)
    predictor = f.make_predictor(p, r, v, mass, Float64(t))
    ed = stagen2_edges(wrapper)
    choice = select_corridor_maneuver(predictor; lower_m=ed.lo, upper_m=ed.hi, margin_m=ed.margin, target_m=ed.target, policy=policy.policy,
        max_dv_mps=max_allowed_dv, deadband_m=policy.deadband_m, monotone_tolerance_m=Float64(policy.config.values["comparison"]["predictor_monotone_tolerance_m"]))
    policy.last_selection = choice; target_pass = policy.periapsis_count + 1
    push!(f.edge_log, (; request_time_s=Float64(t), target_pass, n_meas=ed.n, yhat=ed.yhat, lower_m=ed.lo, upper_m=ed.hi, collapsed=ed.collapsed, delta_v_mps=choice.delta_v_mps, reason=String(choice.reason), m=ed.m, s=ed.s, c=ed.c, t_mult=wrapper.t_mult, w_m=wrapper.w_m, aim=wrapper.aim_ln, target_m=ed.target))
    push!(policy.decisions, (request_time_s=Float64(t), source_orbit=target_pass, target_pass=target_pass, policy=String(policy.policy), delta_v_mps=choice.delta_v_mps,
        predicted_periapsis_altitude_m=choice.predicted_periapsis_altitude_m, target_periapsis_altitude_m=choice.target_periapsis_altitude_m, trigger_reason=String(choice.reason),
        feasible=choice.feasible, predictor_evaluations=choice.evaluations))
    commands[sat_idx] = SimulationModel.CommandTypes.PropulsiveManeuverCommand(valid=true, delta_v_mps=choice.feasible ? abs(choice.delta_v_mps) : 0.0,
        direction_rad=choice.delta_v_mps >= 0.0 ? 0.0 : π, source_orbit=target_pass)
    policy.armed = false
    return nothing
end


"""One Stage N-flight flight (perturbed arm, chance-constrained edge rule). seed fields: as stagen_run_task plus aim (filter constants and quantile multiplier).
The near_passive deadband / config margin is the Stage I value (fixed_margin_m). Plain-data row; never throws."""
function stagen2_run_task(seed)
    t0 = time()
    try
        meta = stagef_meta(seed.config_path); sc = load_stagef(seed.planet; path=seed.config_path); ms0 = stagef_setup(sc)
        epoch = parse_epoch(sc["grid"]["epochs_utc"][seed.epoch_idx]); rot = mars_rot(ms0.planet, epoch)
        cdir = joinpath(seed.cal_root, seed.tag); G = CSV.read(joinpath(cdir, "guess.csv"), DataFrame)[1, :]; N = CSV.read(joinpath(cdir, "nodes.csv"), DataFrame)
        pwf = pw_fit(N.h_m, N.rho_kg_m3, Float64(sc["atmosphere"]["temperature_k"]), "piecewise"); ms = with_fit(ms0, pwf)
        stagef_env!(true, meta)
        model = stagef_gram_model(sc, epoch, stagef_seed(seed.planet * seed.tag)); core_t = model.core.initial_time
        pargs = gram_pargs(ms.planet, Float64(sc["atmosphere_gram"]["entry_interface_km"]))
        dir = joinpath(seed.out, seed.run_id); mkpath(dir); horizon = Float64(seed.horizon_days) * SECONDS_PER_DAY
        state0 = mars_initial_state(ms, rot, Float64(G.seed_alt_m), seed.inc, seed.raan, seed.omega)
        base = gram_args(ms, sc, model, epoch, horizon, dir; state=state0)
        meas = NamedTuple[]
        margin = Float64(meta["stagef"]["fixed_margin_m"]); cfg, hq, hl = stagef_corridor_config(ms, Float64(G.vrel_first_mps), margin); c = cfg.values["comparison"]
        inner = CorridorGuidancePolicy(cfg, :favorable, Symbol(meta["stagef"]["policy"]), Float64(c["max_burn_delta_v_mps"]), Float64(c["max_burn_duration_s"]),
            Float64(c["total_delta_v_budget_mps"]), Float64(c["near_passive_deadband_m"]), Float64(c["target_periapsis_altitude_m"]), deg2rad(Float64(c["planning_lead_deg"])),
            0, false, nothing, NamedTuple[], NamedTuple[], String[])
        st = stagef_est_state(ms, epoch); et0 = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel())
        maker = (p, r, v, mass, t) -> (dv -> stagef_predict_alt(st, r, v, Float64(dv), et0 + t))
        w2m = Float64(sc["band"]["w_cm2_to_w_m2"])
        fpol = StageFPolicy(inner, maker, :adaptive, pwf, Float64(G.vrel_first_mps), Float64(sc["band"]["q_low_w_cm2"]) * w2m, Float64(sc["band"]["q_high_w_cm2"]) * w2m,
            margin, meas, "", 0.0, 0.0, Float64(meta["analysis"]["adapt_min_band_m"]), NamedTuple[])
        pol = StageN2Policy(fpol, Float64(seed.kr), Float64(seed.kq), Float64(seed.km0), Float64(seed.kv0), Float64(seed.t_mult), Float64(seed.w_m), Float64(seed.aim))
        v = cfg.values["vehicle"]; thruster = BaseThrusterModel([Float64(v["thrust_n"])], [0.0], [0.0], [-1.0], [-1.0], [Float64(v["specific_impulse_s"])])
        args = SimulationConfiguration(file_paths=base.file_paths, simulation_settings=base.simulation_settings, mission_configuration=base.mission_configuration,
            environment_model=base.environment_model, dynamics_model=base.dynamics_model, guidance_model=GuidanceModel((pol,), [1.0e12]), navigation_model=base.navigation_model,
            control_model=ControlModel((thruster,), [1.0e12]), initial_time=base.initial_time, integration_tolerances=base.integration_tolerances, solver_config=base.solver_config)
        obs = NamedTuple[]
        extra = (stagef_apsis_callback(obs, meas, args, inner, model, pargs, pwf), _preapoapsis_guidance_callback(inner))
        tr = time(); sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=extra); run_wall = time() - tr
        logpath = joinpath(dir, "gram_density_perturbation_log.csv")
        rows = stagef_table(sol, args, ms, rot, sc, pwf, logpath); isfile(logpath) && rm(logpath; force=true); isfile(logpath * ".summary.toml") && rm(logpath * ".summary.toml"; force=true)
        write_rows_m(joinpath(dir, "passes.csv"), rows)
        ve = Float64(sc["vehicle"]["specific_impulse_s"]) * 9.80665; dv = ve * log(ms.wet_mass_kg / Float64(sol.u[end].sc[1].mass)); burns = length(inner.burn_plans)
        stopped = !isempty(inner.failures); failures = join(inner.failures, " | ")
        write_rows_m(joinpath(dir, "decisions.csv"), fpol.edge_log); write_rows_m(joinpath(dir, "measured.csv"), meas)
        ret = string(sol.retcode)
        return (; planet=seed.planet, tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl="cc", code=seed.code, eps=Float64(seed.eps), w_m=Float64(seed.w_m), t_mult=Float64(seed.t_mult), aim=Float64(seed.aim), status="ok", retcode=ret, passes=length(rows), t_end_days=sol.t[end] / SECONDS_PER_DAY,
            horizon_days=seed.horizon_days, completed=ret == "Success" && !stopped, controller_stopped=stopped, controller_failures=failures, total_dv_mps=dv, burns=burns,
            collapsed_decisions=count(e -> e.collapsed, fpol.edge_log), gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second),
            run_wall_s=run_wall, wall_s=time() - t0, error="")
    catch err
        return (; planet=seed.planet, tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl="cc", code=seed.code, eps=NaN, w_m=NaN, t_mult=NaN, aim=NaN,
            status="run_error", retcode="", passes=0, t_end_days=NaN, horizon_days=seed.horizon_days, completed=false, controller_stopped=false, controller_failures="", total_dv_mps=NaN, burns=0,
            collapsed_decisions=0, gram_core_epoch="", run_wall_s=NaN, wall_s=time() - t0, error=replace(first(sprint(showerror, err), 500), "\n" => " "))
    end
end
