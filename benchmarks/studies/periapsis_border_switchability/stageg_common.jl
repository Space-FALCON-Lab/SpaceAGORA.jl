# Stage G: asymmetric adaptive band. Reuses stagef_common.jl unchanged and adds only the asymmetric edge rule, its policy wrapper and the case function of the F2-style run.
# The thermal-protecting (high-heat) edge keeps its one-sided empirical error quantile e_hi; the low-heat edge takes NO buffer (edge at the nominal predicted Q_low shifted only by the predicted
# density residual yhat) and the corridor never collapses to its center: if the band is narrower than min_band_m the thermal edge stays and the low-heat edge moves up to lo + min_band_m.
# No src/ or examples/ change; Stage F files are not edited.
include(joinpath(@__DIR__, "stagef_common.jl"))

mutable struct StageGPolicy <: SimulationModel.AbstractGuidanceModel
    f::StageFPolicy               # carries inner policy, predictor, fit, measurements and edge log; its mode and e_lo are unused here
end

"""Asymmetric corridor edges at this decision (altitudes): lo = high-Q (thermal) edge with the e_hi buffer, hi = low-Q edge at the predicted nominal Q_low (no buffer)."""
function stageg_edges(pol)   # duck-typed on the StageFPolicy fields used (lets the selftest pass a NamedTuple)
    V3 = pol.vrel_first^3; n = length(pol.meas)
    yhat = stagef_predict_next(pol.pred_name, Float64[m.y for m in pol.meas], Float64[m.lon_rad for m in pol.meas])
    lo = altitude_of_density_fast(pol.fit, 2 * pol.q_high_w_m2 / V3 * exp(-(yhat + pol.e_hi)))
    hi = altitude_of_density_fast(pol.fit, 2 * pol.q_low_w_m2 / V3 * exp(-yhat))
    narrow = !(hi - lo >= pol.min_band_m)
    narrow && (hi = lo + pol.min_band_m)
    return (lo=lo, hi=hi, margin=0.0, target=(lo + hi) / 2, yhat=yhat, collapsed=narrow, n=n)
end

function SimulationModel.GuidanceHooks.calcGuidanceEffect!(wrapper::StageGPolicy, u, p, t::Float64, sat_idx::Int64)
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
    ed = stageg_edges(f)
    choice = select_corridor_maneuver(predictor; lower_m=ed.lo, upper_m=ed.hi, margin_m=ed.margin, target_m=ed.target, policy=policy.policy,
        max_dv_mps=max_allowed_dv, deadband_m=policy.deadband_m, monotone_tolerance_m=Float64(policy.config.values["comparison"]["predictor_monotone_tolerance_m"]))
    policy.last_selection = choice; target_pass = policy.periapsis_count + 1
    push!(f.edge_log, (; request_time_s=Float64(t), target_pass, n_meas=ed.n, yhat=ed.yhat, lower_m=ed.lo, upper_m=ed.hi, collapsed=ed.collapsed, delta_v_mps=choice.delta_v_mps, reason=String(choice.reason)))
    push!(policy.decisions, (request_time_s=Float64(t), source_orbit=target_pass, target_pass=target_pass, policy=String(policy.policy), delta_v_mps=choice.delta_v_mps,
        predicted_periapsis_altitude_m=choice.predicted_periapsis_altitude_m, target_periapsis_altitude_m=choice.target_periapsis_altitude_m, trigger_reason=String(choice.reason),
        feasible=choice.feasible, predictor_evaluations=choice.evaluations))
    commands[sat_idx] = SimulationModel.CommandTypes.PropulsiveManeuverCommand(valid=true, delta_v_mps=choice.feasible ? abs(choice.delta_v_mps) : 0.0,
        direction_rad=choice.delta_v_mps >= 0.0 ? 0.0 : π, source_orbit=target_pass)
    policy.armed = false
    return nothing
end

"""One Stage G run (perturbed arm, asymmetric adaptive controller). Same seed fields and construction as stagef_run_task with ctrl "adaptive"; only the policy wrapper differs.
Plain-data row; never throws."""
function stageg_run_task(seed)
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
            margin, meas, String(seed.adapt_pred), Float64(seed.adapt_e_hi), 0.0, Float64(meta["analysis"]["adapt_min_band_m"]), NamedTuple[])
        pol = StageGPolicy(fpol)
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
        return (; planet=seed.planet, tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl=seed.ctrl, status="ok", retcode=ret, passes=length(rows), t_end_days=sol.t[end] / SECONDS_PER_DAY,
            horizon_days=seed.horizon_days, completed=ret == "Success" && !stopped, controller_stopped=stopped, controller_failures=failures, total_dv_mps=dv, burns=burns,
            collapsed_decisions=count(e -> e.collapsed, fpol.edge_log), gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second),
            run_wall_s=run_wall, wall_s=time() - t0, error="")
    catch err
        return (; planet=seed.planet, tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl=seed.ctrl, status="run_error", retcode="", passes=0, t_end_days=NaN, horizon_days=seed.horizon_days,
            completed=false, controller_stopped=false, controller_failures="", total_dv_mps=NaN, burns=0, collapsed_decisions=0, gram_core_epoch="", run_wall_s=NaN, wall_s=time() - t0,
            error=replace(first(sprint(showerror, err), 500), "\n" => " "))
    end
end
