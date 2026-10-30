# Stage L: the dV / campaign-time / violation frontier. Reuses stagef_common.jl and stageh_common.jl (included, not edited) and adds only (1) the arm code and (2) one case function that flies either
# a fixed-margin arm (any margin; the Stage F :fixed policy) or a member of the continuous family (the Stage H edge rule with an arbitrary SIGNED low-heat offset g and thermal buffer e_hi).
# Family edge rule (stageh_edges, unchanged): thermal edge lo = alt(2 Q_high / V^3 exp(-(yhat + e_hi))), low-heat edge hi = alt(2 Q_low / V^3 exp(-(yhat + g))), g = s * |e_lo| signed:
#   g > 0 widens toward lower heating (Stage H), g = 0 is Stage G, g < 0 narrows (Stage F's sign; s = -1 is Stage F's e_lo). Narrow band: thermal edge kept, hi = lo + min_band_m (Stage G/H anchoring).
include(joinpath(@__DIR__, "stageh_common.jl"))

"""Arm code (run directory prefix); must equal arm_code() in stagel_analysis.py. kind "fixedm" uses margin_m, kind "adaptive_l" uses p_hi and s_frac."""
function stagel_arm_code(kind::AbstractString, margin_m::Real, p_hi::Real, s_frac::Real)
    kind == "fixedm" && return "fm" * lpad(string(round(Int, margin_m)), 5, '0')
    kind == "adaptive_l" && return "fam_p" * lpad(string(round(Int, 100p_hi)), 3, '0') * "_s" * (s_frac < 0 ? "m" : "p") * lpad(string(round(Int, 100abs(s_frac))), 3, '0')
    throw(ArgumentError("unknown kind $kind"))
end

"""One Stage L flight (perturbed arm). seed fields as stageh_run_task plus kind ("fixedm" | "adaptive_l"), margin_m, p_hi, s_frac (labels; the controller reads adapt_e_hi, adapt_e_lo = g directly).
The near_passive deadband / config margin is the Stage I value (fixed_margin_m) for the family and margin_m for fixedm. Plain-data row; never throws."""
function stagel_run_task(seed)
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
        meas = NamedTuple[]; fixed = seed.kind == "fixedm"
        margin = fixed ? Float64(seed.margin_m) : Float64(meta["stagef"]["fixed_margin_m"]); cfg, hq, hl = stagef_corridor_config(ms, Float64(G.vrel_first_mps), margin); c = cfg.values["comparison"]
        inner = CorridorGuidancePolicy(cfg, :favorable, Symbol(meta["stagef"]["policy"]), Float64(c["max_burn_delta_v_mps"]), Float64(c["max_burn_duration_s"]),
            Float64(c["total_delta_v_budget_mps"]), Float64(c["near_passive_deadband_m"]), Float64(c["target_periapsis_altitude_m"]), deg2rad(Float64(c["planning_lead_deg"])),
            0, false, nothing, NamedTuple[], NamedTuple[], String[])
        st = stagef_est_state(ms, epoch); et0 = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel())
        maker = (p, r, v, mass, t) -> (dv -> stagef_predict_alt(st, r, v, Float64(dv), et0 + t))
        w2m = Float64(sc["band"]["w_cm2_to_w_m2"])
        fpol = StageFPolicy(inner, maker, fixed ? :fixed : :adaptive, pwf, Float64(G.vrel_first_mps), Float64(sc["band"]["q_low_w_cm2"]) * w2m, Float64(sc["band"]["q_high_w_cm2"]) * w2m,
            margin, meas, fixed ? "" : String(seed.adapt_pred), fixed ? 0.0 : Float64(seed.adapt_e_hi), fixed ? 0.0 : Float64(seed.adapt_e_lo), Float64(meta["analysis"]["adapt_min_band_m"]), NamedTuple[])
        pol = fixed ? fpol : StageHPolicy(fpol)
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
        return (; planet=seed.planet, tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl=seed.kind, code=seed.code, margin_m=Float64(seed.margin_m), p_hi=Float64(seed.p_hi), s_frac=Float64(seed.s_frac),
            adapt_e_hi=Float64(seed.adapt_e_hi), adapt_e_lo=Float64(seed.adapt_e_lo), status="ok", retcode=ret, passes=length(rows), t_end_days=sol.t[end] / SECONDS_PER_DAY,
            horizon_days=seed.horizon_days, completed=ret == "Success" && !stopped, controller_stopped=stopped, controller_failures=failures, total_dv_mps=dv, burns=burns,
            collapsed_decisions=count(e -> e.collapsed, fpol.edge_log), gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second),
            run_wall_s=run_wall, wall_s=time() - t0, error="")
    catch err
        return (; planet=seed.planet, tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl=seed.kind, code=seed.code, margin_m=NaN, p_hi=NaN, s_frac=NaN, adapt_e_hi=NaN, adapt_e_lo=NaN,
            status="run_error", retcode="", passes=0, t_end_days=NaN, horizon_days=seed.horizon_days, completed=false, controller_stopped=false, controller_failures="", total_dv_mps=NaN, burns=0,
            collapsed_decisions=0, gram_core_epoch="", run_wall_s=NaN, wall_s=time() - t0, error=replace(first(sprint(showerror, err), 500), "\n" => " "))
    end
end
