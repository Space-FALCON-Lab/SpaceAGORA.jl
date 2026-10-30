# Half-lever test I1, GRAM-on twin: one flight of the burn / no-burn twin (tasks/2026-10-15_periapsis-region-border-targeting/halflever_twin_results.md).
# Loaded on the process workers AFTER stageo_common.jl (worker_files order), and by the driver halflever_twin_run.jl. No src/ or examples/ change, no new burn model.
# Flight = the Stage O Odyssey-band flight set-up (stageo_run_task: Odyssey start, calibrated first-pass periapsis, Mars50c 20x20 + Sun, GRAM truth density with the pass-density perturbation
# SPACEAGORA_GRAM_DENSITY_PERTURBATION=pass, reseed per pass, 1.0 s atmosphere step) with the KEEPER REPLACED by the vacuum test's burn machinery (halflever_vacuum.jl): ONE burn of `dv` m/s through the
# engine's apoapsis-centered scheduler at apoapsis j after the start (j = 1: the apoapsis the N2 keepers' first burns use; j = 2: the F-O keeper's), no other burn. dv = 0 is the no-burn twin (same callbacks, same scheduler
# call with a zero command). Both flights stop at the first apoapsis after the `npass`-th periapsis.
# Extra logging (no effect on the dynamics, identical in every flight): a per-accepted-step record below `gate_altitude_m` (state, geodetic altitude, GRAM mean density, held perturbation factor, the drag
# acceleration reconstructed from the engine's CorridorBallisticDrag formula with the perturbed density) -> steps.csv; every periapsis / apoapsis event -> events.csv; the burn request and window -> burn.csv.
# The altitude crossings, osculating elements and cumulative drag delta-v are derived from steps.csv by halflever_twin_analyze.py.
using OrdinaryDiffEq

"""Per-step record callback: below `gate` meters of radius above the equatorial radius every accepted step is logged with the density and the drag acceleration (body-equator inertial frame)."""
function twin_step_callback(rows::Vector{NamedTuple}, args, model, pargs, ms, rot, gate::Real)
    env = args.environment_model; planet = env.planet; et0 = SimulationModel.ephemerides_time_seconds(args.initial_time, env.ephemerides_model)
    Rt = transpose(rot); rp_e = Float64(planet.Rp_e); beta = Float64(ms.beta); mref = Float64(ms.wet_mass_kg); g = Float64(gate)
    condition(u, t, integrator) = norm(u.sc[1].pos) - rp_e < g
    function affect!(integrator)
        t = Float64(integrator.t); s = integrator.u.sc[1]
        lpi = SimulationModel.planet_frame_lpi(planet, et0 + t, env.ephemerides_model); rpp = lpi * SVector{3, Float64}(s.pos)
        alt, lat = geodetic_altitude_latitude(rpp, planet); lon = SimulationModel.SimulationCallbacks.rtolatlong(SVector{3, Float64}(rpp), planet)[3]
        rho = gram_rho(model, pargs, alt, lat, lon, t); r = stagef_pass_r(integrator.p, t)
        vpp = lpi * SVector{3, Float64}(s.vel) - cross(planet.ω, rpp); sp = norm(vpp); bcur = beta * Float64(s.mass) / mref
        a_ii = transpose(lpi) * (-0.5 * rho * r * sp / bcur * vpp)   # CorridorBallisticDrag.wrench / mass, density = GRAM mean x held perturbation factor
        a = Rt * a_ii; p = Rt * SVector{3, Float64}(s.pos); v = Rt * SVector{3, Float64}(s.vel)
        push!(rows, (; t_s=t, rx=p[1], ry=p[2], rz=p[3], vx=v[1], vy=v[2], vz=v[3], mass_kg=Float64(s.mass), alt_m=alt, lat_rad=lat, lon_rad=lon, rho_mean=rho, r_pert=r, adx=a[1], ady=a[2], adz=a[3]))
        return nothing
    end
    return DiscreteCallback(condition, affect!; save_positions=(false, false))
end

"""seed fields: config_path, tag (GRAM seed string tag, e.g. ODY_s1), name (run id), dv (m/s, signed; 0 = no-burn twin), burn_j (1 = the first apoapsis after the start, 2 = the second), npass, horizon_days (safety), out, cal_root, gate_altitude_m. Plain-data row; never throws."""
function twin_task(seed)
    t0 = time()
    try
        meta, sc = stageo_load(seed.config_path, "odyssey"); SO = meta["stageo"]; ms0 = stagef_setup(sc)
        epoch = parse_epoch(sc["grid"]["epochs_utc"][1]); rot = mars_rot(ms0.planet, epoch); Rt = transpose(rot)
        cdir = joinpath(seed.cal_root, "odyssey"); G = CSV.read(joinpath(cdir, "guess.csv"), DataFrame)[1, :]; N = CSV.read(joinpath(cdir, "nodes.csv"), DataFrame)
        pwf = pw_fit(N.h_m, N.rho_kg_m3, Float64(sc["atmosphere"]["temperature_k"]), "piecewise"); ms = with_fit(ms0, pwf)
        stagef_env!(true, meta)
        model = stagef_gram_model(sc, epoch, stagef_seed("mars" * seed.tag)); core_t = model.core.initial_time
        pargs = gram_pargs(ms.planet, Float64(sc["atmosphere_gram"]["entry_interface_km"]))
        dir = joinpath(seed.out, seed.name); mkpath(dir)
        horizon = Float64(seed.horizon_days) * 86400.0   # safety only: the run ends at the first apoapsis after periapsis `npass`
        state0 = mars_initial_state(ms, rot, Float64(G.seed_alt_m), Float64(SO["inclination_deg"]), Float64(SO["raan_deg"]), Float64(SO["omega_deg"]))
        base = gram_args(ms, sc, model, epoch, horizon, dir; state=state0)
        margin = Float64(meta["stagef"]["fixed_margin_m"]); cfg, _, _ = stagef_corridor_config(ms, Float64(G.vrel_first_mps), margin)
        lead = deg2rad(Float64(cfg.values["comparison"]["planning_lead_deg"])); v = cfg.values["vehicle"]
        ctrl = CorridorFixedBurnGuidance(0.0, 1, false)
        thruster = BaseThrusterModel([Float64(v["thrust_n"])], [0.0], [0.0], [-1.0], [-1.0], [Float64(v["specific_impulse_s"])])
        args = SimulationConfiguration(file_paths=base.file_paths, simulation_settings=base.simulation_settings, mission_configuration=base.mission_configuration,
            environment_model=base.environment_model, dynamics_model=base.dynamics_model, guidance_model=GuidanceModel((ctrl,), [1.0e12]), navigation_model=base.navigation_model,
            control_model=ControlModel((thruster,), [1.0e12]), initial_time=base.initial_time, integration_tolerances=base.integration_tolerances,
            solver_config=SolverConfig(solver_mode=base.solver_config.solver_mode, maxiters=Int(SO["solver_maxiters"])))
        mu = args.environment_model.planet.μ; target_nu = π - lead; dv = Float64(seed.dv); bj = Int(seed.burn_j); npass = Int(seed.npass)
        events = NamedTuple[]; burn = NamedTuple[]; steps = NamedTuple[]; count = Ref(0); nper = Ref(0)
        function rec!(integrator, kind)
            integrator.t > 1e-3 || return nothing
            kind == "periapsis" && (nper[] += 1)
            s = integrator.u.sc[1]; r = Rt * SVector{3, Float64}(s.pos); w = Rt * SVector{3, Float64}(s.vel)
            push!(events, (; kind, pass=nper[], t_s=Float64(integrator.t), rx=r[1], ry=r[2], rz=r[3], vx=w[1], vy=w[2], vz=w[3], mass_kg=Float64(s.mass)))
            kind == "apoapsis" && nper[] >= npass && terminate!(integrator)
            return nothing
        end
        apsis = ContinuousCallback((u, t, i) -> dot(u.sc[1].pos, u.sc[1].vel), i -> rec!(i, "periapsis"), i -> rec!(i, "apoapsis"); save_positions=(false, false))
        nu_of(u) = _corridor_element_summary(SVector{3, Float64}(u.sc[1].pos), SVector{3, Float64}(u.sc[1].vel), mu).nu
        function affect!(integrator)
            dot(integrator.u.sc[1].pos, integrator.u.sc[1].vel) > 0.0 || return nothing
            count[] += 1; j = count[]
            d = j == bj ? dv : 0.0; ctrl.delta_v_mps = d; ctrl.source_orbit = j; ctrl.armed = true
            SimulationModel.SimulationCallbacks.schedule_event_driven_thruster_controls!(integrator, 1)
            ctrl.armed = false
            plan = integrator.p.shared_buffers.maneuver_burn_plans[1]; th = integrator.p.args.control_model.control_effectors[1]
            d != 0.0 && !(plan.valid && plan.source_orbit == j && plan.stop_burn_s > plan.start_burn_s) && throw(ErrorException("burn plan rejected at t=$(integrator.t) dv=$d"))
            if j <= bj   # one row per request up to the burn's apoapsis (plan fields of a zero command are whatever the scheduler left: read plan_source_orbit)
                push!(burn, (; j, plan_source_orbit=plan.source_orbit, request_t_s=Float64(integrator.t), nu_request_rad=nu_of(integrator.u), dv_cmd_mps=d, plan_valid=plan.valid, plan_start_s=plan.start_burn_s, plan_stop_s=plan.stop_burn_s,
                    win_start_s=th.start_burn_time[1], win_stop_s=th.stop_burn_time[1], direction_rad=plan.direction_rad, thrust_n=plan.thrust_n, isp_s=plan.isp_s, mass_kg=Float64(integrator.u.sc[1].mass)))
            end
            return nothing
        end
        pre = ContinuousCallback((u, t, i) -> sin(nu_of(u) - target_nu), affect!, nothing; save_positions=(false, false))
        stepcb = twin_step_callback(steps, args, model, pargs, ms, rot, Float64(seed.gate_altitude_m))
        tr = time(); sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=(apsis, pre, stepcb)); run_wall = time() - tr
        logpath = joinpath(dir, "gram_density_perturbation_log.csv")
        rows = stagef_table(sol, args, ms, rot, sc, pwf, logpath); isfile(logpath) && rm(logpath; force=true); isfile(logpath * ".summary.toml") && rm(logpath * ".summary.toml"; force=true)
        write_rows_m(joinpath(dir, "passes.csv"), rows); write_rows_m(joinpath(dir, "events.csv"), events); write_rows_m(joinpath(dir, "steps.csv"), steps); write_rows_m(joinpath(dir, "burn.csv"), burn)
        ve = Float64(sc["vehicle"]["specific_impulse_s"]) * 9.80665; used = ve * log(ms.wet_mass_kg / Float64(sol.u[end].sc[1].mass))
        return (; tag=seed.tag, name=seed.name, dv_cmd_mps=dv, status="ok", retcode=string(sol.retcode), passes=length(rows), n_events=length(events), n_steps=length(steps), n_requests=count[],
            t_end_s=Float64(sol.t[end]), dv_used_mps=used, gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second),
            run_wall_s=run_wall, wall_s=time() - t0, error="")
    catch err
        return (; tag=seed.tag, name=seed.name, dv_cmd_mps=Float64(seed.dv), status="run_error", retcode="", passes=0, n_events=0, n_steps=0, n_requests=0, t_end_s=NaN, dv_used_mps=NaN,
            gram_core_epoch="", run_wall_s=NaN, wall_s=time() - t0, error=replace(first(sprint(showerror, err), 600), "\n" => " "))
    end
end
