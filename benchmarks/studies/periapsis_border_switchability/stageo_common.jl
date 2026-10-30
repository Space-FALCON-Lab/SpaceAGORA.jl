# Stage O: the periapsis-corridor keeper against the flown Mars Odyssey main phase (plan: tasks/2026-10-15_periapsis-region-border-targeting/odyssey_comparison_plan.md, with the 2026-10-24 revision).
# Reuses stagen2_common.jl (included, not edited; it loads the compatibility shim stagen_flight_shim.jl first, then the Stage F/H files, the online Kalman filter and the N2 policy). New code only:
#   (i)   the start orbit from Odyssey's parameters (config [stageo]: Table III apoapsis radius, manifest i / RAAN / omega, epoch at apoapsis, per-band calibrated periapsis), not the navigation state;
#   (ii)  the band view (Odyssey 0.14 to 0.32 or the paper band) applied to the flat per-planet config before any helper reads it;
#   (iii) the stop-at-end-radius callback (first apoapsis at or below ra_target_m) and the per-pass log of rho_peak_kg_m3, rp_m, ra_m (scoring heat rate, plan R1.2);
#   (iv)  the atmosphere-step override of the integrator gate G1;
#   (v)   an optional relative-speed tracker (config track_vrel, default off = N2 exactly as frozen).
# One case function for the fixed (Stage I) and the chance-constrained (N2) keepers. No src/ or examples/ change.
using OrdinaryDiffEq
include(joinpath(@__DIR__, "stagen2_common.jl"))

# ---- pass bookkeeping (pure, selftested) --------------------------------------------------------------------------------------------------------

"""Per-pass bookkeeping. A pass window runs from one apoapsis to the next; its density maximum is sampled along the window and recorded at the closing apoapsis together with that apoapsis radius."""
mutable struct PassBook
    ra_target_m::Float64
    n_peri::Int
    cur_max::Float64
    recs::Vector{NamedTuple}      # (pass, rho_peak_kg_m3, ra_m)
    stopped::Bool
end
PassBook(ra_target_m::Real) = PassBook(Float64(ra_target_m), 0, 0.0, NamedTuple[], false)
book_sample!(pk::PassBook, rho::Real) = (isfinite(rho) && rho > pk.cur_max && (pk.cur_max = Float64(rho)); nothing)
book_periapsis!(pk::PassBook, rho_event::Real) = (pk.n_peri += 1; book_sample!(pk, rho_event); nothing)
"""Close the window at an apoapsis of radius `radius_m`. Returns true when the run must stop (a pass has been flown and the radius is at or below the target)."""
function book_apoapsis!(pk::PassBook, radius_m::Real)::Bool
    pk.n_peri > 0 && (isempty(pk.recs) || pk.recs[end].pass != pk.n_peri) && push!(pk.recs, (; pass=pk.n_peri, rho_peak_kg_m3=pk.cur_max, ra_m=Float64(radius_m)))
    pk.cur_max = 0.0
    stop = pk.n_peri > 0 && radius_m <= pk.ra_target_m
    stop && (pk.stopped = true)
    return stop
end
"""A run that ends mid-pass (horizon) still logs the open pass, with no following apoapsis (ra_m = NaN)."""
book_flush!(pk::PassBook) = (pk.n_peri > 0 && (isempty(pk.recs) || pk.recs[end].pass != pk.n_peri) && push!(pk.recs, (; pass=pk.n_peri, rho_peak_kg_m3=pk.cur_max, ra_m=NaN)); nothing)

"""Inertial vis-viva speed at the periapsis of the osculating orbit through (r, v) (m, m/s), gravitational parameter mu."""
function periapsis_speed(r::AbstractVector, v::AbstractVector, mu::Real)
    rn = sqrt(sum(abs2, r)); vn2 = sum(abs2, v); a = 1 / (2 / rn - vn2 / mu)
    h = (r[2] * v[3] - r[3] * v[2], r[3] * v[1] - r[1] * v[3], r[1] * v[2] - r[2] * v[1]); p = (h[1]^2 + h[2]^2 + h[3]^2) / mu
    e = sqrt(max(0.0, 1 - p / a)); rp = a * (1 - e)
    return sqrt(mu * (2 / rp - 1 / a))
end

# ---- config views ------------------------------------------------------------------------------------------------------------------------------

"""Flat Mars config for one band: [common] + [planets.mars] with the Odyssey start apoapsis radius injected, the band replaced by [stageo.bands.<band>], and (optionally) the atmosphere step overridden.
Returns (meta, sc)."""
function stageo_load(path::AbstractString, band::AbstractString; dt_atm_s::Real=NaN)
    meta = stagef_meta(path); SO = meta["stageo"]; sc = load_stagef("mars"; path=path)
    haskey(SO["bands"], band) || throw(ArgumentError("unknown band $band"))
    sc["orbit"]["apoapsis_radius_m"] = Float64(SO["ra0_m"])
    sc["band"] = merge(Dict{String, Any}(sc["band"]), Dict{String, Any}(SO["bands"][band]))
    isfinite(dt_atm_s) && (sc["integration"]["dt_max_atmosphere_s"] = Float64(dt_atm_s))
    return meta, sc
end

"""Calibration view (as stagei_cal_view, copied rather than included so the Stage F files are not loaded twice): q_low is replaced by q_target^2 / q_high so the helpers' geometric-mean target is exactly q_target."""
function stageo_cal_view(sc)
    v = deepcopy(sc); b = v["band"]; qt = Float64(b["q_target_w_cm2"]); b["q_low_w_cm2"] = qt^2 / Float64(b["q_high_w_cm2"])
    return v
end

# ---- phase cal: one fit-only calibration case per band --------------------------------------------------------------------------------------------

"""Calibration case (a copy of stagei_cal_task on the Odyssey start orbit and the band view); writes <seed.out>/<band>/{guess.csv, nodes.csv, calibration.csv}. Plain-data row; never throws.
The calibrated first-pass periapsis altitude puts pass 1 at the band's target Q (plan R1.1: the sixth orbit element is each keeper's own calibration, not a flown value)."""
function stageo_cal_task(seed)
    t0 = time()
    try
        meta, sc = stageo_load(seed.config_path, seed.band); scc = stageo_cal_view(sc); ms = stagef_setup(sc)
        epoch = parse_epoch(sc["grid"]["epochs_utc"][1]); rot = mars_rot(ms.planet, epoch)
        et0 = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel()); stagef_env!(false, meta)
        model = stagef_gram_model(sc, epoch, stagef_seed("mars" * seed.gram_tag)); pargs = gram_pargs(ms.planet, Float64(sc["atmosphere_gram"]["entry_interface_km"]))
        g = pi_rule_altitude(scc, ms, rot, et0, model, pargs, seed.inc, seed.raan, seed.omega)
        f = gram_exp_fit(sc, model, pargs, g.h_g, g.site)
        temp = Float64(sc["atmosphere"]["temperature_k"]); top = Float64(sc["atmosphere"]["fit_top_altitude_m"])
        pwf = pw_fit(f.pw_h, f.pw_rho, temp, "piecewise"); msc = with_fit(ms, exp_fit(f.rho_g, g.h_g, f.H, top, temp, "cell fit"))
        dir = joinpath(seed.out, seed.band); mkpath(dir)
        seed_alt, recs = gram_calibrate(msc, scc, model, rot, epoch, dir, g.h_g; inc_deg=seed.inc, raan_deg=seed.raan, omega_deg=seed.omega, pwfit=pwf)
        write_rows_m(joinpath(dir, "calibration.csv"), recs)
        write_rows_m(joinpath(dir, "nodes.csv"), [(h_m=a, rho_kg_m3=b) for (a, b) in zip(f.pw_h, f.pw_rho)])
        core_t = model.core.initial_time
        row = (; planet="mars", tag=seed.band, status="ok", h_guess_m=g.h_g, seed_alt_m=seed_alt, vrel_first_mps=recs[end].vrel_mps, q_first_w_cm2=recs[end].q_w_cm2, q_target_w_cm2=recs[end].q_target,
            site_lat_rad=g.site.lat, site_lon_rad=g.site.lon, fit_H_m=f.H, pw_max_abs_lnrho=f.pw_max_abs_lnrho, hloc_at_hg_m=f.hloc_at_hg, n_cal_iter=length(recs),
            gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second), wall_s=time() - t0, error="")
        write_rows_m(joinpath(dir, "guess.csv"), [row]); return row
    catch err
        return (; planet="mars", tag=seed.band, status="cal_error", h_guess_m=NaN, seed_alt_m=NaN, vrel_first_mps=NaN, q_first_w_cm2=NaN, q_target_w_cm2=NaN, site_lat_rad=NaN, site_lon_rad=NaN, fit_H_m=NaN,
            pw_max_abs_lnrho=NaN, hloc_at_hg_m=NaN, n_cal_iter=0, gram_core_epoch="", wall_s=time() - t0, error=replace(first(sprint(showerror, err), 400), "\n" => " "))
    end
end

# ---- callbacks -----------------------------------------------------------------------------------------------------------------------------------

"""A copy of stagef_apsis_callback that also feeds the PassBook: the periapsis event density (mean x perturbation factor) opens the pass maximum, and the apoapsis closes the pass and terminates the run when
its radius is at or below the book's target."""
function stageo_apsis_callback(obs::Vector{NamedTuple}, meas::Vector{NamedTuple}, args, controller, model, pargs, pwfit, pk::PassBook)
    env = args.environment_model; planet = env.planet; et0 = SimulationModel.ephemerides_time_seconds(args.initial_time, env.ephemerides_model); geodetic_of = _geodetic_of(args)
    condition(u, t, integrator) = dot(u.sc[1].pos, u.sc[1].vel)
    function record!(integrator, kind::String)
        integrator.t > 1e-3 || return nothing
        t = Float64(integrator.t); sc1 = integrator.u.sc[1]
        push!(obs, (time_s=t, kind=kind, radius_m=norm(sc1.pos), physical_altitude_m=geodetic_of(t, sc1.pos)[1]))
        if kind == "apoapsis"
            book_apoapsis!(pk, norm(sc1.pos)) && terminate!(integrator)
            return nothing
        end
        controller.periapsis_count += 1
        lpi = SimulationModel.planet_frame_lpi(planet, et0 + t, env.ephemerides_model); rp = lpi * SVector{3, Float64}(sc1.pos); alt, lat = geodetic_altitude_latitude(rp, planet)
        lon = SimulationModel.SimulationCallbacks.rtolatlong(SVector{3, Float64}(rp), planet)[3]
        rho = gram_rho(model, pargs, alt, lat, lon, t); r = stagef_pass_r(integrator.p, t)
        push!(meas, (; t_s=t, alt_m=alt, lat_rad=lat, lon_rad=lon, rho_mean=rho, r=r, y=log(rho * r / corridor_density(pwfit, alt))))
        book_periapsis!(pk, rho * r)
        return nothing
    end
    return ContinuousCallback(condition, i -> record!(i, "periapsis"), i -> record!(i, "apoapsis"); save_positions=(false, false))
end

"""Running maximum of the truth density (GRAM mean x the held perturbation factor) at every accepted step below `gate_alt_m` altitude; the periapsis event density is added by the apsis callback."""
function stageo_peak_callback(args, model, pargs, pk::PassBook, gate_alt_m::Real)
    env = args.environment_model; planet = env.planet; et0 = SimulationModel.ephemerides_time_seconds(args.initial_time, env.ephemerides_model); rp_e = Float64(planet.Rp_e); gate = Float64(gate_alt_m)
    condition(u, t, integrator) = norm(u.sc[1].pos) - rp_e < gate
    function affect!(integrator)
        t = Float64(integrator.t); sc1 = integrator.u.sc[1]
        lpi = SimulationModel.planet_frame_lpi(planet, et0 + t, env.ephemerides_model); rp = lpi * SVector{3, Float64}(sc1.pos); alt, lat = geodetic_altitude_latitude(rp, planet)
        lon = SimulationModel.SimulationCallbacks.rtolatlong(SVector{3, Float64}(rp), planet)[3]
        book_sample!(pk, gram_rho(model, pargs, alt, lat, lon, t) * stagef_pass_r(integrator.p, t))
        return nothing
    end
    return DiscreteCallback(condition, affect!; save_positions=(false, false))
end

"""Guidance wrapper around the Stage I fixed policy or the N2 policy. With `track` it rescales the policy's first-pass relative speed (which converts Q edges to density edges) by the ratio of the current
osculating periapsis speed to the start orbit's (rotation term held; second order at this inclination). With track = false it only delegates, so the flight is the unmodified N2 / Stage I keeper."""
mutable struct StageOPolicy <: SimulationModel.AbstractGuidanceModel
    pol::Any                      # StageN2Policy or StageFPolicy
    f::StageFPolicy
    mu::Float64
    track::Bool
    vp0::Float64                  # periapsis speed of the first call's state (the calibrated start orbit); NaN until then
    vrel0::Float64                # vrel_first at the start
end
function SimulationModel.GuidanceHooks.calcGuidanceEffect!(w::StageOPolicy, u, p, t::Float64, sat_idx::Int64)
    if w.track && sat_idx >= 1
        vp = periapsis_speed(u.sc[sat_idx].pos, u.sc[sat_idx].vel, w.mu)
        isnan(w.vp0) && (w.vp0 = vp)
        w.f.vrel_first = w.vrel0 * vp / w.vp0
    end
    return SimulationModel.GuidanceHooks.calcGuidanceEffect!(w.pol, u, p, t, sat_idx)
end

# ---- phase run: one flight --------------------------------------------------------------------------------------------------------------------------

"""Arm code (run directory prefix): <ARM>_fixed<margin m> for the fixed keepers, <ARM>_<N2 arm code> for the chance-constrained ones; the analysis reads the arm and band from it."""
stageo_arm_code(key::AbstractString, ctrl::AbstractString, margin_m::Real, eps::Real, t_mult::Real, aim::Real, w_m::Real) =
    key * "_" * (ctrl == "fixed" ? "fixed" * lpad(string(round(Int, margin_m)), 4, '0') : stagen2_arm_code(eps, t_mult, aim, w_m))

"""One Stage O flight. seed fields: config_path, band, tag (GRAM seed string tag), inc, raan, omega, arm (key), ctrl ("fixed" | "cc"), code, eps, w_m, aim, t_mult, kr, kq, km0, kv0, run_id, out, cal_root,
horizon_days, dt_atm_s (NaN = config). Plain-data row (the Stage N2 columns plus stop_reason); never throws."""
function stageo_run_task(seed)
    t0 = time()
    try
        meta, sc = stageo_load(seed.config_path, seed.band; dt_atm_s=seed.dt_atm_s); SO = meta["stageo"]; ms0 = stagef_setup(sc)
        epoch = parse_epoch(sc["grid"]["epochs_utc"][1]); rot = mars_rot(ms0.planet, epoch)
        cdir = joinpath(seed.cal_root, seed.band); G = CSV.read(joinpath(cdir, "guess.csv"), DataFrame)[1, :]; N = CSV.read(joinpath(cdir, "nodes.csv"), DataFrame)
        pwf = pw_fit(N.h_m, N.rho_kg_m3, Float64(sc["atmosphere"]["temperature_k"]), "piecewise"); ms = with_fit(ms0, pwf)
        stagef_env!(true, meta)
        model = stagef_gram_model(sc, epoch, stagef_seed("mars" * seed.tag)); core_t = model.core.initial_time
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
        w2m = Float64(sc["band"]["w_cm2_to_w_m2"]); fixed = seed.ctrl == "fixed"
        fpol = StageFPolicy(inner, maker, fixed ? :fixed : :adaptive, pwf, Float64(G.vrel_first_mps), Float64(sc["band"]["q_low_w_cm2"]) * w2m, Float64(sc["band"]["q_high_w_cm2"]) * w2m,
            margin, meas, "", 0.0, 0.0, Float64(meta["analysis"]["adapt_min_band_m"]), NamedTuple[])
        pol = fixed ? fpol : StageN2Policy(fpol, Float64(seed.kr), Float64(seed.kq), Float64(seed.km0), Float64(seed.kv0), Float64(seed.t_mult), Float64(seed.w_m), Float64(seed.aim))
        wrap = StageOPolicy(pol, fpol, Float64(ms.harm.gm_m3s2), Bool(SO["track_vrel"]), NaN, Float64(G.vrel_first_mps))
        v = cfg.values["vehicle"]; thruster = BaseThrusterModel([Float64(v["thrust_n"])], [0.0], [0.0], [-1.0], [-1.0], [Float64(v["specific_impulse_s"])])
        args = SimulationConfiguration(file_paths=base.file_paths, simulation_settings=base.simulation_settings, mission_configuration=base.mission_configuration,
            environment_model=base.environment_model, dynamics_model=base.dynamics_model, guidance_model=GuidanceModel((wrap,), [1.0e12]), navigation_model=base.navigation_model,
            control_model=ControlModel((thruster,), [1.0e12]), initial_time=base.initial_time, integration_tolerances=base.integration_tolerances,
            solver_config=SolverConfig(solver_mode=base.solver_config.solver_mode, maxiters=Int(SO["solver_maxiters"])))   # DiffEq's default 1e5 steps ends a 1 s main-phase flight at day ~82 (gate 1, 2026-10-24)
        obs = NamedTuple[]; pk = PassBook(Float64(SO["ra_target_m"]))
        extra = (stageo_apsis_callback(obs, meas, args, inner, model, pargs, pwf, pk), stageo_peak_callback(args, model, pargs, pk, Float64(SO["peak_gate_altitude_m"])), _preapoapsis_guidance_callback(inner))
        tr = time(); sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=extra); run_wall = time() - tr
        book_flush!(pk); by = Dict(r.pass => r for r in pk.recs)
        logpath = joinpath(dir, "gram_density_perturbation_log.csv")
        rows0 = stagef_table(sol, args, ms, rot, sc, pwf, logpath); isfile(logpath) && rm(logpath; force=true); isfile(logpath * ".summary.toml") && rm(logpath * ".summary.toml"; force=true)
        rows = [merge(r, (; rho_peak_kg_m3=haskey(by, r.pass) ? by[r.pass].rho_peak_kg_m3 : NaN, rp_m=Float64(r.radius_m), ra_m=haskey(by, r.pass) ? by[r.pass].ra_m : NaN)) for r in rows0]
        write_rows_m(joinpath(dir, "passes.csv"), rows)
        ve = Float64(sc["vehicle"]["specific_impulse_s"]) * 9.80665; dv = ve * log(ms.wet_mass_kg / Float64(sol.u[end].sc[1].mass)); burns = length(inner.burn_plans)
        stopped = !isempty(inner.failures); failures = join(inner.failures, " | ")
        write_rows_m(joinpath(dir, "decisions.csv"), fpol.edge_log); write_rows_m(joinpath(dir, "measured.csv"), meas)
        ret = string(sol.retcode)
        stop_reason = pk.stopped ? "end_radius" : stopped ? "controller_stopped" : ret == "Success" ? "horizon" : "retcode_" * ret
        return (; planet="mars", tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl=seed.ctrl, code=seed.code, eps=Float64(seed.eps), w_m=Float64(seed.w_m), t_mult=Float64(seed.t_mult), aim=Float64(seed.aim), status="ok", retcode=ret,
            passes=length(rows), t_end_days=sol.t[end] / SECONDS_PER_DAY, horizon_days=seed.horizon_days, completed=pk.stopped && !stopped, controller_stopped=stopped, controller_failures=failures, total_dv_mps=dv, burns=burns,
            collapsed_decisions=count(e -> e.collapsed, fpol.edge_log), gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second),
            run_wall_s=run_wall, wall_s=time() - t0, error="", stop_reason=stop_reason)
    catch err
        return (; planet="mars", tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl=seed.ctrl, code=seed.code, eps=NaN, w_m=NaN, t_mult=NaN, aim=NaN,
            status="run_error", retcode="", passes=0, t_end_days=NaN, horizon_days=seed.horizon_days, completed=false, controller_stopped=false, controller_failures="", total_dv_mps=NaN, burns=0,
            collapsed_decisions=0, gram_core_epoch="", run_wall_s=NaN, wall_s=time() - t0, error=replace(first(sprint(showerror, err), 500), "\n" => " "), stop_reason="error")
    end
end
