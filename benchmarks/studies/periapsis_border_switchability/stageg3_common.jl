# Stage G3 (WP-G3, gate V-G): free-drift validation flights at other geometries (plan: tasks/2026-10-15_periapsis-region-border-targeting/analytic_framework_and_transfer_plan.md, R4.3, R4.4, R4.10, R5.1).
# Reuses stageo_common.jl (included, not edited; it loads stagen2_common.jl and through it the Stage F/H files, the compatibility shim and the estimator). New code only:
#   (i)   the start orbit per geometry from parameters (config [[stageg3.geometry]]: elements, epoch = base epoch + offset, apoapsis radius and target from [stageo]); a per-geometry fit-only calibration;
#   (ii)  U15 (R4.4, R5.1): the atmosphere-relative speed at the osculating periapsis, periapsis_vrel (below). Where it enters:
#         DECISION: StageGPolicy rescales the keeper's first-pass relative speed (which converts the Q edges to density edges) by V_rel(now) / V_rel(first call), V_rel from the osculating periapsis state and the
#                   planet's rotation vector (Stage O's track_vrel did the same with the inertial vis-viva speed; here the rotation term is included). decision_speed = "first_pass" turns it off.
#         SCORING:  every pass logs vrel_mps (the exact relative speed at the periapsis event, already a column of passes.csv: it is what q_w_cm2 and q_eff_w_cm2 use) and vp_inertial_mps; the analysis scores
#                   Q_s on V_rel and reports the inertial form (Stage O's qs_of) for comparison. The safety stop below uses the same V_rel form.
#   (iii) the zero-burn flight (the keeper policy is never armed: the pre-apoapsis guidance callback is not installed, so no maneuver is ever requested) with the safety terminate
#         (pass-maximum Q_s on V_rel above q_stop, or a periapsis more than n_hloc x H_local below the first-pass periapsis altitude);
#   (iv)  the per-pass log of rho_peak_kg_m3, rp_m, ra_m (as stageo) plus vp_inertial_mps.
# No src/ or examples/ change.
using Dates, TOML, JSON, SHA
include(joinpath(@__DIR__, "stageo_common.jl"))

# ---- U15: periapsis geometry and the atmosphere-relative speed (pure, selftested) ---------------------------------------------------------------

"""Periapsis position and velocity vectors of the osculating orbit through (r, v) (same frame as the inputs; mu the gravitational parameter), plus rp, Vp, the orbit normal and the along-track unit vector at periapsis."""
function periapsis_vectors(r::AbstractVector, v::AbstractVector, mu::Real)
    rn = sqrt(sum(abs2, r)); h = cross(r, v); hn = sqrt(sum(abs2, h))
    ev = ((sum(abs2, v) - mu / rn) .* r .- dot(r, v) .* v) ./ mu; e = sqrt(sum(abs2, ev))
    e > 1.0e-10 || throw(ArgumentError("circular osculating orbit: periapsis undefined"))
    ehat = ev ./ e; hhat = h ./ hn; that = cross(hhat, ehat); rp = hn^2 / mu / (1 + e); vp = hn / rp
    return (; rvec=rp .* ehat, vvec=vp .* that, rp, vp, hhat, that, e)
end

"""Atmosphere-relative speed at the periapsis of the osculating orbit through (r, v): |v_p - omega x r_p| with `omega` the planet's rotation vector in the SAME frame as r, v (the atmosphere co-rotates)."""
periapsis_vrel(r::AbstractVector, v::AbstractVector, mu::Real, omega::AbstractVector) = (pv = periapsis_vectors(r, v, mu); sqrt(sum(abs2, pv.vvec .- cross(omega, pv.rvec))))

"""Along-track component at periapsis of the rotation velocity omega x r_p: equals omega r_p cos(i) (R4.4 identity), i the inclination to the rotation axis."""
along_track_rotation(r::AbstractVector, v::AbstractVector, mu::Real, omega::AbstractVector) = (pv = periapsis_vectors(r, v, mu); dot(cross(omega, pv.rvec), pv.that))

"""Scoring heat rate on the relative speed, W/cm2: Q_s = 0.5 rho_peak V_rel^3 / 1e4 (SI inputs). The inertial form of Stage O (stageo_analysis.qs_of) uses the vis-viva speed instead."""
qs_relative_w_cm2(rho_peak::Real, vrel::Real) = 0.5 * rho_peak * vrel^3 / 1.0e4

# ---- safety stop (pure state machine, selftested) ----------------------------------------------------------------------------------------------

"""Zero-burn safety stop. `enabled` false never fires. q_stop in W/cm2; n_hloc x h_local_m is the allowed periapsis descent below the first periapsis altitude (alt0_m, NaN until the first periapsis)."""
mutable struct SafetyStop
    enabled::Bool
    q_stop_w_cm2::Float64
    n_hloc::Float64
    h_local_m::Float64
    alt0_m::Float64
    vrel_p::Float64               # relative speed at the latest periapsis event
    reason::String                # "" until a stop fires, then "safety_q" or "safety_alt"
end
SafetyStop(enabled::Bool, q_stop::Real, n_hloc::Real, h_local_m::Real) = SafetyStop(enabled, Float64(q_stop), Float64(n_hloc), Float64(h_local_m), NaN, NaN, "")
"""Periapsis event: remember the relative speed; the first periapsis altitude is the reference; returns true (and sets the reason) when the altitude is more than n_hloc x H_local below it."""
function safety_periapsis!(ss::SafetyStop, alt_m::Real, vrel::Real)::Bool
    ss.vrel_p = Float64(vrel); isnan(ss.alt0_m) && (ss.alt0_m = Float64(alt_m))
    fire = ss.enabled && Float64(alt_m) < ss.alt0_m - ss.n_hloc * ss.h_local_m
    fire && isempty(ss.reason) && (ss.reason = "safety_alt")
    return fire
end
"""Closing apoapsis of a pass with density maximum `rho_peak` (kg/m3): returns true (and sets the reason) when the pass-maximum Q_s on V_rel exceeds the threshold."""
function safety_apoapsis!(ss::SafetyStop, rho_peak::Real)::Bool
    fire = ss.enabled && isfinite(ss.vrel_p) && qs_relative_w_cm2(rho_peak, ss.vrel_p) > ss.q_stop_w_cm2
    fire && isempty(ss.reason) && (ss.reason = "safety_q")
    return fire
end

# ---- geometry table and config views ----------------------------------------------------------------------------------------------------------

"""Epoch string of a geometry: the base epoch plus `off_d` days (millisecond resolution, ISO, as parse_epoch reads it)."""
stageg3_epoch_string(base::AbstractString, off_d::Real) = Dates.format(DateTime(base) + Millisecond(round(Int, Float64(off_d) * 86400.0e3)), "yyyy-mm-ddTHH:MM:SS")

const _G3META = Dict{String, Any}()
"""Task directory (env SPACEAGORA_BORDER_TASK_DIR, else found by walking up from this file)."""
function stageg3_task_dir()
    d = get(ENV, "SPACEAGORA_BORDER_TASK_DIR", ""); isempty(d) || return d
    p = @__DIR__
    while true
        c = joinpath(p, "tasks", "2026-10-15_periapsis-region-border-targeting"); isdir(c) && return c
        q = dirname(p); q == p && error("task directory not found (set SPACEAGORA_BORDER_TASK_DIR)"); p = q
    end
end

"""Config dict with the geometry list injected from the final geometry table (config [stageg3.table]: path, sha256, include_g3). The table's sha256 is verified at load: a mismatch is an ERROR.
meta["stageg3"]["geometry"] = the default-set geometries (plus G3 only when include_g3), each a Dict with id, role, inc_deg, omega_deg, raan_deg, epoch_off_d, epoch_utc, field_id; meta["stageg3"]["table_hash"], ["table_start"]."""
function stageg3_meta(path::AbstractString)
    get!(_G3META, abspath(path)) do
        m = TOML.parsefile(path); t = m["stageg3"]["table"]; fp = joinpath(stageg3_task_dir(), String(t["path"])); raw = read(fp); h = bytes2hex(sha256(raw))
        h == String(t["sha256"]) || error("geometry table sha256 mismatch: $fp is $h, config has $(t["sha256"])")
        T = JSON.parse(String(raw)); G = Dict{String, Any}[]
        for g in T["geometries"]
            (g["in_default_set"] || (Bool(t["include_g3"]) && g["role"] == "G3")) || continue
            push!(G, Dict{String, Any}("id" => g["id"], "role" => g["role"], "inc_deg" => Float64(g["inclination_deg"]), "omega_deg" => Float64(g["omega0_deg"]), "raan_deg" => Float64(g["raan0_deg"]),
                "epoch_off_d" => Float64(g["epoch_offset_days"]), "epoch_utc" => g["epoch_utc"], "field_id" => g["id"]))
        end
        m["stageg3"]["geometry"] = G; m["stageg3"]["table_start"] = T["start"]; m["stageg3"]["table_hash"] = h; m
    end
end

"""Geometry record by id (in the loaded set)."""
function stageg3_geometry(meta, gid::AbstractString)
    i = findfirst(g -> String(g["id"]) == gid, meta["stageg3"]["geometry"]); i === nothing && throw(ArgumentError("unknown geometry $gid"))
    return meta["stageg3"]["geometry"][i]
end

"""GRAM core epoch string for an ISO epoch, as the summaries record it (no zero padding, seconds as printed up to the decimal): the GRAM epoch must follow each geometry's epoch."""
function stageg3_gram_prefix(ep::AbstractString)
    d = DateTime(ep); return string(Dates.year(d), "-", Dates.month(d), "-", Dates.day(d), " ", Dates.hour(d), ":", Dates.minute(d), ":", Dates.second(d))
end
"""Runtime guard: error unless the GRAM model's core epoch is the geometry's epoch."""
function stageg3_gram_guard(core_t, ep::AbstractString)
    got = string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second); want = stageg3_gram_prefix(ep)
    startswith(got, want) || error("GRAM core epoch $got does not follow the geometry epoch $ep (expected $want)")
    return got
end

"""Flat Mars config of one geometry: stageo_load (Odyssey band, apoapsis radius, optional atmosphere step) with the epoch replaced by the table's epoch (checked against base epoch + offset, and the table's start
against [stageo]). Returns (meta, sc, geom, epoch_string)."""
function stageg3_load(path::AbstractString, gid::AbstractString; dt_atm_s::Real=NaN)
    meta = stageg3_meta(path); _, sc = stageo_load(path, String(meta["stageg3"]["band"]); dt_atm_s=dt_atm_s); g = stageg3_geometry(meta, gid); ep = String(g["epoch_utc"])
    ep == stageg3_epoch_string(String(meta["stageg3"]["base_epoch_utc"]), Float64(g["epoch_off_d"])) || error("table epoch $ep of $gid differs from base epoch + offset")
    ts = meta["stageg3"]["table_start"]; (Float64(ts["ra0_m"]) == Float64(meta["stageo"]["ra0_m"]) && Float64(ts["ra_target_m"]) == Float64(meta["stageo"]["ra_target_m"]) && Float64(ts["true_anomaly_deg"]) == Float64(sc["orbit"]["true_anomaly_deg"])) || error("table start differs from [stageo] / orbit config")
    sc["grid"]["epochs_utc"] = [ep]
    return meta, sc, g, ep
end

# ---- phase cal: one fit-only calibration case per geometry ----------------------------------------------------------------------------------

"""Calibration case of one geometry (stageo_cal_task on the geometry's elements and epoch); writes <seed.out>/<gid>/{guess.csv, nodes.csv, calibration.csv}. Plain-data row; never throws."""
function stageg3_cal_task(seed)
    t0 = time()
    try
        meta, sc, g, ep = stageg3_load(seed.config_path, seed.gid); scc = stageo_cal_view(sc); ms = stagef_setup(sc)
        epoch = parse_epoch(sc["grid"]["epochs_utc"][1]); rot = mars_rot(ms.planet, epoch)
        et0 = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel()); stagef_env!(false, meta)
        model = stagef_gram_model(sc, epoch, stagef_seed("mars" * seed.gram_tag)); pargs = gram_pargs(ms.planet, Float64(sc["atmosphere_gram"]["entry_interface_km"]))
        inc, raan, omega = Float64(g["inc_deg"]), Float64(g["raan_deg"]), Float64(g["omega_deg"])
        gg = pi_rule_altitude(scc, ms, rot, et0, model, pargs, inc, raan, omega)
        f = gram_exp_fit(sc, model, pargs, gg.h_g, gg.site)
        temp = Float64(sc["atmosphere"]["temperature_k"]); top = Float64(sc["atmosphere"]["fit_top_altitude_m"])
        pwf = pw_fit(f.pw_h, f.pw_rho, temp, "piecewise"); msc = with_fit(ms, exp_fit(f.rho_g, gg.h_g, f.H, top, temp, "cell fit"))
        dir = joinpath(seed.out, seed.gid); mkpath(dir)
        seed_alt, recs = gram_calibrate(msc, scc, model, rot, epoch, dir, gg.h_g; inc_deg=inc, raan_deg=raan, omega_deg=omega, pwfit=pwf)
        write_rows_m(joinpath(dir, "calibration.csv"), recs)
        write_rows_m(joinpath(dir, "nodes.csv"), [(h_m=a, rho_kg_m3=b) for (a, b) in zip(f.pw_h, f.pw_rho)])
        core_t = model.core.initial_time; stageg3_gram_guard(core_t, ep)
        row = (; planet="mars", tag=seed.gid, status="ok", h_guess_m=gg.h_g, seed_alt_m=seed_alt, vrel_first_mps=recs[end].vrel_mps, q_first_w_cm2=recs[end].q_w_cm2, q_target_w_cm2=recs[end].q_target,
            site_lat_rad=gg.site.lat, site_lon_rad=gg.site.lon, fit_H_m=f.H, pw_max_abs_lnrho=f.pw_max_abs_lnrho, hloc_at_hg_m=f.hloc_at_hg, n_cal_iter=length(recs),
            gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second), wall_s=time() - t0, error="", epoch_utc=ep)
        write_rows_m(joinpath(dir, "guess.csv"), [row]); return row
    catch err
        return (; planet="mars", tag=seed.gid, status="cal_error", h_guess_m=NaN, seed_alt_m=NaN, vrel_first_mps=NaN, q_first_w_cm2=NaN, q_target_w_cm2=NaN, site_lat_rad=NaN, site_lon_rad=NaN, fit_H_m=NaN,
            pw_max_abs_lnrho=NaN, hloc_at_hg_m=NaN, n_cal_iter=0, gram_core_epoch="", wall_s=time() - t0, error=replace(first(sprint(showerror, err), 400), "\n" => " "), epoch_utc="")
    end
end

# ---- callbacks --------------------------------------------------------------------------------------------------------------------------------

"""stageo_apsis_callback with (a) the exact relative and inertial periapsis speeds logged per pass (vlog rows (pass, vrel_p_mps, vp_inertial_mps)) and (b) the safety stop: the periapsis altitude test at the
periapsis event and the pass-maximum Q_s test at the closing apoapsis (the pass density maximum is read before the window is closed). Terminates the run on the end radius or a safety stop."""
function stageg3_apsis_callback(obs::Vector{NamedTuple}, meas::Vector{NamedTuple}, vlog::Vector{NamedTuple}, args, controller, model, pargs, pwfit, pk::PassBook, ss::SafetyStop)
    env = args.environment_model; planet = env.planet; et0 = SimulationModel.ephemerides_time_seconds(args.initial_time, env.ephemerides_model); geodetic_of = _geodetic_of(args)
    condition(u, t, integrator) = dot(u.sc[1].pos, u.sc[1].vel)
    function record!(integrator, kind::String)
        integrator.t > 1e-3 || return nothing
        t = Float64(integrator.t); sc1 = integrator.u.sc[1]
        push!(obs, (time_s=t, kind=kind, radius_m=norm(sc1.pos), physical_altitude_m=geodetic_of(t, sc1.pos)[1]))
        if kind == "apoapsis"
            rho_pk = pk.cur_max
            book_apoapsis!(pk, norm(sc1.pos)) && terminate!(integrator)
            !pk.stopped && pk.n_peri > 0 && safety_apoapsis!(ss, rho_pk) && terminate!(integrator)
            return nothing
        end
        controller.periapsis_count += 1
        lpi = SimulationModel.planet_frame_lpi(planet, et0 + t, env.ephemerides_model); rp = lpi * SVector{3, Float64}(sc1.pos); alt, lat = geodetic_altitude_latitude(rp, planet)
        lon = SimulationModel.SimulationCallbacks.rtolatlong(SVector{3, Float64}(rp), planet)[3]
        vrel = norm(lpi * SVector{3, Float64}(sc1.vel) - cross(planet.ω, rp))
        rho = gram_rho(model, pargs, alt, lat, lon, t); r = stagef_pass_r(integrator.p, t)
        push!(meas, (; t_s=t, alt_m=alt, lat_rad=lat, lon_rad=lon, rho_mean=rho, r=r, y=log(rho * r / corridor_density(pwfit, alt))))
        book_periapsis!(pk, rho * r); push!(vlog, (; pass=pk.n_peri, vrel_p_mps=vrel, vp_inertial_mps=norm(sc1.vel)))
        safety_periapsis!(ss, alt, vrel) && terminate!(integrator)
        return nothing
    end
    return ContinuousCallback(condition, i -> record!(i, "periapsis"), i -> record!(i, "apoapsis"); save_positions=(false, false))
end

# ---- guidance wrapper (U15 decision quantity) ---------------------------------------------------------------------------------------------------

"""Guidance wrapper around the Stage I fixed keeper. mode :relative (U15) sets the keeper's first-pass relative speed to V_rel(first call) x V_rel(osculating periapsis now) / V_rel(osculating periapsis at the first call),
the osculating V_rel from periapsis_vrel with the planet's rotation vector in the inertial frame (omega_fn(t) = lpi(t)' * planet.omega); :first_pass only delegates (the keeper as Stage O track_vrel = false).
Documented limit: the speed is that of the osculating periapsis BEFORE the decision's burn (a burn of 0.1 m/s moves the periapsis speed by about 0.01 m/s, far below the 1% effect tracked)."""
mutable struct StageGPolicy <: SimulationModel.AbstractGuidanceModel
    f::StageFPolicy
    mu::Float64
    mode::Symbol
    omega_fn::Any                 # t_s -> rotation vector in the inertial frame (3-vector)
    vrel0::Float64                # first-pass relative speed from the calibration
    vr0::Float64                  # osculating V_rel at the first call; NaN until then
end
function SimulationModel.GuidanceHooks.calcGuidanceEffect!(w::StageGPolicy, u, p, t::Float64, sat_idx::Int64)
    if w.mode === :relative && sat_idx >= 1
        vr = periapsis_vrel(u.sc[sat_idx].pos, u.sc[sat_idx].vel, w.mu, w.omega_fn(t))
        isnan(w.vr0) && (w.vr0 = vr)
        w.f.vrel_first = w.vrel0 * vr / w.vr0
    end
    return SimulationModel.GuidanceHooks.calcGuidanceEffect!(w.f, u, p, t, sat_idx)
end

# ---- phase run: one flight --------------------------------------------------------------------------------------------------------------------------

"""Arm code (run directory component): FO_fixed<margin m> for the Stage I fixed keeper, ZB_zeroburn for the zero-burn flight."""
stageg3_arm_code(keeper::AbstractString, margin_m::Real) = keeper == "fixed" ? "FO_fixed" * lpad(string(round(Int, margin_m)), 4, '0') : keeper == "zero_burn" ? "ZB_zeroburn" : throw(ArgumentError("unknown keeper $keeper"))

"""One gate flight. seed fields: config_path, gid, role, keeper ("fixed" | "zero_burn"), tag (GRAM seed string tag), arm, code, run_id, out, cal_root, horizon_days, dt_atm_s (NaN = config), smoke_force (lower the Q safety
threshold to the config's smoke value), decision_speed ("relative" | "first_pass"). Plain-data row (the Stage O columns plus geometry, keeper and safety columns); never throws."""
function stageg3_run_task(seed)
    t0 = time()
    try
        meta, sc, geom, ep = stageg3_load(seed.config_path, seed.gid; dt_atm_s=seed.dt_atm_s); G3 = meta["stageg3"]; SO = meta["stageo"]; ms0 = stagef_setup(sc)
        inc, raan, omega = Float64(geom["inc_deg"]), Float64(geom["raan_deg"]), Float64(geom["omega_deg"])
        epoch = parse_epoch(ep); rot = mars_rot(ms0.planet, epoch); zero = seed.keeper == "zero_burn"
        cdir = joinpath(seed.cal_root, seed.gid); G = CSV.read(joinpath(cdir, "guess.csv"), DataFrame)[1, :]; N = CSV.read(joinpath(cdir, "nodes.csv"), DataFrame)
        pwf = pw_fit(N.h_m, N.rho_kg_m3, Float64(sc["atmosphere"]["temperature_k"]), "piecewise"); ms = with_fit(ms0, pwf)
        stagef_env!(true, meta)
        model = stagef_gram_model(sc, epoch, stagef_seed("mars" * seed.tag)); core_t = model.core.initial_time; stageg3_gram_guard(core_t, ep)
        pargs = gram_pargs(ms.planet, Float64(sc["atmosphere_gram"]["entry_interface_km"]))
        dir = joinpath(seed.out, seed.run_id); mkpath(dir); horizon = Float64(seed.horizon_days) * SECONDS_PER_DAY
        state0 = mars_initial_state(ms, rot, Float64(G.seed_alt_m), inc, raan, omega)
        base = gram_args(ms, sc, model, epoch, horizon, dir; state=state0)
        meas = NamedTuple[]; vlog = NamedTuple[]
        margin = Float64(G3["keeper_margin_m"]); cfg, hq, hl = stagef_corridor_config(ms, Float64(G.vrel_first_mps), margin); c = cfg.values["comparison"]
        inner = CorridorGuidancePolicy(cfg, :favorable, Symbol(meta["stagef"]["policy"]), Float64(c["max_burn_delta_v_mps"]), Float64(c["max_burn_duration_s"]),
            Float64(c["total_delta_v_budget_mps"]), Float64(c["near_passive_deadband_m"]), Float64(c["target_periapsis_altitude_m"]), deg2rad(Float64(c["planning_lead_deg"])),
            0, false, nothing, NamedTuple[], NamedTuple[], String[])
        # zero-burn flights never plan: no estimator context, and the pre-apoapsis guidance callback is not installed, so the keeper is never armed
        maker = zero ? ((p, r, v, mass, t) -> error("zero-burn flight must never plan a maneuver")) : begin
            st = stagef_est_state(ms, epoch); et0e = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel())
            (p, r, v, mass, t) -> (dv -> stagef_predict_alt(st, r, v, Float64(dv), et0e + t))
        end
        w2m = Float64(sc["band"]["w_cm2_to_w_m2"])
        fpol = StageFPolicy(inner, maker, :fixed, pwf, Float64(G.vrel_first_mps), Float64(sc["band"]["q_low_w_cm2"]) * w2m, Float64(sc["band"]["q_high_w_cm2"]) * w2m,
            margin, meas, "", 0.0, 0.0, Float64(meta["analysis"]["adapt_min_band_m"]), NamedTuple[])
        ds = Symbol(seed.decision_speed); ds in (:relative, :first_pass) || error("decision_speed must be relative or first_pass")
        et0p = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel()); planet = ms.planet
        omega_fn = t -> (SimulationModel.planet_frame_lpi(planet, et0p + t, SpiceEphemeridesModel()))' * SVector{3, Float64}(planet.ω)
        wrap = StageGPolicy(fpol, Float64(ms.harm.gm_m3s2), ds, omega_fn, Float64(G.vrel_first_mps), NaN)
        v = cfg.values["vehicle"]; thruster = BaseThrusterModel([Float64(v["thrust_n"])], [0.0], [0.0], [-1.0], [-1.0], [Float64(v["specific_impulse_s"])])
        args = SimulationConfiguration(file_paths=base.file_paths, simulation_settings=base.simulation_settings, mission_configuration=base.mission_configuration,
            environment_model=base.environment_model, dynamics_model=base.dynamics_model, guidance_model=GuidanceModel((wrap,), [1.0e12]), navigation_model=base.navigation_model,
            control_model=ControlModel((thruster,), [1.0e12]), initial_time=base.initial_time, integration_tolerances=base.integration_tolerances,
            solver_config=SolverConfig(solver_mode=base.solver_config.solver_mode, maxiters=Int(SO["solver_maxiters"])))
        S3 = G3["safety"]; qstop = seed.smoke_force ? Float64(G3["smoke_forced_safety_q_w_cm2"]) : Float64(S3["q_stop_w_cm2"])
        ss = SafetyStop(zero, qstop, Float64(S3["n_hloc"]), Float64(G.hloc_at_hg_m))
        obs = NamedTuple[]; pk = PassBook(Float64(SO["ra_target_m"]))
        extra = (stageg3_apsis_callback(obs, meas, vlog, args, inner, model, pargs, pwf, pk, ss), stageo_peak_callback(args, model, pargs, pk, Float64(SO["peak_gate_altitude_m"])))
        zero || (extra = (extra..., _preapoapsis_guidance_callback(inner)))
        tr = time(); sol = run_simulation(args; return_solution=true, isolate_state=false, extra_callbacks=extra); run_wall = time() - tr
        book_flush!(pk); by = Dict(r.pass => r for r in pk.recs); vl = Dict(r.pass => r for r in vlog)
        logpath = joinpath(dir, "gram_density_perturbation_log.csv")
        rows0 = stagef_table(sol, args, ms, rot, sc, pwf, logpath); isfile(logpath) && rm(logpath; force=true); isfile(logpath * ".summary.toml") && rm(logpath * ".summary.toml"; force=true)
        rows = [merge(r, (; rho_peak_kg_m3=haskey(by, r.pass) ? by[r.pass].rho_peak_kg_m3 : NaN, rp_m=Float64(r.radius_m), ra_m=haskey(by, r.pass) ? by[r.pass].ra_m : NaN,
                           vp_inertial_mps=haskey(vl, r.pass) ? vl[r.pass].vp_inertial_mps : NaN)) for r in rows0]
        write_rows_m(joinpath(dir, "passes.csv"), rows)
        ve = Float64(sc["vehicle"]["specific_impulse_s"]) * 9.80665; dv = ve * log(ms.wet_mass_kg / Float64(sol.u[end].sc[1].mass)); burns = length(inner.burn_plans)
        stopped = !isempty(inner.failures); failures = join(inner.failures, " | ")
        write_rows_m(joinpath(dir, "decisions.csv"), fpol.edge_log); write_rows_m(joinpath(dir, "measured.csv"), meas)
        ret = string(sol.retcode)
        stop_reason = pk.stopped ? "end_radius" : !isempty(ss.reason) ? ss.reason : stopped ? "controller_stopped" : ret == "Success" ? "horizon" : "retcode_" * ret
        return (; planet="mars", tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl=seed.keeper, code=seed.code, status="ok", retcode=ret,
            passes=length(rows), t_end_days=sol.t[end] / SECONDS_PER_DAY, horizon_days=seed.horizon_days, completed=pk.stopped && !stopped, controller_stopped=stopped, controller_failures=failures, total_dv_mps=dv, burns=burns,
            collapsed_decisions=count(e -> e.collapsed, fpol.edge_log), gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second),
            run_wall_s=run_wall, wall_s=time() - t0, error="", stop_reason=stop_reason,
            geometry=seed.gid, role=seed.role, keeper=seed.keeper, epoch_utc=ep, inc_deg=inc, omega_deg=omega, raan_deg=raan, epoch_off_d=Float64(geom["epoch_off_d"]), decision_speed=seed.decision_speed,
            safety_q_w_cm2=qstop, safety_alt0_m=ss.alt0_m, safety_h_local_m=ss.h_local_m, safety_reason=ss.reason)
    catch err
        return (; planet="mars", tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl=seed.keeper, code=seed.code, status="run_error", retcode="", passes=0, t_end_days=NaN, horizon_days=seed.horizon_days,
            completed=false, controller_stopped=false, controller_failures="", total_dv_mps=NaN, burns=0, collapsed_decisions=0, gram_core_epoch="", run_wall_s=NaN, wall_s=time() - t0,
            error=replace(first(sprint(showerror, err), 500), "\n" => " "), stop_reason="error",
            geometry=seed.gid, role=seed.role, keeper=seed.keeper, epoch_utc="", inc_deg=NaN, omega_deg=NaN, raan_deg=NaN, epoch_off_d=NaN, decision_speed=seed.decision_speed,
            safety_q_w_cm2=NaN, safety_alt0_m=NaN, safety_h_local_m=NaN, safety_reason="")
    end
end
