# Venus gate V-G(V) (plan R8.3, R8.4; scoring R4.3, R4.4, R5.1, U19): free-drift validation flights at the Venus gate geometries VG1, VG2, VGref.
# Reuses stageg3_common.jl (INCLUDED, never edited, sha256 2fc55ffe...4bbd): the U15 periapsis_vrel, the safety stop state machine, the apsis callback with the safety tests and the relative-speed log, the U15 guidance wrapper
# StageGPolicy, the GRAM epoch guard, and through it stageo_common.jl / stagen2_common.jl / the Stage F files (planet-generic; Stage L flew them at Venus). New code only:
#   (i)   the Venus config view: [planets.venus.*] plus the carried band and the start orbit from the Venus geometry table (config [stageg3v.table], sha256 verified at load, a mismatch is an ERROR);
#   (ii)  the per-geometry fit-only calibration and the flight case function, the Mars gate's with "venus" for "mars" and [stageg3v] for [stageg3]/[stageo].
# No src/ or examples/ change.
using Dates, TOML, JSON, SHA
include(joinpath(@__DIR__, "stageg3_common.jl"))

const _G3VMETA = Dict{String, Any}()
"""Config dict with the geometry list injected from the Venus geometry table (config [stageg3v.table]: path, sha256). The table's sha256 is verified at load: a mismatch is an ERROR. The table's start block is checked against
[stageg3v.start]. meta["stageg3v"]["geometry"] = geometries (id, role, inc_deg, omega_deg, raan_deg, epoch_off_d, epoch_utc, field_id)."""
function stageg3v_meta(path::AbstractString)
    get!(_G3VMETA, abspath(path)) do
        m = TOML.parsefile(path); t = m["stageg3v"]["table"]; fp = joinpath(stageg3_task_dir(), String(t["path"])); raw = read(fp); h = bytes2hex(sha256(raw))
        h == String(t["sha256"]) || error("geometry table sha256 mismatch: $fp is $h, config has $(t["sha256"])")
        T = JSON.parse(String(raw)); G = Dict{String, Any}[]
        for g in T["geometries"]
            g["in_default_set"] || continue
            push!(G, Dict{String, Any}("id" => g["id"], "role" => g["role"], "inc_deg" => Float64(g["inclination_deg"]), "omega_deg" => Float64(g["omega0_deg"]), "raan_deg" => Float64(g["raan0_deg"]),
                "epoch_off_d" => Float64(g["epoch_offset_days"]), "epoch_utc" => g["epoch_utc"], "field_id" => g["id"]))
        end
        st = m["stageg3v"]["start"]; ts = T["start"]
        (Float64(ts["ra0_m"]) == Float64(st["ra0_m"]) && Float64(ts["ra_target_m"]) == Float64(st["ra_target_m"]) && Float64(ts["true_anomaly_deg"]) == Float64(st["true_anomaly_deg"])) || error("table start differs from [stageg3v.start]")
        m["stageg3v"]["geometry"] = G; m["stageg3v"]["table_start"] = ts; m["stageg3v"]["table_hash"] = h; m
    end
end

"""Geometry record by id."""
function stageg3v_geometry(meta, gid::AbstractString)
    i = findfirst(g -> String(g["id"]) == gid, meta["stageg3v"]["geometry"]); i === nothing && throw(ArgumentError("unknown geometry $gid"))
    return meta["stageg3v"]["geometry"][i]
end

"""Flat Venus config of one geometry: [common] + [planets.venus] with the start apoapsis radius injected, the band replaced by [stageg3v.band], the atmosphere step optionally overridden and the epoch replaced by the table's
(checked against base epoch + offset). Returns (meta, sc, geom, epoch_string)."""
function stageg3v_load(path::AbstractString, gid::AbstractString; dt_atm_s::Real=NaN)
    meta = stageg3v_meta(path); G3 = meta["stageg3v"]; sc = load_stagef("venus"; path=path)
    Float64(sc["orbit"]["apoapsis_radius_m"]) == Float64(G3["start"]["ra0_m"]) || error("[planets.venus.orbit] apoapsis radius differs from [stageg3v.start]")
    sc["orbit"]["apoapsis_radius_m"] = Float64(G3["start"]["ra0_m"])
    sc["band"] = merge(Dict{String, Any}(sc["band"]), Dict{String, Any}(G3["band"]))
    isfinite(dt_atm_s) && (sc["integration"]["dt_max_atmosphere_s"] = Float64(dt_atm_s))
    g = stageg3v_geometry(meta, gid); ep = String(g["epoch_utc"])
    ep == stageg3_epoch_string(String(G3["base_epoch_utc"]), Float64(g["epoch_off_d"])) || error("table epoch $ep of $gid differs from base epoch + offset")
    sc["grid"]["epochs_utc"] = [ep]
    return meta, sc, g, ep
end

# ---- phase cal: one fit-only calibration case per geometry ----------------------------------------------------------------------------------

"""Calibration case of one geometry (the Mars gate's on Venus); writes <seed.out>/<gid>/{guess.csv, nodes.csv, calibration.csv}. Plain-data row; never throws."""
function stageg3v_cal_task(seed)
    t0 = time()
    try
        meta, sc, g, ep = stageg3v_load(seed.config_path, seed.gid); scc = stageo_cal_view(sc); ms = stagef_setup(sc)
        epoch = parse_epoch(sc["grid"]["epochs_utc"][1]); rot = mars_rot(ms.planet, epoch)
        et0 = SimulationModel.ephemerides_time_seconds(epoch, SpiceEphemeridesModel()); stagef_env!(false, meta)
        model = stagef_gram_model(sc, epoch, stagef_seed("venus" * seed.gram_tag)); pargs = gram_pargs(ms.planet, Float64(sc["atmosphere_gram"]["entry_interface_km"]))
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
        row = (; planet="venus", tag=seed.gid, status="ok", h_guess_m=gg.h_g, seed_alt_m=seed_alt, vrel_first_mps=recs[end].vrel_mps, q_first_w_cm2=recs[end].q_w_cm2, q_target_w_cm2=recs[end].q_target,
            site_lat_rad=gg.site.lat, site_lon_rad=gg.site.lon, fit_H_m=f.H, pw_max_abs_lnrho=f.pw_max_abs_lnrho, hloc_at_hg_m=f.hloc_at_hg, n_cal_iter=length(recs),
            gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second), wall_s=time() - t0, error="", epoch_utc=ep)
        write_rows_m(joinpath(dir, "guess.csv"), [row]); return row
    catch err
        return (; planet="venus", tag=seed.gid, status="cal_error", h_guess_m=NaN, seed_alt_m=NaN, vrel_first_mps=NaN, q_first_w_cm2=NaN, q_target_w_cm2=NaN, site_lat_rad=NaN, site_lon_rad=NaN, fit_H_m=NaN,
            pw_max_abs_lnrho=NaN, hloc_at_hg_m=NaN, n_cal_iter=0, gram_core_epoch="", wall_s=time() - t0, error=replace(first(sprint(showerror, err), 400), "\n" => " "), epoch_utc="")
    end
end

# ---- phase run: one flight --------------------------------------------------------------------------------------------------------------------------

"""One Venus gate flight (the Mars gate's stageg3_run_task with the Venus config view). seed fields: config_path, gid, role, keeper ("fixed" | "zero_burn"), tag, arm, code, run_id, out, cal_root, horizon_days, dt_atm_s, smoke_force, decision_speed."""
function stageg3v_run_task(seed)
    t0 = time()
    try
        meta, sc, geom, ep = stageg3v_load(seed.config_path, seed.gid; dt_atm_s=seed.dt_atm_s); G3 = meta["stageg3v"]; SO = G3; ms0 = stagef_setup(sc)
        inc, raan, omega = Float64(geom["inc_deg"]), Float64(geom["raan_deg"]), Float64(geom["omega_deg"])
        epoch = parse_epoch(ep); rot = mars_rot(ms0.planet, epoch); zero = seed.keeper == "zero_burn"
        cdir = joinpath(seed.cal_root, seed.gid); G = CSV.read(joinpath(cdir, "guess.csv"), DataFrame)[1, :]; N = CSV.read(joinpath(cdir, "nodes.csv"), DataFrame)
        pwf = pw_fit(N.h_m, N.rho_kg_m3, Float64(sc["atmosphere"]["temperature_k"]), "piecewise"); ms = with_fit(ms0, pwf)
        stagef_env!(true, meta)
        model = stagef_gram_model(sc, epoch, stagef_seed("venus" * seed.tag)); core_t = model.core.initial_time; stageg3_gram_guard(core_t, ep)
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
        obs = NamedTuple[]; pk = PassBook(Float64(G3["start"]["ra_target_m"]))
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
        return (; planet="venus", tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl=seed.keeper, code=seed.code, status="ok", retcode=ret,
            passes=length(rows), t_end_days=sol.t[end] / SECONDS_PER_DAY, horizon_days=seed.horizon_days, completed=pk.stopped && !stopped, controller_stopped=stopped, controller_failures=failures, total_dv_mps=dv, burns=burns,
            collapsed_decisions=count(e -> e.collapsed, fpol.edge_log), gram_core_epoch=string(core_t.year, "-", core_t.month, "-", core_t.day, " ", core_t.hour, ":", core_t.minute, ":", core_t.second),
            run_wall_s=run_wall, wall_s=time() - t0, error="", stop_reason=stop_reason,
            geometry=seed.gid, role=seed.role, keeper=seed.keeper, epoch_utc=ep, inc_deg=inc, omega_deg=omega, raan_deg=raan, epoch_off_d=Float64(geom["epoch_off_d"]), decision_speed=seed.decision_speed,
            safety_q_w_cm2=qstop, safety_alt0_m=ss.alt0_m, safety_h_local_m=ss.h_local_m, safety_reason=ss.reason)
    catch err
        return (; planet="venus", tag=seed.tag, run_id=seed.run_id, arm=seed.arm, ctrl=seed.keeper, code=seed.code, status="run_error", retcode="", passes=0, t_end_days=NaN, horizon_days=seed.horizon_days,
            completed=false, controller_stopped=false, controller_failures="", total_dv_mps=NaN, burns=0, collapsed_decisions=0, gram_core_epoch="", run_wall_s=NaN, wall_s=time() - t0,
            error=replace(first(sprint(showerror, err), 500), "\n" => " "), stop_reason="error",
            geometry=seed.gid, role=seed.role, keeper=seed.keeper, epoch_utc="", inc_deg=NaN, omega_deg=NaN, raan_deg=NaN, epoch_off_d=NaN, decision_speed=seed.decision_speed,
            safety_q_w_cm2=NaN, safety_alt0_m=NaN, safety_h_local_m=NaN, safety_reason="")
    end
end
