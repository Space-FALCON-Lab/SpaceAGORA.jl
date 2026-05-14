@inline function _extract_pos_vel(x)
    if hasproperty(x, :pos) && hasproperty(x, :vel)
        return SVector{3, Float64}(x[1], x[2], x[3]), SVector{3, Float64}(x[4], x[5], x[6])
    end
    return SVector{3, Float64}(x[1], x[2], x[3]), SVector{3, Float64}(x[4], x[5], x[6])
end

@inline function _extract_mass_kg(x)::Float64
    if hasproperty(x, :mass)
        return Float64(x[7])
    end
    return length(x) >= 7 ? Float64(x[7]) : NaN
end

@inline function _write_density_buffers!(
    p,
    sat_idx::Int,
    rho::Float64,
    T::Float64,
    wind_vec::SVector{3, Float64},
    t::Float64=p.shared_buffers.current_time[],
)
    if sat_idx <= length(p.shared_buffers.densities)
        p.shared_buffers.densities[sat_idx] = rho
    end
    if sat_idx <= length(p.shared_buffers.temperatures)
        p.shared_buffers.temperatures[sat_idx] = T
    end
    if sat_idx <= length(p.shared_buffers.winds)
        p.shared_buffers.winds[sat_idx] = wind_vec
    end
    if sat_idx <= length(p.shared_buffers.density_sample_t)
        p.shared_buffers.density_sample_t[sat_idx] = t
    end
    return nothing
end

@inline function _write_density_time_buffers!(p, num_sats::Int, t::Float64)::Nothing
    times = p.shared_buffers.density_sample_t
    limit = min(num_sats, length(times))
    @inbounds for sat_idx in 1:limit
        times[sat_idx] = t
    end
    return nothing
end

@inline function _density_segment_end_t(p, t::Float64, cache_cfg)::Float64
    segment_end_t = p.shared_buffers.solve_segment_end_time[]
    if isfinite(segment_end_t) && segment_end_t > t
        return segment_end_t
    end
    mission_end_t = Float64(p.args.mission_configuration.mission_time)
    if isfinite(mission_end_t) && mission_end_t > t
        return mission_end_t
    end
    return t + cache_cfg.orbit_horizon_s
end

@inline function _stage_environment_kinematics(x, p, t::Float64)
    engine = _simulation_engine_module()
    pos_ii, vel_ii = _extract_pos_vel(x)
    planet_frame = engine.sample_planet_frame(x, p, 1, t)
    return (
        l_pi=planet_frame.l_pi,
        pos_ii=pos_ii,
        vel_ii=vel_ii,
        pos_pp=planet_frame.pos_pp,
        vel_pp=planet_frame.vel_pp,
        alt=planet_frame.alt_m,
        lat=planet_frame.lat_rad,
        lon=planet_frame.lon_rad,
    )
end

@inline function _buffered_stage_environment_state(x, p, sat_idx::Int, t::Float64)
    kin = _stage_environment_kinematics(x, p, t)
    atmosphere = _simulation_engine_module().sample_buffered_atmosphere(x, p, sat_idx, t)
    return merge(kin, (rho=atmosphere.rho_kg_m3, T=atmosphere.temperature_k, wind=atmosphere.wind_pp))
end

function _density_state_from_kinematics!(
    p,
    sat_idx::Int,
    pos_ii::SVector{3, Float64},
    vel_ii::SVector{3, Float64},
    current_mass_kg::Float64,
    alt::Float64,
    lat::Float64,
    lon::Float64,
    t::Float64,
    density_model,
    cache_cfg,
    stats_enabled::Bool,
    target_include_j2::Bool,
    caches::Vector{Union{Nothing, GramTrackCache}}
)::Tuple{Float64, Float64, SVector{3, Float64}}
    # Vacuum-predicted GRAM density cache: interpolate from a pre-built spline on
    # log(ρ) along the drag-free trajectory.  Only active inside the atmosphere
    # (in_atmosphere flag) to avoid wasteful builds during coast arcs.
    if _vacuum_gram_cache_enabled()
        in_atm = sat_idx <= length(p.shared_buffers.in_atmosphere) &&
                 p.shared_buffers.in_atmosphere[sat_idx]
        if in_atm
            vacuum_cache = _vacuum_gram_cache_for_sat!(p.shared_buffers.vacuum_gram_caches, sat_idx)
            return _query_vacuum_gram_cache!(
                vacuum_cache, density_model, p, pos_ii, vel_ii, alt, t,
                _vacuum_gram_cache_npoints(),
                _vacuum_gram_cache_horizon_s(),
                _vacuum_gram_cache_deviation_m()
            )
        end
    end

    if stats_enabled
        _gram_runtime_stats_update!(s -> begin
            s.density_calls += 1
        end)
    end

    if _gram_track_cache_enabled(cache_cfg, density_model)
        if stats_enabled
            _gram_runtime_stats_update!(s -> begin
                s.cache_enabled_calls += 1
            end)
        end
        cache_horizon_s, cache_alt_tol_m, cache_ang_tol_rad, cache_points = _gram_track_cache_profile(cache_cfg, p, alt)
        cache = _gram_density_cache_for_sat!(caches, sat_idx)
        seg = if stats_enabled
            seg0 = _gram_track_cache_segment(cache, t)
            if seg0 === nothing
                _gram_runtime_stats_update!(s -> begin
                    s.cache_misses += 1
                    s.miss_time_window += 1
                end)
                nothing
            else
                idx0, x0 = seg0
                alt_interp = _lerp(cache.alts[idx0], cache.alts[idx0 + 1], x0)
                lat_interp = _lerp_angle_rad(cache.lats[idx0], cache.lats[idx0 + 1], x0)
                lon_interp = _lerp_angle_rad(cache.lons[idx0], cache.lons[idx0 + 1], x0)
                alt_err = abs(alt - alt_interp)
                lat_err = _angdiff_rad(lat, lat_interp)
                lon_err = _angdiff_rad(lon, lon_interp)
                _gram_runtime_stats_update!(s -> begin
                    s.state_error_samples += 1
                    s.alt_err_abs_max_m = max(s.alt_err_abs_max_m, alt_err)
                    s.lat_err_abs_max_deg = max(s.lat_err_abs_max_deg, rad2deg(lat_err))
                    s.lon_err_abs_max_deg = max(s.lon_err_abs_max_deg, rad2deg(lon_err))
                    s.alt_err_abs_sum_m += alt_err
                    s.lat_err_abs_sum_deg += rad2deg(lat_err)
                    s.lon_err_abs_sum_deg += rad2deg(lon_err)
                end)
                if abs(alt - alt_interp) <= cache_alt_tol_m &&
                   _angdiff_rad(lat, lat_interp) <= cache_ang_tol_rad &&
                   _angdiff_rad(lon, lon_interp) <= cache_ang_tol_rad
                    _gram_runtime_stats_update!(s -> begin
                        s.cache_hits += 1
                    end)
                    seg0
                else
                    _gram_runtime_stats_update!(s -> begin
                        s.cache_misses += 1
                        s.miss_state_tolerance += 1
                    end)
                    nothing
                end
            end
        else
            _gram_track_cache_ready(cache, t, alt, lat, lon, cache_alt_tol_m, cache_ang_tol_rad)
        end
        if seg !== nothing
            idx, x = seg
            return _gram_track_cache_eval(cache, idx, x)
        end
        segment_end_t = _density_segment_end_t(p, t, cache_cfg)
        return _gram_track_cache_refresh!(
            cache,
            density_model,
            p,
            pos_ii,
            vel_ii,
            alt,
            lat,
            lon,
            t,
            cache_horizon_s,
            cache_points,
            cache_alt_tol_m,
            cache_ang_tol_rad,
            cache_cfg.transition_band_m,
            segment_end_t,
            target_include_j2,
            sat_idx,
            current_mass_kg
        )
    end

    if stats_enabled
        _gram_runtime_stats_update!(s -> begin
            s.direct_calls += 1
        end)
    end
    return getDensity(density_model, alt, lat, lon, t, true, p)
end

function _stage_environment_state(x, p, sat_idx::Int, t::Float64; write_buffers::Bool=true)
    kin = _stage_environment_kinematics(x, p, t)
    atmosphere = write_buffers ?
        _simulation_engine_module().sample_buffered_atmosphere(x, p, sat_idx, t) :
        _simulation_engine_module().sample_atmosphere(x, p, sat_idx, t; write_buffers=false)
    return merge(kin, (rho=atmosphere.rho_kg_m3, T=atmosphere.temperature_k, wind=atmosphere.wind_pp))
end

function get_density_callback(num_sats::Int, args::SimulationConfiguration)
    return get_density_callback(num_sats, args.dynamics_model.dynamic_effectors, args)
end

@inline _density_callback_et(p, t::Float64) = p.shared_buffers.et_start[] + t

function get_density_callback(num_sats::Int, effectors::Tuple, args::SimulationConfiguration)
    cache_cfg = _gram_track_cache_config()
    stats_enabled = _gram_runtime_stats_enabled()
    target_include_j2 = _gram_track_cache_target_use_j2() && _uses_j2_gravity_effector(effectors)
    function update_density_sat!(i::Int, p, u, t::Float64)
        density_model = _density_model_for_sat(p, i)
        caches = p.shared_buffers.gram_density_cache
        kin = _stage_environment_kinematics(u.sc[i], p, Float64(t))
        rho, T, wind_vec = _density_state_from_kinematics!(
            p,
            i,
            kin.pos_ii,
            kin.vel_ii,
            _extract_mass_kg(u.sc[i]),
            kin.alt,
            kin.lat,
            kin.lon,
            Float64(t),
            density_model,
            cache_cfg,
            stats_enabled,
            target_include_j2,
            caches,
        )
        _write_density_buffers!(p, i, rho, T, wind_vec, Float64(t))
        return nothing
    end

    condition(u, t, integrator) = true

    function affect!(integrator)
        p = integrator.p
        u = integrator.u
        density_models = p.shared_buffers.density_models
        fallback_density_model = p.args.environment_model.density_model
        decision = _density_callback_thread_decision(args, num_sats)
        use_threads = decision.use_threads
        use_batch = false
        batch_model = nothing
        use_gram_isolated_pool = false
        if _density_batch_enabled(num_sats)
            batch_model = _density_batch_model_for_callback(density_models, fallback_density_model, num_sats)
            if !(batch_model === nothing) && !_gram_track_cache_enabled(cache_cfg, batch_model)
                use_batch = true
            end
        end
        if !use_batch
            gram_pool_model = _gram_isolated_pool_batch_model_for_callback(density_models, fallback_density_model, num_sats)
            if !(gram_pool_model === nothing) && !_gram_track_cache_enabled(cache_cfg, gram_pool_model)
                use_batch = true
                batch_model = gram_pool_model
                use_gram_isolated_pool = true
            end
        elseif batch_model isa EnvironmentModels.GRAMAtmosphereModel
            use_gram_isolated_pool = true
        end
        started_ns = time_ns()

        if use_batch
            # Batch calls are only used when the active density path can evaluate
            # all spacecraft together without per-satellite cache state.
            if stats_enabled
                _gram_runtime_stats_update!(s -> begin
                    s.density_calls += num_sats
                    s.direct_calls += num_sats
                end)
            end
            alts = p.shared_buffers.density_batch_altitudes
            lats = p.shared_buffers.density_batch_latitudes
            lons = p.shared_buffers.density_batch_longitudes
            if use_threads
                ParallelPolicy.threaded_foreach_persistent(:density_callback, num_sats, decision.allotment) do i
                    @inbounds begin
                        kin = _stage_environment_kinematics(u.sc[i], p, Float64(integrator.t))
                        alts[i] = kin.alt
                        lats[i] = kin.lat
                        lons[i] = kin.lon
                    end
                end
            else
                @inbounds for i in 1:num_sats
                    kin = _stage_environment_kinematics(u.sc[i], p, Float64(integrator.t))
                    alts[i] = kin.alt
                    lats[i] = kin.lat
                    lons[i] = kin.lon
                end
            end
            pooled = use_gram_isolated_pool && _gram_isolated_pool_batch_eval!(
                p.shared_buffers.densities,
                p.shared_buffers.temperatures,
                p.shared_buffers.winds,
                batch_model,
                alts,
                lats,
                lons,
                Float64(integrator.t),
                true,
                p;
                allotment_hint=decision.allotment
            )
            if !pooled
                getDensityBatch!(
                    p.shared_buffers.densities,
                    p.shared_buffers.temperatures,
                    p.shared_buffers.winds,
                    batch_model,
                    alts,
                    lats,
                    lons,
                    Float64(integrator.t),
                    true,
                    p
                )
            end
            _write_density_time_buffers!(p, num_sats, Float64(integrator.t))
        elseif use_threads
            ParallelPolicy.threaded_foreach_persistent(:density_callback, num_sats, decision.allotment) do i
                @inbounds update_density_sat!(i, p, u, Float64(integrator.t))
            end
        else
            @inbounds for i in 1:num_sats
                update_density_sat!(i, p, u, Float64(integrator.t))
            end
        end
        if decision.policy_applied
            ParallelPolicy.record_policy_observation!(
                :density_callback;
                mode=decision.mode,
                num_items=num_sats,
                use_threads=use_threads,
                elapsed_ns=(time_ns() - started_ns)
            )
        end
    end

    return DiscreteCallback(condition, affect!, initialize=(cb, u, t, integrator) -> affect!(integrator))
end
