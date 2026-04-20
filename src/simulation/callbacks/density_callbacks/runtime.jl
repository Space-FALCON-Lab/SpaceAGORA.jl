function get_density_callback(num_sats::Int, args::SimulationConfiguration)
    return get_density_callback(num_sats, args.dynamics_model.dynamic_effectors, args)
end

@inline _density_callback_et(p, t::Float64) = p.shared_buffers.et_start[] + t

function get_density_callback(num_sats::Int, effectors::Tuple, args::SimulationConfiguration)
    cache_cfg = _gram_track_cache_config()
    stats_enabled = _gram_runtime_stats_enabled()
    target_include_j2 = _gram_track_cache_target_use_j2() && _uses_j2_gravity_effector(effectors)
    function update_density_sat!(i::Int, p, u, t, segment_end_t::Float64, density_model, caches::Vector{Union{Nothing, GramTrackCache}})
        if stats_enabled
            _gram_runtime_stats_update!(s -> begin
                s.density_calls += 1
            end)
        end
        pos = SVector{3, Float64}(u.sc[i].pos)
        vel = SVector{3, Float64}(u.sc[i].vel)
        et = _density_callback_et(p, Float64(t))
        rp, _ = r_intor_p!(pos, vel, p.args.environment_model.planet, et, p.args.environment_model.ephemerides_model)
        alt, lat, lon = rtolatlong(rp, p.args.environment_model.planet, p.args.environment_model.ephemerides_model)
        ρ, T, wind_vec = if _gram_track_cache_enabled(cache_cfg, density_model)
            if stats_enabled
                _gram_runtime_stats_update!(s -> begin
                    s.cache_enabled_calls += 1
                end)
            end
            cache_horizon_s, cache_alt_tol_m, cache_ang_tol_rad, cache_points = _gram_track_cache_profile(cache_cfg, p, alt)
            cache = _gram_density_cache_for_sat!(caches, i)
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
                _gram_track_cache_eval(cache, idx, x)
            else
                _gram_track_cache_refresh!(
                    cache,
                    density_model,
                    p,
                    pos,
                    vel,
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
                    i,
                    Float64(u.sc[i].mass)
                )
            end
        else
            if stats_enabled
                _gram_runtime_stats_update!(s -> begin
                    s.direct_calls += 1
                end)
            end
            getDensity(density_model, alt, lat, lon, t, true, p)
        end
        p.shared_buffers.densities[i] = ρ
        p.shared_buffers.temperatures[i] = T
        p.shared_buffers.winds[i] = wind_vec
        return nothing
    end

    # Logic for when the callback should trigger
    # Returning true means it runs at every integration step
    condition(u, t, integrator) = true

    # Logic for what happens when it triggers
    function affect!(integrator)
        p = integrator.p
        u = integrator.u
        density_models = p.shared_buffers.density_models
        fallback_density_model = p.args.environment_model.density_model
        gram_density_caches = p.shared_buffers.gram_density_cache
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
        segment_end_t = try
            Float64(integrator.sol.prob.tspan[2])
        catch
            integrator.t + cache_cfg.orbit_horizon_s
        end

        if use_batch
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
                        pos = SVector{3, Float64}(u.sc[i].pos)
                        vel = SVector{3, Float64}(u.sc[i].vel)
                        et = _density_callback_et(p, Float64(integrator.t))
                        rp, _ = r_intor_p!(pos, vel, p.args.environment_model.planet, et, p.args.environment_model.ephemerides_model)
                        alt, lat, lon = rtolatlong(rp, p.args.environment_model.planet, p.args.environment_model.ephemerides_model)
                        alts[i] = alt
                        lats[i] = lat
                        lons[i] = lon
                    end
                end
            else
                @inbounds for i in 1:num_sats
                    pos = SVector{3, Float64}(u.sc[i].pos)
                    vel = SVector{3, Float64}(u.sc[i].vel)
                    et = _density_callback_et(p, Float64(integrator.t))
                    rp, _ = r_intor_p!(pos, vel, p.args.environment_model.planet, et, p.args.environment_model.ephemerides_model)
                    alt, lat, lon = rtolatlong(rp, p.args.environment_model.planet, p.args.environment_model.ephemerides_model)
                    alts[i] = alt
                    lats[i] = lat
                    lons[i] = lon
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
        elseif use_threads
            if isempty(density_models)
                density_model = fallback_density_model
                ParallelPolicy.threaded_foreach_persistent(:density_callback, num_sats, decision.allotment) do i
                    @inbounds update_density_sat!(i, p, u, integrator.t, segment_end_t, density_model, gram_density_caches)
                end
            else
                ParallelPolicy.threaded_foreach_persistent(:density_callback, num_sats, decision.allotment) do i
                    density_model = @inbounds density_models[i]
                    @inbounds update_density_sat!(i, p, u, integrator.t, segment_end_t, density_model, gram_density_caches)
                end
            end
        else
            if isempty(density_models)
                density_model = fallback_density_model
                @inbounds for i in 1:num_sats
                    update_density_sat!(i, p, u, integrator.t, segment_end_t, density_model, gram_density_caches)
                end
            else
                @inbounds for i in 1:num_sats
                    density_model = density_models[i]
                    update_density_sat!(i, p, u, integrator.t, segment_end_t, density_model, gram_density_caches)
                end
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
