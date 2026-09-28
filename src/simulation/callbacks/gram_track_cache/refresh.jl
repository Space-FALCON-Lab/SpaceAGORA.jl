function _gram_track_cache_fill_from_trajectory!(
    cache::GramTrackCache,
    trajectory
)
    n = length(trajectory)
    n >= 2 || throw(ArgumentError("GRAM trajectory must contain at least 2 points, got $n."))
    if length(cache.times) != n
        resize!(cache.times, n)
        resize!(cache.alts, n)
        resize!(cache.lats, n)
        resize!(cache.lons, n)
        resize!(cache.rhos, n)
        resize!(cache.Ts, n)
        resize!(cache.winds, n)
    end

    @inbounds for k in 1:n
        pt = trajectory[k]
        pos = pt.position
        dyn = pt.dynamics
        wnd = pt.winds

        cache.times[k] = Float64(pos.elapsedTime)
        cache.alts[k] = Float64(pos.height) * 1e3
        cache.lats[k] = deg2rad(Float64(pos.latitude))
        lon_rad = deg2rad(Float64(pos.longitude))
        # Normalise to (-π, π] without trig — just a conditional add/subtract.
        if lon_rad > π
            lon_rad -= 2π
        elseif lon_rad ≤ -π
            lon_rad += 2π
        end
        cache.lons[k] = lon_rad
        cache.rhos[k] = Float64(dyn.density)
        cache.Ts[k] = Float64(dyn.temperature)
        cache.winds[k] = SVector{3, Float64}(
            Float64(wnd.perturbedEWWind),
            Float64(wnd.perturbedNSWind),
            Float64(wnd.perturbedVerticalWind)
        )
    end

    return nothing
end

@inline function _gram_kepler_or_linear_target(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet,
    dt::Float64,
    include_j2::Bool
)::Tuple{Float64, Float64, Float64}
    target = _gram_kepler_target(pos, vel, planet, dt; include_j2=include_j2)
    return target === nothing ? _gram_linear_target(pos, vel, planet, dt) : target
end

function _gram_track_cache_refresh!(
    cache::GramTrackCache,
    density_model,
    p,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    alt::Float64,
    lat::Float64,
    lon::Float64,
    t::Float64,
    horizon_s::Float64,
    n_points::Int,
    alt_tol_m::Float64 = 500.0,
    ang_tol_rad::Float64 = deg2rad(1.0),
    transition_band_m::Float64 = 0.0,
    segment_end_t::Float64 = NaN,
    target_include_j2::Bool = true,
    sat_index::Int = 1,
    current_mass_kg::Float64 = NaN
)
    stats_enabled = _gram_runtime_stats_enabled()
    refresh_t0_ns = stats_enabled ? time_ns() : 0
    try
        planet = p.args.environment_model.planet
        ephemerides_model = p.args.environment_model.ephemerides_model
        base_horizon_s = max(1e-3, horizon_s)
        include_j2 = target_include_j2
        EI_m = p.args.environment_model.EI * 1e3
        in_atm_band = alt <= EI_m + max(0.0, transition_band_m)
        mission_is_orbit = p.args.mission_configuration.mission_type == MissionOrbits

        dt_segment = base_horizon_s
        target_alt, target_lat, target_lon = _gram_kepler_or_linear_target(pos, vel, planet, dt_segment, include_j2)

        if _gram_track_cache_periapsis_split_enabled() && in_atm_band
            # Drag passage: build the cache segment to periapsis when possible.
            peri_target = _gram_periapsis_target(pos, vel, planet; include_j2=include_j2)
            if peri_target !== nothing
                dt_segment, target_alt, target_lat, target_lon = peri_target
            else
                # Entry-like/open trajectories: use Allen-Eggers-type endpoint targeting.
                used_entry_model = false
                if _gram_entry_target_mode() == :allen_eggers
                    dt_entry = if isfinite(segment_end_t)
                        max(1e-3, segment_end_t - t)
                    else
                        base_horizon_s
                    end
                    if isfinite(dt_entry) && dt_entry > 1e-6
                        mass_kg = _gram_entry_mass_kg(p, sat_index, current_mass_kg)
                        area_m2 = _gram_entry_reference_area_m2(p, sat_index)
                        entry_target = _gram_entry_target_allen_eggers(
                            pos,
                            vel,
                            planet,
                            dt_entry,
                            mass_kg,
                            area_m2;
                            cd=_gram_entry_target_cd()
                        )
                        if entry_target !== nothing
                            dt_segment = dt_entry
                            target_alt, target_lat, target_lon = entry_target
                            used_entry_model = true
                        end
                    end
                end

                # If entry model is disabled/unavailable, use solver endpoint propagation.
                if !used_entry_model && isfinite(segment_end_t)
                    dt_to_end = segment_end_t - t
                    if dt_to_end > 1e-6
                        dt_segment = max(1e-3, dt_to_end)
                        target_alt, target_lat, target_lon = _gram_kepler_or_linear_target(pos, vel, planet, dt_segment, include_j2)
                    end
                end
            end
        elseif mission_is_orbit
            # Orbit mission: build one segment over a full orbital period.
            orbit_target = _gram_orbit_period_target(pos, vel, planet; include_j2=include_j2)
            if orbit_target !== nothing
                dt_segment, target_alt, target_lat, target_lon = orbit_target
            end
        else
            # Time mission: propagate to the remaining solver time span.
            used_tspan_endpoint = false
            if isfinite(segment_end_t)
                dt_to_end = segment_end_t - t
                if dt_to_end > 1e-6
                    dt_segment = max(1e-3, dt_to_end)
                    target_alt, target_lat, target_lon = _gram_kepler_or_linear_target(pos, vel, planet, dt_segment, include_j2)
                    used_tspan_endpoint = true
                end
            end
            if !used_tspan_endpoint
                orbit_target = _gram_orbit_period_target(pos, vel, planet; include_j2=include_j2)
                if orbit_target !== nothing
                    dt_segment, target_alt, target_lat, target_lon = orbit_target
                end
            end
        end

        base_n = max(2, n_points)
        dt_per_sample = max(1e-3, base_horizon_s / max(1, base_n - 1))
        n_from_time = Int(ceil(dt_segment / dt_per_sample)) + 1
        expected_length_m, radius_m = _gram_expected_track_length_m(
            planet,
            alt,
            lat,
            lon,
            target_alt,
            target_lat,
            target_lon
        )
        target_spacing_m = _gram_track_cache_target_spacing_m(alt_tol_m, ang_tol_rad, radius_m)
        n_from_length = Int(ceil(expected_length_m / target_spacing_m)) + 1
        n = max(base_n, n_from_time, n_from_length)
        n = min(n, _gram_track_cache_max_npos())

        if stats_enabled
            _gram_runtime_stats_update!(s -> begin
                s.refresh_calls += 1
                s.refresh_points_total += n
                s.refresh_points_max = max(s.refresh_points_max, n)
            end)
        end

        if length(cache.times) != n
            resize!(cache.times, n)
            resize!(cache.alts, n)
            resize!(cache.lats, n)
            resize!(cache.lons, n)
            resize!(cache.rhos, n)
            resize!(cache.Ts, n)
            resize!(cache.winds, n)
        end

        gram_traj_built = false
        if _gram_track_trajectory_supported(density_model)
            gram_driver = getproperty(density_model, :gram)
            gram_atmosphere = getproperty(density_model, :gram_atmosphere)
            denom = max(1, n - 1)
            Δalt_km = (target_alt - alt) * 1e-3 / denom
            Δlat_deg = rad2deg(_angle_delta_rad(lat, target_lat)) / denom
            Δlon_deg = rad2deg(_angle_delta_rad(lon, target_lon)) / denom
            Δt = dt_segment / denom

            gen_track = () -> Base.invokelatest(
                getproperty(gram_driver, :generate_trajectory),
                gram_atmosphere;
                initial_height=alt * 1e-3,
                initial_latitude=rad2deg(lat),
                initial_longitude=rad2deg(lon),
                initial_elapsed_time=t,
                delta_height=Δalt_km,
                delta_latitude=Δlat_deg,
                delta_longitude=Δlon_deg,
                delta_elapsed_time=Δt,
                n_points=n,
                update_initial_perturbations=true
            )

            trajectory = if EnvironmentModels._gram_use_global_lock()
                lock(tracked_lock(:gram_cache)) do
                    gen_track()
                end
            else
                gen_track()
            end
            _gram_track_cache_fill_from_trajectory!(cache, trajectory)
            gram_traj_built = true
        end

        if !gram_traj_built
            @inbounds for k in 1:n
                x = (k - 1) / (n - 1)
                t_k = t + x * dt_segment
                et_k = p.shared_buffers.et_start[] + t_k
                pos_k = pos + vel * (t_k - t)
                rp_k, _ = r_intor_p!(pos_k, vel, planet, et_k, ephemerides_model)
                alt_k, lat_k, lon_k = rtolatlong(rp_k, planet, ephemerides_model)
                cache.times[k] = t_k
                cache.alts[k] = alt_k
                cache.lats[k] = lat_k
                cache.lons[k] = lon_k
            end
            if !_gram_isolated_pool_batch_eval!(
                cache.rhos,
                cache.Ts,
                cache.winds,
                density_model,
                cache.alts,
                cache.lats,
                cache.lons,
                cache.times,
                true,
                p
            )
                getDensityBatch!(
                    cache.rhos,
                    cache.Ts,
                    cache.winds,
                    density_model,
                    cache.alts,
                    cache.lats,
                    cache.lons,
                    cache.times,
                    true,
                    p
                )
            end
        end
        cache.valid = true
        cache.t0 = cache.times[1]
        cache.t1 = cache.times[end]
        cache.index_hint = 1
        if stats_enabled
            _gram_runtime_stats_update!(s -> begin
                s.refresh_elapsed_s += (time_ns() - refresh_t0_ns) * 1e-9
            end)
        end
        return cache.rhos[1], cache.Ts[1], cache.winds[1]
    catch err
        cache.valid = false
        if stats_enabled
            _gram_runtime_stats_update!(s -> begin
                s.refresh_failures += 1
                s.refresh_elapsed_s += (time_ns() - refresh_t0_ns) * 1e-9
            end)
        end
        if !_gram_track_cache_warning_emitted[]
            _gram_track_cache_warning_emitted[] = true
            @warn "GRAM track cache refresh failed; falling back to direct GRAM sampling." exception=err
        end
        return getDensity(density_model, alt, lat, lon, t, true, p)
    end
end
