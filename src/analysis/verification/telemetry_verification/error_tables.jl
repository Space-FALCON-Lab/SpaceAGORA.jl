function _orbit_rows_errors(
    cfg::OrbitEventsScenarioConfig,
    args::SimulationConfiguration,
    results_df::DataFrame,
    max_points::Int;
    bias_by_event::Dict{String, Float64}=Dict{String, Float64}()
)
    extrema = _extract_extrema_series(results_df, args.environment_model.planet, cfg.orbit_altitude_mode)
    tele_peri = _load_telemetry_curve(cfg.telemetry_peri_path, max_points)
    tele_apo = _load_telemetry_curve(cfg.telemetry_apo_path, max_points)
    peri_bias = get(bias_by_event, "peri", 0.0)
    apo_bias = get(bias_by_event, "apo", 0.0)
    peri_step = length(tele_peri.orbit) >= 2 ? median(diff(tele_peri.orbit)) : 1.0
    apo_step = length(tele_apo.orbit) >= 2 ? median(diff(tele_apo.orbit)) : 1.0
    peri_sim_axis = tele_peri.orbit[1] .+ peri_step .* collect(0:(length(extrema.peri.altitude)-1))
    apo_sim_axis = tele_apo.orbit[1] .+ apo_step .* collect(0:(length(extrema.apo.altitude)-1))
    peri_summary, peri_errors = _compare_orbit_curve(
        cfg.name,
        "peri",
        tele_peri.orbit,
        tele_peri.altitude,
        extrema.peri.altitude .+ peri_bias;
        sim_axis=peri_sim_axis
    )
    apo_summary, apo_errors = _compare_orbit_curve(
        cfg.name,
        "apo",
        tele_apo.orbit,
        tele_apo.altitude,
        extrema.apo.altitude .+ apo_bias;
        sim_axis=apo_sim_axis
    )
    return [peri_summary, apo_summary], [peri_errors, apo_errors]
end

@inline function _telemetry_altitude_km(
    pos_m::SVector{3, Float64},
    vel_mps::SVector{3, Float64},
    planet,
    altitude_mode::Symbol
)::Float64
    if altitude_mode == :vacuum
        return (norm(pos_m) - planet.Rp_e) * 1e-3
    elseif altitude_mode == :oblate
        pos_p, _ = r_intor_p!(pos_m, vel_mps, planet)
        return rtolatlong(pos_p, planet)[1] * 1e-3
    end
    throw(ArgumentError("Unsupported altitude_mode='$altitude_mode' in telemetry altitude calculation"))
end

function _time_aligned_rows_errors(
    cfg::TimeAlignedScenarioConfig,
    args::SimulationConfiguration,
    results_df::DataFrame,
    telemetry;
    bias_by_event::Dict{String, Float64}=Dict{String, Float64}()
)
    if cfg.comparison_mode == :orbit_events
        tele_speed_kmps = sqrt.(telemetry.vx_kmps .^ 2 .+ telemetry.vy_kmps .^ 2 .+ telemetry.vz_kmps .^ 2)
        tele_extrema = _extract_extrema_from_time_aligned_telemetry(
            telemetry.time_s,
            telemetry.altitude_km,
            cfg.extrema_min_separation_s;
            speed_kmps=tele_speed_kmps
        )
        sim_extrema = _extract_extrema_series(results_df, args.environment_model.planet, cfg.orbit_altitude_mode)
        peri_bias = get(bias_by_event, "peri", 0.0)
        apo_bias = get(bias_by_event, "apo", 0.0)
        peri_speed_bias = get(bias_by_event, "peri_speed", 0.0)
        apo_speed_bias = get(bias_by_event, "apo_speed", 0.0)

        peri_summary, peri_errors = _compare_orbit_curve(
            cfg.name,
            "peri",
            tele_extrema.peri.orbit,
            tele_extrema.peri.altitude,
            sim_extrema.peri.altitude .+ peri_bias
        )
        apo_summary, apo_errors = _compare_orbit_curve(
            cfg.name,
            "apo",
            tele_extrema.apo.orbit,
            tele_extrema.apo.altitude,
            sim_extrema.apo.altitude .+ apo_bias
        )

        peri_speed_summary, peri_speed_errors = _compare_orbit_curve(
            cfg.name,
            "peri_speed",
            tele_extrema.peri.orbit,
            tele_extrema.peri.speed_kmps,
            sim_extrema.peri.speed_kmps .+ peri_speed_bias
        )
        apo_speed_summary, apo_speed_errors = _compare_orbit_curve(
            cfg.name,
            "apo_speed",
            tele_extrema.apo.orbit,
            tele_extrema.apo.speed_kmps,
            sim_extrema.apo.speed_kmps .+ apo_speed_bias
        )

        return [peri_summary, apo_summary, peri_speed_summary, apo_speed_summary], [peri_errors, apo_errors, peri_speed_errors, apo_speed_errors]
    end

    sim_time = _to_float_vector(results_df.time, "sim-time")
    sim_x_m = _require_column(results_df, ["sc1_pos_1", "sc1_position_1"], "sim-position-x")
    sim_y_m = _require_column(results_df, ["sc1_pos_2", "sc1_position_2"], "sim-position-y")
    sim_z_m = _require_column(results_df, ["sc1_pos_3", "sc1_position_3"], "sim-position-z")
    sim_vx_mps = _require_column(results_df, ["sc1_vel_1", "sc1_velocity_1"], "sim-velocity-x")
    sim_vy_mps = _require_column(results_df, ["sc1_vel_2", "sc1_velocity_2"], "sim-velocity-y")
    sim_vz_mps = _require_column(results_df, ["sc1_vel_3", "sc1_velocity_3"], "sim-velocity-z")

    if cfg.comparison_frame == :planet_fixed
        et0 = _initial_time_et(args.initial_time)
        @inbounds for i in eachindex(sim_time)
            pos_fixed_m, vel_fixed_mps = _j2000_to_planet_fixed_state(
                cfg.planet_name,
                et0 + sim_time[i],
                SVector{3, Float64}(sim_x_m[i], sim_y_m[i], sim_z_m[i]),
                SVector{3, Float64}(sim_vx_mps[i], sim_vy_mps[i], sim_vz_mps[i])
            )
            sim_x_m[i] = pos_fixed_m[1]
            sim_y_m[i] = pos_fixed_m[2]
            sim_z_m[i] = pos_fixed_m[3]
            sim_vx_mps[i] = vel_fixed_mps[1]
            sim_vy_mps[i] = vel_fixed_mps[2]
            sim_vz_mps[i] = vel_fixed_mps[3]
        end
    end

    sim_x_km = sim_x_m .* 1e-3 .+ get(bias_by_event, "state_x_time", 0.0)
    sim_y_km = sim_y_m .* 1e-3 .+ get(bias_by_event, "state_y_time", 0.0)
    sim_z_km = sim_z_m .* 1e-3 .+ get(bias_by_event, "state_z_time", 0.0)
    sim_altitude_km = Vector{Float64}(undef, length(sim_time))
    @inbounds for i in eachindex(sim_time)
        sim_altitude_km[i] = _telemetry_altitude_km(
            SVector{3, Float64}(sim_x_m[i], sim_y_m[i], sim_z_m[i]),
            SVector{3, Float64}(sim_vx_mps[i], sim_vy_mps[i], sim_vz_mps[i]),
            args.environment_model.planet,
            cfg.orbit_altitude_mode
        )
    end
    sim_altitude_km .+= get(bias_by_event, "altitude_time", 0.0)

    altitude_summary, altitude_errors = _compare_time_series(
        cfg.name,
        "altitude_time",
        telemetry.time_s,
        telemetry.altitude_km,
        sim_time,
        sim_altitude_km
    )
    x_summary, x_errors = _compare_time_series(
        cfg.name,
        "state_x_time",
        telemetry.time_s,
        telemetry.x_km,
        sim_time,
        sim_x_km
    )
    y_summary, y_errors = _compare_time_series(
        cfg.name,
        "state_y_time",
        telemetry.time_s,
        telemetry.y_km,
        sim_time,
        sim_y_km
    )
    z_summary, z_errors = _compare_time_series(
        cfg.name,
        "state_z_time",
        telemetry.time_s,
        telemetry.z_km,
        sim_time,
        sim_z_km
    )
    summaries = [altitude_summary, x_summary, y_summary, z_summary]
    errors = [altitude_errors, x_errors, y_errors, z_errors]

    has_velocity_truth = !any(isnothing, (cfg.telemetry_vx_col, cfg.telemetry_vy_col, cfg.telemetry_vz_col))
    if has_velocity_truth
        sim_vx_kmps = sim_vx_mps .* 1e-3 .+ get(bias_by_event, "state_vx_time", 0.0)
        sim_vy_kmps = sim_vy_mps .* 1e-3 .+ get(bias_by_event, "state_vy_time", 0.0)
        sim_vz_kmps = sim_vz_mps .* 1e-3 .+ get(bias_by_event, "state_vz_time", 0.0)

        vx_summary, vx_errors = _compare_time_series(
            cfg.name, "state_vx_time", telemetry.time_s, telemetry.vx_kmps, sim_time, sim_vx_kmps
        )
        vy_summary, vy_errors = _compare_time_series(
            cfg.name, "state_vy_time", telemetry.time_s, telemetry.vy_kmps, sim_time, sim_vy_kmps
        )
        vz_summary, vz_errors = _compare_time_series(
            cfg.name, "state_vz_time", telemetry.time_s, telemetry.vz_kmps, sim_time, sim_vz_kmps
        )
        append!(summaries, [vx_summary, vy_summary, vz_summary])
        append!(errors, [vx_errors, vy_errors, vz_errors])
    end

    return summaries, errors
end

@inline _tolerances_for(cfg::OrbitEventsScenarioConfig, profile::Symbol) = profile == :quick ? cfg.tolerances_quick : cfg.tolerances_full
