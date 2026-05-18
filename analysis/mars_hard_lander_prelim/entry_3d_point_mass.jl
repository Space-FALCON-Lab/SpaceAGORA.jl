struct Trajectory3DResult
    summary::NamedTuple
    trajectory::DataFrame
end

@inline function _safe_unit(v::SVector{3, Float64}, fallback::SVector{3, Float64})
    v_norm = norm(v)
    return v_norm > 1e-12 ? v / v_norm : fallback
end

@inline function _enu_basis(r_vec::SVector{3, Float64})
    up = _safe_unit(r_vec, SVector{3, Float64}(1.0, 0.0, 0.0))
    z_axis = SVector{3, Float64}(0.0, 0.0, 1.0)
    east = _safe_unit(cross(z_axis, up), SVector{3, Float64}(0.0, 1.0, 0.0))
    north = _safe_unit(cross(up, east), SVector{3, Float64}(0.0, 0.0, 1.0))
    return east, north, up
end

@inline function _lat_lon(r_vec::SVector{3, Float64})
    r_norm = norm(r_vec)
    lat = asin(clamp(r_vec[3] / max(r_norm, eps()), -1.0, 1.0))
    lon = atan(r_vec[2], r_vec[1])
    return lat, lon
end

@inline function _wind_frame(v_rel_enu::SVector{3, Float64})
    xw = _safe_unit(-v_rel_enu, SVector{3, Float64}(-1.0, 0.0, 0.0))
    up_enu = SVector{3, Float64}(0.0, 0.0, 1.0)
    zw = _safe_unit(up_enu - dot(up_enu, xw) * xw, SVector{3, Float64}(0.0, 0.0, 1.0))
    yw = _safe_unit(cross(xw, zw), SVector{3, Float64}(0.0, 1.0, 0.0))
    return xw, yw, zw
end

function _lateral_command(config::PrelimConfig, deployed::Bool, crossrange_km::Float64, crossrange_rate_mps::Float64)
    deployed || return 0.0
    return clamp(
        -(crossrange_km / config.lateral_guidance_crossrange_scale_km) -
        (crossrange_rate_mps / config.lateral_guidance_crossrange_rate_scale_mps),
        -1.0,
        1.0,
    )
end

function _sample_3d_state(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    aero_case::CalibratedAeroCase,
    policy::ControlPolicy,
    state::SVector{6, Float64},
    time_s::Float64,
    density_scale::Float64,
    use_winds::Bool,
    panel_command_override::Union{Nothing, DifferentialPanelCommand}=nothing,
)
    r_vec = SVector{3, Float64}(state[1], state[2], state[3])
    v_vec = SVector{3, Float64}(state[4], state[5], state[6])
    r_norm = norm(r_vec)
    h_m = max(r_norm - config.planet.Rp_m, 0.0)
    east, north, up = _enu_basis(r_vec)
    lat, lon = _lat_lon(r_vec)

    v_e = dot(v_vec, east)
    v_n = dot(v_vec, north)
    v_u = dot(v_vec, up)
    rho_kg_m3, temperature_k, wind_enu = atmosphere_state(
        adapter,
        h_m,
        lat,
        lon,
        time_s;
        include_wind=use_winds,
    )
    rho_kg_m3 *= density_scale

    v_rel_enu = SVector{3, Float64}(v_e, v_n, v_u) - wind_enu
    v_rel_mag = max(norm(v_rel_enu), 1.0)
    mach = mach_number(adapter, v_rel_mag, temperature_k)
    deployed = deployed_active(policy, h_m)
    trim = trim_active(policy, h_m)
    crossrange_km = config.planet.Rp_m * lat / 1e3
    lateral_cmd = panel_command_override === nothing ?
        _lateral_command(config, deployed, crossrange_km, v_n) :
        NaN
    active_panel_command = (deployed && trim) ? panel_command_override : nothing
    loads = aerodynamic_loads_3d(
        config,
        aero_case,
        deployed,
        mach,
        config.planet.γ,
        panel_command_override === nothing ? lateral_cmd : 0.0,
        active_panel_command,
        trim,
    )
    q_pa = 0.5 * rho_kg_m3 * v_rel_mag^2

    xw, yw, zw = _wind_frame(v_rel_enu)
    aero_enu = (q_pa / aero_case.mass_kg) * (
        loads.CDA_m2 * xw +
        loads.CLA_m2 * zw +
        loads.CSA_m2 * yw
    )
    aero_cart = aero_enu[1] * east + aero_enu[2] * north + aero_enu[3] * up
    gravity_cart = -(config.planet.μ / r_norm^3) * r_vec
    acc_cart = gravity_cart + aero_cart

    drag_accel_mps2 = q_pa * loads.CDA_m2 / aero_case.mass_kg
    lift_accel_mps2 = q_pa * abs(loads.CLA_m2) / aero_case.mass_kg
    side_accel_mps2 = q_pa * abs(loads.CSA_m2) / aero_case.mass_kg
    total_aero_g = sqrt(drag_accel_mps2^2 + lift_accel_mps2^2 + side_accel_mps2^2) / _earth_g()
    heating_proxy = sqrt(max(rho_kg_m3, 0.0)) * v_rel_mag^3
    kn = NaN
    horizontal_speed_mps = sqrt(v_e^2 + v_n^2)

    return (
        state=state,
        derivative=SVector{6, Float64}(v_vec[1], v_vec[2], v_vec[3], acc_cart[1], acc_cart[2], acc_cart[3]),
        altitude_km=h_m / 1e3,
        downrange_km=config.planet.Rp_m * lon / 1e3,
        crossrange_km=crossrange_km,
        heading_deg=rad2deg(atan(v_e, v_n)),
        velocity_mps=norm(v_vec),
        flight_path_angle_deg=rad2deg(atan(v_u, max(horizontal_speed_mps, 1.0))),
        rho_kg_m3=rho_kg_m3,
        temperature_k=temperature_k,
        q_pa=q_pa,
        drag_accel_mps2=drag_accel_mps2,
        lift_accel_mps2=lift_accel_mps2,
        side_accel_mps2=side_accel_mps2,
        total_aero_g=total_aero_g,
        heating_proxy=heating_proxy,
        mach=mach,
        beta_eff_kg_m2=loads.beta_eff_kg_m2,
        wind_east_mps=wind_enu[1],
        wind_north_mps=wind_enu[2],
        wind_up_mps=wind_enu[3],
        lateral_command=lateral_cmd,
        deployed=deployed,
        state_label=String(loads.state_label),
        kn=kn,
    )
end

function _rk4_step(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    aero_case::CalibratedAeroCase,
    policy::ControlPolicy,
    state::SVector{6, Float64},
    sample,
    time_s::Float64,
    dt_s::Float64,
    density_scale::Float64,
    use_winds::Bool,
    panel_command_override::Union{Nothing, DifferentialPanelCommand}=nothing,
)
    k1 = sample.derivative
    k2 = _sample_3d_state(
        config,
        adapter,
        aero_case,
        policy,
        state + 0.5 * dt_s * k1,
        time_s + 0.5 * dt_s,
        density_scale,
        use_winds,
        panel_command_override,
    ).derivative
    k3 = _sample_3d_state(
        config,
        adapter,
        aero_case,
        policy,
        state + 0.5 * dt_s * k2,
        time_s + 0.5 * dt_s,
        density_scale,
        use_winds,
        panel_command_override,
    ).derivative
    k4 = _sample_3d_state(
        config,
        adapter,
        aero_case,
        policy,
        state + dt_s * k3,
        time_s + dt_s,
        density_scale,
        use_winds,
        panel_command_override,
    ).derivative
    return state + (dt_s / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
end

function _push_trajectory_sample!(
    time_s_vec::Vector{Float64},
    altitude_km::Vector{Float64},
    downrange_km::Vector{Float64},
    crossrange_km::Vector{Float64},
    heading_deg_vec::Vector{Float64},
    velocity_mps::Vector{Float64},
    flight_path_angle_deg::Vector{Float64},
    rho_kg_m3::Vector{Float64},
    temperature_k::Vector{Float64},
    q_pa::Vector{Float64},
    drag_accel_mps2::Vector{Float64},
    lift_accel_mps2::Vector{Float64},
    side_accel_mps2::Vector{Float64},
    total_aero_g::Vector{Float64},
    heating_proxy::Vector{Float64},
    mach_values::Vector{Float64},
    beta_eff_kg_m2::Vector{Float64},
    wind_east_mps::Vector{Float64},
    wind_north_mps::Vector{Float64},
    wind_up_mps::Vector{Float64},
    lateral_command_vec::Vector{Float64},
    deployed_vec::Vector{Bool},
    state_label_vec::Vector{String},
    kn_values::Vector{Float64},
    time_s::Float64,
    sample,
)
    push!(time_s_vec, time_s)
    push!(altitude_km, sample.altitude_km)
    push!(downrange_km, sample.downrange_km)
    push!(crossrange_km, sample.crossrange_km)
    push!(heading_deg_vec, sample.heading_deg)
    push!(velocity_mps, sample.velocity_mps)
    push!(flight_path_angle_deg, sample.flight_path_angle_deg)
    push!(rho_kg_m3, sample.rho_kg_m3)
    push!(temperature_k, sample.temperature_k)
    push!(q_pa, sample.q_pa)
    push!(drag_accel_mps2, sample.drag_accel_mps2)
    push!(lift_accel_mps2, sample.lift_accel_mps2)
    push!(side_accel_mps2, sample.side_accel_mps2)
    push!(total_aero_g, sample.total_aero_g)
    push!(heating_proxy, sample.heating_proxy)
    push!(mach_values, sample.mach)
    push!(beta_eff_kg_m2, sample.beta_eff_kg_m2)
    push!(wind_east_mps, sample.wind_east_mps)
    push!(wind_north_mps, sample.wind_north_mps)
    push!(wind_up_mps, sample.wind_up_mps)
    push!(lateral_command_vec, sample.lateral_command)
    push!(deployed_vec, sample.deployed)
    push!(state_label_vec, sample.state_label)
    push!(kn_values, sample.kn)
    return nothing
end

function solve_entry_trajectory_3d(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    aero_case::CalibratedAeroCase,
    policy::ControlPolicy;
    save_trajectory::Bool=true,
    h0_m::Real=config.h0_m,
    V0_mps::Real=config.V0_mps,
    gamma0_deg::Real=config.gamma0_deg,
    heading_deg::Real=90.0,
    density_scale::Real=1.0,
    use_winds::Bool=config.monte_carlo_use_winds,
    panel_command_override::Union{Nothing, DifferentialPanelCommand}=nothing,
)::Trajectory3DResult
    r0 = SVector{3, Float64}(config.planet.Rp_m + Float64(h0_m), 0.0, 0.0)
    east0 = SVector{3, Float64}(0.0, 1.0, 0.0)
    north0 = SVector{3, Float64}(0.0, 0.0, 1.0)
    up0 = SVector{3, Float64}(1.0, 0.0, 0.0)
    γ0 = deg2rad(Float64(gamma0_deg))
    ψ0 = deg2rad(Float64(heading_deg))
    v_h = Float64(V0_mps) * cos(γ0)
    v_local0 = SVector{3, Float64}(v_h * sin(ψ0), v_h * cos(ψ0), Float64(V0_mps) * sin(γ0))
    v0 = v_local0[1] * east0 + v_local0[2] * north0 + v_local0[3] * up0
    current_state = SVector{6, Float64}(r0[1], r0[2], r0[3], v0[1], v0[2], v0[3])
    current_time_s = 0.0
    density_scale_float = Float64(density_scale)
    integration_dt_s = config.saveat_s

    time_s_vec = Float64[]
    altitude_km = Float64[]
    downrange_km = Float64[]
    crossrange_km = Float64[]
    heading_deg_vec = Float64[]
    velocity_mps = Float64[]
    flight_path_angle_deg = Float64[]
    rho_kg_m3 = Float64[]
    temperature_k = Float64[]
    q_pa = Float64[]
    drag_accel_mps2 = Float64[]
    lift_accel_mps2 = Float64[]
    side_accel_mps2 = Float64[]
    total_aero_g = Float64[]
    heating_proxy = Float64[]
    mach_values = Float64[]
    beta_eff_kg_m2 = Float64[]
    wind_east_mps = Float64[]
    wind_north_mps = Float64[]
    wind_up_mps = Float64[]
    lateral_command_vec = Float64[]
    deployed_vec = Bool[]
    state_label_vec = String[]
    kn_values = Float64[]

    current_sample = _sample_3d_state(
        config,
        adapter,
        aero_case,
        policy,
        current_state,
        current_time_s,
        density_scale_float,
        use_winds,
        panel_command_override,
    )
    peak_dynamic_pressure_pa = current_sample.q_pa
    peak_drag_accel_mps2 = current_sample.drag_accel_mps2
    peak_side_accel_mps2 = current_sample.side_accel_mps2
    peak_total_decel_earth_g = current_sample.total_aero_g
    peak_heating_proxy = current_sample.heating_proxy
    final_sample = current_sample
    final_time_s = current_time_s

    if save_trajectory
        _push_trajectory_sample!(
            time_s_vec,
            altitude_km,
            downrange_km,
            crossrange_km,
            heading_deg_vec,
            velocity_mps,
            flight_path_angle_deg,
            rho_kg_m3,
            temperature_k,
            q_pa,
            drag_accel_mps2,
            lift_accel_mps2,
            side_accel_mps2,
            total_aero_g,
            heating_proxy,
            mach_values,
            beta_eff_kg_m2,
            wind_east_mps,
            wind_north_mps,
            wind_up_mps,
            lateral_command_vec,
            deployed_vec,
            state_label_vec,
            kn_values,
            current_time_s,
            current_sample,
        )
    end

    integration_stopped_at_surface = false
    while current_time_s < config.max_time_s
        current_sample.altitude_km <= 0.0 && break
        dt_s = min(integration_dt_s, config.max_time_s - current_time_s)
        next_state = _rk4_step(
            config,
            adapter,
            aero_case,
            policy,
            current_state,
            current_sample,
            current_time_s,
            dt_s,
            density_scale_float,
            use_winds,
            panel_command_override,
        )
        next_time_s = current_time_s + dt_s
        next_sample = _sample_3d_state(
            config,
            adapter,
            aero_case,
            policy,
            next_state,
            next_time_s,
            density_scale_float,
            use_winds,
            panel_command_override,
        )

        if next_sample.altitude_km <= 0.0
            h_prev_m = current_sample.altitude_km * 1e3
            h_next_m = next_sample.altitude_km * 1e3
            frac = clamp(h_prev_m / max(h_prev_m - h_next_m, eps()), 0.0, 1.0)
            impact_time_s = current_time_s + frac * dt_s
            impact_state = current_state + frac * (next_state - current_state)
            impact_sample = _sample_3d_state(
                config,
                adapter,
                aero_case,
                policy,
                impact_state,
                impact_time_s,
                density_scale_float,
                use_winds,
                panel_command_override,
            )
            peak_dynamic_pressure_pa = max(peak_dynamic_pressure_pa, impact_sample.q_pa)
            peak_drag_accel_mps2 = max(peak_drag_accel_mps2, impact_sample.drag_accel_mps2)
            peak_side_accel_mps2 = max(peak_side_accel_mps2, impact_sample.side_accel_mps2)
            peak_total_decel_earth_g = max(peak_total_decel_earth_g, impact_sample.total_aero_g)
            peak_heating_proxy = max(peak_heating_proxy, impact_sample.heating_proxy)
            final_sample = impact_sample
            final_time_s = impact_time_s
            if save_trajectory
                _push_trajectory_sample!(
                    time_s_vec,
                    altitude_km,
                    downrange_km,
                    crossrange_km,
                    heading_deg_vec,
                    velocity_mps,
                    flight_path_angle_deg,
                    rho_kg_m3,
                    temperature_k,
                    q_pa,
                    drag_accel_mps2,
                    lift_accel_mps2,
                    side_accel_mps2,
                    total_aero_g,
                    heating_proxy,
                    mach_values,
                    beta_eff_kg_m2,
                    wind_east_mps,
                    wind_north_mps,
                    wind_up_mps,
                    lateral_command_vec,
                    deployed_vec,
                    state_label_vec,
                    kn_values,
                    impact_time_s,
                    impact_sample,
                )
            end
            integration_stopped_at_surface = true
            break
        end

        current_state = next_state
        current_time_s = next_time_s
        current_sample = next_sample
        peak_dynamic_pressure_pa = max(peak_dynamic_pressure_pa, current_sample.q_pa)
        peak_drag_accel_mps2 = max(peak_drag_accel_mps2, current_sample.drag_accel_mps2)
        peak_side_accel_mps2 = max(peak_side_accel_mps2, current_sample.side_accel_mps2)
        peak_total_decel_earth_g = max(peak_total_decel_earth_g, current_sample.total_aero_g)
        peak_heating_proxy = max(peak_heating_proxy, current_sample.heating_proxy)
        final_sample = current_sample
        final_time_s = current_time_s
        if save_trajectory
            _push_trajectory_sample!(
                time_s_vec,
                altitude_km,
                downrange_km,
                crossrange_km,
                heading_deg_vec,
                velocity_mps,
                flight_path_angle_deg,
                rho_kg_m3,
                temperature_k,
                q_pa,
                drag_accel_mps2,
                lift_accel_mps2,
                side_accel_mps2,
                total_aero_g,
                heating_proxy,
                mach_values,
                beta_eff_kg_m2,
                wind_east_mps,
                wind_north_mps,
                wind_up_mps,
                lateral_command_vec,
                deployed_vec,
                state_label_vec,
                kn_values,
                current_time_s,
                current_sample,
            )
        end
    end

    trajectory_df = save_trajectory ? DataFrame(
        time_s=time_s_vec,
        altitude_km=altitude_km,
        downrange_km=downrange_km,
        crossrange_km=crossrange_km,
        heading_deg=heading_deg_vec,
        velocity_mps=velocity_mps,
        flight_path_angle_deg=flight_path_angle_deg,
        rho_kg_m3=rho_kg_m3,
        temperature_k=temperature_k,
        q_pa=q_pa,
        drag_accel_mps2=drag_accel_mps2,
        lift_accel_mps2=lift_accel_mps2,
        side_accel_mps2=side_accel_mps2,
        total_aero_g=total_aero_g,
        heating_proxy=heating_proxy,
        mach=mach_values,
        beta_eff_kg_m2=beta_eff_kg_m2,
        wind_east_mps=wind_east_mps,
        wind_north_mps=wind_north_mps,
        wind_up_mps=wind_up_mps,
        lateral_command=lateral_command_vec,
        kn=kn_values,
        deployed=deployed_vec,
        state_label=state_label_vec,
    ) : DataFrame()

    impact_velocity_mps = final_sample.velocity_mps
    impact_g_load_proxy = impact_velocity_mps^2 / (2.0 * config.impact_stop_distance_m * _earth_g())
    summary = (
        atmosphere_path=atmosphere_label(adapter),
        impact_downrange_km=final_sample.downrange_km,
        impact_crossrange_km=final_sample.crossrange_km,
        flight_time_s=final_time_s,
        impact_velocity_mps=impact_velocity_mps,
        impact_flight_path_angle_deg=final_sample.flight_path_angle_deg,
        impact_heading_deg=final_sample.heading_deg,
        peak_dynamic_pressure_pa=peak_dynamic_pressure_pa,
        peak_drag_accel_mps2=peak_drag_accel_mps2,
        peak_side_accel_mps2=peak_side_accel_mps2,
        peak_total_decel_earth_g=peak_total_decel_earth_g,
        impact_g_load_proxy_earth_g=impact_g_load_proxy,
        max_heating_proxy=peak_heating_proxy,
        integration_stopped_at_surface=integration_stopped_at_surface || abs(final_sample.altitude_km) <= 0.25,
    )
    return Trajectory3DResult(summary, save_trajectory ? trajectory_df : DataFrame())
end
