struct ControlPolicy
    name::String
    kind::Symbol
    h_switch_m::Float64
    h_jettison_m::Float64
end

struct EntryProblem
    config::PrelimConfig
    adapter::MarsAtmosphereAdapter
    aero_case::CalibratedAeroCase
    policy::ControlPolicy
    density_scale::Float64
end

struct TrajectoryResult
    summary::NamedTuple
    trajectory::DataFrame
end

@inline function deployed_active(policy::ControlPolicy, h_m::Float64)
    if policy.kind === :body_only
        return false
    elseif policy.kind === :fixed_deployed || policy.kind === :fixed_deployed_trim
        return h_m > policy.h_jettison_m
    elseif policy.kind === :bang_bang || policy.kind === :bang_bang_trim
        h_m <= policy.h_jettison_m && return false
        return h_m <= policy.h_switch_m
    elseif policy.kind === :reverse_bang_bang || policy.kind === :reverse_bang_bang_trim
        h_m <= policy.h_jettison_m && return false
        return h_m > policy.h_switch_m
    end
    throw(ArgumentError("Unsupported control policy kind $(repr(policy.kind))."))
end

@inline function trim_active(policy::ControlPolicy, h_m::Float64)
    if policy.kind === :fixed_deployed_trim || policy.kind === :bang_bang_trim || policy.kind === :reverse_bang_bang_trim
        return deployed_active(policy, h_m)
    end
    return false
end

@inline function _earth_g()
    return 9.80665
end

function _entry_dynamics!(du, u, p::EntryProblem, t)
    r_m, theta_rad, V_mps, gamma_rad = u
    h_m = max(r_m - p.config.planet.Rp_m, 0.0)
    rho_kg_m3, temperature_k = density_temperature(p.adapter, h_m, theta_rad, t)
    rho_kg_m3 *= p.density_scale
    mach = mach_number(p.adapter, max(V_mps, 1.0), temperature_k)
    deployed = deployed_active(p.policy, h_m)
    trim = trim_active(p.policy, h_m)
    loads = aerodynamic_loads(p.aero_case, deployed, mach, p.config.planet.γ, trim)
    q_pa = 0.5 * rho_kg_m3 * V_mps^2
    drag_accel_mps2 = q_pa * loads.CDA_m2 / p.aero_case.mass_kg
    lift_accel_mps2 = q_pa * loads.CLA_m2 / p.aero_case.mass_kg
    gravity_mps2 = p.config.planet.μ / r_m^2

    du[1] = V_mps * sin(gamma_rad)
    du[2] = V_mps * cos(gamma_rad) / r_m
    du[3] = -drag_accel_mps2 - gravity_mps2 * sin(gamma_rad)
    du[4] = lift_accel_mps2 / max(V_mps, 1.0) + (V_mps / r_m - gravity_mps2 / max(V_mps, 1.0)) * cos(gamma_rad)
    return nothing
end

function _surface_callback(config::PrelimConfig)
    condition(u, _, _) = u[1] - config.planet.Rp_m
    affect!(integrator) = terminate!(integrator)
    return ContinuousCallback(condition, affect!; rootfind=true)
end

function _interp_descending(x_query::Float64, x_desc::AbstractVector{<:Real}, y_desc::AbstractVector{<:Real})
    n = length(x_desc)
    n == length(y_desc) || throw(ArgumentError("Interpolation vectors must have the same length."))
    x_query >= x_desc[1] && return Float64(y_desc[1])
    x_query <= x_desc[end] && return Float64(y_desc[end])
    @inbounds for i in 1:(n - 1)
        x_hi = Float64(x_desc[i])
        x_lo = Float64(x_desc[i + 1])
        if x_hi >= x_query >= x_lo
            y_hi = Float64(y_desc[i])
            y_lo = Float64(y_desc[i + 1])
            w = (x_query - x_lo) / max(x_hi - x_lo, eps())
            return y_lo + w * (y_hi - y_lo)
        end
    end
    return Float64(y_desc[end])
end

function solve_entry_trajectory(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    aero_case::CalibratedAeroCase,
    policy::ControlPolicy;
    save_trajectory::Bool=true,
    h0_m::Real=config.h0_m,
    V0_mps::Real=config.V0_mps,
    gamma0_deg::Real=config.gamma0_deg,
    theta0_rad::Real=config.theta0_rad,
    density_scale::Real=1.0,
)::TrajectoryResult
    u0 = [
        config.planet.Rp_m + Float64(h0_m),
        Float64(theta0_rad),
        Float64(V0_mps),
        deg2rad(Float64(gamma0_deg)),
    ]
    problem = ODEProblem(
        _entry_dynamics!,
        u0,
        (0.0, config.max_time_s),
        EntryProblem(config, adapter, aero_case, policy, Float64(density_scale)),
    )
    solution = solve(
        problem,
        Tsit5();
        abstol=config.solver_abstol,
        reltol=config.solver_reltol,
        callback=_surface_callback(config),
        saveat=config.saveat_s,
        save_start=true,
    )

    time_s = Float64.(solution.t)
    altitudes_km = Float64[]
    downrange_km = Float64[]
    velocities_mps = Float64[]
    gamma_deg = Float64[]
    rho_kg_m3 = Float64[]
    temperature_k = Float64[]
    q_pa = Float64[]
    drag_accel_mps2 = Float64[]
    lift_accel_mps2 = Float64[]
    total_aero_g = Float64[]
    heating_proxy = Float64[]
    mach_values = Float64[]
    beta_eff_kg_m2 = Float64[]
    cda_m2 = Float64[]
    cla_m2 = Float64[]
    kn_values = Float64[]
    deployed_state = Bool[]
    state_label = String[]

    for (idx, state) in enumerate(solution.u)
        r_m = Float64(state[1])
        theta_rad = Float64(state[2])
        V_mps = Float64(state[3])
        gamma_rad = Float64(state[4])
        h_m = max(r_m - config.planet.Rp_m, 0.0)
        rho_now, temperature_now = density_temperature(adapter, h_m, theta_rad, time_s[idx])
        rho_now *= Float64(density_scale)
        mach_now = mach_number(adapter, max(V_mps, 1.0), temperature_now)
        deployed_now = deployed_active(policy, h_m)
        trim_now = trim_active(policy, h_m)
        loads = aerodynamic_loads(aero_case, deployed_now, mach_now, config.planet.γ, trim_now)
        q_now = 0.5 * rho_now * V_mps^2
        drag_now = q_now * loads.CDA_m2 / aero_case.mass_kg
        lift_now = q_now * loads.CLA_m2 / aero_case.mass_kg
        total_now = sqrt(drag_now^2 + lift_now^2) / _earth_g()
        heat_now = sqrt(max(rho_now, 0.0)) * V_mps^3
        kn_now = NaN

        push!(altitudes_km, h_m / 1e3)
        push!(downrange_km, config.planet.Rp_m * theta_rad / 1e3)
        push!(velocities_mps, V_mps)
        push!(gamma_deg, rad2deg(gamma_rad))
        push!(rho_kg_m3, rho_now)
        push!(temperature_k, temperature_now)
        push!(q_pa, q_now)
        push!(drag_accel_mps2, drag_now)
        push!(lift_accel_mps2, lift_now)
        push!(total_aero_g, total_now)
        push!(heating_proxy, heat_now)
        push!(mach_values, mach_now)
        push!(beta_eff_kg_m2, loads.beta_eff_kg_m2)
        push!(cda_m2, loads.CDA_m2)
        push!(cla_m2, loads.CLA_m2)
        push!(kn_values, kn_now)
        push!(deployed_state, deployed_now)
        push!(state_label, String(loads.state_label))
    end

    trajectory_df = DataFrame(
        time_s=time_s,
        altitude_km=altitudes_km,
        downrange_km=downrange_km,
        velocity_mps=velocities_mps,
        flight_path_angle_deg=gamma_deg,
        rho_kg_m3=rho_kg_m3,
        temperature_k=temperature_k,
        q_pa=q_pa,
        drag_accel_mps2=drag_accel_mps2,
        lift_accel_mps2=lift_accel_mps2,
        total_aero_g=total_aero_g,
        heating_proxy=heating_proxy,
        mach=mach_values,
        beta_eff_kg_m2=beta_eff_kg_m2,
        CDA_m2=cda_m2,
        CLA_m2=cla_m2,
        kn=kn_values,
        deployed=deployed_state,
        state_label=state_label,
    )

    impact_velocity_mps = velocities_mps[end]
    impact_g_load_proxy = impact_velocity_mps^2 / (2.0 * config.impact_stop_distance_m * _earth_g())
    summary = (
        atmosphere_path=atmosphere_label(adapter),
        policy=policy.name,
        h_switch_m=policy.h_switch_m,
        h_jettison_m=policy.h_jettison_m,
        impact_downrange_km=downrange_km[end],
        flight_time_s=time_s[end],
        impact_velocity_mps=impact_velocity_mps,
        impact_flight_path_angle_deg=gamma_deg[end],
        peak_dynamic_pressure_pa=maximum(q_pa),
        peak_drag_accel_mps2=maximum(drag_accel_mps2),
        peak_total_decel_earth_g=maximum(total_aero_g),
        max_heating_proxy=maximum(heating_proxy),
        impact_g_load_proxy_earth_g=impact_g_load_proxy,
        integration_stopped_at_surface=abs(altitudes_km[end]) <= 0.25,
    )

    return TrajectoryResult(summary, save_trajectory ? trajectory_df : DataFrame())
end
