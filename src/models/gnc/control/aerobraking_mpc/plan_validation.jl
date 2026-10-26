"""
    interpolate_mpc_plan(times_s, areas_m2, time_s)

Linearly interpolate an MPC exposed-area plan in physical time. This is the
same interpolation convention used by `AerobrakingMPCControlModel` at runtime.
"""
function interpolate_mpc_plan(times_s, areas_m2, time_s::Real)
    length(times_s) == length(areas_m2) || throw(DimensionMismatch(
        "MPC plan time and area vectors must have equal length."))
    isempty(times_s) && throw(ArgumentError("MPC plan cannot be empty."))
    t = Float64(time_s)
    t <= times_s[1] && return Float64(areas_m2[1])
    t >= times_s[end] && return Float64(areas_m2[end])
    index = clamp(searchsortedlast(times_s, t), 1, length(times_s) - 1)
    fraction = (t - times_s[index]) / (times_s[index + 1] - times_s[index])
    return muladd(fraction, areas_m2[index + 1] - areas_m2[index], areas_m2[index])
end

"""Linearly sample a scalar history on a monotonically increasing time grid."""
function interpolate_mpc_history(times_s, values, time_s::Real)
    return interpolate_mpc_plan(times_s, values, time_s)
end

"""
    evaluate_cartesian_mpc_outputs(position, velocity, elapsed_time_s,
        area_m2, params, config; density)

Evaluate the physical output convention used by the KS aerobraking MPC from
an accepted Cartesian state.
"""
function evaluate_cartesian_mpc_outputs(
    position_ii_m,
    velocity_ii_m,
    elapsed_time_s::Real,
    area_m2::Real,
    params::AerobrakingMPCParams,
    config::AerobrakingMPCConfig;
    density::Function,
)
    position = SVector{3, Float64}(position_ii_m)
    velocity = SVector{3, Float64}(velocity_ii_m)
    radius = norm(position)
    altitude = radius - params.Re
    rho = max(0.0, Float64(applicable(density, altitude, elapsed_time_s,
        position) ? density(altitude, elapsed_time_s, position) :
        density(altitude, elapsed_time_s)))
    relative_velocity = velocity - ks_rotation_cross_matrix(params) * position
    speed = norm(relative_velocity)
    area = Float64(area_m2)
    panel_fraction = commanded_area_fraction(config, area)
    drag = 0.5 * rho * speed^2 * config.drag_coefficient * area
    heat_rate = 0.5 * rho * speed^3 * panel_fraction / 1.0e4
    energy = 0.5 * dot(velocity, velocity) - params.μ / radius
    return (
        altitude_m=altitude,
        density_kg_m3=rho,
        relative_speed_m_s=speed,
        drag_n=drag,
        heat_rate_w_cm2=heat_rate,
        specific_energy_j_kg=energy,
        specific_energy_mj_kg=energy / 1.0e6,
    )
end

"""
    evaluate_ks_mpc_outputs(state, area_m2, params, config; density)

Evaluate the four physical outputs used by the KS aerobraking MPC from a
nonlinear KS state: altitude, drag force, panel kinetic-energy-flux heat rate,
and two-body specific orbital energy. The heat-rate output uses W/cm², matching
the public `AerobrakingMPCConfig` convention.
"""
function evaluate_ks_mpc_outputs(
    state,
    area_m2::Real,
    params::AerobrakingMPCParams,
    config::AerobrakingMPCConfig;
    density::Function,
)
    cartesian = ks_state_to_cartesian(state)
    return evaluate_cartesian_mpc_outputs(
        cartesian.position_ii_m,
        cartesian.velocity_ii_m,
        state[10],
        area_m2,
        params,
        config;
        density=density,
    )
end

"""Trapezoidally integrate a W/cm² heat-rate history to J/cm²."""
function cumulative_mpc_heat_load(heat_rate_w_cm2, times_s)
    length(heat_rate_w_cm2) == length(times_s) || throw(DimensionMismatch(
        "Heat-rate and time vectors must have equal length."))
    result = zeros(Float64, length(times_s))
    for index in 2:length(times_s)
        result[index] = result[index - 1] + 0.5 *
            (heat_rate_w_cm2[index - 1] + heat_rate_w_cm2[index]) *
            (times_s[index] - times_s[index - 1])
    end
    return result
end

"""
    propagate_ks_mpc_plan(reference, problem, area_plan, config, params; kwargs...)

Independently propagate an untouched MPC area plan through the reusable
nonlinear KS dynamics. This is a validation utility; the operational
SpaceAGORA callback remains `AerobrakingMPCControlModel`.
"""
function propagate_ks_mpc_plan(
    reference,
    problem::AerobrakingMPCProblem,
    area_plan,
    config::AerobrakingMPCConfig,
    params::AerobrakingMPCParams;
    density::Function,
    delta_s::Real=reference.delta_s,
    max_steps::Integer=20_000,
    cutoff_altitude_m::Real=reference.cutoff_altitude_m,
)
    length(area_plan) == problem.N || throw(DimensionMismatch(
        "Area plan has $(length(area_plan)) values but the MPC problem has $(problem.N) nodes."))
    all(diff(problem.t) .> 0.0) || throw(ArgumentError(
        "MPC problem times must be strictly increasing."))
    step_s = Float64(delta_s)
    step_s > 0.0 || throw(ArgumentError("KS fictitious-time step must be positive."))

    states = Vector{Vector{Float64}}()
    areas = Float64[]
    state = collect(reference.states[1, :])
    push!(states, copy(state))
    push!(areas, interpolate_mpc_plan(problem.t, area_plan, state[10]))
    exited = false
    for _ in 1:Int(max_steps)
        area = interpolate_mpc_plan(problem.t, area_plan, state[10])
        output = evaluate_ks_mpc_outputs(
            state, area, params, config; density=density)
        next_state = ks_rk4_step(
            state,
            params,
            area,
            step_s;
            density_kg_m3=density,
            drag_coefficient=config.drag_coefficient,
            mass_kg=config.mass_kg,
            use_drag=true,
        )
        next_output = evaluate_ks_mpc_outputs(
            next_state,
            interpolate_mpc_plan(problem.t, area_plan, next_state[10]),
            params,
            config;
            density=density,
        )
        next_cartesian = ks_state_to_cartesian(next_state)
        next_radial_velocity = dot(
            next_cartesian.position_ii_m, next_cartesian.velocity_ii_m) /
            norm(next_cartesian.position_ii_m)
        if output.altitude_m <= cutoff_altitude_m &&
                next_output.altitude_m > cutoff_altitude_m &&
                next_radial_velocity > 0.0
            fraction = (cutoff_altitude_m - output.altitude_m) /
                (next_output.altitude_m - output.altitude_m)
            state = state .+ fraction .* (next_state .- state)
            push!(states, copy(state))
            push!(areas, interpolate_mpc_plan(problem.t, area_plan, state[10]))
            exited = true
            break
        end
        all(isfinite, next_state) || throw(ErrorException(
            "Nonlinear KS plan propagation became non-finite."))
        state = next_state
        push!(states, copy(state))
        push!(areas, interpolate_mpc_plan(problem.t, area_plan, state[10]))
    end
    exited || throw(ErrorException(
        "Nonlinear KS plan propagation did not reach the outbound cutoff altitude."))

    state_matrix = reduce(vcat, transpose.(states))
    absolute_times = state_matrix[:, 10]
    outputs = [evaluate_ks_mpc_outputs(
        view(state_matrix, index, :), areas[index], params, config;
        density=density) for index in axes(state_matrix, 1)]
    altitude = getproperty.(outputs, :altitude_m)
    rho = getproperty.(outputs, :density_kg_m3)
    drag = getproperty.(outputs, :drag_n)
    heat_rate = getproperty.(outputs, :heat_rate_w_cm2)
    energy = getproperty.(outputs, :specific_energy_mj_kg)
    heat_load = cumulative_mpc_heat_load(heat_rate, absolute_times)
    slew = vcat(NaN, diff(areas) ./ diff(absolute_times))
    return (
        time_s=absolute_times .- first(absolute_times),
        absolute_time_s=absolute_times,
        states=state_matrix,
        altitude_m=altitude,
        density_kg_m3=rho,
        area_m2=areas,
        heat_rate_w_cm2=heat_rate,
        heat_load_j_cm2=heat_load,
        drag_n=drag,
        slew_m2_s=slew,
        energy_mj_kg=energy,
    )
end
