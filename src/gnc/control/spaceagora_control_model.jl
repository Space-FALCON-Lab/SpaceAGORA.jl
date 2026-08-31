#=
"""
    Runtime SpaceAGORA control hook.

    This hook solves or reuses the MPC area plan, applies active limits, and
    rotates the selected solar-panel links during run_simulation.
"""
=#
Base.@kwdef mutable struct AerobrakingMPCControlModel <: SimulationModel.AbstractControlEffectorModel
    config::AerobrakingMPCConfig
    state::AerobrakingMPCState = AerobrakingMPCState()
    controlled_panel_links::Tuple{Vararg{Int}}
    control_dt_s::Float64
    min_alpha_rad::Float64
    max_alpha_rad::Float64
    solve_interval_s::Float64
    build_reference_on_tick::Bool
    qp_max_nodes::Union{Nothing, Int}
    reference::AerobrakingMPCReferenceConfig
    fallback_area_m2::Float64
end

function _mpc_control_sat_state(u, i::Int)
    return hasproperty(u, :sc) ? u.sc[i] : u
end

function _mpc_control_pos_vel(sc)
    pos = hasproperty(sc, :pos) ? SVector{3, Float64}(sc.pos) : SVector{3, Float64}(sc[1], sc[2], sc[3])
    vel = hasproperty(sc, :vel) ? SVector{3, Float64}(sc.vel) : SVector{3, Float64}(sc[4], sc[5], sc[6])
    return pos, vel
end

function _mpc_apply_panel_alpha!(
    model::AerobrakingMPCControlModel,
    spacecraft,
    alpha::Float64,
)
    links = spacecraft.links
    for idx in model.controlled_panel_links
        1 <= idx <= length(links) || throw(ArgumentError("Controlled panel link $(idx) is out of bounds."))
        link = links[idx]
        link.root && throw(ArgumentError("Controlled panel link $(idx) is a root link and cannot be rotated."))
        axis = SVector{3, Float64}(abs.(link.r))
        SimulationModel.Kinematics.rotate_link(link, axis, pi / 2 - alpha)
        link.α = alpha
    end
    return nothing
end

function _mpc_desired_area(model::AerobrakingMPCControlModel)
    if model.config.mode isa MaxEnergyDepletionMode
        return model.config.bus_reference_area_m2 + model.config.controllable_area_m2
    end
    return model.fallback_area_m2
end

function _mpc_interpolated_plan_area(model::AerobrakingMPCControlModel, t::Float64)
    times = model.state.plan_time_s
    areas = model.state.plan_area_m2
    if !(isfinite(model.state.plan_epoch_s) && length(times) == length(areas) && !isempty(times))
        return _mpc_desired_area(model)
    end
    τ = max(0.0, t - model.state.plan_epoch_s)
    τ <= times[1] && return areas[1]
    τ >= times[end] && return areas[end]
    idx = searchsortedlast(times, τ)
    idx = clamp(idx, 1, length(times) - 1)
    t0 = times[idx]
    t1 = times[idx + 1]
    t1 <= t0 && return areas[idx]
    λ = (τ - t0) / (t1 - t0)
    return (1.0 - λ) * areas[idx] + λ * areas[idx + 1]
end

function _mpc_relative_speed(pos, vel, params::AerobrakingMPCParams)
    Ωx = @SMatrix [0.0 -params.Ω 0.0; params.Ω 0.0 0.0; 0.0 0.0 0.0]
    return norm(vel - Ωx * pos)
end

function _mpc_runtime_atmosphere(p, altitude_m::Float64, t::Float64)
    model = p.args.environment_model.density_model
    wind = p.args.environment_model.wind
    rho, temperature, _ = SimulationModel.getDensity(model, altitude_m, 0.0, 0.0, t, wind, p)
    return max(0.0, Float64(rho)), max(Float64(temperature), eps(Float64))
end

function _mpc_heat_rate_w_cm2(p, speed::Float64, density::Float64, temperature::Float64, alpha::Float64)
    planet = p.args.environment_model.planet
    sound_speed = sqrt(max(0.0, planet.γ * planet.R * temperature))
    if !(isfinite(speed) && isfinite(sound_speed) && speed > 0.0 && sound_speed > 0.0 && density > 0.0)
        return 0.0
    end
    mach = speed / sound_speed
    S = sqrt(0.5 * planet.γ) * mach
    qdot = SimulationModel.getHeatRate(
        p.args.environment_model.thermal_model,
        S,
        temperature,
        density,
        speed,
        alpha,
    )
    return isfinite(qdot) ? max(0.0, Float64(qdot)) : 0.0
end

function _mpc_alpha_for_heat_rate_limit(
    model::AerobrakingMPCControlModel,
    p,
    speed::Float64,
    density::Float64,
    temperature::Float64,
    desired_alpha::Float64,
    limit_w_cm2::Float64,
)
    min_alpha = model.min_alpha_rad
    q_min = _mpc_heat_rate_w_cm2(p, speed, density, temperature, min_alpha)
    q_max = _mpc_heat_rate_w_cm2(p, speed, density, temperature, desired_alpha)
    q_max <= limit_w_cm2 && return desired_alpha
    q_min >= limit_w_cm2 && return min_alpha
    lo = min_alpha
    hi = desired_alpha
    for _ in 1:48
        mid = 0.5 * (lo + hi)
        q_mid = _mpc_heat_rate_w_cm2(p, speed, density, temperature, mid)
        if q_mid <= limit_w_cm2
            lo = mid
        else
            hi = mid
        end
    end
    return lo
end

function _mpc_area_command_from_constraints!(
    model::AerobrakingMPCControlModel,
    p,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    t::Float64,
    desired_area_m2::Union{Nothing, Float64}=nothing,
)
    config = model.config
    state = model.state
    params = mpc_params_from_spaceagora(p.args)
    min_area = config.bus_reference_area_m2
    max_area = config.bus_reference_area_m2 + config.controllable_area_m2
    desired_area = desired_area_m2 === nothing ? _mpc_desired_area(model) : clamp(desired_area_m2, min_area, max_area)
    area_limit = max_area
    active_limit = :none

    altitude_m = norm(pos) - params.Re
    density, temperature = _mpc_runtime_atmosphere(p, altitude_m, t)
    speed = _mpc_relative_speed(pos, vel, params)
    dynamic_pressure = 0.5 * density * speed^2
    desired_alpha = alpha_from_commanded_area(
        config,
        desired_area;
        min_alpha_rad=model.min_alpha_rad,
        max_alpha_rad=model.max_alpha_rad,
    )

    if config.use_constraints && config.use_drag_constraint && dynamic_pressure > eps(Float64)
        drag_area_limit = config.drag_max_n / (dynamic_pressure * config.drag_coefficient)
        if drag_area_limit < area_limit
            area_limit = drag_area_limit
            active_limit = :drag
        end
    end

    if config.use_constraints && config.use_qdot_constraint
        alpha_limit = _mpc_alpha_for_heat_rate_limit(
            model,
            p,
            speed,
            density,
            temperature,
            desired_alpha,
            config.qdot_max_w_cm2,
        )
        heat_area_limit = commanded_area_from_alpha(
            config,
            alpha_limit;
            min_alpha_rad=model.min_alpha_rad,
            max_alpha_rad=model.max_alpha_rad,
        )
        if heat_area_limit < area_limit
            area_limit = heat_area_limit
            active_limit = :heat_rate
        end
    end

    dt = isfinite(state.last_command_time_s) ? max(0.0, t - state.last_command_time_s) : 0.0
    if dt > 0.0
        previous_heat_rate = isfinite(state.last_heat_rate_w_cm2) ? max(0.0, state.last_heat_rate_w_cm2) : 0.0
        state.estimated_heat_load_j_cm2 += previous_heat_rate * dt
    end
    if config.use_constraints && config.use_heat_load_constraint && dt > 0.0
        remaining_load = config.heat_load_max_j_cm2 - state.estimated_heat_load_j_cm2
        if remaining_load <= 0.0
            area_limit = min(area_limit, min_area)
            active_limit = :heat_load
        else
            allowable_rate = remaining_load / dt
            alpha_limit = _mpc_alpha_for_heat_rate_limit(
                model,
                p,
                speed,
                density,
                temperature,
                desired_alpha,
                allowable_rate,
            )
            heat_load_area_limit = commanded_area_from_alpha(
                config,
                alpha_limit;
                min_alpha_rad=model.min_alpha_rad,
                max_alpha_rad=model.max_alpha_rad,
            )
            if heat_load_area_limit < area_limit
                area_limit = heat_load_area_limit
                active_limit = :heat_load
            end
        end
    end

    commanded_area = clamp(min(desired_area, area_limit), min_area, max_area)
    if config.use_constraints && config.use_slew_constraint &&
            isfinite(state.held_commanded_area_m2) && dt > 0.0
        slew = config.area_slew_max_m2_s * dt
        commanded_area = clamp(
            commanded_area,
            state.held_commanded_area_m2 - slew,
            state.held_commanded_area_m2 + slew,
        )
        commanded_area = clamp(commanded_area, min_area, max_area)
        if abs(commanded_area - min(desired_area, area_limit)) > 1.0e-9
            active_limit = :slew
        end
    end

    alpha = alpha_from_commanded_area(
        config,
        commanded_area;
        min_alpha_rad=model.min_alpha_rad,
        max_alpha_rad=model.max_alpha_rad,
    )
    heat_rate = _mpc_heat_rate_w_cm2(p, speed, density, temperature, alpha)
    drag = dynamic_pressure * config.drag_coefficient * commanded_area

    state.held_commanded_area_m2 = commanded_area
    state.held_alpha_rad = alpha
    state.last_command_time_s = t
    state.last_heat_rate_w_cm2 = heat_rate
    state.last_drag_n = drag
    state.last_density_kg_m3 = density
    state.active_limit = active_limit
    state.selected_mode = objective_kind(config.mode)
    state.solver_status = active_limit == :none ? :area_command : :constraint_limited_area_command
    return commanded_area
end

function _mpc_update_plan!(
    model::AerobrakingMPCControlModel,
    u,
    p,
    t::Float64,
    sat_idx::Int,
)
    params = mpc_params_from_spaceagora(p.args)
    sc_state = _mpc_control_sat_state(u, sat_idx)
    pos, vel = _mpc_control_pos_vel(sc_state)
    if t - model.state.last_solve_time_s < model.solve_interval_s
        desired_area = _mpc_interpolated_plan_area(model, t)
        _mpc_area_command_from_constraints!(model, p, pos, vel, t, desired_area)
        return nothing
    end
    density = density_function_from_spaceagora(p.args)
    model.state.last_solve_time_s = t
    if !model.build_reference_on_tick
        _mpc_area_command_from_constraints!(model, p, pos, vel, t)
        model.state.plan_epoch_s = NaN
        model.state.plan_time_s = Float64[]
        model.state.plan_area_m2 = Float64[]
        model.state.predicted_terminal_energy = NaN
        model.state.last_solution = nothing
        return nothing
    end
    try
        nominal_area = _mpc_desired_area(model)
        ref = build_reference_drag_pass(
            params,
            pos,
            vel;
            config=model.config,
            reference=model.reference,
            nominal_area_m2=nominal_area,
            density=density,
        )
        problem = build_mpc_problem(
            ref,
            params,
            model.config;
            density=density,
            max_nodes=model.qp_max_nodes,
        )
        solution = solve_mpc_qp(problem, model.config)
        model.state.last_solution = solution
        if solution.ok && !isempty(solution.commanded_area_m2)
            plan_t = collect(problem.t .- first(problem.t))
            plan_area = collect(solution.commanded_area_m2)
            area = _mpc_area_command_from_constraints!(
                model,
                p,
                pos,
                vel,
                t,
                first(plan_area),
            )
            model.state.plan_epoch_s = t
            model.state.plan_time_s = plan_t
            model.state.plan_area_m2 = plan_area
            model.state.predicted_terminal_energy = solution.terminal_energy
            model.state.solver_status = solution.solver_status
        else
            area = _mpc_area_command_from_constraints!(model, p, pos, vel, t)
            model.state.plan_epoch_s = t
            model.state.plan_time_s = collect(problem.t .- first(problem.t))
            model.state.plan_area_m2 = fill(area, length(problem.t))
            model.state.predicted_terminal_energy = solution.terminal_energy
            model.state.solver_status = solution.solver_status
        end
    catch err
        _mpc_area_command_from_constraints!(model, p, pos, vel, t)
        model.state.solver_status = :reference_or_qp_failed
        model.state.plan_epoch_s = NaN
        model.state.plan_time_s = Float64[]
        model.state.plan_area_m2 = Float64[]
        model.state.last_solution = nothing
    end
    return nothing
end

function SimulationModel.calcControlEffect!(
    model::AerobrakingMPCControlModel,
    u,
    p,
    t::Float64,
    sat_idx::Int64,
)
    _mpc_update_plan!(model, u, p, Float64(t), Int(sat_idx))
    area = isfinite(model.state.held_commanded_area_m2) ?
        model.state.held_commanded_area_m2 :
        _mpc_desired_area(model)
    alpha = alpha_from_commanded_area(
        model.config,
        area;
        min_alpha_rad=model.min_alpha_rad,
        max_alpha_rad=model.max_alpha_rad,
    )
    model.state.held_commanded_area_m2 = area
    model.state.held_alpha_rad = alpha
    spacecraft = p.args.dynamics_model.spacecraft[Int(sat_idx)]
    _mpc_apply_panel_alpha!(model, spacecraft, alpha)
    return nothing
end

function SimulationModel.calcControlForceTorque(
    model::AerobrakingMPCControlModel,
    u::AbstractVector,
    p,
    i::Int64,
    t::Float64,
)
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function _mpc_save_control(model::AerobrakingMPCControlModel)
    area = isfinite(model.state.held_commanded_area_m2) ?
        model.state.held_commanded_area_m2 :
        _mpc_desired_area(model)
    alpha = isfinite(model.state.held_alpha_rad) ?
        model.state.held_alpha_rad :
        alpha_from_commanded_area(
            model.config,
            area;
            min_alpha_rad=model.min_alpha_rad,
            max_alpha_rad=model.max_alpha_rad,
        )
    return area, alpha
end

function _mpc_limit_code(limit::Symbol)::Float64
    limit == :none && return 0.0
    limit == :drag && return 1.0
    limit == :heat_rate && return 2.0
    limit == :heat_load && return 3.0
    limit == :slew && return 4.0
    return -1.0
end

function mpc_control_save_fields(model::AerobrakingMPCControlModel)
    return (
        SimulationModel.SaveField(
            :mpc_commanded_area_m2,
            (u, t, integrator) -> [
                _mpc_save_control(integrator.p.args.control_model.control_effectors[1])[1]
                for _ in eachindex(integrator.p.args.dynamics_model.spacecraft)
            ];
            per_satellite=true,
            column_prefix="mpc_commanded_area_m2",
        ),
        SimulationModel.SaveField(
            :mpc_alpha_rad,
            (u, t, integrator) -> [
                _mpc_save_control(integrator.p.args.control_model.control_effectors[1])[2]
                for _ in eachindex(integrator.p.args.dynamics_model.spacecraft)
            ];
            per_satellite=true,
            column_prefix="mpc_alpha_rad",
        ),
        SimulationModel.SaveField(
            :mpc_active_limit_code,
            (u, t, integrator) -> [
                _mpc_limit_code(integrator.p.args.control_model.control_effectors[1].state.active_limit)
                for _ in eachindex(integrator.p.args.dynamics_model.spacecraft)
            ];
            per_satellite=true,
            column_prefix="mpc_active_limit_code",
        ),
        SimulationModel.SaveField(
            :mpc_heat_rate_w_cm2,
            (u, t, integrator) -> [
                integrator.p.args.control_model.control_effectors[1].state.last_heat_rate_w_cm2
                for _ in eachindex(integrator.p.args.dynamics_model.spacecraft)
            ];
            per_satellite=true,
            column_prefix="mpc_heat_rate_w_cm2",
        ),
        SimulationModel.SaveField(
            :mpc_drag_n,
            (u, t, integrator) -> [
                integrator.p.args.control_model.control_effectors[1].state.last_drag_n
                for _ in eachindex(integrator.p.args.dynamics_model.spacecraft)
            ];
            per_satellite=true,
            column_prefix="mpc_drag_n",
        ),
        SimulationModel.SaveField(
            :mpc_heat_load_j_cm2,
            (u, t, integrator) -> [
                integrator.p.args.control_model.control_effectors[1].state.estimated_heat_load_j_cm2
                for _ in eachindex(integrator.p.args.dynamics_model.spacecraft)
            ];
            per_satellite=true,
            column_prefix="mpc_heat_load_j_cm2",
        ),
    )
end
