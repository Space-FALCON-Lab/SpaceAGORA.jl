# Runtime control hook for MPC planning, constraint enforcement, and panel actuation.
"""
Aerobraking MPC control effector for one spacecraft.

At its control ticks it propagates a nominal KS drag pass from the current
state, solves the condensed QP ([`solve_mpc_qp`](@ref)) and follows the
resulting exposed-area plan with the panels in `controlled_panel_links`,
clamping each command to the runtime drag, heat-rate, heat-load and slew
limits. With `solve_trigger_altitude_m` it solves once per pass, when the
spacecraft descends through that altitude; otherwise at most every
`solve_interval_s` seconds, except while it climbs out inside the reference
interface (`reference.h_cut_m`). The prediction runs in a non-rotating frame
aligned with the planet-fixed frame at the solve time, so J2 and the
co-rotating atmosphere act about the planet's pole, and the atmosphere is
sampled at the run's time and planet-fixed position. The heat-load estimate
restarts at each atmosphere entry. Declare it in a `ControlModel` or on the
spacecraft.
"""
Base.@kwdef mutable struct AerobrakingMPCControlModel <: AbstractControlEffectorModel
    config::AerobrakingMPCConfig
    state::AerobrakingMPCState = AerobrakingMPCState()
    spacecraft_index::Int
    controlled_panel_links::Tuple{Vararg{Int}}
    min_alpha_rad::Float64
    max_alpha_rad::Float64
    solve_interval_s::Float64
    build_reference_on_tick::Bool
    qp_max_nodes::Union{Nothing, Int}
    reference::AerobrakingMPCReferenceConfig
    fallback_area_m2::Float64
    prediction_latitude_rad::Float64
    prediction_longitude_rad::Float64
    prediction_wind::Bool
    solve_trigger_altitude_m::Union{Nothing, Float64} = nothing
    heat_rate_model::Symbol = :spaceagora_thermal
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
    t::Float64,
)
    # Shared panel actuation: rotates hinge-less panels directly and moves
    # mounted hinges within their limits, as the EDG controller does.
    _apply_solar_panel_aoa!(
        SolarPanelAngleOfAttackControlModel(model.controlled_panel_links),
        spacecraft, alpha, t)
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

# Rotation from the inertial (J2000) frame to the planet-fixed frame at run time
# `t`, as the engine uses for the atmosphere; identity outside a run.
function _mpc_planet_frame_rotation(p, t::Float64)::SMatrix{3, 3, Float64, 9}
    hasproperty(p, :shared_buffers) || return SMatrix{3, 3, Float64, 9}(I)
    environment = p.args.environment_model
    return ephemerides_frame_lpi(environment.ephemerides_model, environment.planet,
        p.shared_buffers.et_start[] + t)
end

# The planet-fixed position and planet-relative velocity of an inertial state.
function _mpc_planet_fixed_state(p, pos::SVector{3, Float64},
    vel::SVector{3, Float64}, t::Float64)
    rotation = _mpc_planet_frame_rotation(p, t)
    omega = _planet_rotation_rate_rad_s(p.args.environment_model.planet)
    pos_pf = rotation * pos
    vel_pf = rotation * vel - cross(SVector{3, Float64}(0.0, 0.0, omega), pos_pf)
    return pos_pf, vel_pf, rotation
end

function _mpc_runtime_atmosphere(model::AerobrakingMPCControlModel, p,
    pos::SVector{3, Float64}, vel::SVector{3, Float64}, t::Float64)
    density_model = p.args.environment_model.density_model
    wind = p.args.environment_model.wind
    pos_pf, vel_pf, _ = _mpc_planet_fixed_state(p, pos, vel, t)
    altitude_m, latitude_rad, longitude_rad =
        rtolatlong(pos_pf, p.args.environment_model.planet)
    rho, temperature, _ = getDensity(
        density_model,
        altitude_m,
        latitude_rad,
        longitude_rad,
        t,
        model.prediction_wind && wind,
        p,
    )
    return max(0.0, Float64(rho)), max(Float64(temperature), eps(Float64)), norm(vel_pf)
end

function _mpc_heat_rate_w_cm2(model::AerobrakingMPCControlModel, p,
    speed::Float64, density::Float64, temperature::Float64, alpha::Float64)
    if model.heat_rate_model === :kinetic_energy_flux
        panel_fraction = sin(clamp(alpha, 0.0, pi / 2))
        return max(0.0, 0.5 * density * speed^3 * panel_fraction / 1.0e4)
    end
    model.heat_rate_model === :spaceagora_thermal || throw(ArgumentError(
        "Unsupported MPC runtime heat-rate model $(repr(model.heat_rate_model)). " *
        "Use :spaceagora_thermal or :kinetic_energy_flux."))
    planet = p.args.environment_model.planet
    sound_speed = sqrt(max(0.0, planet.γ * planet.R * temperature))
    if !(isfinite(speed) && isfinite(sound_speed) && speed > 0.0 && sound_speed > 0.0 && density > 0.0)
        return 0.0
    end
    mach = speed / sound_speed
    S = sqrt(0.5 * planet.γ) * mach
    get_heat_rate = getproperty(parentmodule(@__MODULE__), :getHeatRate)
    qdot = get_heat_rate(
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
    q_min = _mpc_heat_rate_w_cm2(
        model, p, speed, density, temperature, min_alpha)
    q_max = _mpc_heat_rate_w_cm2(
        model, p, speed, density, temperature, desired_alpha)
    q_max <= limit_w_cm2 && return desired_alpha
    q_min >= limit_w_cm2 && return min_alpha
    lo = min_alpha
    hi = desired_alpha
    for _ in 1:48
        mid = 0.5 * (lo + hi)
        q_mid = _mpc_heat_rate_w_cm2(
            model, p, speed, density, temperature, mid)
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

    density, temperature, speed = _mpc_runtime_atmosphere(model, p, pos, vel, t)
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
    heat_rate = _mpc_heat_rate_w_cm2(
        model, p, speed, density, temperature, alpha)
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
    if model.solve_trigger_altitude_m !== nothing
        altitude_m = norm(pos) - params.Re
        radial_velocity_m_s = dot(pos, vel) / norm(pos)
        if altitude_m > model.solve_trigger_altitude_m
            model.state.solve_armed = true
            desired_area = _mpc_interpolated_plan_area(model, t)
            _mpc_area_command_from_constraints!(model, p, pos, vel, t, desired_area)
            return nothing
        elseif !(model.state.solve_armed && radial_velocity_m_s < 0.0)
            desired_area = _mpc_interpolated_plan_area(model, t)
            _mpc_area_command_from_constraints!(model, p, pos, vel, t, desired_area)
            return nothing
        end
    end
    # A reference from inside the interface on the way out would coast to the next
    # pass; keep following the current plan until the spacecraft has left.
    outbound_inside = model.solve_trigger_altitude_m === nothing &&
        norm(pos) - params.Re <= model.reference.h_cut_m && dot(pos, vel) >= 0.0
    if outbound_inside || t - model.state.last_solve_time_s < model.solve_interval_s
        desired_area = _mpc_interpolated_plan_area(model, t)
        _mpc_area_command_from_constraints!(model, p, pos, vel, t, desired_area)
        return nothing
    end
    density = density_function_from_spaceagora(
        p;
        latitude=model.prediction_latitude_rad,
        longitude=model.prediction_longitude_rad,
        wind=model.prediction_wind,
        epoch_s=t,
    )
    model.state.last_solve_time_s = t
    model.state.solve_armed = false
    model.state.solve_count += 1
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
        # Prediction frame: inertial axes turned onto the planet-fixed axes at `t`.
        rotation = _mpc_planet_frame_rotation(p, t)
        ref = build_reference_drag_pass(
            params,
            rotation * pos,
            rotation * vel;
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
        model.state.last_error = nothing
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
            model.state.plan_epoch_s = t + first(problem.t)
            model.state.plan_time_s = plan_t
            model.state.plan_area_m2 = plan_area
            model.state.predicted_terminal_energy = solution.terminal_energy
            model.state.solver_status = solution.solver_status
        else
            area = _mpc_area_command_from_constraints!(model, p, pos, vel, t)
            model.state.plan_epoch_s = t + first(problem.t)
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
        model.state.last_error = sprint(showerror, err)
    end
    return nothing
end

function update_control!(
    model::AerobrakingMPCControlModel,
    u,
    p,
    t::Float64,
    sat_idx::Int64,
)
    Int(sat_idx) == model.spacecraft_index || return nothing
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
    _mpc_apply_panel_alpha!(model, spacecraft, alpha, Float64(t))
    return nothing
end

function control_wrench(
    model::AerobrakingMPCControlModel,
    x::StateSample,
    p,
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

"""
    mpc_control_save_fields(model) -> Tuple{SaveField...}

Per-spacecraft result columns for an [`AerobrakingMPCControlModel`](@ref): the
commanded area and panel angle, the limiting constraint code (0 none, 1 drag,
2 heat rate, 3 heat load, 4 slew), and the runtime heat rate, drag and heat load.
"""
function mpc_control_save_fields(model::AerobrakingMPCControlModel)
    SaveField = getproperty(parentmodule(@__MODULE__), :SaveField)
    return (
        SaveField(
            :mpc_commanded_area_m2,
            (u, t, integrator) -> [
                i == model.spacecraft_index ? _mpc_save_control(model)[1] : NaN
                for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
            ];
            per_satellite=true,
            column_prefix="mpc_commanded_area_m2",
        ),
        SaveField(
            :mpc_alpha_rad,
            (u, t, integrator) -> [
                i == model.spacecraft_index ? _mpc_save_control(model)[2] : NaN
                for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
            ];
            per_satellite=true,
            column_prefix="mpc_alpha_rad",
        ),
        SaveField(
            :mpc_active_limit_code,
            (u, t, integrator) -> [
                i == model.spacecraft_index ? _mpc_limit_code(model.state.active_limit) : NaN
                for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
            ];
            per_satellite=true,
            column_prefix="mpc_active_limit_code",
        ),
        SaveField(
            :mpc_heat_rate_w_cm2,
            (u, t, integrator) -> [
                i == model.spacecraft_index ? model.state.last_heat_rate_w_cm2 : NaN
                for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
            ];
            per_satellite=true,
            column_prefix="mpc_heat_rate_w_cm2",
        ),
        SaveField(
            :mpc_drag_n,
            (u, t, integrator) -> [
                i == model.spacecraft_index ? model.state.last_drag_n : NaN
                for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
            ];
            per_satellite=true,
            column_prefix="mpc_drag_n",
        ),
        SaveField(
            :mpc_heat_load_j_cm2,
            (u, t, integrator) -> [
                i == model.spacecraft_index ? model.state.estimated_heat_load_j_cm2 : NaN
                for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
            ];
            per_satellite=true,
            column_prefix="mpc_heat_load_j_cm2",
        ),
    )
end

# GNC runtime contract: the MPC acts on one spacecraft, holds an exposed-area and
# panel-angle command between ticks, and can be declared on a spacecraft.
GNCRuntime.spacecraft_scope(m::AerobrakingMPCControlModel) = BoundSpacecraft((m.spacecraft_index,))
GNCRuntime.held_command(m::AerobrakingMPCControlModel, spacecraft_index) =
    spacecraft_index == m.spacecraft_index ?
        (commanded_area_m2=m.state.held_commanded_area_m2, alpha_rad=m.state.held_alpha_rad) : nothing
GNCRuntime.held_command_time(m::AerobrakingMPCControlModel, spacecraft_index) =
    spacecraft_index == m.spacecraft_index ? m.state.last_command_time_s : NaN

function SimulationLifecycle.bind_spacecraft(m::AerobrakingMPCControlModel, sat_idx::Int)
    m.spacecraft_index == sat_idx && return m
    b = SimulationLifecycle._shallow_copy(m)
    b.spacecraft_index = sat_idx
    b.state = deepcopy(m.state)   # its own plan and held command, not the template's
    return b
end

# The heat-load estimate covers one pass: it restarts when the spacecraft enters
# the atmosphere, as the QP's heat-load budget does.
SimulationLifecycle.requires_atmosphere_events(::AerobrakingMPCControlModel) = true
function SimulationLifecycle.atmosphere_transition!(
    m::AerobrakingMPCControlModel, u, p, t, spacecraft_index, inside)
    (inside && spacecraft_index == m.spacecraft_index) || return nothing
    m.state.estimated_heat_load_j_cm2 = 0.0
    m.state.last_command_time_s = Float64(t)
    return nothing
end
