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

# The area commanded without a plan, also the nominal area of the reference pass:
# fully exposed panels for maximum depletion, `fallback_area_m2` for targeting.
function _mpc_desired_area(model::AerobrakingMPCControlModel)
    if model.config.mode isa MaxEnergyDepletionMode
        return model.config.bus_reference_area_m2 + model.config.controllable_area_m2
    end
    return model.fallback_area_m2
end

# The held plan's area at run time `t`; plan times are relative to `plan_epoch_s`.
function _mpc_interpolated_plan_area(model::AerobrakingMPCControlModel, t::Float64)
    times = model.state.plan_time_s
    areas = model.state.plan_area_m2
    if !(isfinite(model.state.plan_epoch_s) && length(times) == length(areas) && !isempty(times))
        return _mpc_desired_area(model)
    end
    return interpolate_mpc_plan(times, areas, t - model.state.plan_epoch_s)
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

# Density, temperature and planet-relative speed from the simulation's own
# atmosphere at the current state, for the runtime limits.
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

# Runtime panel heat rate (W/cm²) at angle of attack `alpha`: the configured
# SpaceAGORA thermal model (`:spaceagora_thermal`, default) or the kinetic-energy
# flux the QP predicts (`:kinetic_energy_flux`).
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

# The largest angle of attack in [min_alpha_rad, desired_alpha] whose runtime heat
# rate stays within `limit_w_cm2`, by bisection (the heat rate grows with alpha).
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

# The exposed area at the angle of attack `_mpc_alpha_for_heat_rate_limit` returns.
function _mpc_area_for_heat_rate_limit(model::AerobrakingMPCControlModel, p,
    speed::Float64, density::Float64, temperature::Float64, desired_alpha::Float64,
    limit_w_cm2::Float64)
    alpha_limit = _mpc_alpha_for_heat_rate_limit(
        model, p, speed, density, temperature, desired_alpha, limit_w_cm2)
    return commanded_area_from_alpha(model.config, alpha_limit;
        min_alpha_rad=model.min_alpha_rad, max_alpha_rad=model.max_alpha_rad)
end

# Runtime safety filter, applied at every control tick: clamp the plan's (or the
# default) area to the drag, heat-rate, remaining heat-load and slew limits
# evaluated with the simulation's atmosphere and thermal model, then hold it.
# Records the command, its heat rate and drag, and the limit that bound it.
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
        heat_area_limit = _mpc_area_for_heat_rate_limit(
            model, p, speed, density, temperature, desired_alpha, config.qdot_max_w_cm2)
        if heat_area_limit < area_limit
            area_limit = heat_area_limit
            active_limit = :heat_rate
        end
    end

    # Accumulate this pass's heat load with the heat rate held since the last tick.
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
            # Spend at most the remaining budget over the next control interval.
            heat_load_area_limit = _mpc_area_for_heat_rate_limit(
                model, p, speed, density, temperature, desired_alpha, remaining_load / dt)
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
    state.active_limit = active_limit
    state.solver_status = active_limit == :none ? :area_command : :constraint_limited_area_command
    return commanded_area
end

function _mpc_clear_plan!(state::AerobrakingMPCState)
    state.plan_epoch_s = NaN
    state.plan_time_s = Float64[]
    state.plan_area_m2 = Float64[]
    state.predicted_terminal_energy = NaN
    state.last_solution = nothing
    return state
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
    altitude_m = norm(pos) - params.Re
    descending = dot(pos, vel) < 0.0
    # With a trigger altitude, solve once per pass: the trigger re-arms above that
    # altitude and fires at the first tick below it while descending. Without one,
    # solve every `solve_interval_s`, except while climbing out inside the
    # reference interface, where a new reference would coast to the next pass.
    trigger = model.solve_trigger_altitude_m
    trigger !== nothing && altitude_m > trigger && (model.state.solve_armed = true)
    solve_now = trigger === nothing ?
        (descending || altitude_m > model.reference.h_cut_m) :
        (altitude_m <= trigger && model.state.solve_armed && descending)
    if !(solve_now && t - model.state.last_solve_time_s >= model.solve_interval_s)
        _mpc_area_command_from_constraints!(
            model, p, pos, vel, t, _mpc_interpolated_plan_area(model, t))
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
        # Runtime limits only, on the default area.
        _mpc_area_command_from_constraints!(model, p, pos, vel, t)
        _mpc_clear_plan!(model.state)
        return nothing
    end
    try
        # Prediction frame: inertial axes turned onto the planet-fixed axes at `t`.
        rotation = _mpc_planet_frame_rotation(p, t)
        ref = build_reference_drag_pass(
            params,
            rotation * pos,
            rotation * vel;
            config=model.config,
            reference=model.reference,
            nominal_area_m2=_mpc_desired_area(model),
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
        # Follow the solved plan, or hold the runtime-limited default area over the
        # pass when OSQP did not solve. Plan times start at the interface crossing.
        area = solution.ok ?
            _mpc_area_command_from_constraints!(model, p, pos, vel, t, first(solution.commanded_area_m2)) :
            _mpc_area_command_from_constraints!(model, p, pos, vel, t)
        state = model.state
        state.plan_epoch_s = t + first(problem.t)
        state.plan_time_s = collect(problem.t .- first(problem.t))
        state.plan_area_m2 = solution.ok ? collect(solution.commanded_area_m2) : fill(area, problem.N)
        state.predicted_terminal_energy = solution.terminal_energy
        state.solver_status = solution.solver_status
        state.last_solution = solution
        state.last_error = nothing
    catch err
        _mpc_area_command_from_constraints!(model, p, pos, vel, t)
        _mpc_clear_plan!(model.state)
        model.state.solver_status = :reference_or_qp_failed
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

# Result columns read the run's own controller: a run deep-copies its GNC models
# (unless `isolate_state=false`) and binds those declared on a spacecraft, so the
# model passed to the save-field functions only selects the columns.
function _mpc_running_effector(integrator, spacecraft_index::Int, ::Type{T}) where {T}
    for effector in integrator.p.args.control_model.control_effectors
        effector isa T && spacecraft_index in GNCRuntime.spacecraft_scope(effector).indices &&
            return effector
    end
    return nothing
end

# The controller whose command is applied: the model itself, or a campaign's active phase.
_mpc_active_control(model::AerobrakingMPCControlModel) = model

# One per-spacecraft column `name` of `value(effector)` for the running effector of
# type `T`, NaN for spacecraft it does not act on.
function _mpc_save_field(name::Symbol, ::Type{T}, value) where {T}
    SaveField = getproperty(parentmodule(@__MODULE__), :SaveField)
    return SaveField(
        name,
        (u, t, integrator) -> [
            begin
                effector = _mpc_running_effector(integrator, i, T)
                effector === nothing ? NaN : Float64(value(effector))
            end
            for i in eachindex(integrator.p.args.dynamics_model.spacecraft)
        ];
        per_satellite=true,
        column_prefix=String(name),
    )
end

# The command and runtime-limit columns shared by the single and campaign controllers.
function _mpc_control_columns(::Type{T}) where {T}
    active(effector) = _mpc_active_control(effector)
    return (
        _mpc_save_field(:mpc_commanded_area_m2, T, e -> _mpc_save_control(active(e))[1]),
        _mpc_save_field(:mpc_alpha_rad, T, e -> _mpc_save_control(active(e))[2]),
        _mpc_save_field(:mpc_active_limit_code, T, e -> _mpc_limit_code(active(e).state.active_limit)),
        _mpc_save_field(:mpc_heat_rate_w_cm2, T, e -> active(e).state.last_heat_rate_w_cm2),
        _mpc_save_field(:mpc_drag_n, T, e -> active(e).state.last_drag_n),
        _mpc_save_field(:mpc_heat_load_j_cm2, T, e -> active(e).state.estimated_heat_load_j_cm2),
    )
end

"""
    mpc_control_save_fields(model) -> Tuple{SaveField...}

Per-spacecraft result columns for the [`AerobrakingMPCControlModel`](@ref)s of a
run: the commanded area and panel angle, the limiting constraint code (0 none,
1 drag, 2 heat rate, 3 heat load, 4 slew), and the runtime heat rate (W/cm²),
drag (N) and heat load of the pass (J/cm²). The values come from the run's own
controllers, so they hold with `isolate_state` and for a controller declared on
a spacecraft.
"""
mpc_control_save_fields(::AerobrakingMPCControlModel) = _mpc_control_columns(AerobrakingMPCControlModel)

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
