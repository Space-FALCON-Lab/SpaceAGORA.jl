using ..FrameTransforms
using ComponentArrays
using Logging

const _STANDARD_GRAVITY_MPS2 = 9.80665

@inline function _control_effector_log_enabled(p)::Bool
    if get(ENV, "SPACEAGORA_DEBUG_CONTROL", "0") == "1"
        return true
    end
    return try
        hasproperty(p, :shared_buffers) &&
        hasproperty(p.shared_buffers, :debug_control) &&
        Bool(p.shared_buffers.debug_control[])
    catch
        false
    end
end

@inline _control_effector_strict_exceptions() = get(ENV, "SPACEAGORA_STRICT_CONTROL_EXCEPTIONS", "0") == "1"

@inline function _trace_bool_enabled(raw::AbstractString)::Bool
    token = lowercase(strip(raw))
    return token in ("1", "true", "yes", "on")
end

@inline function _maneuver_trace_enabled()::Bool
    return _trace_bool_enabled(get(ENV, "SPACEAGORA_TRACE_MANEUVERS", "0")) ||
           !isempty(strip(get(ENV, "SPACEAGORA_MANEUVER_TRACE_CSV", "")))
end

@inline function _maneuver_trace_path()::String
    raw = strip(get(ENV, "SPACEAGORA_MANEUVER_TRACE_CSV", ""))
    return isempty(raw) ? "/tmp/spaceagora_maneuver_trace.csv" : raw
end

@inline function _safe_orbit_counter(p, i::Int64)::Int64
    return try
        Int64(p.orbit_counter[i])
    catch
        Int64(-1)
    end
end

@inline function _guidance_maneuver_command(p, i::Int64)
    if !hasproperty(p, :shared_buffers) || !hasproperty(p.shared_buffers, :maneuver_commands)
        return nothing
    end

    commands = p.shared_buffers.maneuver_commands
    if i < 1 || i > length(commands)
        return nothing
    end

    command = commands[i]
    return command.valid ? command : nothing
end

@inline function _burn_plan_buffer(p)
    if !hasproperty(p, :shared_buffers) || !hasproperty(p.shared_buffers, :maneuver_burn_plans)
        return nothing
    end
    return p.shared_buffers.maneuver_burn_plans
end

@inline function _active_burn_plan(p, i::Int64)
    plans = _burn_plan_buffer(p)
    if plans === nothing || i < 1 || i > length(plans)
        return nothing
    end
    plan = plans[i]
    return plan.valid ? plan : nothing
end

@inline function _set_burn_plan!(p, i::Int64, plan::PropulsiveBurnPlan)
    plans = _burn_plan_buffer(p)
    if plans === nothing || i < 1 || i > length(plans)
        return nothing
    end
    plans[i] = plan
    return nothing
end

@inline function _clear_burn_plan!(p, i::Int64)
    plans = _burn_plan_buffer(p)
    if plans === nothing || i < 1 || i > length(plans)
        return nothing
    end
    plans[i] = PropulsiveBurnPlan()
    return nothing
end

@inline function _commanded_maneuver(controlModel::BaseThrusterModel, p, i::Int64)
    command = _guidance_maneuver_command(p, i)
    if command !== nothing
        delta_v_mps = Float64(command.delta_v_mps)
        direction_rad = Float64(command.direction_rad)
        if isfinite(delta_v_mps) && delta_v_mps < 0.0
            delta_v_mps = abs(delta_v_mps)
            direction_rad = π
        end
        controlModel.Δv[i] = delta_v_mps
        controlModel.direction[i] = direction_rad
        _note_schedule!(controlModel, i)
        return (
            delta_v_mps=delta_v_mps,
            direction_rad=direction_rad,
            source_orbit=command.source_orbit,
        )
    end

    if i < 1 || i > length(controlModel.Δv) || i > length(controlModel.direction)
        return nothing
    end

    delta_v_cmd = Float64(controlModel.Δv[i])
    direction_rad = Float64(controlModel.direction[i])
    if isfinite(delta_v_cmd)
        if delta_v_cmd > 0.0
            direction_rad = 0.0
        elseif delta_v_cmd < 0.0
            delta_v_cmd = abs(delta_v_cmd)
            direction_rad = π
        end
        controlModel.Δv[i] = delta_v_cmd
        controlModel.direction[i] = direction_rad
        _note_schedule!(controlModel, i)
    end

    return (
        delta_v_mps=delta_v_cmd,
        direction_rad=direction_rad,
        source_orbit=Int64(-1),
    )
end

@inline function _model_burn_window(controlModel::BaseThrusterModel, i::Int64)
    if i < 1 || i > length(controlModel.start_burn_time)
        return NaN, NaN
    end
    return controlModel.start_burn_time[i], controlModel.stop_burn_time[i]
end

@inline function _effective_burn_window(controlModel::BaseThrusterModel, p, i::Int64)
    start_time, stop_time = _model_burn_window(controlModel, i)
    if isfinite(start_time) && isfinite(stop_time) && stop_time > start_time
        return start_time, stop_time
    end
    plan = _active_burn_plan(p, i)
    if plan !== nothing
        return plan.start_burn_s, plan.stop_burn_s
    end
    return start_time, stop_time
end

@inline function _effective_direction_rad(controlModel::BaseThrusterModel, p, i::Int64)::Float64
    plan = _active_burn_plan(p, i)
    if plan !== nothing
        return plan.direction_rad
    end
    if i < 1 || i > length(controlModel.direction)
        return NaN
    end
    return Float64(controlModel.direction[i])
end

@inline function _effective_thrust_isp(controlModel::BaseThrusterModel, p, i::Int64)
    plan = _active_burn_plan(p, i)
    if plan !== nothing
        return plan.thrust_n, plan.isp_s
    end
    if i < 1 || i > length(controlModel.thrust) || i > length(controlModel.Isp)
        return NaN, NaN
    end
    return Float64(controlModel.thrust[i]), Float64(controlModel.Isp[i])
end

@inline function _available_propellant_kg(p, i::Int64, current_mass_kg::Float64)::Union{Nothing, Float64}
    if !hasproperty(p, :args) || !hasproperty(p.args, :dynamics_model)
        return nothing
    end
    spacecraft = p.args.dynamics_model.spacecraft
    if i < 1 || i > length(spacecraft)
        return nothing
    end
    sc = spacecraft[i]
    if !isfinite(sc.prop_mass) || sc.prop_mass <= 0.0
        return nothing
    end
    return max(0.0, current_mass_kg - sc.dry_mass)
end

# The trusted size check of a constant-thrust burn of `delta_v_mps` from
# `mass_kg`: the propellant it takes, its impulse and its duration, with the
# reason `:accepted`, or the reason it cannot be planned (the other values are
# then NaN, except for `:insufficient_propellant`, which keeps them).
# `available_propellant_kg === nothing` means the vehicle declares no
# propellant limit. The burn scheduler checks every burn it plans with it, and
# the burn-plan family (`BurnPlanFamily`) validates proposed plans with it.
function _burn_plan_sizing(
    mass_kg::Float64,
    delta_v_mps::Float64,
    direction_rad::Float64,
    thrust_n::Float64,
    isp_s::Float64,
    available_propellant_kg::Union{Nothing, Float64}
)
    refused(reason::Symbol) = (reason=reason, propellant_required_kg=NaN, commanded_impulse_n_s=NaN, burn_duration_s=NaN)
    (isfinite(mass_kg) && mass_kg > 0.0) || return refused(:invalid_mass)
    (isfinite(delta_v_mps) && delta_v_mps > 0.0) || return refused(:invalid_delta_v)
    isfinite(direction_rad) || return refused(:invalid_direction)
    (isfinite(thrust_n) && thrust_n > 0.0) || return refused(:invalid_thrust)
    (isfinite(isp_s) && isp_s > 0.0) || return refused(:invalid_isp)

    exhaust_velocity_mps = isp_s * _STANDARD_GRAVITY_MPS2
    mass_fraction = exp(-delta_v_mps / exhaust_velocity_mps)
    (isfinite(mass_fraction) && 0.0 < mass_fraction < 1.0) || return refused(:invalid_mass_fraction)

    propellant_required_kg = mass_kg * (1.0 - mass_fraction)
    commanded_impulse_n_s = propellant_required_kg * exhaust_velocity_mps
    burn_duration_s = commanded_impulse_n_s / thrust_n
    if !(isfinite(propellant_required_kg) && propellant_required_kg > 0.0 &&
         isfinite(commanded_impulse_n_s) && commanded_impulse_n_s > 0.0 &&
         isfinite(burn_duration_s) && burn_duration_s > 0.0)
        return refused(:invalid_burn_size)
    end

    affordable = available_propellant_kg === nothing || !(propellant_required_kg > available_propellant_kg + 1e-9)
    return (reason=affordable ? :accepted : :insufficient_propellant, propellant_required_kg=propellant_required_kg,
            commanded_impulse_n_s=commanded_impulse_n_s, burn_duration_s=burn_duration_s)
end

# The scheduler's plan for `maneuver` (not yet timed), or `nothing`, with the
# reason from `_burn_plan_sizing`.
function _sized_burn_plan(controlModel::BaseThrusterModel, p, i::Int64, mass_kg::Float64, maneuver)
    thrust_n, isp_s = _effective_thrust_isp(controlModel, p, i)
    delta_v_mps = Float64(maneuver.delta_v_mps)
    direction_rad = Float64(maneuver.direction_rad)
    sizing = _burn_plan_sizing(mass_kg, delta_v_mps, direction_rad, thrust_n, isp_s,
                               _available_propellant_kg(p, i, mass_kg))
    sizing.reason === :accepted || return (nothing, sizing.reason)
    plan = PropulsiveBurnPlan(
        valid=true,
        delta_v_mps=delta_v_mps,
        direction_rad=direction_rad,
        source_orbit=Int64(maneuver.source_orbit),
        thrust_n=thrust_n,
        isp_s=isp_s,
        commanded_impulse_n_s=sizing.commanded_impulse_n_s,
        propellant_required_kg=sizing.propellant_required_kg,
    )
    return plan, sizing.reason
end

function _validated_burn_plan(
    controlModel::BaseThrusterModel,
    p,
    i::Int64,
    mass_kg::Float64,
    maneuver
)::Union{Nothing, PropulsiveBurnPlan}
    return first(_sized_burn_plan(controlModel, p, i, mass_kg, maneuver))
end

# Keep a refused plan's reason, one record per refusal, for `burn_plan_rejections`.
function _record_burn_rejection!(controlModel::BaseThrusterModel, i::Int64, t::Float64, reason::Symbol, maneuver)
    1 <= i <= length(controlModel.plan_rejections) || return nothing
    push!(controlModel.plan_rejections[i], (time_s=t, rejection=reason,
        delta_v_mps=Float64(maneuver.delta_v_mps), direction_rad=Float64(maneuver.direction_rad),
        source_orbit=Int64(maneuver.source_orbit)))
    return nothing
end

"""
    burn_plan_rejections(model::BaseThrusterModel, spacecraft_index) -> Vector{NamedTuple}

The burns the scheduler refused to plan for the spacecraft at `spacecraft_index`
(its position in the run, from 1) in the run, oldest first, as copies. Each
record has `time_s` (elapsed s), the `rejection` reason, and the commanded
`delta_v_mps`, `direction_rad` and `source_orbit`. The reasons are those of the
burn-plan family's validation (`validate_burn_plan`): `:invalid_mass`,
`:invalid_delta_v`, `:invalid_direction`, `:invalid_thrust`, `:invalid_isp`,
`:invalid_mass_fraction`, `:invalid_burn_size` and `:insufficient_propellant`.
A command of exactly zero Δv asks for no burn and leaves no record. A refused
command leaves any burn already planned in place. The records start empty at
each run's start; after `run_simulation` with the default `isolate_state=true`,
read them from the configuration that ran (`solution.prob.p.args`).
"""
function burn_plan_rejections(model::BaseThrusterModel, spacecraft_index::Integer)
    1 <= spacecraft_index <= length(model.plan_rejections) || throw(ArgumentError(
        "BaseThrusterModel has no spacecraft $(spacecraft_index); it serves $(length(model.plan_rejections))."))
    return copy(model.plan_rejections[spacecraft_index])
end

# One spacecraft's schedule entry, as the configured and written records hold it.
@inline _schedule_entry(model::BaseThrusterModel, i::Int) =
    (model.Δv[i], model.direction[i], model.start_burn_time[i], model.stop_burn_time[i])

# Records the scheduler's write of spacecraft `i`'s schedule.
@inline function _note_schedule!(model::BaseThrusterModel, i::Int)
    1 <= i <= length(model.written_schedule) && (model.written_schedule[i] = _schedule_entry(model, i))
    return nothing
end

# Each run starts from the configured schedule (G-8): an entry the scheduler wrote
# in the previous run and nobody changed since returns to the value that run
# started from; any other entry is the configuration as it stands. The burn state,
# the mounted levels and the refusal records start empty.
function SimulationLifecycle.initialize_control!(model::BaseThrusterModel, u, p, t)
    n = length(model.thrust)
    first_run = isempty(model.configured_schedule)
    first_run && resize!(model.configured_schedule, n)
    for i in 1:n
        current = _schedule_entry(model, i)
        written = model.written_schedule[i]
        start = first_run || written === nothing ? current :
            ntuple(k -> isequal(current[k], written[k]) ? model.configured_schedule[i][k] : current[k], 4)
        model.Δv[i], model.direction[i], model.start_burn_time[i], model.stop_burn_time[i] = start
        model.configured_schedule[i] = start
        model.written_schedule[i] = nothing
    end
    fill!(model.burn_windows, (NaN, NaN, 0.0))
    fill!(model.burning, false)
    foreach(empty!, model.mounted_levels)
    foreach(empty!, model.plan_rejections)
    return nothing
end

# The window of a burn of `burn_time_s` centred on the next apoapsis after `t`,
# from the osculating elements `a`, `e` and the wrapped true anomaly
# `ν_wrapped` about a body of gravitational parameter `μ`, or `nothing` when the
# timing is not finite. The scheduler and `ApoapsisBurnPlanner` share it.
function _apoapsis_burn_window(t::Float64, a, e, ν_wrapped::Float64, μ, burn_time_s::Float64)
    ψ = 2*atan(sqrt((1-e)/(1+e))*tan(ν_wrapped/2))
    M = ψ - e*sin(ψ)
    n = sqrt(μ/a^3)
    if !isfinite(M) || !isfinite(n) || n <= 0.0
        return nothing
    end
    time_to_apoapsis = (π - M)/n
    apoapsis_time = t + time_to_apoapsis
    if !isfinite(apoapsis_time)
        return nothing
    end
    return apoapsis_time - burn_time_s / 2, apoapsis_time + burn_time_s / 2
end

@inline function _constant_thrust_burn_duration_s(
    mass_kg::Float64,
    delta_v_mps::Float64,
    thrust_n::Float64,
    isp_s::Float64
)::Float64
    if !(isfinite(mass_kg) && mass_kg > 0.0)
        return NaN
    end
    if !(isfinite(delta_v_mps) && delta_v_mps >= 0.0)
        return NaN
    end
    if delta_v_mps == 0.0
        return 0.0
    end
    if !(isfinite(thrust_n) && thrust_n > 0.0)
        return NaN
    end
    if !(isfinite(isp_s) && isp_s > 0.0)
        return NaN
    end

    g0 = 9.80665
    exhaust_velocity = isp_s * g0
    mass_fraction_used = 1.0 - exp(-delta_v_mps / exhaust_velocity)
    return (mass_kg * exhaust_velocity / thrust_n) * mass_fraction_used
end

function _trace_maneuver_event!(
    event::String,
    controlModel::BaseThrusterModel,
    p,
    i::Int64,
    t::Float64;
    start_burn_s::Float64=NaN,
    stop_burn_s::Float64=NaN,
    alt_m::Float64=NaN,
    e::Float64=NaN,
    nu_rad::Float64=NaN,
    a_m::Float64=NaN
)
    _maneuver_trace_enabled() || return nothing
    path = _maneuver_trace_path()
    orbit_counter = _safe_orbit_counter(p, i)
    dv_cmd = i <= length(controlModel.Δv) ? Float64(controlModel.Δv[i]) : NaN
    direction = i <= length(controlModel.direction) ? Float64(controlModel.direction[i]) : NaN
    thrust_n = i <= length(controlModel.thrust) ? Float64(controlModel.thrust[i]) : NaN
    burn_duration_s = (isfinite(start_burn_s) && isfinite(stop_burn_s)) ? (stop_burn_s - start_burn_s) : NaN

    # Each row is one write to a file opened with the system's append flag, so
    # runs tracing to the same file at once add whole rows and never overwrite
    # each other's (`open(path, "a")` seeks to the end only once, when opened).
    # Runs that create the file at the same moment may each write the header.
    row = join((event, t, i, orbit_counter, dv_cmd, direction, thrust_n, start_burn_s, stop_burn_s,
                burn_duration_s, alt_m, e, nu_rad, a_m), ",") * "\n"
    mkpath(dirname(path))
    if !isfile(path)
        row = "event,t_s,spacecraft_idx,orbit_counter,dv_cmd_mps,direction_rad,thrust_n,start_burn_s,stop_burn_s,burn_duration_s,alt_m,e,nu_rad,a_m\n" * row
    end
    file = Base.Filesystem.open(path, Base.Filesystem.JL_O_WRONLY | Base.Filesystem.JL_O_CREAT | Base.Filesystem.JL_O_APPEND, 0o644)
    try
        write(file, row)
    finally
        close(file)
    end
    return nothing
end

@inline function _control_effector_exception_fallback(p, spacecraft_idx::Int, err, bt)
    if _control_effector_log_enabled(p)
        @warn "update_control! orbital-element conversion failed for spacecraft index $(spacecraft_idx); burn scheduling skipped." exception=(err, bt)
    end
    if _control_effector_strict_exceptions()
        throw(err)
    end
    return nothing
end

"""
    calcControlMassFlowRate(model, u, p, spacecraft_index, t) -> Float64

Deprecated, use [`control_propellant_rate`](@ref)`(model, x::StateSample, p, t)`;
this name will be removed in SpaceAGORA 0.5.0. **The sign changes:** this hook
returns the spacecraft's mass derivative in kg/s, negative while propellant is
used, and `control_propellant_rate` returns the propellant used, positive.

Until then the dynamics call an effector's own method of this hook when its
type has no `control_propellant_rate`, and `run_simulation` warns once per type
before integration starts. The arguments are those of the deprecated
`calcControlForceTorque`. An effector with neither method uses no propellant.

A direct call on an effector that implements `control_propellant_rate` returns
that hook's result negated, and warns once per type that it is deprecated
(under `julia --depwarn=error` it throws instead). On any other effector it
returns `0.0`.
"""
function calcControlMassFlowRate end

"""
    calcReactionWheelTorque(model, u, p, spacecraft_index, t) -> Union{Nothing, AbstractVector}

Deprecated, use [`control_wheel_torque`](@ref)`(model, x::StateSample, p, t)`;
this name will be removed in SpaceAGORA 0.5.0. The torque, its frame and how
the dynamics apply it are unchanged.

Until then the dynamics call an effector's own method of this hook when its
type has no `control_wheel_torque`, and `run_simulation` warns once per type
before integration starts. The arguments are those of the deprecated
`calcControlForceTorque`. An effector with neither method drives no wheel.

A direct call on an effector that implements `control_wheel_torque` returns
that hook's result, and warns once per type that it is deprecated (under
`julia --depwarn=error` it throws instead). On any other effector it returns
`nothing`.
"""
function calcReactionWheelTorque end

"""
    control_thruster_levels(controlModel::BaseThrusterModel, i)

The commanded throttle of the maneuver engine: 1 while spacecraft `i`'s burn
is on, inside its `[start, stop)` window, and 0 otherwise. It follows the
burn state that `refresh_burn!` keeps, which switches exactly at the window's
ends, as the thrust does. A maneuver burn drives one engine, so the vector
covers the first thruster only and the vehicle's remaining thrusters stay
idle.
"""
function control_thruster_levels(controlModel::BaseThrusterModel, i::Int)
    (1 <= i <= length(controlModel.burning)) || return nothing
    controlModel.actuation === :mounted && return nothing   # the mounted thrusters' own levels report it
    return [controlModel.burning[i] ? 1.0 : 0.0]
end

# Internal: the firing levels an effector reports for spacecraft `i` at time `t`,
# for the saved `thruster_level` field. A sample is saved after the callbacks of
# the step it falls in, so held levels can belong to the end of that step; an
# effector that knows its levels over time answers for the sample's own time.
# The default is the held levels.
_thruster_levels_at(effector, i::Int, t) = control_thruster_levels(effector, i)

# The burn scheduler's level at `t` from the recorded window: 1 inside
# `[start, stop)` and 0 outside, as the burn itself switches.
function _thruster_levels_at(controlModel::BaseThrusterModel, i::Int, t)
    (1 <= i <= length(controlModel.burn_windows)) || return nothing
    controlModel.actuation === :mounted && return nothing   # the mounted thrusters' own levels report it
    start_s, stop_s, _ = controlModel.burn_windows[i]
    return [isfinite(start_s) && isfinite(stop_s) && start_s <= t < stop_s ? 1.0 : 0.0]
end

"""
    control_wrench(model::BaseThrusterModel, x::StateSample, p, t)

The maneuver engine's force and torque on the spacecraft at
`x.spacecraft_index` (see the generic `control_wrench` for the arguments,
frames and units). During the effective burn window, thrust is parallel or
antiparallel to the inertial velocity according to the configured direction;
torque is zero. Outside the burn window or for invalid thrust inputs, both
vectors are zero. At the window's ends, where the run stops, the thrust
follows the burn's state, which the run switches between integration steps
(see `refresh_burn!`): a step that ends at the start coasts, a step that ends
at the stop burns.
"""
function control_wrench(controlModel::BaseThrusterModel, x::StateSample, p::ODEParams, t::Float64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    # Calculate the control force and torque based on the thruster model and current state
    i = x.spacecraft_index
    if i < 1 || i > length(controlModel.thrust)
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    start_time, stop_time = _effective_burn_window(controlModel, p, i)
    if _burn_active(controlModel, i, t, start_time, stop_time)
        thrust_mag, _ = _effective_thrust_isp(controlModel, p, i)
        vel_vec = x.vel_ii
        vel_mag = norm(vel_vec)
        direction_rad = _effective_direction_rad(controlModel, p, i)
        if vel_mag == 0.0 || !isfinite(vel_mag) || !isfinite(thrust_mag) || thrust_mag <= 0.0 || !isfinite(direction_rad)
            return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
        end
        thrust_dir = normalize(vel_vec) # Prograde direction
        thrust_dir *= cos(direction_rad) >= 0.0 ? 1.0 : -1.0 # 0 -> prograde, π -> retrograde
        force = thrust_mag * thrust_dir
        # For simplicity, assume torque is zero for now (can be updated later to include offset thrusters, gimbaled thrusters, etc.)
        torque = SVector{3, Float64}(0.0, 0.0, 0.0)
    else
        force = SVector{3, Float64}(0.0, 0.0, 0.0)
        torque = SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    return force, torque
end

# The propellant the maneuver engine uses, in kg/s: `norm(force) / (Isp * g0)`
# during an active valid burn, 0.0 otherwise.
function control_propellant_rate(controlModel::BaseThrusterModel, x::StateSample, p::ODEParams, t::Float64)::Float64
    i = x.spacecraft_index
    if i < 1 || i > length(controlModel.thrust)
        return 0.0
    end

    start_time, stop_time = _effective_burn_window(controlModel, p, i)
    if !_burn_active(controlModel, i, t, start_time, stop_time)
        return 0.0
    end

    _, Isp = _effective_thrust_isp(controlModel, p, i)
    if !isfinite(Isp) || Isp <= 0.0
        return 0.0
    end

    force, _ = control_wrench(controlModel, x, p, t)
    applied_thrust = norm(force)
    if !isfinite(applied_thrust) || applied_thrust <= 0.0
        return 0.0
    end

    return propellant_mass_rate(applied_thrust, Isp)
end

"""
    update_control!(model::BaseThrusterModel, u::ComponentVector, p, t, i)

For `BaseThrusterModel`, `u` is the full simulation state (`u.sc[i]` selects
the spacecraft), `p` holds the simulation parameters, and `t` is elapsed time
in seconds. This method schedules eligible maneuvers around apoapsis, updates
the model's burn times and the stored burn plan, preserves an active burn's
schedule, and clears completed burns. It does not directly compute or store
force and torque; [`control_wrench`](@ref) evaluates the scheduled
thrust during dynamics evaluation.
"""
function update_control!(controlModel::BaseThrusterModel, u::ComponentVector, p::ODEParams, t::Float64, i::Int64)
    # Calculate the control effect (force and torque) based on the control model and current state, and store it in the shared buffers for use in the dynamics calculations
    if i < 1 || i > length(controlModel.start_burn_time)
        return
    end
    # The window this spacecraft's schedule was last recorded with (by
    # `refresh_burn!`, which runs after every scheduling call), to tell a new
    # schedule from an update in the trace.
    recorded_start, recorded_stop, _ = controlModel.burn_windows[i]

    active_plan = _active_burn_plan(p, i)
    model_start_time, model_stop_time = _model_burn_window(controlModel, i)
    start_time, stop_time = if isfinite(model_start_time) && isfinite(model_stop_time) && model_stop_time > model_start_time
        model_start_time, model_stop_time
    elseif active_plan === nothing
        model_start_time, model_stop_time
    else
        active_plan.start_burn_s, active_plan.stop_burn_s
    end
    if isfinite(start_time) && isfinite(stop_time) && stop_time > start_time
        in_burn_window = t >= start_time - 1e-9 && t <= stop_time + 1e-9
        # Lock the schedule only after ignition begins and until burn stop.
        if in_burn_window
            return
        end
        # Burn completed: clear the schedule so future campaign maneuvers can be planned.
        if t > stop_time + 1e-9
            _trace_maneuver_event!("schedule_clear", controlModel, p, i, t; start_burn_s=start_time, stop_burn_s=stop_time)
            controlModel.start_burn_time[i] = -1.0
            controlModel.stop_burn_time[i] = -1.0
            _note_schedule!(controlModel, i)
            _clear_burn_plan!(p, i)
            recorded_start, recorded_stop = NaN, NaN
        end
    end
    pos = SVector{3, Float64}(u.sc[i].pos)
    vel = SVector{3, Float64}(u.sc[i].vel)
    mass = u.sc[i].mass
    if !isfinite(mass) || mass <= 0.0
        return
    end
    maneuver = _commanded_maneuver(controlModel, p, i)
    maneuver === nothing && return

    # Calculate the current orbital elements from the state vector
    oe = try
        planet = p.args.environment_model.planet
        rvtoorbitalelement(pos, vel, planet)
    catch err
        _control_effector_exception_fallback(p, i, err, catch_backtrace())
        return
    end
    a, e, _, _, _, ν = oe
    if !isfinite(a) || !isfinite(e) || !isfinite(ν) || a <= 0.0 || e < 0.0 || e >= 1.0
        return
    end
    
    # Don't start maneuver until spacecraft has exited the atmosphere and is pre-apoapsis.
    # For near-circular orbits, ν from OE conversion can be arbitrary; allow scheduling.
    alt = norm(pos) - p.args.environment_model.planet.Rp_e
    circular_e_tol = 1e-8
    ν_wrapped = _wrap_2pi(Float64(ν))
    pre_apoapsis = e <= circular_e_tol || ν_wrapped < π - 1e-12
    if alt >= p.args.environment_model.EI * 1000 - 1e-6 && pre_apoapsis
        plan, reason = _sized_burn_plan(controlModel, p, i, mass, maneuver)
        if plan === nothing
            # A zero Δv asks for no burn; any other refusal keeps its reason.
            maneuver.delta_v_mps == 0.0 || _record_burn_rejection!(controlModel, i, t, reason, maneuver)
            return
        end
        burn_time = _constant_thrust_burn_duration_s(mass, plan.delta_v_mps, plan.thrust_n, plan.isp_s)
        if !isfinite(burn_time) || burn_time < 0.0
            return
        end

        # Set the burn symmetric about the estimated time of the next apoapsis.
        window = _apoapsis_burn_window(t, a, e, ν_wrapped, p.args.environment_model.planet.μ, burn_time)
        window === nothing && return
        start_burn_time, stop_burn_time = window
        controlModel.start_burn_time[i] = start_burn_time
        controlModel.stop_burn_time[i] = stop_burn_time
        _note_schedule!(controlModel, i)
        new_plan = PropulsiveBurnPlan(
            valid=true,
            delta_v_mps=plan.delta_v_mps,
            direction_rad=plan.direction_rad,
            source_orbit=plan.source_orbit,
            start_burn_s=start_burn_time,
            stop_burn_s=stop_burn_time,
            thrust_n=plan.thrust_n,
            isp_s=plan.isp_s,
            commanded_impulse_n_s=plan.commanded_impulse_n_s,
            propellant_required_kg=plan.propellant_required_kg
        )
        _set_burn_plan!(p, i, new_plan)

        same_window = isfinite(recorded_start) && isfinite(recorded_stop) &&
            abs(recorded_start - start_burn_time) <= 1e-9 &&
            abs(recorded_stop - stop_burn_time) <= 1e-9
        if !same_window
            event = (isfinite(recorded_start) && isfinite(recorded_stop)) ? "schedule_update" : "schedule_set"
            _trace_maneuver_event!(
                event,
                controlModel,
                p,
                i,
                t;
                start_burn_s=start_burn_time,
                stop_burn_s=stop_burn_time,
                alt_m=alt,
                e=e,
                nu_rad=ν_wrapped,
                a_m=a
            )
        end
    end
end

# ---- Mounted actuation --------------------------------------------------------
#
# With `actuation = :mounted` the scheduler fires a mounted engine and points it
# with the mounted reaction wheels (see the `BaseThrusterModel` docstring).

actuation_mode(model::BaseThrusterModel) = model.actuation

function actuator_levels(model::BaseThrusterModel, i::Int)
    (model.actuation === :mounted && 1 <= i <= length(model.mounted_levels)) || return nothing
    levels = model.mounted_levels[i]
    isempty(levels) && return nothing
    return levels
end

# The window that `refresh_burn!` records for spacecraft `i`, as (start, stop,
# thrust), or (NaN, NaN, 0.0) when no valid burn is set.
@inline function _recordable_burn_window(model::BaseThrusterModel, p, i::Int)
    start_s, stop_s = _effective_burn_window(model, p, i)
    thrust_n, _ = _effective_thrust_isp(model, p, i)
    direction_rad = _effective_direction_rad(model, p, i)
    valid = isfinite(start_s) && isfinite(stop_s) && stop_s > start_s &&
            isfinite(thrust_n) && thrust_n > 0.0 && isfinite(direction_rad)
    return valid ? (Float64(start_s), Float64(stop_s), Float64(thrust_n)) : (NaN, NaN, 0.0)
end

# Whether spacecraft `i`'s effective burn differs from the one `refresh_burn!`
# last recorded. The scheduler records every burn it plans in the same pass, so
# this finds a plan installed in the burn-plan buffer between its passes, by an
# `AlgorithmGuidance` of the burn-plan family.
@inline _burn_window_changed(model::BaseThrusterModel, p, i::Int)::Bool =
    !isequal(_recordable_burn_window(model, p, i), model.burn_windows[i])

_burn_window_pending(model::BaseThrusterModel, p)::Bool =
    any(i -> _burn_window_changed(model, p, i), eachindex(model.burn_windows))

"""
    refresh_burn!(model::BaseThrusterModel, p, i, t) -> Nothing

Record spacecraft `i`'s effective burn window and thrust and set whether the
burn is on at time `t` (inside `[start, stop)`). Called when burns are
scheduled, at the window's ends, where the solver stops, and, in `:mounted`
actuation, at each control tick, where it also sizes the level vector to the
mounted actuators and sets the engine's level. Between those calls the burn's
state holds, so no integration step sees it change. The maneuver trace records
`burn_start` and `burn_end` here, when the burn switches.
"""
function refresh_burn!(model::BaseThrusterModel, p, i::Int, t::Real)
    (1 <= i <= length(model.burning)) || return nothing
    was_burning = model.burning[i]
    previous_start, previous_stop, _ = model.burn_windows[i]
    window = _recordable_burn_window(model, p, i)
    start_s, stop_s, thrust_n = window
    valid = isfinite(start_s)
    model.burn_windows[i] = window
    burning = valid && start_s <= t < stop_s
    model.burning[i] = burning
    model.actuation === :mounted && _refresh_mounted_engine!(model, p, i, burning ? Float64(thrust_n) : 0.0)
    if burning != was_burning
        # A burn that ends with its schedule cleared reports the window it ran.
        window_start, window_stop = burning || valid ? (Float64(start_s), Float64(stop_s)) : (previous_start, previous_stop)
        _trace_maneuver_event!(burning ? "burn_start" : "burn_end", model, p, Int64(i), Float64(t);
            start_burn_s=window_start, stop_burn_s=window_stop)
    end
    return nothing
end

# Inside the window the burn is on and outside it off; at an end, where the
# solver stops, it follows the flag `refresh_burn!` keeps, so the step that
# ends there and the step after it each see one thrust.
@inline _burn_active(model::BaseThrusterModel, i::Int, t, start_s, stop_s)::Bool =
    (start_s < t < stop_s) || ((t == start_s || t == stop_s) && model.burning[i])

# `:mounted`: the level vector sized to the mounted actuators, the engine at `engine_level`.
function _refresh_mounted_engine!(model::BaseThrusterModel, p, i::Int, engine_level::Float64)
    spacecraft = p.args.dynamics_model.spacecraft[i]
    n_thrusters = mounted_thruster_count(spacecraft)
    model.engine <= n_thrusters || throw(ArgumentError(
        "Mounted burns: engine $(model.engine) on spacecraft $(i), which carries $(n_thrusters) mounted thrusters. " *
        "Mount the engine with add_actuator!, or set `engine` to its mounting position."))
    levels = model.mounted_levels[i]
    n = mounted_actuator_count(spacecraft)
    if length(levels) != n
        resize!(levels, n)
        fill!(levels, 0.0)
    end
    levels[model.engine] = engine_level
    return nothing
end

"""
    burn_edge(model::BaseThrusterModel, t) -> Bool

Whether `t` is an end of a spacecraft's burn window: an instant the solver
stops at for the burn to switch.
"""
burn_edge(model::BaseThrusterModel, t::Real)::Bool = any(w -> t == w[1] || t == w[2], model.burn_windows)

# The engine (the `engine`-th mounted thruster) and its link.
function _mounted_engine(model::BaseThrusterModel, spacecraft)
    k = 0
    for link in spacecraft.links, thruster in link.thrusters
        k += 1
        k == model.engine && return link, thruster
    end
    throw(ArgumentError("Mounted burns: no mounted thruster at position $(model.engine)."))
end

# The body torque that turns the engine's thrust direction toward the burn
# direction (prograde or retrograde) and follows the velocity's turn along the
# orbit: J (ωn² θ ê - 2 ζ ωn (ω - ω_target)). `nothing` without a direction.
function _burn_pointing_torque(model::BaseThrusterModel, spacecraft, state, p, i::Int)
    pos = SVector{3, Float64}(state.pos[1], state.pos[2], state.pos[3])
    vel = SVector{3, Float64}(state.vel[1], state.vel[2], state.vel[3])
    speed = norm(vel)
    (isfinite(speed) && speed > 0.0) || return nothing
    direction_rad = _effective_direction_rad(model, p, i)
    isfinite(direction_rad) || return nothing
    v_hat = vel / speed
    target_ii = cos(direction_rad) >= 0.0 ? v_hat : -v_hat
    q = SVector{4, Float64}(state.q[1], state.q[2], state.q[3], state.q[4])
    ω = SVector{3, Float64}(state.ω[1], state.ω[2], state.ω[3])
    to_body = rot(q)                                       # inertial to body
    target_b = to_body * target_ii
    link, engine = _mounted_engine(model, spacecraft)
    d = thrust_direction(engine)                           # turned by its gimbals, if any
    nd = norm(d)
    (isfinite(nd) && nd > 0.0) || return nothing
    axis_b = actuator_context(link).link_to_body * (d / nd)
    e = cross(axis_b, target_b)
    s = norm(e)
    c = dot(axis_b, target_b)
    θ = atan(s, c)
    e_hat = if s > 1.0e-12
        e / s
    elseif c < 0.0                                          # pointing the opposite way: turn about any perpendicular axis
        perpendicular = cross(axis_b, abs(axis_b[1]) < 0.9 ? SVector(1.0, 0.0, 0.0) : SVector(0.0, 1.0, 0.0))
        perpendicular / norm(perpendicular)
    else
        SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    r = norm(pos)
    accel = -(p.args.environment_model.planet.μ / r^3) * pos
    v_hat_rate = (accel - dot(accel, v_hat) * v_hat) / speed
    ω_target_b = to_body * cross(v_hat, v_hat_rate)
    ωn = model.pointing_frequency_rad_s
    ζ = model.pointing_damping
    return SMatrix{3, 3, Float64, 9}(spacecraft.inertia_tensor) * (ωn^2 * θ * e_hat - 2.0 * ζ * ωn * (ω - ω_target_b))
end

"""
    mounted_burn_tick!(model::BaseThrusterModel, u, p, t, i) -> Nothing

The control tick of a scheduler in `:mounted` actuation for spacecraft `i`:
refresh the burn state, then, while a burn is planned or running, set the
mounted wheels' levels to the least-norm split of the pointing torque (within
their torque limits); the wheels idle otherwise.
"""
function mounted_burn_tick!(model::BaseThrusterModel, u, p, t::Float64, i::Int)
    model.actuation === :mounted || return nothing
    (1 <= i <= length(model.mounted_levels)) || return nothing
    refresh_burn!(model, p, i, t)
    spacecraft = p.args.dynamics_model.spacecraft[i]
    levels = model.mounted_levels[i]
    n_thrusters = mounted_thruster_count(spacecraft)
    n_wheels = mounted_wheel_count(spacecraft)
    for j in 1:n_wheels
        levels[n_thrusters + j] = 0.0
    end
    _, stop_s, _ = model.burn_windows[i]
    (n_wheels > 0 && isfinite(stop_s) && t <= stop_s) || return nothing
    state = u.sc[i]
    hasproperty(state, :q) || throw(ArgumentError(
        "Mounted burns point the engine with the wheels and need the attitude state: set orientation_sim=true (spacecraft $(i))."))
    torque = _burn_pointing_torque(model, spacecraft, state, p, i)
    torque === nothing && return nothing
    assembly = spacecraft.root.rw_assembly
    first_mounted = assembly.n_wheels - n_wheels
    axes = assembly.J_rw[:, (first_mounted + 1):assembly.n_wheels]
    upper = [wheel.max_torque_nm for wheel in assembly.wheels]
    wheel_levels = allocate_actuators(axes, torque, -upper, upper)
    for j in 1:n_wheels
        levels[n_thrusters + j] = wheel_levels[j]
    end
    return nothing
end
