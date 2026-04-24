include(joinpath(@__DIR__, "..", "..", "core", "interfaces", "reference_system.jl"))
using ComponentArrays
using Logging

const _MANEUVER_TRACE_LOCK = ReentrantLock()
const _MANEUVER_TRACE_LAST_WINDOW = Dict{Tuple{UInt64, Int64}, Tuple{Float64, Float64}}()
const _MANEUVER_TRACE_BURN_ACTIVE = Dict{Tuple{UInt64, Int64}, Bool}()
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

@inline function _maneuver_trace_key(controlModel::BaseThrusterModel, i::Int64)::Tuple{UInt64, Int64}
    return (UInt64(objectid(controlModel.start_burn_time)), i)
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

function _validated_burn_plan(
    controlModel::BaseThrusterModel,
    p,
    i::Int64,
    mass_kg::Float64,
    maneuver
)::Union{Nothing, PropulsiveBurnPlan}
    thrust_n, isp_s = _effective_thrust_isp(controlModel, p, i)
    delta_v_mps = Float64(maneuver.delta_v_mps)
    direction_rad = Float64(maneuver.direction_rad)

    if !(isfinite(mass_kg) && mass_kg > 0.0 &&
         isfinite(delta_v_mps) && delta_v_mps > 0.0 &&
         isfinite(direction_rad) &&
         isfinite(thrust_n) && thrust_n > 0.0 &&
         isfinite(isp_s) && isp_s > 0.0)
        return nothing
    end

    exhaust_velocity_mps = isp_s * _STANDARD_GRAVITY_MPS2
    mass_fraction = exp(-delta_v_mps / exhaust_velocity_mps)
    if !(isfinite(mass_fraction) && 0.0 < mass_fraction < 1.0)
        return nothing
    end

    propellant_required_kg = mass_kg * (1.0 - mass_fraction)
    commanded_impulse_n_s = propellant_required_kg * exhaust_velocity_mps
    burn_duration_s = commanded_impulse_n_s / thrust_n
    if !(isfinite(propellant_required_kg) && propellant_required_kg > 0.0 &&
         isfinite(commanded_impulse_n_s) && commanded_impulse_n_s > 0.0 &&
         isfinite(burn_duration_s) && burn_duration_s > 0.0)
        return nothing
    end

    available_propellant_kg = _available_propellant_kg(p, i, mass_kg)
    if available_propellant_kg !== nothing && propellant_required_kg > available_propellant_kg + 1e-9
        return nothing
    end

    return PropulsiveBurnPlan(
        valid=true,
        delta_v_mps=delta_v_mps,
        direction_rad=direction_rad,
        source_orbit=Int64(maneuver.source_orbit),
        thrust_n=thrust_n,
        isp_s=isp_s,
        commanded_impulse_n_s=commanded_impulse_n_s,
        propellant_required_kg=propellant_required_kg,
    )
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

    lock(_MANEUVER_TRACE_LOCK) do
        mkpath(dirname(path))
        new_file = !isfile(path)
        open(path, "a") do io
            if new_file
                println(io, "event,t_s,spacecraft_idx,orbit_counter,dv_cmd_mps,direction_rad,thrust_n,start_burn_s,stop_burn_s,burn_duration_s,alt_m,e,nu_rad,a_m")
            end
            println(io,
                string(event), ",",
                string(t), ",",
                string(i), ",",
                string(orbit_counter), ",",
                string(dv_cmd), ",",
                string(direction), ",",
                string(thrust_n), ",",
                string(start_burn_s), ",",
                string(stop_burn_s), ",",
                string(burn_duration_s), ",",
                string(alt_m), ",",
                string(e), ",",
                string(nu_rad), ",",
                string(a_m)
            )
        end
    end
    return nothing
end

@inline function _control_effector_exception_fallback(p, spacecraft_idx::Int, err, bt)
    if _control_effector_log_enabled(p)
        @warn "calcControlEffect! orbital-element conversion failed for spacecraft index $(spacecraft_idx); burn scheduling skipped." exception=(err, bt)
    end
    if _control_effector_strict_exceptions()
        throw(err)
    end
    return nothing
end

function calcControlMassFlowRate(controlModel::AbstractControlEffectorModel, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)::Float64
    return 0.0
end

function calcControlMassFlowRate(controlModel, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)::Float64
    return 0.0
end

"""
calcControlForceTorque(controlModel::BaseThrusterModel, x::AbstractVector, p::ODEParams, i::Int64, t::Float64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}

Calculate the control force and torque based on the thruster model and current state, called in the dynamics loop to get the current thruster force
- `controlModel`: The thruster model containing thrust magnitudes, directions, burn times, and specific impulses for each thruster
- `x`: The current state vector of the spacecraft
- `p`: The ODE parameters containing simulation configuration and other relevant data
- `i`: The index of the spacecraft for which to calculate the control force/torque
- `t`: The current time in the simulation

Returns a tuple containing the total control force and torque as 3D vectors
"""
function calcControlForceTorque(controlModel::BaseThrusterModel, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    # Calculate the control force and torque based on the thruster model and current state
    if i < 1 || i > length(controlModel.thrust)
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    start_time, stop_time = _effective_burn_window(controlModel, p, i)
    if t >= start_time && t <= stop_time
        thrust_mag, _ = _effective_thrust_isp(controlModel, p, i)
        vel_vec = SVector{3, Float64}(u.vel)
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

function calcControlMassFlowRate(controlModel::BaseThrusterModel, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)::Float64
    if i < 1 || i > length(controlModel.thrust)
        return 0.0
    end

    start_time, stop_time = _effective_burn_window(controlModel, p, i)
    if !(t >= start_time && t <= stop_time)
        return 0.0
    end

    _, Isp = _effective_thrust_isp(controlModel, p, i)
    if !isfinite(Isp) || Isp <= 0.0
        return 0.0
    end

    force, _ = calcControlForceTorque(controlModel, u, p, i, t)
    applied_thrust = norm(force)
    if !isfinite(applied_thrust) || applied_thrust <= 0.0
        return 0.0
    end

    return -applied_thrust / (Isp * _STANDARD_GRAVITY_MPS2)
end

"""
calcControlEffect!(controlModel::BaseThrusterModel, u::ComponentVector, p::ODEParams, t::Float64, i::Int64)
Calculate the control effect (force and torque) based on the control model and current state, and store it in the shared buffers for use in the dynamics calculations

Args
- `controlModel`: The thruster model containing thrust magnitudes, directions, burn times, and specific impulses for each thruster
- `u`: The current state vector of the spacecraft as a ComponentVector
- `p`: The ODE parameters containing simulation configuration and other relevant data
- `t`: The current time in the simulation
- `i`: The index of the spacecraft for which to calculate the control effect

Returns
- Updates the control force and torque in the shared buffers for the specified spacecraft index
"""
function calcControlEffect!(controlModel::BaseThrusterModel, u::ComponentVector, p::ODEParams, t::Float64, i::Int64)
    # Calculate the control effect (force and torque) based on the control model and current state, and store it in the shared buffers for use in the dynamics calculations
    if i < 1 || i > length(controlModel.start_burn_time)
        return
    end
    trace_key = _maneuver_trace_key(controlModel, i)

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
        was_in_burn_window = get(_MANEUVER_TRACE_BURN_ACTIVE, trace_key, false)
        if in_burn_window && !was_in_burn_window
            _MANEUVER_TRACE_BURN_ACTIVE[trace_key] = true
            _trace_maneuver_event!("burn_start", controlModel, p, i, t; start_burn_s=start_time, stop_burn_s=stop_time)
        elseif !in_burn_window && was_in_burn_window && t > stop_time + 1e-9
            _MANEUVER_TRACE_BURN_ACTIVE[trace_key] = false
            _trace_maneuver_event!("burn_end", controlModel, p, i, t; start_burn_s=start_time, stop_burn_s=stop_time)
        end

        # Lock the schedule only after ignition begins and until burn stop.
        if in_burn_window
            return
        end
        # Burn completed: clear the schedule so future campaign maneuvers can be planned.
        if t > stop_time + 1e-9
            _trace_maneuver_event!("schedule_clear", controlModel, p, i, t; start_burn_s=start_time, stop_burn_s=stop_time)
            controlModel.start_burn_time[i] = -1.0
            controlModel.stop_burn_time[i] = -1.0
            _clear_burn_plan!(p, i)
            _MANEUVER_TRACE_BURN_ACTIVE[trace_key] = false
            pop!(_MANEUVER_TRACE_LAST_WINDOW, trace_key, nothing)
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
        plan = _validated_burn_plan(controlModel, p, i, mass, maneuver)
        if plan === nothing
            return
        end
        burn_time = _constant_thrust_burn_duration_s(mass, plan.delta_v_mps, plan.thrust_n, plan.isp_s)
        if !isfinite(burn_time) || burn_time < 0.0
            return
        end

        # Estimate the time of apoapsis and set the burn to be symmetric about that time
        ψ = 2*atan(sqrt((1-e)/(1+e))*tan(ν_wrapped/2))

        M = ψ - e*sin(ψ)
        n = sqrt(p.args.environment_model.planet.μ/a^3)
        if !isfinite(M) || !isfinite(n) || n <= 0.0
            return
        end
        time_to_apoapsis = (π - M)/n
        apoapsis_time = t + time_to_apoapsis
        if !isfinite(apoapsis_time)
            return
        end
        # Calculate start/end time as symmetric about the apoapsis time
        start_burn_time = apoapsis_time - burn_time / 2
        stop_burn_time = apoapsis_time + burn_time / 2
        controlModel.start_burn_time[i] = start_burn_time
        controlModel.stop_burn_time[i] = stop_burn_time
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

        prev_window = get(_MANEUVER_TRACE_LAST_WINDOW, trace_key, (NaN, NaN))
        same_window = isfinite(prev_window[1]) && isfinite(prev_window[2]) &&
            abs(prev_window[1] - start_burn_time) <= 1e-9 &&
            abs(prev_window[2] - stop_burn_time) <= 1e-9
        if !same_window
            event = (isfinite(prev_window[1]) && isfinite(prev_window[2])) ? "schedule_update" : "schedule_set"
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
            _MANEUVER_TRACE_LAST_WINDOW[trace_key] = (start_burn_time, stop_burn_time)
        end
    end
end
