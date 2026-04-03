include(joinpath(@__DIR__, "..", "..", "core", "interfaces", "reference_system.jl"))
using ComponentArrays
using Logging

const _MANEUVER_TRACE_LOCK = ReentrantLock()
const _MANEUVER_TRACE_LAST_WINDOW = Dict{Tuple{UInt64, Int64}, Tuple{Float64, Float64}}()
const _MANEUVER_TRACE_BURN_ACTIVE = Dict{Tuple{UInt64, Int64}, Bool}()

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

@inline function _commanded_maneuver!(controlModel::BaseThrusterModel, p, i::Int64)
    command = _guidance_maneuver_command(p, i)
    if command !== nothing
        controlModel.Δv[i] = command.delta_v_mps
        controlModel.direction[i] = command.direction_rad
        return (
            delta_v_mps=command.delta_v_mps,
            direction_rad=command.direction_rad,
            source_orbit=command.source_orbit,
        )
    end

    if i < 1 || i > length(controlModel.Δv) || i > length(controlModel.direction)
        return nothing
    end

    return (
        delta_v_mps=Float64(controlModel.Δv[i]),
        direction_rad=Float64(controlModel.direction[i]),
        source_orbit=Int64(-1),
    )
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
    if i < 1 || i > length(controlModel.start_burn_time)
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    start_time = controlModel.start_burn_time[i]
    stop_time = controlModel.stop_burn_time[i]
    if t >= start_time && t <= stop_time
        thrust_mag = controlModel.thrust[i]
        vel_vec = SVector{3, Float64}(u.vel)
        vel_mag = norm(vel_vec)
        if vel_mag == 0.0 || !isfinite(vel_mag)
            return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
        end
        thrust_dir = normalize(vel_vec) # Prograde direction
        thrust_dir *= cos(controlModel.direction[i]) >= 0.0 ? 1.0 : -1.0 # 0 -> prograde, π -> retrograde
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
    if i < 1 || i > length(controlModel.start_burn_time)
        return 0.0
    end

    start_time = controlModel.start_burn_time[i]
    stop_time = controlModel.stop_burn_time[i]
    if !(t >= start_time && t <= stop_time)
        return 0.0
    end

    Isp = controlModel.Isp[i]
    if !isfinite(Isp) || Isp <= 0.0
        return 0.0
    end

    force, _ = calcControlForceTorque(controlModel, u, p, i, t)
    applied_thrust = norm(force)
    if !isfinite(applied_thrust) || applied_thrust <= 0.0
        return 0.0
    end

    g0 = 9.80665 # Standard gravity [m/s^2]
    return -applied_thrust / (Isp * g0)
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

    # Default behavior: track/recenter the burn window until ignition,
    # then lock it once the burn has started.
    start_time = controlModel.start_burn_time[i]
    stop_time = controlModel.stop_burn_time[i]
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
    maneuver = _commanded_maneuver!(controlModel, p, i)
    maneuver === nothing && return

    # Calculate the current orbital elements from the state vector
    # State is in J2000; convert to PCI (body-equatorial) for rvtoorbitalelement
    oe = try
        planet = p.args.environment_model.planet
        pos_pci = planet.J2000_to_pci * pos
        vel_pci = planet.J2000_to_pci * vel
        rvtoorbitalelement(pos_pci, vel_pci, planet)
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
    pre_apoapsis = e <= circular_e_tol || ν < π - 1e-12
    if alt >= p.args.environment_model.EI * 1000 - 1e-6 && pre_apoapsis
        # Calculate the burn time required to achieve the desired Δv based on the current mass and thrust
        Δv = maneuver.delta_v_mps
        if !isfinite(Δv) || Δv <= 0.0
            return
        end
        thrust_mag = controlModel.thrust[i]
        # Calculate the burn duration for a constant-thrust impulsive approximation.
        if thrust_mag <= 0.0 || !isfinite(thrust_mag)
            return
        end
        burn_time = mass * Δv / thrust_mag
        if !isfinite(burn_time) || burn_time < 0.0
            return
        end
        # println("e: $e, a: $a, Δv: $Δv, burn_time: $burn_time")
        # Estimate the time of apoapsis and set the burn to be symmetric about that time
        ψ = 2*atan(sqrt((1-e)/(1+e))*tan(ν/2))

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
        start_burn_time = apoapsis_time - burn_time/2
        stop_burn_time = apoapsis_time + burn_time/2
        
        # Update the start/end time fields in the control model for the current spacecraft
        controlModel.start_burn_time[i] = start_burn_time
        controlModel.stop_burn_time[i] = stop_burn_time

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
                nu_rad=ν,
                a_m=a
            )
            _MANEUVER_TRACE_LAST_WINDOW[trace_key] = (start_burn_time, stop_burn_time)
        end
    end
end
