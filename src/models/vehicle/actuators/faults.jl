"""
    ActuatorFaults

Failure modes of mounted actuators: an [`ActuatorFault`](@ref) names a
mounted actuator, a kind and a time window, and the engine applies it to that
actuator's level between the window's ends. Mount faults with
[`add_fault!`](@ref).
"""
module ActuatorFaults

using StaticArrays
using ..AbstractTypes: AbstractActuator
using ..Components: PanelHinge, ThrusterGimbal, ActuatorFault, FaultSet
using ..SpacecraftModels: Link, SpacecraftModel
using ..ActuatorMounting: mounted_actuators
using ..Actuators: actuator_limits

export add_fault!, has_faults, faulted_levels, refresh_faults!, fault_edge, resolve_faults!

"""
    add_fault!(model::SpacecraftModel, fault::ActuatorFault) -> fault

Give `model` a fault on one of its mounted actuators (the record mounted with
`add_actuator!`, matched by identity).
"""
function add_fault!(model::SpacecraftModel, fault::ActuatorFault)
    any(entry -> entry.actuator === fault.actuator, mounted_actuators(model)) || throw(ArgumentError(
        "ActuatorFault names an actuator that is not mounted on this spacecraft; mount it with add_actuator! first."))
    push!(model.root.faults.faults, fault)
    return fault
end

"""Whether `model` carries any actuator fault."""
has_faults(model::SpacecraftModel)::Bool = !isempty(model.root.faults.faults)

"""
    resolve_faults!(model) -> model

Record each fault's position in the mounting order, as the level vectors use
it. Called at the start of a run.
"""
function resolve_faults!(model::SpacecraftModel)
    mounted = mounted_actuators(model)
    for f in model.root.faults.faults
        k = findfirst(entry -> entry.actuator === f.actuator, mounted)
        k === nothing && throw(ArgumentError("ActuatorFault names an actuator no longer mounted on spacecraft $(model.id)."))
        f.index = k
    end
    return model
end

"""
    refresh_faults!(model, levels, t) -> Bool

Set which faults are active at time `t` (inside `[start_s, stop_s)`), capturing
a stuck actuator's level at onset from `levels` (the command it held, or a
hinge's or gimbal's current angle). Returns whether any fault switched. Called at the
start of a run and at the faults' window ends, where the solver stops.
"""
function refresh_faults!(model::SpacecraftModel, levels, t::Real)::Bool
    changed = false
    for f in model.root.faults.faults
        active = f.start_s <= t < f.stop_s
        if active && !f.active
            f.held_level = _onset_level(f, levels)
        end
        active == f.active && continue
        f.active = active
        changed = true
    end
    return changed
end

function _onset_level(f::ActuatorFault, levels)::Float64
    isnan(f.value) || return f.value
    f.actuator isa Union{PanelHinge, ThrusterGimbal} && return f.actuator.angle_rad
    (levels === nothing || f.index < 1 || f.index > length(levels)) && return 0.0
    return Float64(levels[f.index])
end

"""Whether `t` is an end of one of `model`'s fault windows."""
fault_edge(model::SpacecraftModel, t::Real)::Bool =
    any(f -> t == f.start_s || t == f.stop_s, model.root.faults.faults)

"""
    faulted_levels(model, levels) -> AbstractVector

`levels` (one per mounted actuator) as the active faults leave them: a failed
actuator at zero (a failed hinge or gimbal holds its angle), a stuck one at its held
level, a degraded one scaled with its limits, a biased one offset. `levels`
itself when no fault is active.
"""
function faulted_levels(model::SpacecraftModel, levels::AbstractVector)
    fs = model.root.faults
    any(f -> f.active, fs.faults) || return levels
    out = fs.scratch
    resize!(out, length(levels))
    copyto!(out, levels)
    for f in fs.faults
        (f.active && 1 <= f.index <= length(out)) || continue
        out[f.index] = _apply_fault(f, out[f.index])
    end
    return out
end

function _apply_fault(f::ActuatorFault, level::Float64)::Float64
    if f.kind === :failed
        return f.actuator isa Union{PanelHinge, ThrusterGimbal} ? f.held_level : 0.0
    elseif f.kind === :stuck
        return f.held_level
    elseif f.kind === :degraded
        lo, hi = actuator_limits(f.actuator)
        return clamp(level * f.value, f.value * Float64(lo), f.value * Float64(hi))
    else   # :biased
        return level + f.value
    end
end

end # module ActuatorFaults
