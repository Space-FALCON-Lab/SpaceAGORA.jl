# The actuator-fault record (failure modes; see `ActuatorFaults` in faults.jl).
# Included by `Components` (src/models/vehicle/spacecraft/components.jl), which
# owns it.

const _FAULT_KINDS = (:failed, :stuck, :degraded, :biased)

"""
    ActuatorFault(actuator; kind, start_s=0.0, stop_s=Inf, value=NaN)

A failure of one mounted actuator (the record mounted with `add_actuator!`)
between `start_s` and `stop_s` (seconds into the run):

- `:failed`: it applies level zero; a failed hinge or gimbal holds its angle.
- `:stuck`: it holds `value`, or, when `value` is NaN, the level it was
  commanded to (a hinge or gimbal: the angle it had) when the fault began.
- `:degraded`: its level and its limits scale by `value`, from 0 to 1.
- `:biased`: `value` is added to its command before saturation.

Faults act on the mounted path only: the wrench, the propellant and power,
hinge moves, the `actuator_level` and thruster-level columns and the plumes
see them; the controller sees its own commands. On a thruster with firing
dynamics a fault acts on the thrust it delivers, from the fault's start
rather than from the next control tick: failed, none; stuck, the held level;
degraded, the fraction `value`; biased, offset by `value` within its limits.
Mount with `add_fault!`.
"""
mutable struct ActuatorFault
    actuator::AbstractActuator
    kind::Symbol
    start_s::Float64
    stop_s::Float64
    value::Float64
    # Run state: active between the window's ends (switched between steps), the
    # level a stuck actuator holds, its position in the mounting order.
    active::Bool
    held_level::Float64
    index::Int

    function ActuatorFault(actuator::AbstractActuator, kind::Symbol, start_s::Real, stop_s::Real, value::Real)
        kind in _FAULT_KINDS || throw(ArgumentError("ActuatorFault kind must be one of $(_FAULT_KINDS), got $(repr(kind))."))
        (isfinite(start_s) && !isnan(stop_s) && stop_s > start_s) ||
            throw(ArgumentError("ActuatorFault needs a finite start_s before stop_s, got $(start_s) and $(stop_s)."))
        if kind === :degraded
            (0.0 <= value <= 1.0) || throw(ArgumentError("A :degraded fault needs a value from 0 to 1, got $(value)."))
        elseif kind === :biased
            isfinite(value) || throw(ArgumentError("A :biased fault needs a finite value, got $(value)."))
        elseif kind === :stuck
            (isnan(value) || isfinite(value)) || throw(ArgumentError("A :stuck fault needs a finite value or NaN, got $(value)."))
        end
        return new(actuator, kind, Float64(start_s), Float64(stop_s), Float64(value), false, NaN, 0)
    end
end

ActuatorFault(actuator::AbstractActuator; kind::Symbol, start_s::Real=0.0, stop_s::Real=Inf, value::Real=NaN) =
    ActuatorFault(actuator, kind, start_s, stop_s, value)

"""The faults of a spacecraft (kept on its root link) and a scratch level vector."""
mutable struct FaultSet
    faults::Vector{ActuatorFault}
    scratch::Vector{Float64}
end
FaultSet() = FaultSet(ActuatorFault[], Float64[])
