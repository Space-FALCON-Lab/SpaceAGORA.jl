"""
The actuator contract: what every actuator kind implements, what a control
effector implements to command the actuators mounted on a spacecraft, and the
allocation of a wrench command over a set of actuators.

Each actuator kind lives beside this file, in one file per kind (the thruster
in `thruster.jl`). Mounting actuators on links and evaluating a mounted set
are in `mounting.jl`, which needs the spacecraft types defined after this
module.
"""
module Actuators

using LinearAlgebra
using StaticArrays
using ..AbstractTypes: AbstractActuator
using ..EffectorSampling: EffectorEnvironmentRequirements

export ActuatorContext, ActuatorRequirements, STANDARD_GRAVITY_MPS2, propellant_mass_rate
export actuator_wrench, actuator_limits, actuator_saturate, actuator_consumption
export actuator_fraction, actuator_unit_wrench, actuator_requirements
export actuation_mode, actuator_levels
export allocate_actuators

"""Standard gravity used to turn specific impulse into exhaust velocity, m/s²."""
const STANDARD_GRAVITY_MPS2 = 9.80665

"""
    propellant_mass_rate(thrust_n, isp_s) -> Float64

The propellant a rocket engine uses per second (kg/s) at `thrust_n` newtons
with a specific impulse of `isp_s` seconds: `thrust_n / (isp_s * g0)`. Every
thrust model of the package uses this one formula; callers handle zero or
invalid inputs.
"""
@inline propellant_mass_rate(thrust_n::Real, isp_s::Real)::Float64 = Float64(thrust_n) / (Float64(isp_s) * STANDARD_GRAVITY_MPS2)

const _ZERO3 = SVector{3, Float64}(0.0, 0.0, 0.0)
const _IDENTITY3 = SMatrix{3, 3, Float64, 9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
const _IDENTITY_Q = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)

"""
    ActuatorContext(; link_r_body, link_to_body, mass_kg, q_ib, ω_body, magnetic_field_body, environment=nothing)

What an actuator kind may use when it is evaluated:

- `link_r_body`: the mounting link's centre of mass in the spacecraft body
  frame, relative to the body reference point, metres;
- `link_to_body`: the rotation from the link frame to the body frame;
- `mass_kg`: the vehicle's current mass;
- `q_ib`: the attitude quaternion, inertial to body, scalar last;
- `ω_body`: the body angular velocity, rad/s;
- `magnetic_field_body`: the local magnetic field in the body frame, T, when
  a mounted actuator asks for it (see [`actuator_requirements`](@ref)); NaN
  otherwise;
- `environment`: the sampled environment (an `EnvironmentSample`: planet
  frame, atmosphere, Sun) when a mounted actuator asks for any of it;
  `nothing` otherwise.

The link pose is the link's configured geometry. The attitude and the rate are
the identity and zero when the run does not propagate attitude.
"""
struct ActuatorContext{E}
    link_r_body::SVector{3, Float64}
    link_to_body::SMatrix{3, 3, Float64, 9}
    mass_kg::Float64
    q_ib::SVector{4, Float64}
    ω_body::SVector{3, Float64}
    magnetic_field_body::SVector{3, Float64}
    environment::E
end

ActuatorContext(link_r_body, link_to_body, mass_kg, q_ib, ω_body, magnetic_field_body) =
    ActuatorContext(link_r_body, link_to_body, mass_kg, q_ib, ω_body, magnetic_field_body, nothing)

ActuatorContext(;
    link_r_body=_ZERO3,
    link_to_body=_IDENTITY3,
    mass_kg::Real=NaN,
    q_ib=_IDENTITY_Q,
    ω_body=_ZERO3,
    magnetic_field_body=SVector{3, Float64}(NaN, NaN, NaN),
    environment=nothing,
) = ActuatorContext(SVector{3, Float64}(link_r_body), SMatrix{3, 3, Float64, 9}(link_to_body),
                    Float64(mass_kg), SVector{4, Float64}(q_ib), SVector{3, Float64}(ω_body),
                    SVector{3, Float64}(magnetic_field_body), environment)

"""
    ActuatorRequirements(; planet_frame=false, atmosphere=false, solar=false, magnetic_field=false)

What an actuator kind needs sampled before it is evaluated: the planet-fixed
frame, the atmosphere (density, temperature, wind), the Sun's position, the
local magnetic field. The engine samples the union over the actuators mounted
on a spacecraft, once per evaluation, and passes it in the
[`ActuatorContext`](@ref): the field as `magnetic_field_body`, the rest as
`environment`.
"""
struct ActuatorRequirements
    planet_frame::Bool
    atmosphere::Bool
    solar::Bool
    magnetic_field::Bool
end

ActuatorRequirements(; planet_frame::Bool=false, atmosphere::Bool=false, solar::Bool=false, magnetic_field::Bool=false) =
    ActuatorRequirements(planet_frame, atmosphere, solar, magnetic_field)

Base.:|(a::ActuatorRequirements, b::ActuatorRequirements) = ActuatorRequirements(
    a.planet_frame | b.planet_frame, a.atmosphere | b.atmosphere, a.solar | b.solar, a.magnetic_field | b.magnetic_field)

"""Whether the requirements ask for any sampled environment besides the magnetic field."""
needs_environment(r::ActuatorRequirements)::Bool = r.planet_frame | r.atmosphere | r.solar

"""The same requirements as the effector sampler takes them (no third bodies)."""
effector_requirements(r::ActuatorRequirements) =
    EffectorEnvironmentRequirements(planet_frame=r.planet_frame, atmosphere=r.atmosphere, solar=r.solar)

"""
    actuator_requirements(actuator) -> ActuatorRequirements

Optional: what the actuator needs sampled to be evaluated. Default: nothing.
The torque rod asks for the magnetic field.
"""
@inline actuator_requirements(::Any) = ActuatorRequirements()

"""
    actuator_wrench(actuator, level, ctx::ActuatorContext) -> (force_body, torque_body)

Required for every actuator kind: the force (N) and the torque (N m) the
actuator puts on the vehicle at the commanded `level`, both in the spacecraft
body frame, the torque about the body reference point. `level` is in the
kind's unit (newtons for a thruster) and has already been passed through
[`actuator_saturate`](@ref).
"""
function actuator_wrench end

"""
    actuator_limits(actuator) -> (lower, upper)

Required for every actuator kind: the bounds of its level, in the kind's unit.
An infinite bound means unlimited.
"""
function actuator_limits end

"""
    actuator_saturate(actuator, level) -> Float64

The level the hardware can produce. Default: `level` clamped to
[`actuator_limits`](@ref); a non-finite level gives the idle level, zero
clamped to the limits.
"""
@inline function actuator_saturate(actuator, level)::Float64
    lower, upper = actuator_limits(actuator)
    lo = Float64(lower)
    hi = Float64(upper)
    x = Float64(level)
    isfinite(x) || return clamp(0.0, lo, hi)
    return clamp(x, lo, hi)
end

"""
    actuator_consumption(actuator, level, ctx::ActuatorContext) -> (mass_rate_kg_s, power_w)

Optional: the propellant the actuator uses per second at `level` (kg/s,
positive while consuming) and the electrical power it draws (W). Default: no
consumption.
"""
@inline actuator_consumption(::Any, ::Any, ::Any) = (mass_rate_kg_s = 0.0, power_w = 0.0)

"""
    actuator_fraction(actuator, level) -> Float64

Optional: `level` as a fraction of the actuator's capability, from 0 to 1. It
fills the `thruster_level` result columns and the viewer's plumes. Default: the
magnitude of `level` over the larger magnitude of the two limits; zero when
that capability is zero or not finite.
"""
@inline function actuator_fraction(actuator, level)::Float64
    lower, upper = actuator_limits(actuator)
    capability = max(abs(Float64(lower)), abs(Float64(upper)))
    (isfinite(capability) && capability > 0.0) || return 0.0
    fraction = abs(Float64(level)) / capability
    return isfinite(fraction) ? clamp(fraction, 0.0, 1.0) : 0.0
end

"""
    actuator_unit_wrench(actuator, ctx::ActuatorContext) -> (force_body, torque_body)

Optional: the wrench at level one, which is the actuator's column in an
allocation matrix. Default: `actuator_wrench(actuator, 1.0, ctx)`, right for a
kind whose wrench is linear in its level.
"""
@inline actuator_unit_wrench(actuator, ctx) = actuator_wrench(actuator, 1.0, ctx)

"""
    actuation_mode(controller) -> Symbol

How a control effector acts on the vehicle.

- `:ideal`, the default: the effector's own `control_wrench`,
  `control_propellant_rate` and `control_wheel_torque` give the wrench and the
  mass rate.
- `:mounted`: the engine derives the wrench, the propellant mass rate and the
  thruster firing levels from the actuators mounted on the spacecraft, at the
  levels [`actuator_levels`](@ref) returns. A reaction-wheel torque the
  effector reports through `control_wheel_torque` is still applied.
"""
@inline actuation_mode(::Any) = :ideal

"""
    actuator_levels(controller, i) -> Union{Nothing, AbstractVector{<:Real}}

For a control effector in `:mounted` mode: the level of every actuator mounted
on spacecraft `i` (its position in the run, not its ID), in mounting order
(see `mounted_actuators`), held between control ticks. Return `nothing` for a
spacecraft the effector does not drive.

The engine calls this at every right-hand-side evaluation, and in the power,
fault and result readers, for the spacecraft the controller's `spacecraft_scope`
names (for every spacecraft in SpaceAGORA 0.6 and earlier), so compute the levels
in `update_control!` and return stored values here.

The defaults of this and the other optional actuator methods take untyped
arguments, so a method written for one kind or one controller, typed or not,
is always the one that runs.
"""
@inline actuator_levels(::Any, ::Any) = nothing

@inline function _bounds(bound::Real, n::Int)::Vector{Float64}
    return fill(Float64(bound), n)
end

@inline function _bounds(bound::AbstractVector{<:Real}, n::Int)::Vector{Float64}
    length(bound) == n || throw(DimensionMismatch("allocate_actuators: $(length(bound)) bounds for $(n) actuators."))
    return Float64.(bound)
end

"""
    allocate_actuators(B, command, lower, upper; method=:least_norm, minimum_fraction=0.0) -> Vector{Float64}

Actuator levels that produce `command` through the allocation matrix `B`, one
column per actuator holding its wrench at level one (see `allocation_matrix`),
within `lower` and `upper` (scalars, or one bound per actuator).

- `:least_norm`: the least-norm solution `pinv(B) * command`, each level
  clipped to its bounds; a non-finite level gives zero. `command` has one entry
  per row of `B`.
- `:axis_clamp`: six actuators whose force columns push along +x, -x, +y, -y,
  +z, -z, in that order. Each body axis of the force `command` (its first
  three entries) is served by the jet on its side, up to that jet's upper
  bound; jets below `minimum_fraction` of their upper bound are turned off.
  Every level is then clamped to its bounds, so a positive lower bound keeps
  an idle jet at that bound.
"""
function allocate_actuators(
    B::AbstractMatrix{<:Real},
    command::AbstractVector{<:Real},
    lower,
    upper;
    method::Symbol=:least_norm,
    minimum_fraction::Real=0.0,
)::Vector{Float64}
    n = size(B, 2)
    lo = _bounds(lower, n)
    hi = _bounds(upper, n)
    if method === :least_norm
        size(B, 1) == length(command) ||
            throw(DimensionMismatch("allocate_actuators: a command of length $(length(command)) for $(size(B, 1)) matrix rows."))
        return _least_norm_levels(pinv(Matrix{Float64}(B)), command, lo, hi)
    elseif method === :axis_clamp
        _check_axis_pairs(B)
        isfinite(minimum_fraction) && 0 <= minimum_fraction < 1 ||
            throw(ArgumentError("allocate_actuators: minimum_fraction must be finite in [0, 1)."))
        length(command) >= 3 || throw(DimensionMismatch("allocate_actuators: :axis_clamp needs a force command of three entries."))
        return _axis_clamp_levels(command, lo, hi, Float64(minimum_fraction))
    end
    throw(ArgumentError("allocate_actuators: method must be :least_norm or :axis_clamp, got $(repr(method))."))
end

# Same arithmetic, in the same order, as the descent controller's jet levels:
# u_j = P[j, 1] c_1 + P[j, 2] c_2 + ..., then clipped.
function _least_norm_levels(P::AbstractMatrix{Float64}, command::AbstractVector{<:Real},
                            lo::Vector{Float64}, hi::Vector{Float64})::Vector{Float64}
    n = size(P, 1)
    levels = zeros(Float64, n)
    @inbounds for j in 1:n
        u = P[j, 1] * Float64(command[1])
        for k in 2:size(P, 2)
            u += P[j, k] * Float64(command[k])
        end
        levels[j] = isfinite(u) ? clamp(u, lo[j], hi[j]) : 0.0
    end
    return levels
end

function _check_axis_pairs(B::AbstractMatrix{<:Real})
    size(B, 2) == 6 && size(B, 1) >= 3 ||
        throw(ArgumentError("allocate_actuators: :axis_clamp needs six actuator columns with at least three force rows, got $(size(B))."))
    @inbounds for j in 1:6
        axis = (j + 1) ÷ 2
        sign = isodd(j) ? 1.0 : -1.0
        f = SVector{3, Float64}(B[1, j], B[2, j], B[3, j])
        n = norm(f)
        (isfinite(n) && n > 0.0 && abs(f[axis] / n - sign) <= 1e-9) ||
            throw(ArgumentError("allocate_actuators: :axis_clamp needs columns along +x, -x, +y, -y, +z, -z; column $(j) is $(f)."))
    end
    return nothing
end

# Same arithmetic as the RPO six-axis allocator, then each level clamped to
# both bounds (with the allocator's zero lower bounds this changes nothing).
function _axis_clamp_levels(command::AbstractVector{<:Real}, lo::Vector{Float64}, hi::Vector{Float64}, minimum_fraction::Float64)::Vector{Float64}
    f = zeros(Float64, 6)
    desired = SVector{3, Float64}(command[1], command[2], command[3])
    @inbounds for axis in 1:3
        pos_idx = 2 * axis - 1
        neg_idx = 2 * axis
        if desired[axis] >= 0.0
            f[pos_idx] = min(desired[axis], hi[pos_idx])
        else
            f[neg_idx] = min(-desired[axis], hi[neg_idx])
        end
    end
    for j in 1:6
        f[j] <= minimum_fraction * hi[j] && (f[j] = 0.0)
        f[j] = clamp(f[j], lo[j], hi[j])
    end
    return f
end

end # module Actuators
