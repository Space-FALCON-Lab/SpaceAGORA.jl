# Magnetic dipoles: the passive `Magnet` and the torque-rod actuator kind.
# Included by `Components` (src/models/vehicle/spacecraft/components.jl), which
# owns and exports both records; the actuator contract is in `actuators.jl`.

"""
    Magnet(; m, location)

A permanent magnetic dipole fixed to a link: `m` in the body frame, A m².
The `MagneticTorqueRodModel` effector applies the torque of every link's
`magnets`; a `Magnet` is not commanded.
"""
@kwdef mutable struct Magnet
    m::MVector{3, Float64} = MVector{3, Float64}(zeros(3)) # Magnetic dipole moment in the body frame, A*m^2
    location::MVector{3, Float64} = MVector{3, Float64}(zeros(3)) # Location in the link frame, relative to the CoM of the link, m
end

# ---- The torque-rod actuator kind --------------------------------------------

"""
    Magnetorquer(; axis, max_dipole_am2=Inf, power_at_max_dipole_w=0.0)

One magnetic torque rod mounted on a spacecraft link with `add_actuator!`. Its
level is the dipole it is commanded to, in ampere square metres, along `axis`
(in the link frame); limits ±`max_dipole_am2`, unlimited by default. The torque
on the vehicle is the dipole crossed with the local magnetic field, which the
engine takes from the run's `MagneticTorqueRodModel` field source (a tilted
dipole or IGRF, Earth only), or the dipole when the run has none.

Power: a rod's dipole is proportional to its current, so its resistive power
grows with the square of the dipole, `power_at_max_dipole_w * (level /
max_dipole_am2)^2`. A nonzero `power_at_max_dipole_w` needs a finite
`max_dipole_am2`.

`Magnet` above is a permanent dipole, not an actuator: the
`MagneticTorqueRodModel` effector applies it whether or not anything commands
the vehicle.
"""
struct Magnetorquer <: AbstractActuator
    axis::SVector{3, Float64}
    max_dipole_am2::Float64
    power_at_max_dipole_w::Float64

    function Magnetorquer(axis, max_dipole_am2::Real, power_at_max_dipole_w::Real)
        a = SVector{3, Float64}(axis)
        n = norm(a)
        (isfinite(n) && n > 0.0) || throw(ArgumentError("Magnetorquer axis must be finite and nonzero, got $(a)."))
        (max_dipole_am2 > 0.0 && !isnan(max_dipole_am2)) || throw(ArgumentError("Magnetorquer max_dipole_am2 must be positive, got $(max_dipole_am2)."))
        (isfinite(power_at_max_dipole_w) && power_at_max_dipole_w >= 0.0) ||
            throw(ArgumentError("Magnetorquer power_at_max_dipole_w must be finite and nonnegative, got $(power_at_max_dipole_w)."))
        (power_at_max_dipole_w == 0.0 || isfinite(max_dipole_am2)) ||
            throw(ArgumentError("Magnetorquer power_at_max_dipole_w needs a finite max_dipole_am2."))
        return new(a / n, Float64(max_dipole_am2), Float64(power_at_max_dipole_w))
    end
end

Magnetorquer(axis, max_dipole_am2::Real) = Magnetorquer(axis, max_dipole_am2, 0.0)
Magnetorquer(; axis, max_dipole_am2::Real=Inf, power_at_max_dipole_w::Real=0.0) =
    Magnetorquer(axis, max_dipole_am2, power_at_max_dipole_w)

actuator_limits(rod::Magnetorquer) = (-rod.max_dipole_am2, rod.max_dipole_am2)
actuator_requirements(::Magnetorquer) = ActuatorRequirements(magnetic_field=true)

"""
    actuator_consumption(rod::Magnetorquer, level, ctx) -> (mass_rate_kg_s, power_w)

No propellant; the resistive power `power_at_max_dipole_w * (level / max_dipole_am2)^2`.
"""
@inline function actuator_consumption(rod::Magnetorquer, level::Real, ::ActuatorContext)
    rod.power_at_max_dipole_w == 0.0 && return (mass_rate_kg_s = 0.0, power_w = 0.0)
    fraction = Float64(level) / rod.max_dipole_am2
    return (mass_rate_kg_s = 0.0, power_w = rod.power_at_max_dipole_w * fraction * fraction)
end

"""
    actuator_wrench(rod::Magnetorquer, level, ctx) -> (force_body, torque_body)

The torque `m x B` of a dipole of `level` A m² along the rod's axis in the
body-frame field `ctx.magnetic_field_body` (T); no force. No field (NaN) gives
no torque.
"""
@inline function actuator_wrench(rod::Magnetorquer, level::Real, ctx::ActuatorContext)
    zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
    field = ctx.magnetic_field_body
    all(isfinite, field) || return (zero3, zero3)
    dipole_body = ctx.link_to_body * (Float64(level) * rod.axis)
    return zero3, cross(dipole_body, field)
end
