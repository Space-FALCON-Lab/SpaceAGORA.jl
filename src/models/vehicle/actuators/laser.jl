# The laser-terminal actuator kind. Included by `Components`
# (src/models/vehicle/spacecraft/components.jl), which owns and exports the
# record; the actuator contract it implements is in `actuators.jl`.

const _SPEED_OF_LIGHT_MPS = 299_792_458.0

# ---- Record ---------------------------------------------------------------

"""
    LaserTerminal(; boresight, max_power_w, location=zero, magnification=1.0, range_m=200e3, efficiency=1.0)

One optical terminal mounted on a spacecraft link with `add_actuator!`. Its
level is the optical power it transmits, in watts, from zero to
`max_power_w`, along `boresight` (in the link frame) from `location` (link
frame, metres). The beam carries momentum away, so the vehicle recoils with a
force of `level / c` opposite the boresight; the terminal draws
`level / efficiency` watts of electrical power.

`magnification` and `range_m` describe an inter-satellite link between two
terminals of this design: the open-cavity magnification factor and the
maximum range of `LaserMomentumLink(terminal)`, the link kind whose force
pushes linked spacecraft apart.
"""
struct LaserTerminal <: AbstractActuator
    boresight::SVector{3, Float64}
    location::SVector{3, Float64}
    max_power_w::Float64
    magnification::Float64
    range_m::Float64
    efficiency::Float64

    function LaserTerminal(boresight, location, max_power_w::Real, magnification::Real, range_m::Real, efficiency::Real)
        axis = SVector{3, Float64}(boresight)
        n = norm(axis)
        (isfinite(n) && n > 0.0) || throw(ArgumentError("LaserTerminal boresight must be finite and nonzero, got $(axis)."))
        (isfinite(max_power_w) && max_power_w > 0.0) || throw(ArgumentError("LaserTerminal max_power_w must be finite and positive, got $(max_power_w)."))
        (isfinite(magnification) && magnification > 0.0) || throw(ArgumentError("LaserTerminal magnification must be finite and positive, got $(magnification)."))
        (isfinite(range_m) && range_m > 0.0) || throw(ArgumentError("LaserTerminal range_m must be finite and positive, got $(range_m)."))
        (isfinite(efficiency) && 0.0 < efficiency <= 1.0) || throw(ArgumentError("LaserTerminal efficiency must be in (0, 1], got $(efficiency)."))
        return new(axis / n, SVector{3, Float64}(location), Float64(max_power_w), Float64(magnification),
                   Float64(range_m), Float64(efficiency))
    end
end

LaserTerminal(; boresight, max_power_w::Real, location=SVector{3, Float64}(0.0, 0.0, 0.0),
              magnification::Real=1.0, range_m::Real=200e3, efficiency::Real=1.0) =
    LaserTerminal(boresight, location, max_power_w, magnification, range_m, efficiency)

# ---- Methods of the actuator contract --------------------------------------

actuator_limits(terminal::LaserTerminal) = (0.0, terminal.max_power_w)

"""
    actuator_wrench(terminal::LaserTerminal, level, ctx) -> (force_body, torque_body)

The recoil of transmitting `level` watts: a force of `level / c` opposite the
boresight, applied at the terminal's location, in the body frame.
"""
@inline function actuator_wrench(terminal::LaserTerminal, level::Real, ctx::ActuatorContext)
    force_body = ctx.link_to_body * ((-Float64(level) / _SPEED_OF_LIGHT_MPS) * terminal.boresight)
    arm_body = ctx.link_to_body * terminal.location + ctx.link_r_body
    return force_body, cross(arm_body, force_body)
end

"""
    actuator_consumption(terminal::LaserTerminal, level, ctx) -> (mass_rate_kg_s, power_w)

No propellant; `level / efficiency` watts of electrical power.
"""
@inline actuator_consumption(terminal::LaserTerminal, level::Real, ::ActuatorContext) =
    (mass_rate_kg_s = 0.0, power_w = Float64(level) / terminal.efficiency)
