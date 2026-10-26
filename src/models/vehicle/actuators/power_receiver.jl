# The power-receiver terminal kind. Included by `Components`
# (src/models/vehicle/spacecraft/components.jl), which owns and exports the
# record; the actuator contract it implements is in `actuators.jl`.

# ---- Record ---------------------------------------------------------------

"""
    PowerReceiver(; normal, radius_m, efficiency, max_power_w=Inf, location=zero)

The receiving end of a power beam, mounted on a spacecraft link with
`add_actuator!`. It is interlink terminal hardware: a `PowerBeamLink` can end on
it, and the power it delivers is generation in its spacecraft's power budget.

- `radius_m`: the radius (m) of its circular aperture, which sets how much of a
  beam it captures.
- `efficiency`: the fraction of the captured optical power it turns into
  electrical power, in (0, 1].
- `max_power_w`: its converter's limit (W) on the power it delivers.
- `normal`: the aperture's facing, and `location`, its position, both in the
  link frame. They are kept for the incidence-angle and pointing models to
  come (I-7); the delivered power does not use them yet.

It applies no force or torque and draws no power: its level is always 0. It
mounts through the generic actuator path, so it takes its place among the
actuators of user-defined kinds in the mounting order.
"""
struct PowerReceiver <: AbstractActuator
    normal::SVector{3, Float64}
    location::SVector{3, Float64}
    radius_m::Float64
    efficiency::Float64
    max_power_w::Float64

    function PowerReceiver(normal, location, radius_m::Real, efficiency::Real, max_power_w::Real)
        axis = SVector{3, Float64}(normal)
        n = norm(axis)
        (isfinite(n) && n > 0.0) || throw(ArgumentError("PowerReceiver normal must be finite and nonzero, got $(axis)."))
        (isfinite(radius_m) && radius_m > 0.0) || throw(ArgumentError("PowerReceiver radius_m must be finite and positive, got $(radius_m)."))
        (isfinite(efficiency) && 0.0 < efficiency <= 1.0) || throw(ArgumentError("PowerReceiver efficiency must be in (0, 1], got $(efficiency)."))
        (max_power_w > 0.0 && !isnan(max_power_w)) || throw(ArgumentError("PowerReceiver max_power_w must be positive, got $(max_power_w)."))
        return new(axis / n, SVector{3, Float64}(location), Float64(radius_m), Float64(efficiency), Float64(max_power_w))
    end
end

PowerReceiver(; normal, radius_m::Real, efficiency::Real, max_power_w::Real=Inf,
              location=SVector{3, Float64}(0.0, 0.0, 0.0)) =
    PowerReceiver(normal, location, radius_m, efficiency, max_power_w)

# ---- Methods of the actuator contract --------------------------------------

actuator_limits(::PowerReceiver) = (0.0, 0.0)

@inline actuator_wrench(::PowerReceiver, ::Real, ::ActuatorContext) =
    (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
