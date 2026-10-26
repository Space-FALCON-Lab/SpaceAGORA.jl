# The radio-terminal kind. Included by `Components`
# (src/models/vehicle/spacecraft/components.jl), which owns and exports the
# record; the actuator contract it implements is in `actuators.jl`.

# ---- Record ---------------------------------------------------------------

"""
    RadioTerminal(; boresight, max_power_w, gain_dbi, efficiency, noise_temperature_k=500.0, max_rate_bps=Inf, location=zero)

A radio mounted on a spacecraft link with `add_actuator!`. It is interlink
terminal hardware: a `DataLink` can start or end on it. Its level is the radio
power it transmits, in watts, from zero to `max_power_w`.

- `gain_dbi`: its antenna gain (dBi), the same for sending and receiving.
- `efficiency`: the fraction of its electrical draw it transmits, in (0, 1].
  It draws `level / efficiency` watts while it transmits, and nothing while it
  receives (a receiver's electronics are a constant `PowerLoad`).
- `noise_temperature_k`: its system noise temperature (K) when it receives.
- `max_rate_bps`: the most data it sends or receives, in bits per second.
- `boresight`: the antenna's facing, and `location`, its position, both in the
  link frame. They are kept for the pointing models to come (I-7); the rate
  does not use them yet.

It applies no force or torque. The battery's shedding treats it as a laser
terminal (tier 1, see `terminal_shed_tier`). It mounts through the generic
actuator path, so it takes its place among the actuators of user-defined kinds
in the mounting order.
"""
struct RadioTerminal <: AbstractActuator
    boresight::SVector{3, Float64}
    location::SVector{3, Float64}
    max_power_w::Float64
    gain_dbi::Float64
    efficiency::Float64
    noise_temperature_k::Float64
    max_rate_bps::Float64

    function RadioTerminal(boresight, location, max_power_w::Real, gain_dbi::Real, efficiency::Real,
                           noise_temperature_k::Real, max_rate_bps::Real)
        axis = SVector{3, Float64}(boresight)
        n = norm(axis)
        (isfinite(n) && n > 0.0) || throw(ArgumentError("RadioTerminal boresight must be finite and nonzero, got $(axis)."))
        (isfinite(max_power_w) && max_power_w > 0.0) || throw(ArgumentError("RadioTerminal max_power_w must be finite and positive, got $(max_power_w)."))
        isfinite(gain_dbi) || throw(ArgumentError("RadioTerminal gain_dbi must be finite, got $(gain_dbi)."))
        (isfinite(efficiency) && 0.0 < efficiency <= 1.0) || throw(ArgumentError("RadioTerminal efficiency must be in (0, 1], got $(efficiency)."))
        (isfinite(noise_temperature_k) && noise_temperature_k > 0.0) || throw(ArgumentError(
            "RadioTerminal noise_temperature_k must be finite and positive, got $(noise_temperature_k)."))
        (max_rate_bps > 0.0 && !isnan(max_rate_bps)) || throw(ArgumentError("RadioTerminal max_rate_bps must be positive, got $(max_rate_bps)."))
        return new(axis / n, SVector{3, Float64}(location), Float64(max_power_w), Float64(gain_dbi), Float64(efficiency),
                   Float64(noise_temperature_k), Float64(max_rate_bps))
    end
end

RadioTerminal(; boresight, max_power_w::Real, gain_dbi::Real, efficiency::Real, noise_temperature_k::Real=500.0,
              max_rate_bps::Real=Inf, location=SVector{3, Float64}(0.0, 0.0, 0.0)) =
    RadioTerminal(boresight, location, max_power_w, gain_dbi, efficiency, noise_temperature_k, max_rate_bps)

# ---- Methods of the actuator contract --------------------------------------

actuator_limits(terminal::RadioTerminal) = (0.0, terminal.max_power_w)

@inline actuator_wrench(::RadioTerminal, ::Real, ::ActuatorContext) =
    (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))

"""
    actuator_consumption(terminal::RadioTerminal, level, ctx) -> (mass_rate_kg_s, power_w)

No propellant; `level / efficiency` watts of electrical power.
"""
@inline actuator_consumption(terminal::RadioTerminal, level::Real, ::ActuatorContext) =
    (mass_rate_kg_s = 0.0, power_w = Float64(level) / terminal.efficiency)
