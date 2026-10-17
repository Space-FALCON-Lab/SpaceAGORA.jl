# The control-moment-gyro actuator kind. Included by `Components`
# (src/models/vehicle/spacecraft/components.jl) before the reaction-wheel
# assembly, which holds the mounted CMGs; the actuator contract it implements
# is in `actuators.jl`.

# ---- Record ------------------------------------------------------------------

"""
    ControlMomentGyro(; rotor_momentum_nms, gimbal_axis, spin_axis, max_gimbal_rate_rad_s=Inf, gimbal_angle_rad=0.0,
                        min_gimbal_angle_rad=-Inf, max_gimbal_angle_rad=Inf)

A single-gimbal control moment gyro on the root link: a rotor spinning at a
constant momentum `rotor_momentum_nms` (N m s) whose spin axis turns about
`gimbal_axis` (body frame). At gimbal angle δ the rotor's momentum is
`h = rotor_momentum_nms (cos δ s + sin δ (g × s))`, with `s = spin_axis`
(perpendicular to `g = gimbal_axis`) at δ = 0. Its level is the gimbal rate
δ̇ in rad/s, within `±max_gimbal_rate_rad_s`. The gimbal stops at
`min_gimbal_angle_rad` and `max_gimbal_angle_rad`: there it does not move
further, whatever the command (`cmg_angle_saturate`).

The gimbal angle is integrated with the vehicle's state (`cmg_angles`,
present only when a CMG is mounted), starting from `gimbal_angle_rad`. The
vehicle feels the torque `-dh/dt = -δ̇ (g × h)` and the rotor's momentum in
the gyroscopic term of its attitude equation. Mount with `add_actuator!` on
the root link; `cmg_steering` turns a desired torque into gimbal rates.
"""
mutable struct ControlMomentGyro <: AbstractActuator
    rotor_momentum_nms::Float64
    gimbal_axis::SVector{3, Float64}
    spin_axis::SVector{3, Float64}
    max_gimbal_rate_rad_s::Float64
    gimbal_angle_rad::Float64
    min_gimbal_angle_rad::Float64
    max_gimbal_angle_rad::Float64

    function ControlMomentGyro(rotor_momentum_nms::Real, gimbal_axis, spin_axis, max_gimbal_rate_rad_s::Real, gimbal_angle_rad::Real,
                               min_gimbal_angle_rad::Real, max_gimbal_angle_rad::Real)
        (isfinite(rotor_momentum_nms) && rotor_momentum_nms > 0.0) ||
            throw(ArgumentError("ControlMomentGyro rotor_momentum_nms must be finite and positive, got $(rotor_momentum_nms)."))
        g = SVector{3, Float64}(gimbal_axis)
        s = SVector{3, Float64}(spin_axis)
        ng, ns = norm(g), norm(s)
        (isfinite(ng) && ng > 0.0 && isfinite(ns) && ns > 0.0) ||
            throw(ArgumentError("ControlMomentGyro gimbal_axis and spin_axis must be finite and nonzero."))
        g, s = g / ng, s / ns
        abs(dot(g, s)) <= 1e-9 || throw(ArgumentError("ControlMomentGyro spin_axis must be perpendicular to gimbal_axis."))
        (max_gimbal_rate_rad_s > 0.0 && !isnan(max_gimbal_rate_rad_s)) ||
            throw(ArgumentError("ControlMomentGyro max_gimbal_rate_rad_s must be positive, got $(max_gimbal_rate_rad_s)."))
        isfinite(gimbal_angle_rad) || throw(ArgumentError("ControlMomentGyro gimbal_angle_rad must be finite, got $(gimbal_angle_rad)."))
        (isnan(min_gimbal_angle_rad) || isnan(max_gimbal_angle_rad) || min_gimbal_angle_rad > max_gimbal_angle_rad) &&
            throw(ArgumentError("ControlMomentGyro needs min_gimbal_angle_rad <= max_gimbal_angle_rad, got $(min_gimbal_angle_rad) and $(max_gimbal_angle_rad)."))
        min_gimbal_angle_rad <= gimbal_angle_rad <= max_gimbal_angle_rad ||
            throw(ArgumentError("ControlMomentGyro gimbal_angle_rad $(gimbal_angle_rad) is outside its stops [$(min_gimbal_angle_rad), $(max_gimbal_angle_rad)]."))
        return new(Float64(rotor_momentum_nms), g, s, Float64(max_gimbal_rate_rad_s), Float64(gimbal_angle_rad),
                   Float64(min_gimbal_angle_rad), Float64(max_gimbal_angle_rad))
    end
end

ControlMomentGyro(rotor_momentum_nms::Real, gimbal_axis, spin_axis, max_gimbal_rate_rad_s::Real, gimbal_angle_rad::Real) =
    ControlMomentGyro(rotor_momentum_nms, gimbal_axis, spin_axis, max_gimbal_rate_rad_s, gimbal_angle_rad, -Inf, Inf)

ControlMomentGyro(; rotor_momentum_nms::Real, gimbal_axis, spin_axis, max_gimbal_rate_rad_s::Real=Inf, gimbal_angle_rad::Real=0.0,
                  min_gimbal_angle_rad::Real=-Inf, max_gimbal_angle_rad::Real=Inf) =
    ControlMomentGyro(rotor_momentum_nms, gimbal_axis, spin_axis, max_gimbal_rate_rad_s, gimbal_angle_rad,
                      min_gimbal_angle_rad, max_gimbal_angle_rad)

# ---- Methods of the actuator contract --------------------------------------

actuator_limits(cmg::ControlMomentGyro) = (-cmg.max_gimbal_rate_rad_s, cmg.max_gimbal_rate_rad_s)

"""
    actuator_wrench(cmg::ControlMomentGyro, level, ctx) -> (force_body, torque_body)

None here: a CMG's torque depends on its gimbal angle, a state; the engine
applies `-δ̇ (g × h)` in the attitude equation.
"""
@inline actuator_wrench(::ControlMomentGyro, ::Real, ::ActuatorContext) =
    (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))

"""
    cmg_momentum(cmg::ControlMomentGyro, angle) -> SVector{3}

The rotor's momentum in the body frame at gimbal angle `angle`, N m s.
"""
@inline cmg_momentum(cmg::ControlMomentGyro, angle::Real)::SVector{3, Float64} =
    cmg.rotor_momentum_nms * (cos(angle) * cmg.spin_axis + sin(angle) * cross(cmg.gimbal_axis, cmg.spin_axis))

"""
    cmg_torque_axis(cmg::ControlMomentGyro, angle) -> SVector{3}

`g × h` at gimbal angle `angle`: the CMG puts the torque `-δ̇ (g × h)` on the
vehicle at gimbal rate δ̇ (the column of the CMG's torque Jacobian, negated).
"""
@inline cmg_torque_axis(cmg::ControlMomentGyro, angle::Real)::SVector{3, Float64} =
    cross(cmg.gimbal_axis, cmg_momentum(cmg, angle))

"""
    cmg_angle_saturate(cmg::ControlMomentGyro, rate, angle) -> Float64

The gimbal rate the CMG applies at gimbal angle `angle` when commanded `rate`:
`rate`, or zero when the gimbal is at or beyond a stop and `rate` would drive
it further. The attitude equation applies it to the gimbal angle and to the
torque alike.
"""
@inline function cmg_angle_saturate(cmg::ControlMomentGyro, rate::Real, angle::Real)::Float64
    r = Float64(rate)
    ((r > 0.0 && angle >= cmg.max_gimbal_angle_rad) || (r < 0.0 && angle <= cmg.min_gimbal_angle_rad)) && return 0.0
    return r
end
