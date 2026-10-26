# The thruster-gimbal actuator kind. Included by `Components`
# (src/models/vehicle/spacecraft/components.jl), which owns and exports the
# record; the actuator contract it implements is in `actuators.jl`.

# ---- Record and state --------------------------------------------------------

# The gimbal's run state: its angle, the time of its latest tick and the angle before it.
const _GimbalRunState = NamedTuple{(:angle_rad, :last_update_s, :previous_angle_rad), Tuple{Float64, Float64, Float64}}

"""
    ThrusterGimbal(thruster; axis, min_angle_rad=-Inf, max_angle_rad=Inf, max_rate_rad_s=Inf, angle_rad=0.0)

One axis of a gimbal that turns the thrust of `thruster` (a `Thruster`
mounted on the same link) about `axis`, in the link frame. Its level is the
commanded gimbal angle in radians: the thrust direction turns by the
gimbal's current angle about `axis` (right-handed), about the thruster's
location, which stays fixed. The angle stays within
[`min_angle_rad`, `max_angle_rad`] and moves at most `max_rate_rad_s` between
control ticks; all three default to unlimited. `angle_rad` holds the current
angle, `previous_angle_rad` the angle before its latest tick (at
`last_update_s`). As for a `PanelHinge`, the angles a run starts from are kept
apart from them: each run with `isolate_state=false` starts from the previous
run's start, unless they were changed after that run's last tick.

A two-axis gimbal is two `ThrusterGimbal`s on the same thruster; their
rotations compose in mounting order (the first mounted turns the direction
first). Mount with `add_actuator!` on the thruster's link.

Like a hinge, a gimbal changes geometry, not the wrench by itself: the engine
applies its level at the control tick of the effector that commands it,
between accepted integration steps, and the thruster's wrench, propellant use
and allocation column follow the turned direction from then on.
"""
mutable struct ThrusterGimbal <: AbstractActuator
    thruster::Thruster
    axis::SVector{3, Float64}
    min_angle_rad::Float64
    max_angle_rad::Float64
    max_rate_rad_s::Float64
    angle_rad::Float64
    last_update_s::Float64
    previous_angle_rad::Float64
    run_start::_GimbalRunState                  # the run state a run starts from, first the configured one
    run_written::Union{Nothing, _GimbalRunState} # the run's last tick's write of it (`ActuatorMounting._start_run_state!`)

    function ThrusterGimbal(thruster::Thruster, axis, min_angle_rad::Real, max_angle_rad::Real, max_rate_rad_s::Real, angle_rad::Real)
        a = SVector{3, Float64}(axis)
        n = norm(a)
        (isfinite(n) && n > 0.0) || throw(ArgumentError("ThrusterGimbal axis must be finite and nonzero, got $(a)."))
        (min_angle_rad <= max_angle_rad && !isnan(min_angle_rad) && !isnan(max_angle_rad)) ||
            throw(ArgumentError("ThrusterGimbal needs min_angle_rad <= max_angle_rad, got $(min_angle_rad) and $(max_angle_rad)."))
        (max_rate_rad_s > 0.0 && !isnan(max_rate_rad_s)) ||
            throw(ArgumentError("ThrusterGimbal max_rate_rad_s must be positive, got $(max_rate_rad_s)."))
        (isfinite(angle_rad) && min_angle_rad <= angle_rad <= max_angle_rad) ||
            throw(ArgumentError("ThrusterGimbal angle_rad must be finite and within its limits, got $(angle_rad)."))
        return new(thruster, a / n, Float64(min_angle_rad), Float64(max_angle_rad), Float64(max_rate_rad_s), Float64(angle_rad), NaN,
                   Float64(angle_rad), _GimbalRunState((angle_rad, NaN, angle_rad)), nothing)
    end
end

ThrusterGimbal(thruster::Thruster; axis, min_angle_rad::Real=-Inf, max_angle_rad::Real=Inf, max_rate_rad_s::Real=Inf, angle_rad::Real=0.0) =
    ThrusterGimbal(thruster, axis, min_angle_rad, max_angle_rad, max_rate_rad_s, angle_rad)

# ---- Methods of the actuator contract --------------------------------------

actuator_limits(gimbal::ThrusterGimbal) = (gimbal.min_angle_rad, gimbal.max_angle_rad)

"""
    actuator_wrench(gimbal::ThrusterGimbal, level, ctx) -> (force_body, torque_body)

None: a gimbal turns its thruster's thrust, whose wrench carries the turn.
"""
@inline actuator_wrench(::ThrusterGimbal, ::Real, ::ActuatorContext) =
    (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))

"""
    gimbal_step(gimbal::ThrusterGimbal, level, t) -> Float64

The angle the gimbal reaches at a control tick at time `t` when commanded to
`level`: the level within the angle limits, then within `max_rate_rad_s`
times the time since the gimbal's previous tick (a rate-limited gimbal holds
its angle on its first tick, which only records the time), as `hinge_step`.
"""
gimbal_step(gimbal::ThrusterGimbal, level::Real, t::Real)::Float64 =
    _rate_limited_angle(actuator_saturate(gimbal, level), gimbal.angle_rad, gimbal.last_update_s, gimbal.max_rate_rad_s, t)

# The rotation by `angle` about the unit vector `axis` (Rodrigues).
@inline function _axis_rotation(axis::SVector{3, Float64}, angle::Float64)::SMatrix{3, 3, Float64, 9}
    K = SMatrix{3, 3, Float64, 9}(0.0, axis[3], -axis[2], -axis[3], 0.0, axis[1], axis[2], -axis[1], 0.0)
    return SMatrix{3, 3, Float64, 9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0) + sin(angle) * K + (1.0 - cos(angle)) * (K * K)
end

# The same rotation at the angles `angle_of(gimbal)` gives, for a time other
# than the current one (the saved thrust directions).
function _gimbal_rotation(thruster::Thruster, gimbals, angle_of)::SMatrix{3, 3, Float64, 9}
    R = SMatrix{3, 3, Float64, 9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    for g in gimbals
        g.thruster === thruster || continue
        R = _axis_rotation(g.axis, Float64(angle_of(g))) * R
    end
    return R
end

# The rotation of `thruster`'s thrust by the gimbals among `gimbals` that turn
# it, at their current angles, in mounting order (the first turns first).
function _gimbal_rotation(thruster::Thruster, gimbals)::SMatrix{3, 3, Float64, 9}
    R = SMatrix{3, 3, Float64, 9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    for g in gimbals
        g.thruster === thruster || continue
        R = _axis_rotation(g.axis, g.angle_rad) * R
    end
    return R
end
