# The hinge actuator kind. Included by `Components`
# (src/models/vehicle/spacecraft/components.jl), which owns and exports the
# record; the actuator contract it implements is in `actuators.jl`.

# ---- Record and state --------------------------------------------------------

# The hinge's run state: its angle and the time of its latest tick.
const _HingeRunState = NamedTuple{(:angle_rad, :last_update_s), Tuple{Float64, Float64}}

"""
    PanelHinge(; axis, min_angle_rad=-Inf, max_angle_rad=Inf, max_rate_rad_s=Inf, angle_rad=0.0)

A hinge that turns one non-root link (a solar panel, a boom) about `axis`, in
the body frame, mounted on that link with `add_actuator!`. Its level is the
commanded hinge angle in radians; the link's attitude is the rotation by the
hinge's current angle about `axis`. The angle stays within
[`min_angle_rad`, `max_angle_rad`] and moves at most `max_rate_rad_s` between
control ticks; all three default to unlimited. `angle_rad` holds the current
angle (the state the engine advances at each control tick).

The angle a run starts from is kept apart from it: each run with
`isolate_state=false` starts from the angle the previous run started from, unless
`angle_rad` was changed after that run's last tick, and then from the new angle.

A hinge changes the vehicle's geometry, not its wrench: the engine applies its
level at the control tick of the effector that commands it, between accepted
integration steps, and the aerodynamic and radiation models see the new link
attitude from then on.
"""
mutable struct PanelHinge <: AbstractActuator
    axis::SVector{3, Float64}
    min_angle_rad::Float64
    max_angle_rad::Float64
    max_rate_rad_s::Float64
    angle_rad::Float64
    last_update_s::Float64
    run_start::_HingeRunState                  # the run state a run starts from, first the configured one
    run_written::Union{Nothing, _HingeRunState} # the run's last tick's write of it (`ActuatorMounting._start_run_state!`)

    function PanelHinge(axis, min_angle_rad::Real, max_angle_rad::Real, max_rate_rad_s::Real, angle_rad::Real)
        a = SVector{3, Float64}(axis)
        n = norm(a)
        (isfinite(n) && n > 0.0) || throw(ArgumentError("PanelHinge axis must be finite and nonzero, got $(a)."))
        (min_angle_rad <= max_angle_rad && !isnan(min_angle_rad) && !isnan(max_angle_rad)) ||
            throw(ArgumentError("PanelHinge needs min_angle_rad <= max_angle_rad, got $(min_angle_rad) and $(max_angle_rad)."))
        (max_rate_rad_s > 0.0 && !isnan(max_rate_rad_s)) || throw(ArgumentError("PanelHinge max_rate_rad_s must be positive, got $(max_rate_rad_s)."))
        (isfinite(angle_rad) && min_angle_rad <= angle_rad <= max_angle_rad) ||
            throw(ArgumentError("PanelHinge angle_rad must be finite and within its limits, got $(angle_rad)."))
        return new(a / n, Float64(min_angle_rad), Float64(max_angle_rad), Float64(max_rate_rad_s), Float64(angle_rad), NaN,
                   _HingeRunState((angle_rad, NaN)), nothing)
    end
end

PanelHinge(; axis, min_angle_rad::Real=-Inf, max_angle_rad::Real=Inf, max_rate_rad_s::Real=Inf, angle_rad::Real=0.0) =
    PanelHinge(axis, min_angle_rad, max_angle_rad, max_rate_rad_s, angle_rad)

# ---- Methods of the actuator contract --------------------------------------

actuator_limits(hinge::PanelHinge) = (hinge.min_angle_rad, hinge.max_angle_rad)

"""
    actuator_wrench(hinge::PanelHinge, level, ctx) -> (force_body, torque_body)

None: a hinge moves geometry. Its level is applied at control ticks.
"""
@inline actuator_wrench(::PanelHinge, ::Real, ::ActuatorContext) =
    (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))

"""
    hinge_step(hinge::PanelHinge, level, t) -> Float64

The angle the hinge reaches at a control tick at time `t` when commanded to
`level`: the level within the angle limits, then within `max_rate_rad_s`
times the time since the hinge's previous tick. A rate-limited hinge holds
its angle on its first tick, which only records the time.
"""
hinge_step(hinge::PanelHinge, level::Real, t::Real)::Float64 =
    _rate_limited_angle(actuator_saturate(hinge, level), hinge.angle_rad, hinge.last_update_s, hinge.max_rate_rad_s, t)

# The angle reached from `angle` toward `target` (already within the angle
# limits) at time `t`: within `max_rate` times the time since `last_update_s`,
# and no move on the first tick of a rate-limited joint (NaN `last_update_s`).
@inline function _rate_limited_angle(target, angle::Float64, last_update_s::Float64, max_rate::Float64, t::Real)::Float64
    if isfinite(max_rate)
        step = isfinite(last_update_s) ? max_rate * max(0.0, Float64(t) - last_update_s) : 0.0
        target = clamp(target, angle - step, angle + step)
    end
    return target
end
