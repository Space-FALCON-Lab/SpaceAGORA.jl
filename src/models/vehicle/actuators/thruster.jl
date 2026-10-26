# The thruster actuator kind. Included by `Components`
# (src/models/vehicle/spacecraft/components.jl), which owns and exports the
# record, so serialized spacecraft keep their type identity; the actuator
# contract it implements is in `actuators.jl`.

# ---- Record ---------------------------------------------------------------

# The firing state of a thruster with firing dynamics: its ramp factor and thrust.
const _FiringRunState = NamedTuple{(:κ, :thrust), Tuple{Float64, Float64}}

"""
    Thruster(; max_thrust, location, direction, Isp, ...)

One jet mounted on a spacecraft link. `direction` is the direction of the
thrust it puts on the vehicle (not the direction of its exhaust), in the link
frame; `location` is the nozzle position in the link frame, relative to the
link's centre of mass, metres; `max_thrust` is its rating in newtons and
`Isp` its specific impulse in seconds. The remaining fields describe its
firing dynamics and hold its state during a run.

A `Thruster` is an [`AbstractActuator`](@ref) whose level is a thrust in
newtons, from zero to `max_thrust`. Mount it with `add_actuator!`.

With `firing_dynamics=true`, a mounted thruster does not deliver its command
at once: at each control tick of the effector commanding it, the commanded
thrust becomes an on-time within the control period, subject to
`min_firing_time` and the Schmitt trigger (`level_on`, `level_off`), and the
thrust ramps with `cutoff_frequency`; the vehicle feels the resulting average
thrust until the next tick. Off by default, the thruster acts at its
saturated command.
"""
@kwdef mutable struct Thruster <: AbstractActuator
    max_thrust::Float64 = 1.0 # Maximum thrust, N
    location::MVector{3, Float64} = MVector{3, Float64}(zeros(3)) # Location in the link frame, relative to the CoM of the link, m
    direction::MVector{3, Float64} = MVector{3, Float64}(zeros(3)) # Unit vector direction of thrust in the link frame, n/d
    Isp::Float64 = 0.0 # Specific impulse of the thruster, s
    cutoff_frequency::Float64 = 1.0 # Governing parameter for thruster ramp-up/ramp-down, rad/s.
    min_firing_time::Float64 = 0.0 # Minimum time for the thruster to be firing, s.
    level_on::Float64 = 0.75 # Upper threshold for the Schmitt trigger, nd
    level_off::Float64 = 0.25 # Lower threshold for the Schmitt trigger, nd
    κ::Float64 = 0.0 # Thrust factor, nd.
    thrust::Float64 = 0.0 # Current thrust magnitude, to be updated during simulation.
    stop_firing_time::Float64 = 0.0 # Time at which the thruster should stop firing, s.
    firing_dynamics::Bool = false # In :mounted actuation, deliver the command through the minimum on-time, Schmitt trigger and thrust ramp above
    gimballed::Bool = false # A ThrusterGimbal turns this thruster's thrust (set when one is mounted)
    gimbal_rotation::SMatrix{3, 3, Float64, 9} = SMatrix{3, 3, Float64, 9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0) # Its gimbals' current rotation, link frame
    run_start::_FiringRunState = _FiringRunState((κ, thrust)) # The firing state a run starts from, first the configured one
    run_written::Union{Nothing, _FiringRunState} = nothing # The run's last tick's write of it (`ActuatorMounting._start_run_state!`)
end

# ---- State ----------------------------------------------------------------
#
# The mutable fields `thrust`, `κ` and `stop_firing_time` above hold the firing
# state; `run_simulation(isolate_state=true)` copies them per run. With
# `isolate_state=false` each run starts `κ` and `thrust` from the previous run's
# start, unless they were changed after that run's last control tick
# (`run_start`, `run_written`; `ActuatorMounting._start_run_state!`). The firing
# dynamics that update them (Schmitt trigger, impulse integration, ramp
# factor) are in `thruster_hooks.jl`, which needs the link types. A checkpoint
# saves `κ` and `thrust` of a mounted thruster with firing dynamics, and a
# resumed run restores them (the engine's callback snapshot,
# `resume_checkpoint.jl`); the run never writes `stop_firing_time`.

# ---- Methods of the actuator contract --------------------------------------

actuator_limits(thruster::Thruster) = (0.0, Float64(thruster.max_thrust))

"""
    thrust_direction(thruster::Thruster) -> SVector{3}

The direction of the thruster's thrust in its link frame: its `direction`,
turned by its gimbals' current angles when a `ThrusterGimbal` is mounted on
it. Not normalized.
"""
@inline thrust_direction(thruster::Thruster)::SVector{3, Float64} =
    thruster.gimballed ? thruster.gimbal_rotation * SVector{3, Float64}(thruster.direction) : SVector{3, Float64}(thruster.direction)

"""
    actuator_wrench(thruster::Thruster, level, ctx) -> (force_body, torque_body)

A thrust of `level` newtons along the thruster's `direction`, applied at its
`location`: the force in the body frame and its torque about the body
reference point. A zero or non-finite direction gives no wrench.
"""
@inline function actuator_wrench(thruster::Thruster, level::Real, ctx::ActuatorContext)
    direction = thrust_direction(thruster)
    n = norm(direction)
    (isfinite(n) && n > 0.0) || return (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
    force_body = ctx.link_to_body * ((Float64(level) / n) * direction)
    arm_body = ctx.link_to_body * SVector{3, Float64}(thruster.location) + ctx.link_r_body
    return force_body, cross(arm_body, force_body)
end

"""
    actuator_consumption(thruster::Thruster, level, ctx) -> (mass_rate_kg_s, power_w)

Propellant use at a thrust of `level` newtons, `level / (Isp g0)`; zero when
`Isp` is not positive. No electrical power.
"""
@inline function actuator_consumption(thruster::Thruster, level::Real, ::ActuatorContext)
    isp = Float64(thruster.Isp)
    thrust = Float64(level)
    rate = (isfinite(isp) && isp > 0.0 && isfinite(thrust) && thrust > 0.0) ?
        propellant_mass_rate(thrust, isp) : 0.0
    return (mass_rate_kg_s = rate, power_w = 0.0)
end

# ---- Builders ---------------------------------------------------------------

"""
    six_axis_thrusters(model::SixAxisThrusterModel) -> Vector{Thruster}

The six jets of a six-axis thruster model as mountable thrusters, in the
model's column order: each one's thrust direction, location, rating, specific
impulse and minimum firing time. Mount them on the vehicle with
`add_actuator!` to drive them in `:mounted` actuation.
"""
function six_axis_thrusters(model::SixAxisThrusterModel)::Vector{Thruster}
    return [Thruster(
        max_thrust=model.max_thrust_n[j],
        location=MVector{3, Float64}(model.locations_body[:, j]),
        direction=MVector{3, Float64}(model.directions_body[:, j]),
        Isp=model.isp_s[j],
        min_firing_time=model.min_firing_time_s,
    ) for j in 1:6]
end
