module ExampleActuatorTemplate

# A new actuator kind and a controller that commands it. The kind says what the
# hardware does at a commanded level; the engine derives the force, torque and
# propellant use on the vehicle from the actuators mounted on it.

using LinearAlgebra
using StaticArrays
using SpaceAGORA

# 1. The hardware record: one cold-gas nozzle pushing along a fixed axis of
#    the link it is mounted on. The level is its thrust in newtons.
struct ExampleColdGasNozzle <: SpaceAGORA.AbstractActuator
    axis::SVector{3, Float64}          # thrust direction in the link frame (unit vector)
    location::SVector{3, Float64}      # nozzle position in the link frame, m
    max_thrust_n::Float64
    isp_s::Float64
end

# 2. Required: the level bounds and the wrench at a (saturated) level, both in
#    the spacecraft body frame, the torque about the body reference point.
SpaceAGORA.actuator_limits(n::ExampleColdGasNozzle) = (0.0, n.max_thrust_n)

function SpaceAGORA.actuator_wrench(n::ExampleColdGasNozzle, level, ctx::SpaceAGORA.ActuatorContext)
    force_body = ctx.link_to_body * (level * n.axis)
    arm_body = ctx.link_to_body * n.location + ctx.link_r_body
    return force_body, cross(arm_body, force_body)
end

# 3. Optional: propellant use (kg/s, positive) and power (W) at a level.
SpaceAGORA.actuator_consumption(n::ExampleColdGasNozzle, level, ctx) =
    (mass_rate_kg_s = level / (n.isp_s * 9.80665), power_w = 0.0)

# 4. A controller in `:mounted` actuation: it keeps one level per mounted
#    actuator of its spacecraft, updated at each control tick, and starts each
#    run from `initial_levels`. It writes no force, mass-flow or firing-level
#    hook; the engine derives them.
mutable struct ExampleNozzleController <: SpaceAGORA.AbstractControlEffectorModel
    spacecraft_index::Int
    levels::Vector{Float64}
    initial_levels::Vector{Float64}     # the levels each run starts from
end

# Configured with one vector: the levels each run starts from, copied into both
# fields, so the levels a run writes never change the configured start.
ExampleNozzleController(spacecraft_index::Integer, levels::AbstractVector) =
    ExampleNozzleController(spacecraft_index, collect(Float64, levels), collect(Float64, levels))

SpaceAGORA.actuation_mode(::ExampleNozzleController) = :mounted
SpaceAGORA.actuator_levels(c::ExampleNozzleController, spacecraft_index::Int) =
    spacecraft_index == c.spacecraft_index ? c.levels : nothing

function SpaceAGORA.update_control!(c::ExampleNozzleController, u, p, t::Float64, spacecraft_index::Int)
    # Compute and store the levels here (they are held until the next tick).
    return nothing
end

# At each run's start: the levels go back to `initial_levels`, so two runs of one
# configuration with `isolate_state=false` start equal (the conformance check
# `initialize_resets`).
function SpaceAGORA.initialize_control!(c::ExampleNozzleController, u, p, t)
    copy!(c.levels, c.initial_levels)
    return nothing
end

# The held levels are run state: saved at a checkpoint as plain data and put back
# when a run resumes from the checkpoint file, so the resumed run holds the same
# levels until its next tick. The conformance check `checkpoint_round_trip` tests
# this; it is required from SpaceAGORA 0.5.0.
SpaceAGORA.checkpoint_state(c::ExampleNozzleController) = (levels=copy(c.levels),)
function SpaceAGORA.restore_checkpoint_state!(c::ExampleNozzleController, data)
    copy!(c.levels, data.levels)
    return nothing
end

# Mounting, for a spacecraft `sc` and one of its links:
#     add_actuator!(sc, sc.root, ExampleColdGasNozzle(SVector(1.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 0.5, 60.0))
# and the controller goes in ControlModel(control_effectors=(ExampleNozzleController(1, [0.1]),), ...);
# every run starts at level 0.1 until the controller's first tick.

end # module ExampleActuatorTemplate
