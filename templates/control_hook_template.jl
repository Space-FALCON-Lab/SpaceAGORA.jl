# A control effector in the feedback profile ("Writing a GNC Algorithm" in the
# docs): the tick runs the control law once per control period and holds the
# command; the output hooks return the held command at every evaluation of the
# equations of motion and never run the law. Replace `control_law` with yours.
module ExampleControlHookTemplate

using LinearAlgebra
using StaticArrays
using SpaceAGORA

Base.@kwdef mutable struct ExampleControlHook <: SpaceAGORA.AbstractControlEffectorModel
    spacecraft_index::Int = 1           # the spacecraft it acts on, its position in the run
    period_s::Float64 = 1.0             # the tick period the law assumes, s
    max_force_n::Float64 = 1.0          # force limit on each inertial axis, N
    isp_s::Float64 = 300.0              # specific impulse of the thrust behind the force, s
    # Run state: the held command and the time of the tick that set it.
    force_ii::SVector{3, Float64} = zero(SVector{3, Float64})
    torque_body::SVector{3, Float64} = zero(SVector{3, Float64})
    tick_s::Float64 = NaN
end

# Declarations, each with a default the engine uses when you leave it out.
SpaceAGORA.spacecraft_scope(model::ExampleControlHook) = BoundSpacecraft((model.spacecraft_index,))
SpaceAGORA.gnc_period_s(model::ExampleControlHook) = model.period_s
SpaceAGORA.skip_inactive(::ExampleControlHook) = true
SpaceAGORA.held_command(model::ExampleControlHook, spacecraft_index) =
    (force_ii=model.force_ii, torque_body=model.torque_body)
SpaceAGORA.held_command_time(model::ExampleControlHook, spacecraft_index) = model.tick_s

# Called instead of the tick while the spacecraft is inactive: the off command.
function SpaceAGORA.release_command!(model::ExampleControlHook, t, spacecraft_index)
    model.force_ii = zero(SVector{3, Float64})
    model.torque_body = zero(SVector{3, Float64})
    model.tick_s = t
    return nothing
end

# Your control law: the command from the spacecraft's state at the tick, as an
# inertial force (N) and a bus-frame torque (N·m). This one commands nothing.
control_law(model::ExampleControlHook, state, t) = (zero(SVector{3, Float64}), zero(SVector{3, Float64}))

# The tick, once per control period: run the law, limit it, hold the command.
# `u.sc[spacecraft_index]` has `pos` (m), `vel` (m/s) and `mass` (kg), and with
# attitude propagation `q` and `ω`.
function SpaceAGORA.update_control!(model::ExampleControlHook, u, p, t, spacecraft_index)
    force_ii, torque_body = control_law(model, u.sc[spacecraft_index], t)
    model.force_ii = clamp.(force_ii, -model.max_force_n, model.max_force_n)
    model.torque_body = torque_body
    model.tick_s = t
    return nothing
end

# The held command: inertial force (N) and bus-frame torque (N·m). `x` is the
# spacecraft's state sample; its position in the run is `x.spacecraft_index`.
# The output hooks are called for every spacecraft, whatever the scope.
function SpaceAGORA.control_wrench(model::ExampleControlHook, x::SpaceAGORA.StateSample, p, t)
    x.spacecraft_index == model.spacecraft_index ||
        return zero(SVector{3, Float64}), zero(SVector{3, Float64})
    return model.force_ii, model.torque_body
end

# The propellant used in kg/s, positive while propellant is used. An effector
# that uses no propellant can leave this method out.
function SpaceAGORA.control_propellant_rate(model::ExampleControlHook, x::SpaceAGORA.StateSample, p, t)
    x.spacecraft_index == model.spacecraft_index || return 0.0
    return propellant_mass_rate(norm(model.force_ii), model.isp_s)
end

# At each run's start: fresh run state, so that two runs of one configuration
# with `isolate_state=false` start equal.
function SpaceAGORA.initialize_control!(model::ExampleControlHook, u, p, t)
    model.force_ii = zero(SVector{3, Float64})
    model.torque_body = zero(SVector{3, Float64})
    model.tick_s = NaN
    return nothing
end

# At a checkpoint: the run state as plain data, put back when a run resumes from
# the checkpoint file (after the reset above), so the resumed run holds the same
# command until its next tick as the uninterrupted run. The conformance check
# `checkpoint_round_trip` tests this; it is required from SpaceAGORA 0.5.0.
SpaceAGORA.checkpoint_state(model::ExampleControlHook) =
    (force_ii=model.force_ii, torque_body=model.torque_body, tick_s=model.tick_s)
function SpaceAGORA.restore_checkpoint_state!(model::ExampleControlHook, data)
    model.force_ii, model.torque_body, model.tick_s = data.force_ii, data.torque_body, data.tick_s
    return nothing
end

# Declared on a spacecraft (`SpacecraftModel(...; control=ControlModel(...))`):
# a copy for that spacecraft that shares the configuration, with fresh run state.
SpaceAGORA.bind_spacecraft(model::ExampleControlHook, spacecraft_index::Int) = ExampleControlHook(;
    spacecraft_index, period_s=model.period_s, max_force_n=model.max_force_n, isp_s=model.isp_s)

end # module ExampleControlHookTemplate
