"""
Mounting actuators on spacecraft links and evaluating the mounted set.

The mounting order is the order every per-actuator vector follows (levels,
allocation columns): the built-in kinds first, in a fixed order, then the
actuators of user-defined kinds.

1. every thruster, link by link in `model.links` order and within a link in
   `link.thrusters` order (the order of the `thruster_level` result columns
   and of the viewer's plumes);
2. every mounted reaction wheel, in mounting order (wheels mount on the root
   link);
3. every torque rod, link by link in `link.magnetorquers` order;
4. every laser terminal, link by link in `link.laser_terminals` order;
5. every hinge, link by link (one per hinged link);
6. every thruster gimbal axis, link by link in `link.gimbals` order;
7. every control moment gyro, in mounting order (CMGs mount on the root
   link);
8. every actuator of a user-defined kind, link by link in
   `link.custom_actuators` order.

Built-in kinds added later take their place before the user-defined kinds.
"""
module ActuatorMounting

using LinearAlgebra
using StaticArrays
using ..AbstractTypes: AbstractActuator
using ..QuaternionMath: rot
using ..Actuators
using ..Components: Thruster, ReactionWheel, wheel_momentum_saturate, wheel_power, Magnetorquer, LaserTerminal, PanelHinge, hinge_step
using ..Components: ThrusterGimbal, gimbal_step, _gimbal_rotation
using ..Components: ControlMomentGyro, cmg_torque_axis, cmg_angle_saturate
using ..SpacecraftModels: Link, SpacecraftModel
using ..Assembly: add_thruster!

export add_actuator!, mounted_actuators, mounted_actuator_count, mounted_thruster_count, mounted_wheel_count, mounted_magnetorquer_count, mounted_laser_count, mounted_hinge_count, apply_hinge_levels!, move_hinge!
export mounted_gimbal_count, apply_gimbal_levels!, move_gimbal!
export mounted_cmg_count, cmg_steering
export actuator_context, allocation_matrix, mounted_wrench, mounted_thruster_fractions
export mounted_requirements, applied_actuator_levels, mounted_power

const _ZERO3 = SVector{3, Float64}(0.0, 0.0, 0.0)
const _IDENTITY3 = SMatrix{3, 3, Float64, 9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
const _IDENTITY_Q = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
const _NAN3 = SVector{3, Float64}(NaN, NaN, NaN)

# Actuators mount only on the spacecraft's own links: the dynamics traverse
# `model.links`, so an actuator on any other link would be silently ignored.
@inline function _check_own_link(model::SpacecraftModel, link::Link)
    any(l -> l === link, model.links) || throw(ArgumentError(
        "add_actuator!: the link is not one of spacecraft $(model.id)'s links; mount actuators on a link in `model.links`."))
    return nothing
end

"""
    add_actuator!(model::SpacecraftModel, link::Link, actuator) -> actuator

Mount `actuator` on `link`, which must be one of `model.links` (any other
link is refused, since the dynamics would never see it). A [`Thruster`](@ref) joins
`link.thrusters` (its direction is normalized and the link's thruster torque
Jacobian gains its column). A [`ReactionWheel`](@ref) mounts on the root link
only and joins its wheel assembly, whose momentum the engine integrates. A
[`Magnetorquer`](@ref) joins `link.magnetorquers` and a [`LaserTerminal`](@ref)
`link.laser_terminals`. A [`PanelHinge`](@ref) mounts on a non-root link, one per
link, and sets that link's attitude to its current angle. A [`ThrusterGimbal`](@ref)
mounts on the link of the thruster it turns and joins `link.gimbals`; the
thruster's thrust turns by its current angle. A [`ControlMomentGyro`](@ref)
mounts on the root link only and joins its wheel assembly's `cmgs`; its
gimbal angle becomes a state the engine integrates. An actuator of any other [`AbstractActuator`](@ref) kind joins
`link.custom_actuators`. Returns the actuator.
"""
function add_actuator!(model::SpacecraftModel, link::Link, thruster::Thruster)
    _check_own_link(model, link)
    add_thruster!(model, link, thruster)
    return thruster
end

function add_actuator!(model::SpacecraftModel, link::Link, wheel::ReactionWheel)
    link === model.root || throw(ArgumentError(
        "Reaction wheels mount on the spacecraft's root link, whose wheel assembly the engine integrates."))
    _check_own_link(model, link)
    assembly = link.rw_assembly
    assembly.J_rw = hcat(assembly.J_rw, Vector{Float64}(wheel.spin_axis))
    assembly.J_rw_pinv = Matrix{Float64}(pinv(assembly.J_rw))
    assembly.n_wheels += 1
    push!(assembly.h_wheels, 0.0)
    push!(assembly.h_dot_wheels, 0.0)
    push!(assembly.wheels, wheel)
    push!(assembly.mounted_rates, 0.0)
    return wheel
end

function add_actuator!(model::SpacecraftModel, link::Link, rod::Magnetorquer)
    _check_own_link(model, link)
    push!(link.magnetorquers, rod)
    return rod
end

function add_actuator!(model::SpacecraftModel, link::Link, terminal::LaserTerminal)
    _check_own_link(model, link)
    push!(link.laser_terminals, terminal)
    return terminal
end

function add_actuator!(model::SpacecraftModel, link::Link, hinge::PanelHinge)
    _check_own_link(model, link)
    link.root && throw(ArgumentError("A hinge turns a non-root link; the root link carries the body frame."))
    isempty(link.hinges) || throw(ArgumentError("This link already has a hinge; mount at most one per link."))
    push!(link.hinges, hinge)
    _set_hinge_attitude!(link, hinge, hinge.angle_rad)
    return hinge
end

function add_actuator!(model::SpacecraftModel, link::Link, gimbal::ThrusterGimbal)
    _check_own_link(model, link)
    any(t -> t === gimbal.thruster, link.thrusters) || throw(ArgumentError(
        "A gimbal mounts on the link of the thruster it turns; mount the thruster on this link first."))
    any(g -> g === gimbal, link.gimbals) && throw(ArgumentError("This gimbal is already mounted."))
    push!(link.gimbals, gimbal)
    _update_gimbal_rotation!(link, gimbal.thruster)
    return gimbal
end

function add_actuator!(model::SpacecraftModel, link::Link, cmg::ControlMomentGyro)
    link === model.root || throw(ArgumentError(
        "Control moment gyros mount on the spacecraft's root link, whose gimbal angles the engine integrates."))
    _check_own_link(model, link)
    assembly = link.rw_assembly
    any(c -> c === cmg, assembly.cmgs) && throw(ArgumentError("This control moment gyro is already mounted."))
    push!(assembly.cmgs, cmg)
    push!(assembly.cmg_rates, 0.0)
    return cmg
end

# The thruster's thrust turned by its gimbals on `link` at their current angles.
function _update_gimbal_rotation!(link::Link, thruster::Thruster)
    thruster.gimbal_rotation = _gimbal_rotation(thruster, link.gimbals)
    thruster.gimballed = true
    return nothing
end

# The link's attitude as the rotation by `angle` about the hinge axis: the
# arithmetic of `Kinematics.rotate_link(link, axis, angle)`.
@inline function _set_hinge_attitude!(link::Link, hinge::PanelHinge, angle::Float64)
    link.q .= SVector{4, Float64}(hinge.axis[1] * sin(angle / 2), hinge.axis[2] * sin(angle / 2),
                                  hinge.axis[3] * sin(angle / 2), cos(angle / 2))
    return nothing
end

function add_actuator!(model::SpacecraftModel, link::Link, actuator::AbstractActuator)
    _check_own_link(model, link)
    push!(link.custom_actuators, actuator)
    return actuator
end

"""
    mounted_thruster_count(model::SpacecraftModel) -> Int

The number of thrusters mounted on `model`'s links.
"""
function mounted_thruster_count(model::SpacecraftModel)::Int
    n = 0
    for link in model.links
        n += length(link.thrusters)
    end
    return n
end

"""
    mounted_wheel_count(model::SpacecraftModel) -> Int

The number of reaction wheels mounted with [`add_actuator!`](@ref). Wheels
given to the root link's assembly directly (its `J_rw` keyword) are not
actuators of the mounted set.
"""
mounted_wheel_count(model::SpacecraftModel)::Int = length(model.root.rw_assembly.wheels)

"""
    mounted_magnetorquer_count(model::SpacecraftModel) -> Int

The number of torque rods mounted on `model`'s links.
"""
function mounted_magnetorquer_count(model::SpacecraftModel)::Int
    n = 0
    for link in model.links
        n += length(link.magnetorquers)
    end
    return n
end

"""
    mounted_laser_count(model::SpacecraftModel) -> Int

The number of laser terminals mounted on `model`'s links.
"""
function mounted_laser_count(model::SpacecraftModel)::Int
    n = 0
    for link in model.links
        n += length(link.laser_terminals)
    end
    return n
end

"""
    mounted_hinge_count(model::SpacecraftModel) -> Int

The number of hinges mounted on `model`'s links.
"""
function mounted_hinge_count(model::SpacecraftModel)::Int
    n = 0
    for link in model.links
        n += length(link.hinges)
    end
    return n
end

"""
    apply_hinge_levels!(model::SpacecraftModel, levels, t) -> model

Move every hinge mounted on `model` toward its level in `levels` (mounting
order) at a control tick at time `t`, within its angle and rate limits (see
`hinge_step`), and set its link's attitude to the angle reached. The engine
calls this after the control tick of an effector in `:mounted` actuation.
"""
function apply_hinge_levels!(model::SpacecraftModel, levels::AbstractVector, t::Real)
    mounted_hinge_count(model) == 0 && return model
    length(levels) == mounted_actuator_count(model) || throw(ArgumentError(
        "Mounted actuation: $(length(levels)) actuator levels for the $(mounted_actuator_count(model)) actuators mounted on spacecraft $(model.id)."))
    k = mounted_thruster_count(model) + mounted_wheel_count(model) + mounted_magnetorquer_count(model) + mounted_laser_count(model)
    for link in model.links, hinge in link.hinges
        k += 1
        move_hinge!(link, hinge, levels[k], t)
    end
    return model
end

"""
    move_hinge!(link::Link, hinge::PanelHinge, level, t) -> Float64

Move `hinge`, mounted on `link`, toward the angle `level` at a control tick at
time `t`, within its angle and rate limits (see `hinge_step`), set the link's
attitude to the angle reached and record the tick. Returns the angle reached.
"""
function move_hinge!(link::Link, hinge::PanelHinge, level::Real, t::Real)::Float64
    angle = hinge_step(hinge, level, t)
    _set_hinge_attitude!(link, hinge, angle)
    hinge.angle_rad = angle
    hinge.last_update_s = Float64(t)
    return angle
end

"""
    mounted_gimbal_count(model::SpacecraftModel) -> Int

The number of thruster gimbal axes mounted on `model`'s links.
"""
function mounted_gimbal_count(model::SpacecraftModel)::Int
    n = 0
    for link in model.links
        n += length(link.gimbals)
    end
    return n
end

"""
    apply_gimbal_levels!(model::SpacecraftModel, levels, t) -> model

Move every thruster gimbal mounted on `model` toward its level in `levels`
(mounting order) at a control tick at time `t`, within its angle and rate
limits (see `gimbal_step`), and turn its thruster's thrust to the angles
reached. The engine calls this after the control tick of an effector in
`:mounted` actuation.
"""
function apply_gimbal_levels!(model::SpacecraftModel, levels::AbstractVector, t::Real)
    mounted_gimbal_count(model) == 0 && return model
    length(levels) == mounted_actuator_count(model) || throw(ArgumentError(
        "Mounted actuation: $(length(levels)) actuator levels for the $(mounted_actuator_count(model)) actuators mounted on spacecraft $(model.id)."))
    k = mounted_thruster_count(model) + mounted_wheel_count(model) + mounted_magnetorquer_count(model) +
        mounted_laser_count(model) + mounted_hinge_count(model)
    for link in model.links, gimbal in link.gimbals
        k += 1
        move_gimbal!(link, gimbal, levels[k], t)
    end
    return model
end

"""
    move_gimbal!(link::Link, gimbal::ThrusterGimbal, level, t) -> Float64

Move `gimbal`, mounted on `link`, toward the angle `level` at a control tick
at time `t`, within its angle and rate limits, turn its thruster's thrust and
record the tick. Returns the angle reached.
"""
function move_gimbal!(link::Link, gimbal::ThrusterGimbal, level::Real, t::Real)::Float64
    angle = gimbal_step(gimbal, level, t)
    gimbal.angle_rad = angle
    gimbal.last_update_s = Float64(t)
    _update_gimbal_rotation!(link, gimbal.thruster)
    return angle
end

"""
    mounted_cmg_count(model::SpacecraftModel) -> Int

The number of control moment gyros mounted on `model` (on its root link).
"""
mounted_cmg_count(model::SpacecraftModel)::Int = length(model.root.rw_assembly.cmgs)

# The number of mounted actuators before the first CMG in the mounting order.
_cmg_level_offset(model::SpacecraftModel)::Int =
    mounted_thruster_count(model) + mounted_wheel_count(model) + mounted_magnetorquer_count(model) +
    mounted_laser_count(model) + mounted_hinge_count(model) + mounted_gimbal_count(model)

"""
    cmg_steering(model::SpacecraftModel, angles, torque_body; λ0=0.01, μ=10.0) -> Vector{Float64}

Gimbal rates for `model`'s mounted control moment gyros, at gimbal angles
`angles` (mounting order), that put `torque_body` (N m, body frame) on the
vehicle: the singularity-robust inverse of the CMGs' torque Jacobian
(Bedrossian): with `A = [g_j × h_j]` normalized by the largest rotor
momentum, `δ̇ = -Aᵀ (A Aᵀ + λ I)⁻¹ τ` and `λ = λ0 exp(-μ det(A Aᵀ))`, which
stays bounded near a singular gimbal configuration at the cost of a torque
error (`λ0 = 0` gives the pseudo-inverse). A CMG at a gimbal stop whose rate
would drive it further is left out and the others are solved again, so it is
commanded zero; rates beyond a CMG's limit are then scaled down together,
keeping the torque's direction.
"""
function cmg_steering(model::SpacecraftModel, angles::AbstractVector, torque_body::AbstractVector; λ0::Real=0.01, μ::Real=10.0)::Vector{Float64}
    cmgs = model.root.rw_assembly.cmgs
    n = length(cmgs)
    length(angles) == n || throw(ArgumentError("cmg_steering: $(length(angles)) gimbal angles for $(n) mounted CMGs."))
    n == 0 && return Float64[]
    href = maximum(c -> c.rotor_momentum_nms, cmgs)
    stopped = falses(n)
    rates = zeros(Float64, n)
    for _ in 0:n
        A = zeros(Float64, 3, n)
        for j in 1:n
            stopped[j] || (A[:, j] .= cmg_torque_axis(cmgs[j], angles[j]) ./ href)
        end
        M = A * transpose(A)
        λ = Float64(λ0) * exp(-Float64(μ) * det(M))
        # pinv: the inverse for λ > 0, the least-squares minimum-norm solution for λ = 0 at a singular A.
        rates = -(transpose(A) * (pinv(M + λ * Matrix{Float64}(I, 3, 3)) * Vector{Float64}(torque_body))) ./ href
        # A CMG at a stop that the solution would drive further is left out.
        blocked = [j for j in 1:n if !stopped[j] && rates[j] != 0.0 && cmg_angle_saturate(cmgs[j], rates[j], angles[j]) == 0.0]
        isempty(blocked) && break
        stopped[blocked] .= true
    end
    scale = 1.0
    for j in 1:n
        limit = cmgs[j].max_gimbal_rate_rad_s
        isfinite(limit) && (scale = max(scale, abs(rates[j]) / limit))
    end
    return rates ./ scale
end

"""
    mounted_actuator_count(model::SpacecraftModel) -> Int

The number of actuators mounted on `model`, of every kind.
"""
function mounted_actuator_count(model::SpacecraftModel)::Int
    n = mounted_thruster_count(model) + mounted_wheel_count(model) + mounted_magnetorquer_count(model) + mounted_laser_count(model) +
        mounted_hinge_count(model) + mounted_gimbal_count(model) + mounted_cmg_count(model)
    for link in model.links
        n += length(link.custom_actuators)
    end
    return n
end

"""
    mounted_actuators(model::SpacecraftModel) -> Vector{NamedTuple{(:link, :actuator)}}

Every actuator mounted on `model`, with the link it is mounted on, in mounting
order.
"""
function mounted_actuators(model::SpacecraftModel)
    out = NamedTuple{(:link, :actuator), Tuple{Link, AbstractActuator}}[]
    for link in model.links, thruster in link.thrusters
        push!(out, (link=link, actuator=thruster))
    end
    for wheel in model.root.rw_assembly.wheels
        push!(out, (link=model.root, actuator=wheel))
    end
    for link in model.links, rod in link.magnetorquers
        push!(out, (link=link, actuator=rod))
    end
    for link in model.links, terminal in link.laser_terminals
        push!(out, (link=link, actuator=terminal))
    end
    for link in model.links, hinge in link.hinges
        push!(out, (link=link, actuator=hinge))
    end
    for link in model.links, gimbal in link.gimbals
        push!(out, (link=link, actuator=gimbal))
    end
    for cmg in model.root.rw_assembly.cmgs
        push!(out, (link=model.root, actuator=cmg))
    end
    for link in model.links, actuator in link.custom_actuators
        push!(out, (link=link, actuator=actuator))
    end
    return out
end

"""
    actuator_context(link::Link; mass_kg=NaN, q_ib=identity, ω_body=zero, magnetic_field_body=NaN, environment=nothing) -> ActuatorContext

The evaluation context of an actuator mounted on `link`: the link's configured
pose in the body frame (its `r` and, for a non-root link, its attitude `q`),
the vehicle's mass, attitude, rate and local magnetic field, and the sampled
environment the mounted actuators asked for.
"""
actuator_context(link::Link; mass_kg::Real=NaN, q_ib=_IDENTITY_Q, ω_body=_ZERO3, magnetic_field_body=_NAN3, environment=nothing) =
    actuator_context(link, Float64(mass_kg), SVector{4, Float64}(q_ib), SVector{3, Float64}(ω_body), SVector{3, Float64}(magnetic_field_body), environment)

@inline function actuator_context(link::Link, mass_kg::Float64, q_ib::SVector{4, Float64}, ω_body::SVector{3, Float64},
                                  magnetic_field_body::SVector{3, Float64}=_NAN3, environment=nothing)
    # Same link-to-body map as the thruster torque Jacobian (Kinematics.rotate_to_body).
    link_to_body = link.root ? _IDENTITY3 : SMatrix{3, 3, Float64, 9}(rot(SVector{4, Float64}(link.q))')
    return ActuatorContext(SVector{3, Float64}(link.r), link_to_body, mass_kg, q_ib, ω_body, magnetic_field_body, environment)
end

"""
    mounted_requirements(model::SpacecraftModel) -> ActuatorRequirements

The union of [`actuator_requirements`](@ref) over every actuator mounted on
`model`: what the engine samples once per evaluation for them.
"""
function mounted_requirements(model::SpacecraftModel)::ActuatorRequirements
    req = mounted_magnetorquer_count(model) > 0 ? ActuatorRequirements(magnetic_field=true) : ActuatorRequirements()
    for link in model.links, actuator in link.custom_actuators
        req = req | (actuator_requirements(actuator)::ActuatorRequirements)
    end
    return req
end

# The thrust the `k`-th mounted thruster delivers: its saturated command, or,
# with firing dynamics, the average thrust its last control tick scheduled
# (`thrust`) as its active faults leave it. Faults act on that thrust, not on
# the command the tick scheduled from, so they switch at their own times
# rather than at the next tick: a failed thruster delivers nothing, a stuck
# one its held level, a degraded one the fraction `value`, a biased one the
# offset `value` within its limits. (Without firing dynamics the faults act
# on the level, through `ActuatorFaults.faulted_levels`.)
@inline _delivered_thrust(model::SpacecraftModel, k::Int, thruster::Thruster, level)::Float64 =
    thruster.firing_dynamics ? _firing_thrust(model, k, thruster) : actuator_saturate(thruster, level)

function _firing_thrust(model::SpacecraftModel, k::Int, thruster::Thruster)::Float64
    thrust = Float64(thruster.thrust)
    for f in model.root.faults.faults
        (f.active && f.index == k) || continue
        if f.kind === :failed
            return 0.0
        elseif f.kind === :stuck
            return actuator_saturate(thruster, f.held_level)
        elseif f.kind === :degraded
            thrust *= f.value
        else   # :biased
            lo, hi = actuator_limits(thruster)
            thrust = clamp(thrust + f.value, Float64(lo), Float64(hi))
        end
    end
    return thrust
end

"""
    allocation_matrix(model::SpacecraftModel) -> Matrix{Float64}

The 6 × n allocation matrix of the actuators mounted on `model`, in mounting
order: column `k` stacks the force (rows 1 to 3, N) and the torque about the
body reference point (rows 4 to 6, N m) of actuator `k` at level one, both in
the body frame (see [`actuator_unit_wrench`](@ref)). Use it with
[`allocate_actuators`](@ref).
"""
function allocation_matrix(model::SpacecraftModel)::Matrix{Float64}
    mounted = mounted_actuators(model)
    B = zeros(Float64, 6, length(mounted))
    for (k, entry) in enumerate(mounted)
        force, torque = actuator_unit_wrench(entry.actuator, actuator_context(entry.link))
        B[1:3, k] .= force
        B[4:6, k] .= torque
    end
    return B
end

@inline function _add_actuator(actuator, x::Float64, ctx::ActuatorContext,
                               force::SVector{3, Float64}, torque::SVector{3, Float64}, mass_rate::Float64)
    f, tq = actuator_wrench(actuator, x, ctx)
    use = actuator_consumption(actuator, x, ctx)
    return force + SVector{3, Float64}(f), torque + SVector{3, Float64}(tq), mass_rate + Float64(use.mass_rate_kg_s)
end

# Function barrier for user-defined kinds held in an abstractly typed list.
@noinline _add_custom_actuator(actuator::AbstractActuator, level, ctx, force, torque, mass_rate) =
    _add_actuator(actuator, actuator_saturate(actuator, level), ctx, force, torque, mass_rate)

"""
    mounted_wrench(model::SpacecraftModel, levels; mass_kg=NaN, q_ib=identity, ω_body=zero, wheel_momentum=nothing,
                   magnetic_field_body=NaN, wheel_rates=nothing, environment=nothing)
        -> (force_body, torque_body, mass_rate_kg_s, wheel_torque_body)

The wrench of every actuator mounted on `model` at its level in `levels` (one
per mounted actuator, in mounting order), each level first saturated: the
total force and torque in the body frame, the torque about the body reference
point; the total propellant mass rate (kg/s, positive while consuming); and the
part of the torque the reaction wheels produce, which the engine integrates
into their momentum. `wheel_momentum` holds the root wheel assembly's stored
momenta (the integrated `h_wheels` state); with it, a wheel at its momentum
limit produces no torque that would store more. Torque rods act in
`magnetic_field_body`, the local field in the body frame (T). When
`wheel_rates` is given (one entry per mounted wheel), each mounted wheel's own
momentum rate, minus its saturated level, is subtracted into it: the engine
integrates the mounted wheels from these, wheel by wheel. `environment` is the
sampled environment the mounted actuators asked for, passed to each one's
context. A thruster with firing dynamics acts at the thrust its last control
tick scheduled rather than at its level.
"""
mounted_wrench(model::SpacecraftModel, levels::AbstractVector; mass_kg::Real=NaN, q_ib=_IDENTITY_Q, ω_body=_ZERO3, wheel_momentum=nothing,
               magnetic_field_body=_NAN3, wheel_rates=nothing, environment=nothing) =
    mounted_wrench(model, levels, Float64(mass_kg), SVector{4, Float64}(q_ib), SVector{3, Float64}(ω_body), wheel_momentum,
                   SVector{3, Float64}(magnetic_field_body), wheel_rates, environment)

function mounted_wrench(model::SpacecraftModel, levels::AbstractVector, mass_kg::Float64,
                        q_ib::SVector{4, Float64}, ω_body::SVector{3, Float64}, wheel_momentum=nothing,
                        magnetic_field_body::SVector{3, Float64}=_NAN3, wheel_rates=nothing, environment=nothing)
    n = mounted_actuator_count(model)
    length(levels) == n || throw(ArgumentError(
        "Mounted actuation: $(length(levels)) actuator levels for the $(n) actuators mounted on spacecraft $(model.id)."))
    force = _ZERO3
    torque = _ZERO3
    wheel_torque = _ZERO3
    mass_rate = 0.0
    k = 0
    @inbounds for link in model.links
        isempty(link.thrusters) && continue
        ctx = actuator_context(link, mass_kg, q_ib, ω_body, magnetic_field_body, environment)
        for thruster in link.thrusters
            k += 1
            force, torque, mass_rate = _add_actuator(thruster, _delivered_thrust(model, k, thruster, levels[k]), ctx, force, torque, mass_rate)
        end
    end
    assembly = model.root.rw_assembly
    if !isempty(assembly.wheels)
        ctx = actuator_context(model.root, mass_kg, q_ib, ω_body, magnetic_field_body, environment)
        first_mounted = assembly.n_wheels - length(assembly.wheels)   # wheels given directly come first
        @inbounds for (j, wheel) in enumerate(assembly.wheels)
            k += 1
            x = actuator_saturate(wheel, levels[k])
            if wheel_momentum !== nothing
                x = wheel_momentum_saturate(wheel, x, Float64(wheel_momentum[first_mounted + j]))
            end
            _, tq = actuator_wrench(wheel, x, ctx)
            torque = torque + tq
            wheel_torque = wheel_torque + tq
            wheel_rates === nothing || (wheel_rates[j] -= x)    # the wheel stores the opposite of what it gives the body
        end
    end
    @inbounds for link in model.links
        isempty(link.magnetorquers) && continue
        ctx = actuator_context(link, mass_kg, q_ib, ω_body, magnetic_field_body, environment)
        for rod in link.magnetorquers
            k += 1
            force, torque, mass_rate = _add_actuator(rod, actuator_saturate(rod, levels[k]), ctx, force, torque, mass_rate)
        end
    end
    @inbounds for link in model.links
        isempty(link.laser_terminals) && continue
        ctx = actuator_context(link, mass_kg, q_ib, ω_body, magnetic_field_body, environment)
        for terminal in link.laser_terminals
            k += 1
            force, torque, mass_rate = _add_actuator(terminal, actuator_saturate(terminal, levels[k]), ctx, force, torque, mass_rate)
        end
    end
    # Hinges and gimbals move geometry at control ticks; CMGs act through their
    # gimbal-angle state in the attitude equation. No wrench here.
    k += mounted_hinge_count(model) + mounted_gimbal_count(model) + mounted_cmg_count(model)
    @inbounds for link in model.links
        isempty(link.custom_actuators) && continue
        ctx = actuator_context(link, mass_kg, q_ib, ω_body, magnetic_field_body, environment)
        for actuator in link.custom_actuators
            k += 1
            # The list is abstractly typed: assert the result so the loop stays concretely typed.
            force, torque, mass_rate = _add_custom_actuator(actuator, levels[k], ctx, force, torque, mass_rate)::Tuple{SVector{3, Float64}, SVector{3, Float64}, Float64}
        end
    end
    return force, torque, mass_rate, wheel_torque
end

"""
    mounted_thruster_fractions(model::SpacecraftModel, levels) -> Vector{Float64}

The firing level, from 0 to 1, of each thruster mounted on `model`, from the
first entries of `levels` (mounting order): each saturated level as a fraction
of the thruster's rating.
"""
function mounted_thruster_fractions(model::SpacecraftModel, levels::AbstractVector)::Vector{Float64}
    out = Float64[]
    k = 0
    for link in model.links, thruster in link.thrusters
        k += 1
        x = k <= length(levels) ? _delivered_thrust(model, k, thruster, levels[k]) : 0.0
        push!(out, actuator_fraction(thruster, x))
    end
    return out
end

"""
    applied_actuator_levels(model::SpacecraftModel, levels; wheel_momentum=nothing, cmg_angles=nothing) -> Vector{Float64}

The level each mounted actuator actually applies, in mounting order: its
command saturated to its limits; a thruster with firing dynamics at the thrust
its last tick scheduled; a wheel at its momentum limit at zero when the
command would store more (given `wheel_momentum`); a hinge at its current
angle; a control moment gyro at zero at a gimbal stop when the command would
drive it further (given `cmg_angles`). These fill the `actuator_level` result
columns.
"""
function applied_actuator_levels(model::SpacecraftModel, levels::AbstractVector; wheel_momentum=nothing, cmg_angles=nothing)::Vector{Float64}
    out = zeros(Float64, mounted_actuator_count(model))
    length(levels) == length(out) || return out
    k = 0
    for link in model.links, thruster in link.thrusters
        k += 1
        out[k] = _delivered_thrust(model, k, thruster, levels[k])
    end
    assembly = model.root.rw_assembly
    first_mounted = assembly.n_wheels - length(assembly.wheels)
    for (j, wheel) in enumerate(assembly.wheels)
        k += 1
        x = actuator_saturate(wheel, levels[k])
        wheel_momentum === nothing || (x = wheel_momentum_saturate(wheel, x, Float64(wheel_momentum[first_mounted + j])))
        out[k] = x
    end
    for link in model.links, rod in link.magnetorquers
        k += 1
        out[k] = actuator_saturate(rod, levels[k])
    end
    for link in model.links, terminal in link.laser_terminals
        k += 1
        out[k] = actuator_saturate(terminal, levels[k])
    end
    for link in model.links, hinge in link.hinges
        k += 1
        out[k] = hinge.angle_rad
    end
    for link in model.links, gimbal in link.gimbals
        k += 1
        out[k] = gimbal.angle_rad
    end
    for (j, cmg) in enumerate(model.root.rw_assembly.cmgs)
        k += 1
        x = actuator_saturate(cmg, levels[k])
        cmg_angles === nothing || (x = cmg_angle_saturate(cmg, x, Float64(cmg_angles[j])))
        out[k] = x
    end
    for link in model.links, actuator in link.custom_actuators
        k += 1
        out[k] = actuator_saturate(actuator, levels[k])
    end
    return out
end

"""
    mounted_power(model::SpacecraftModel, levels; wheel_momentum=nothing) -> Float64

The electrical power (W) the mounted actuators draw at their applied levels
(see [`applied_actuator_levels`](@ref)): the sum of each one's
`actuator_consumption(...).power_w`, evaluated without a sampled environment,
except the reaction wheels, whose power `wheel_power` computes from their
stored momenta in `wheel_momentum` (the root assembly's integrated `h_wheels`)
when it is given. It fills the `actuator_power` result column.
"""
function mounted_power(model::SpacecraftModel, levels::AbstractVector; wheel_momentum=nothing)::Float64
    applied = applied_actuator_levels(model, levels; wheel_momentum=wheel_momentum)
    assembly = model.root.rw_assembly
    first_mounted = assembly.n_wheels - length(assembly.wheels)
    total = 0.0
    j = 0
    for (k, entry) in enumerate(mounted_actuators(model))
        if entry.actuator isa ReactionWheel
            j += 1
            h = wheel_momentum === nothing ? NaN : Float64(wheel_momentum[first_mounted + j])
            total += wheel_power(entry.actuator, applied[k], h)
        else
            total += Float64(actuator_consumption(entry.actuator, applied[k], actuator_context(entry.link)).power_w)
        end
    end
    return total
end

end # module ActuatorMounting
