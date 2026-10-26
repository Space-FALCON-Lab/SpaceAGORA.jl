# Descent control for the Apollo-style landing: throttles the descent engine
# to the guidance's thrust command (with a slew limit), points it with an RCS
# attitude controller that tracks the commanded attitude, and applies the
# thrust along the vehicle's actual engine axis, so pointing errors show up
# in the trajectory. Torques are direct body torques bounded by the RCS
# authority; propellant leaves through both the engine and the RCS.

"""
    ApolloDescentControlConfig(; ...)

Descent engine (`dps_isp_s`, `thrust_slew_n_s`), RCS authority per body axis
(`rcs_torque_limit_nm`, `rcs_isp_s`, `rcs_moment_arm_m`), the attitude
controller (`rate_gain` rad/s of body rate per rad of error, capped at
`rate_limit_rad_s`, and `rate_bandwidth` 1/s for the rate loop), and the
touchdown height: the radial clearance of the vehicle reference point at which
the footpads meet the ground.

`engine_trim = true` (with `actuation = :mounted`) trims the descent engine:
at each control cycle the gimbals mounted on it (`ThrusterGimbal`) are
commanded to the angles that put the thrust line through the centre of mass
(the body origin), as the LM's gimbal drive did; their own angle and rate
limits apply. Default false: the gimbals hold their angles.
"""
Base.@kwdef struct ApolloDescentControlConfig
    dps_isp_s::Float64 = 311.0
    thrust_slew_n_s::Float64 = 22_000.0
    rcs_torque_limit_nm::SVector{3, Float64} = SVector{3, Float64}(1_500.0, 1_500.0, 1_500.0)
    rcs_isp_s::Float64 = 290.0
    rcs_moment_arm_m::Float64 = 1.7
    rate_gain::Float64 = 0.5
    rate_limit_rad_s::Float64 = 0.15
    rate_bandwidth::Float64 = 1.2
    touchdown_height_m::Float64 = 2.5
    engine_trim::Bool = false
end

"""
    DescentThrusterLayout

The vehicle's thrusters as the descent control drives them, in the order the
visualization scene lists them (links in order, each link's `thrusters` in
order): the descent engine is the highest-rated one, every other thruster is
an RCS jet. `torque_to_levels` is the least-norm map from a commanded body
torque to the jets' firing levels, `pinv(A)` of the 3 x n matrix whose column
j is the torque a jet produces at full thrust about the spacecraft reference
point. Each jet's location and direction are taken to the spacecraft body
frame through its link's configured pose, as the mounted wrench does;
articulated/gimbaled jets need a separate allocator. In `:ideal` actuation the dynamics apply the bounded
body torque directly; in `:mounted` actuation the jets' allocated levels act.
"""
struct DescentThrusterLayout
    engine::Int
    engine_max_thrust_n::Float64
    jets::Vector{Int}
    torque_arms_nm::Matrix{Float64}      # 3 x n_jets, column j at full thrust
    torque_to_levels::Matrix{Float64}    # n_jets x 3
end

"Per-spacecraft actuator state of the descent control effector."
mutable struct ApolloDescentControlState
    thrust_n::Vector{Float64}
    torque_nm::Vector{SVector{3, Float64}}
    attitude_error_rad::Vector{Float64}
    last_update_s::Vector{Float64}
    # Firing level (0 to 1) of every thruster, in scene order, and the layout it
    # is allocated over (built on the first control cycle from the spacecraft).
    thruster_level::Vector{Vector{Float64}}
    thruster_layout::Vector{Union{Nothing, DescentThrusterLayout}}
    # The mounted actuators' levels in `:mounted` actuation, in mounting order:
    # the thrusters' (the same levels as thrusts in newtons) first, then every
    # other mounted actuator holding (see `_hold_other_mounted_levels!`).
    thruster_force_n::Vector{Vector{Float64}}
end
function ApolloDescentControlState(n::Integer)
    m = Int(n)
    m > 0 || throw(ArgumentError("ApolloDescentControlState needs at least one spacecraft"))
    return ApolloDescentControlState(
        zeros(m),
        fill(SVector{3, Float64}(0.0, 0.0, 0.0), m),
        fill(NaN, m),
        fill(NaN, m),
        [Float64[] for _ in 1:m],
        Union{Nothing, DescentThrusterLayout}[nothing for _ in 1:m],
        [Float64[] for _ in 1:m],
    )
end

# Each run starts the selected spacecraft's actuator entries as a new state holds
# them (G-8); the level vectors are emptied in place, so a handle to them stays the
# spacecraft's, and the layout is rebuilt at the first control cycle.
function _reset_descent_actuators!(act::ApolloDescentControlState, indices)
    for i in indices
        act.thrust_n[i] = 0.0
        act.torque_nm[i] = SVector{3, Float64}(0.0, 0.0, 0.0)
        act.attitude_error_rad[i] = NaN
        act.last_update_s[i] = NaN
        empty!(act.thruster_level[i])
        act.thruster_layout[i] = nothing
        empty!(act.thruster_force_n[i])
    end
    return act
end

"""
    ApolloDescentControlModel(config, guidance_config, state, terrain=NoTerrainModel(); actuation=:ideal)

Control effector paired with `ApolloDescentGuidanceModel` through the
shared [`ApolloDescentState`](@ref). Carries the terrain model so the engine
registers the touchdown event on it. When a configuration has both models,
`run_simulation` refuses it before the run, with an `ArgumentError`, if any
spacecraft is guided or controlled without its counterpart on the same state.
A guidance (or control) model of another type in the configuration may write (or
fly) the state, so the check then skips spacecraft without Apollo guidance (or
control).

`actuation=:ideal` applies the slew-limited engine thrust along body -z and
the clamped attitude torque directly, with the configured specific impulses.
`actuation=:mounted` drives the vehicle's mounted thrusters instead: the
engine at the slewed thrust and the RCS jets at their least-norm allocation
of the clamped torque, so the torque, the force and the propellant use follow
from the jets' geometry, clipping and specific impulses.
"""
struct ApolloDescentControlModel{T <: AbstractTerrainModel, N} <: AbstractControlEffectorModel
    config::ApolloDescentControlConfig
    guidance::ApolloDescentConfig
    state::ApolloDescentState
    terrain::T
    actuators::ApolloDescentControlState
    spacecraft_indices::NTuple{N, Int}
    actuation::Symbol
end
function ApolloDescentControlModel(config::ApolloDescentControlConfig, guidance::ApolloDescentConfig,
                                   state::ApolloDescentState, terrain::AbstractTerrainModel=NoTerrainModel();
                                   spacecraft_indices=eachindex(state.phase), actuation::Symbol=:ideal)
    _validate_descent_config(guidance, terrain)
    for name in (:dps_isp_s, :rcs_isp_s, :rcs_moment_arm_m, :rate_limit_rad_s)
        value = getfield(config, name)
        isfinite(value) && value > 0 || throw(ArgumentError("$name must be finite and positive"))
    end
    for name in (:thrust_slew_n_s, :rate_gain, :rate_bandwidth, :touchdown_height_m)
        value = getfield(config, name)
        isfinite(value) && value >= 0 || throw(ArgumentError("$name must be finite and nonnegative"))
    end
    all(x -> isfinite(x) && x >= 0, config.rcs_torque_limit_nm) ||
        throw(ArgumentError("RCS torque limits must be finite and nonnegative"))
    actuation in (:ideal, :mounted) || throw(ArgumentError("actuation must be :ideal or :mounted, got $(repr(actuation))"))
    config.engine_trim && actuation !== :mounted &&
        throw(ArgumentError("engine_trim drives the descent engine's mounted gimbals; it needs actuation=:mounted"))
    return ApolloDescentControlModel(config, guidance, state, terrain,
        ApolloDescentControlState(length(state.phase)), _descent_indices(state, spacecraft_indices), actuation)
end

import ..GuidanceHooks

# The guidance hands its commands to this control through the shared state.
# When a configuration has both models, a spacecraft controlled on a state that
# no Apollo guidance guides it on never gets a command (the engine stays off),
# and one guided without control is never flown. A configuration with only one
# side is not checked. Nor is a direction where the configuration also has a
# guidance (or control) model of another type: it may write (or fly) the
# state, and which spacecraft it serves is not declared.
function _apollo_unpaired_spacecraft(args)
    guided = IdDict{ApolloDescentState, Set{Int}}()
    controlled = IdDict{ApolloDescentState, Set{Int}}()
    other_guidance = false
    other_control = false
    for g in args.guidance_model.guidance_effectors
        if g isa GuidanceHooks.ApolloDescentGuidanceModel
            union!(get!(Set{Int}, guided, g.state), g.spacecraft_indices)
        else
            other_guidance = true
        end
    end
    for c in args.control_model.control_effectors
        if c isa ApolloDescentControlModel
            union!(get!(Set{Int}, controlled, c.state), c.spacecraft_indices)
        else
            other_control = true
        end
    end
    (isempty(guided) || isempty(controlled)) && return (unguided=Int[], uncontrolled=Int[])
    unpaired(from, to) = sort!(unique!([i for (state, ids) in from for i in setdiff(ids, get(to, state, Set{Int}()))]))
    return (unguided=other_guidance ? Int[] : unpaired(controlled, guided),
            uncontrolled=other_control ? Int[] : unpaired(guided, controlled))
end

function SimulationLifecycle.initialize_control!(c::ApolloDescentControlModel, u, p, t)
    _reset_descent_actuators!(c.actuators, c.spacecraft_indices)
    return nothing
end

# Its spacecraft's actuator entries go into checkpoint files, by field (the shared
# descent state is the guidance's).
function SimulationLifecycle.checkpoint_state(c::ApolloDescentControlModel)
    act, ids = c.actuators, collect(c.spacecraft_indices)
    return NamedTuple{fieldnames(ApolloDescentControlState)}(
        map(name -> deepcopy(getfield(act, name)[ids]), fieldnames(ApolloDescentControlState)))
end
function SimulationLifecycle.restore_checkpoint_state!(c::ApolloDescentControlModel, data)
    for name in fieldnames(ApolloDescentControlState), (k, i) in enumerate(c.spacecraft_indices)
        getfield(c.actuators, name)[i] = getfield(data, name)[k]
    end
    return nothing
end

function SimulationLifecycle.preflight_control(c::ApolloDescentControlModel, args; isolate_state)
    pairing = _apollo_unpaired_spacecraft(args)
    isempty(pairing.unguided) && isempty(pairing.uncontrolled) && return nothing
    problems = String[]
    isempty(pairing.unguided) || push!(problems, "spacecraft $(join(pairing.unguided, ", ")) " *
        "have an ApolloDescentControlModel but no ApolloDescentGuidanceModel on the same ApolloDescentState, " *
        "so their descent engine never receives a command")
    isempty(pairing.uncontrolled) || push!(problems, "spacecraft $(join(pairing.uncontrolled, ", ")) " *
        "have an ApolloDescentGuidanceModel but no ApolloDescentControlModel on the same ApolloDescentState, " *
        "so no control flies their guidance commands")
    throw(ArgumentError(
        "Apollo descent guidance and control are not paired: $(join(problems, "; ")). Build both models " *
        "from one ApolloDescentState with the same spacecraft_indices."))
end

const _DESCENT_G0 = 9.80665

"""
    attitude_error_vector(q, q_cmd) -> SVector{3}

Rotation vector (body axes, radians) from the commanded attitude to the
current one, negated: the direction the body must turn. Both quaternions are
scalar-last body-to-inertial rotations.
"""
@inline function attitude_error_vector(q::SVector{4, Float64}, q_cmd::SVector{4, Float64})::SVector{3, Float64}
    all(isfinite, q) && all(isfinite, q_cmd) && norm(q) > eps(Float64) && norm(q_cmd) > eps(Float64) ||
        throw(ArgumentError("attitude quaternions must be finite and nonzero"))
    current, target = q / norm(q), q_cmd / norm(q_cmd)
    relative = quat_mult(SVector(-current[1], -current[2], -current[3], current[4]), target)
    relative[4] < 0 && (relative = -relative)
    vector = SVector(relative[1], relative[2], relative[3])
    n = norm(vector)
    return n > eps(Float64) ? vector * (2 * atan(n, relative[4]) / n) : 2 * vector
end

"""
    descent_thruster_layout(model) -> Union{Nothing, DescentThrusterLayout}

The vehicle's thrusters in scene order, split into the descent engine (the
highest-rated one) and the RCS jets, with the least-norm map from a commanded
body torque to the jets' firing levels. `nothing` when the spacecraft carries
no thrusters.

A `Thruster`'s `direction` is its thrust direction (the descent engine's is
body -z), so a jet at full thrust puts a force `max_thrust * direction` on the
vehicle and a torque `location x force` about the spacecraft reference point.
"""
function descent_thruster_layout(model)::Union{Nothing, DescentThrusterLayout}
    locations = SVector{3, Float64}[]
    directions = SVector{3, Float64}[]
    max_thrust = Float64[]
    for link in model.links
        # The link-to-body map of the mounted wrench (ActuatorMounting.actuator_context).
        # A root link at the body origin needs none; skipping it keeps that
        # arithmetic, and so every existing layout, bit for bit unchanged.
        offset = SVector{3, Float64}(link.r)
        identity_pose = link.root && iszero(offset)
        to_body = link.root ? nothing : SMatrix{3, 3, Float64, 9}(rot(SVector{4, Float64}(link.q))')
        for thruster in link.thrusters
            location = SVector{3, Float64}(thruster.location)
            direction = SVector{3, Float64}(thruster.direction)
            if !identity_pose
                location = (to_body === nothing ? location : to_body * location) + offset
                direction = to_body === nothing ? direction : to_body * direction
            end
            push!(locations, location)
            push!(directions, direction)
            push!(max_thrust, Float64(thruster.max_thrust))
        end
    end
    isempty(max_thrust) && return nothing
    engine = argmax(max_thrust)
    jets = [k for k in eachindex(max_thrust) if k != engine]
    arms = zeros(3, length(jets))
    for (j, k) in enumerate(jets)
        d = directions[k]
        n = norm(d)
        n > eps(Float64) || continue
        τ = cross(locations[k], max_thrust[k] * (d / n))
        arms[1, j] = τ[1]; arms[2, j] = τ[2]; arms[3, j] = τ[3]
    end
    # Least-norm allocation: u = A' (A A')^-1 τ, through `pinv` so a jet set
    # that cannot reach all three axes still gives the best-fit torque.
    return DescentThrusterLayout(engine, max_thrust[engine], jets, arms, isempty(jets) ? zeros(0, 3) : pinv(arms))
end

"""
    descent_thruster_levels!(levels, layout, thrust_n, torque_nm) -> levels

Firing levels (0 to 1) of every thruster: the descent engine's actual thrust
over its rating, and the RCS jets from the least-norm allocation of the
commanded body torque over their torque arms, each clipped into 0 to 1.
"""
function descent_thruster_levels!(levels::Vector{Float64}, layout::DescentThrusterLayout, thrust_n::Float64, torque_nm::SVector{3, Float64})
    fill!(levels, 0.0)
    if layout.engine_max_thrust_n > 0.0 && isfinite(thrust_n)
        levels[layout.engine] = clamp(thrust_n / layout.engine_max_thrust_n, 0.0, 1.0)
    end
    isempty(layout.jets) && return levels
    for (j, k) in enumerate(layout.jets)
        u = layout.torque_to_levels[j, 1] * torque_nm[1] +
            layout.torque_to_levels[j, 2] * torque_nm[2] +
            layout.torque_to_levels[j, 3] * torque_nm[3]
        levels[k] = isfinite(u) ? clamp(u, 0.0, 1.0) : 0.0
    end
    return levels
end

# Refresh the firing levels for the viewer's plumes. Called once per control
# cycle so the hook below only reads stored state.
@inline function _update_descent_thruster_levels!(model::ApolloDescentControlModel, p::ODEParams, i::Int)
    act = model.actuators
    layout = act.thruster_layout[i]
    spacecraft = p.args.dynamics_model.spacecraft[i]
    if layout === nothing
        layout = descent_thruster_layout(spacecraft)
        layout === nothing && return nothing
        act.thruster_layout[i] = layout
        act.thruster_level[i] = zeros(length(layout.jets) + 1)
        # One mounted level per mounted actuator; the thrusters come first.
        act.thruster_force_n[i] = zeros(max(mounted_actuator_count(spacecraft), length(layout.jets) + 1))
    end
    levels = descent_thruster_levels!(act.thruster_level[i], layout, act.thrust_n[i], act.torque_nm[i])
    forces = act.thruster_force_n[i]
    k = 0
    for link in spacecraft.links, thruster in link.thrusters
        k += 1
        k <= length(forces) && (forces[k] = levels[k] * thruster.max_thrust)
    end
    _hold_other_mounted_levels!(forces, spacecraft, k)
    model.config.engine_trim && _trim_engine_gimbals!(forces, spacecraft, layout.engine)
    return nothing
end

# Engine trim: the descent engine's gimbals at the angles that point its thrust
# along the line from the nozzle through the body origin (the centre of mass),
# so the engine puts no torque on the vehicle. The angles are solved one gimbal
# at a time (each turns the thrust about its own axis, in mounting order), and
# the engine's apply step limits them as it does any gimbal command.
function _trim_engine_gimbals!(forces::Vector{Float64}, spacecraft, engine::Int)
    k = 0
    target = nothing
    gimbals = Tuple{Int, ThrusterGimbal}[]
    engine_thruster = nothing
    for link in spacecraft.links, thruster in link.thrusters
        k += 1
        k == engine || continue
        engine_thruster = thruster
        ctx = actuator_context(link)
        nozzle = ctx.link_to_body * SVector{3, Float64}(thruster.location) + ctx.link_r_body
        d_body = ctx.link_to_body * SVector{3, Float64}(thruster.direction)
        n = norm(nozzle)
        n > 1e-9 || return forces                    # at the centre of mass: any direction is trimmed
        line = dot(d_body, nozzle) < 0.0 ? -nozzle / n : nozzle / n
        target = ctx.link_to_body' * line
    end
    (engine_thruster === nothing || target === nothing) && return forces
    for (m, entry) in enumerate(mounted_actuators(spacecraft))
        a = entry.actuator
        a isa ThrusterGimbal && a.thruster === engine_thruster && m <= length(forces) && push!(gimbals, (m, a))
    end
    isempty(gimbals) && return forces
    d = SVector{3, Float64}(engine_thruster.direction)
    d = d / norm(d)
    angles = Float64[g.angle_rad for (_, g) in gimbals]
    for _ in 1:(length(gimbals) == 1 ? 1 : 12)
        for j in eachindex(gimbals)
            # The thrust as the earlier gimbals turn it, and the target as the
            # later ones would see it (their turns undone, last first).
            v = d
            for m in 1:(j - 1)
                v = _turn_about(gimbals[m][2].axis, angles[m], v)
            end
            w = target
            for m in length(gimbals):-1:(j + 1)
                w = _turn_about(gimbals[m][2].axis, -angles[m], w)
            end
            angles[j] = _best_axis_angle(gimbals[j][2].axis, v, w, angles[j])
        end
    end
    for (j, (m, _)) in enumerate(gimbals)
        forces[m] = angles[j]
    end
    return forces
end

# `v` turned by `angle` about the unit vector `axis` (Rodrigues), as a gimbal turns its thrust.
@inline _turn_about(axis::SVector{3, Float64}, angle::Float64, v::SVector{3, Float64})::SVector{3, Float64} =
    v * cos(angle) + cross(axis, v) * sin(angle) + axis * (dot(axis, v) * (1.0 - cos(angle)))

# The angle about the unit vector `axis` that turns `v` closest to `w`; the
# current angle when either lies along the axis.
@inline function _best_axis_angle(axis::SVector{3, Float64}, v::SVector{3, Float64}, w::SVector{3, Float64}, current::Float64)::Float64
    vp = v - dot(axis, v) * axis
    wp = w - dot(axis, w) * axis
    (norm(vp) > 1e-12 && norm(wp) > 1e-12) || return current
    return atan(dot(axis, cross(vp, wp)), dot(vp, wp))
end

# The controller commands only the thrusters; every other mounted actuator
# holds: a hinge or a thruster gimbal at its current angle, any other kind at
# level zero. Refreshed at each control cycle.
function _hold_other_mounted_levels!(forces::Vector{Float64}, spacecraft, n_thrusters::Int)
    length(forces) > n_thrusters || return forces
    for (k, entry) in enumerate(mounted_actuators(spacecraft))
        k <= n_thrusters && continue
        k > length(forces) && break
        a = entry.actuator
        forces[k] = a isa Union{PanelHinge, ThrusterGimbal} ? a.angle_rad : 0.0
    end
    return forces
end

# The thrusters off; the other mounted actuators keep holding.
@inline function _zero_descent_thrusters!(act::ApolloDescentControlState, i::Int)
    fill!(act.thruster_level[i], 0.0)
    forces = act.thruster_force_n[i]
    @inbounds for k in 1:min(length(act.thruster_level[i]), length(forces))
        forces[k] = 0.0
    end
    return nothing
end

function control_thruster_levels(model::ApolloDescentControlModel, i::Int)
    i in model.spacecraft_indices || return nothing
    return model.actuators.thruster_level[i]
end

actuation_mode(model::ApolloDescentControlModel) = model.actuation

# The mounted thrusters' levels in newtons; none until the first control cycle
# has built the layout.
function actuator_levels(model::ApolloDescentControlModel, i::Int)
    i in model.spacecraft_indices || return nothing
    forces = model.actuators.thruster_force_n[i]
    return isempty(forces) ? nothing : forces
end

@inline function _descent_state_view(u, i::Int)
    return hasproperty(u, :sc) ? u.sc[i] : u
end

# The engine does not tick an inactive spacecraft (skip_inactive): at each control
# tick after touchdown or impact it releases the command instead, engine and jets
# off; the other mounted actuators keep holding.
function GNCRuntime.release_command!(model::ApolloDescentControlModel, t, i)
    i in model.spacecraft_indices || return nothing
    model.actuators.thrust_n[i] = 0.0
    model.actuators.torque_nm[i] = SVector{3, Float64}(0.0, 0.0, 0.0)
    _zero_descent_thrusters!(model.actuators, i)
    return nothing
end

function update_control!(model::ApolloDescentControlModel, u, p::ODEParams, t::Float64, i::Int64)
    n = length(model.actuators.thrust_n)
    i in model.spacecraft_indices || return nothing
    sc = _descent_state_view(u, i)
    hasproperty(sc, :q) && hasproperty(sc, :ω) || throw(ArgumentError("ApolloDescentControlModel needs orientation_sim=true (attitude state q and ω)."))
    cfg = model.config
    act = model.actuators
    # nothing to track until the guidance has issued its first command
    if isnan(model.state.phase_start_s[1, i])
        act.thrust_n[i] = 0.0
        act.torque_nm[i] = SVector{3, Float64}(0.0, 0.0, 0.0)
        act.last_update_s[i] = t
        _update_descent_thruster_levels!(model, p, i)
        return nothing
    end
    dt = isnan(act.last_update_s[i]) ? 0.0 : max(0.0, t - act.last_update_s[i])
    act.last_update_s[i] = t
    # engine: slew toward the guidance command
    cmd = model.state.thrust_cmd_n[i]
    if dt > 0.0
        step = cfg.thrust_slew_n_s * dt
        act.thrust_n[i] = clamp(cmd, act.thrust_n[i] - step, act.thrust_n[i] + step)
    end
    # A first call at zero elapsed time retains the initial actuator thrust;
    # later calls always respect the configured slew limit.
    # attitude: rate command from the error, torque from the rate error, bounded by the RCS
    q = SVector{4, Float64}(sc.q[1], sc.q[2], sc.q[3], sc.q[4])
    ω = SVector{3, Float64}(sc.ω[1], sc.ω[2], sc.ω[3])
    e = attitude_error_vector(q, model.state.attitude_cmd[i])
    act.attitude_error_rad[i] = norm(e)
    ω_des = cfg.rate_gain * e
    nd = norm(ω_des)
    if nd > cfg.rate_limit_rad_s
        ω_des = ω_des * (cfg.rate_limit_rad_s / nd)
    end
    ω_des = ω_des + model.state.rate_cmd[i]
    inertia = p.args.dynamics_model.spacecraft[i].inertia_tensor
    τ = SVector{3, Float64}(inertia * (cfg.rate_bandwidth * (ω_des - ω)))
    lim = cfg.rcs_torque_limit_nm
    act.torque_nm[i] = SVector{3, Float64}(clamp(τ[1], -lim[1], lim[1]), clamp(τ[2], -lim[2], lim[2]), clamp(τ[3], -lim[3], lim[3]))
    _update_descent_thruster_levels!(model, p, i)
    return nothing
end

function control_wrench(model::ApolloDescentControlModel, x::StateSample, p::ODEParams, t::Float64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    i = x.spacecraft_index
    n = length(model.actuators.thrust_n)
    (i in model.spacecraft_indices && p.is_active[i]) || return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    thrust = model.actuators.thrust_n[i]
    q = x.q_ib::SVector{4, Float64}              # the descent controls attitude, so the run propagates it
    A = rot(q)                                   # rows: body axes in inertial coordinates
    engine_axis_i = -SVector{3, Float64}(A[3, 1], A[3, 2], A[3, 3])   # thrust pushes along body -z
    return thrust * engine_axis_i, model.actuators.torque_nm[i]
end

function control_propellant_rate(model::ApolloDescentControlModel, x::StateSample, p::ODEParams, t::Float64)::Float64
    i = x.spacecraft_index
    n = length(model.actuators.thrust_n)
    (i in model.spacecraft_indices && p.is_active[i]) || return 0.0
    cfg = model.config
    dps = propellant_mass_rate(model.actuators.thrust_n[i], cfg.dps_isp_s)
    τ = model.actuators.torque_nm[i]
    rcs = (abs(τ[1]) + abs(τ[2]) + abs(τ[3])) / (cfg.rcs_moment_arm_m * cfg.rcs_isp_s * _DESCENT_G0)
    return dps + rcs
end

function touchdown_spec(model::ApolloDescentControlModel, spacecraft_index::Int)
    spacecraft_index in model.spacecraft_indices || return nothing
    frame = model.state.site_frame
    radius = model.guidance.reference_radius_m
    function on_touchdown(t::Float64, r_p::SVector{3, Float64}, v_p::SVector{3, Float64}, i::Int)
        st = model.state
        st.touchdown_s[i] = t
        st.phase[i] = :landed
        st.phase_start_s[4, i] = t
        st.thrust_cmd_n[i] = 0.0; st.throttle[i] = 0.0
        model.actuators.thrust_n[i] = 0.0
        model.actuators.torque_nm[i] = SVector{3, Float64}(0.0, 0.0, 0.0)
        _zero_descent_thrusters!(model.actuators, i)
        f = frame[i]
        if f !== nothing
            d = r_p - f.origin_p
            v_s = SVector{3, Float64}(dot(v_p, f.x), dot(v_p, f.y), dot(v_p, f.z))
            st.touchdown_v_mps[i] = v_s
            st.touchdown_miss_m[i] = hypot(dot(d, f.x), dot(d, f.y))
        else
            st.touchdown_v_mps[i] = v_p
        end
        return nothing
    end
    return (terrain=model.terrain, reference_radius_m=radius, height_m=model.config.touchdown_height_m, on_touchdown=on_touchdown)
end
