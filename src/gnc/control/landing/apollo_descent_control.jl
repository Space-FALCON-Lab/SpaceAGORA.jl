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
touchdown height: the radar altitude of the vehicle reference point at which
the footpads meet the ground.
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
end

"""
    DescentThrusterLayout

The vehicle's thrusters as the descent control drives them, in the order the
visualization scene lists them (links in order, each link's `thrusters` in
order): the descent engine is the highest-rated one, every other thruster is
an RCS jet. `torque_to_levels` is the least-norm map from a commanded body
torque to the jets' firing levels, `pinv(A)` of the 3 x n matrix whose column
j is the torque a jet produces at full thrust about the spacecraft reference
point.
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
end
function ApolloDescentControlState(n::Integer)
    m = Int(n)
    return ApolloDescentControlState(
        zeros(m),
        fill(SVector{3, Float64}(0.0, 0.0, 0.0), m),
        fill(NaN, m),
        fill(NaN, m),
        [Float64[] for _ in 1:m],
        Union{Nothing, DescentThrusterLayout}[nothing for _ in 1:m],
    )
end

"""
    ApolloDescentControlModel(config, guidance_config, state, terrain=NoTerrainModel())

Control effector paired with [`ApolloDescentGuidanceModel`](@ref) through the
shared [`ApolloDescentState`](@ref). Carries the terrain model so the engine
registers the touchdown event on it.
"""
struct ApolloDescentControlModel{T <: AbstractTerrainModel} <: AbstractControlEffectorModel
    config::ApolloDescentControlConfig
    guidance::ApolloDescentConfig
    state::ApolloDescentState
    terrain::T
    actuators::ApolloDescentControlState
end
function ApolloDescentControlModel(config::ApolloDescentControlConfig, guidance::ApolloDescentConfig, state::ApolloDescentState, terrain::AbstractTerrainModel=NoTerrainModel())
    return ApolloDescentControlModel(config, guidance, state, terrain, ApolloDescentControlState(length(state.phase)))
end

const _DESCENT_G0 = 9.80665

"""
    attitude_error_vector(q, q_cmd) -> SVector{3}

Rotation vector (body axes, radians) from the commanded attitude to the
current one, negated: the direction the body must turn. Both quaternions are
scalar-last body-to-inertial rotations.
"""
@inline function attitude_error_vector(q::SVector{4, Float64}, q_cmd::SVector{4, Float64})::SVector{3, Float64}
    A = rot(q); Ac = rot(q_cmd)
    Re = A * Ac'                       # ≈ I - δ× for the body rotated by δ (body axes) from the command
    s = SVector{3, Float64}(Re[3, 2] - Re[2, 3], Re[1, 3] - Re[3, 1], Re[2, 1] - Re[1, 2]) / 2
    c = clamp((tr(Re) - 1) / 2, -1.0, 1.0)
    θ = acos(c)
    n = norm(s)
    return n > 1e-12 ? s * (θ / n) : s   # s = -δ sin θ / ... scaled to the angle: the turn toward the command
end

"""
    descent_thruster_layout(model) -> Union{Nothing, DescentThrusterLayout}

The vehicle's thrusters in scene order, split into the descent engine (the
highest-rated one) and the RCS jets, with the least-norm map from a commanded
body torque to the jets' firing levels. `nothing` when the spacecraft carries
no thrusters.

A `Thruster`'s `direction` is the exhaust direction (the descent engine's is
body +z while its thrust pushes along body -z), so a jet at full thrust puts a
force `-max_thrust * direction` on the vehicle and a torque
`location x force` about the spacecraft reference point.
"""
function descent_thruster_layout(model)::Union{Nothing, DescentThrusterLayout}
    locations = SVector{3, Float64}[]
    directions = SVector{3, Float64}[]
    max_thrust = Float64[]
    for link in model.links, thruster in link.thrusters
        push!(locations, SVector{3, Float64}(thruster.location))
        push!(directions, SVector{3, Float64}(thruster.direction))
        push!(max_thrust, Float64(thruster.max_thrust))
    end
    isempty(max_thrust) && return nothing
    engine = argmax(max_thrust)
    jets = [k for k in eachindex(max_thrust) if k != engine]
    arms = zeros(3, length(jets))
    for (j, k) in enumerate(jets)
        d = directions[k]
        n = norm(d)
        n > eps(Float64) || continue
        τ = cross(locations[k], -max_thrust[k] * (d / n))
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
    @inbounds for (j, k) in enumerate(layout.jets)
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
    if layout === nothing
        layout = descent_thruster_layout(p.args.dynamics_model.spacecraft[i])
        layout === nothing && return nothing
        act.thruster_layout[i] = layout
        act.thruster_level[i] = zeros(length(layout.jets) + 1)
    end
    descent_thruster_levels!(act.thruster_level[i], layout, act.thrust_n[i], act.torque_nm[i])
    return nothing
end

function control_thruster_levels(model::ApolloDescentControlModel, i::Int)
    (1 <= i <= length(model.actuators.thruster_level)) || return nothing
    return model.actuators.thruster_level[i]
end

@inline function _descent_state_view(u, i::Int)
    return hasproperty(u, :sc) ? u.sc[i] : u
end

function calcControlEffect!(model::ApolloDescentControlModel, u, p::ODEParams, t::Float64, i::Int64)
    n = length(model.actuators.thrust_n)
    1 <= i <= n || return nothing
    sc = _descent_state_view(u, i)
    hasproperty(sc, :q) && hasproperty(sc, :ω) || throw(ArgumentError("ApolloDescentControlModel needs orientation_sim=true (attitude state q and ω)."))
    cfg = model.config
    act = model.actuators
    # nothing to track until the guidance has issued its first command
    if isnan(model.state.phase_start_s[1, i])
        act.thrust_n[i] = 0.0
        act.torque_nm[i] = SVector{3, Float64}(0.0, 0.0, 0.0)
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
    else
        act.thrust_n[i] = cmd
    end
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

function calcControlForceTorque(model::ApolloDescentControlModel, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    n = length(model.actuators.thrust_n)
    (1 <= i <= n) || return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    thrust = model.actuators.thrust_n[i]
    q = SVector{4, Float64}(u.q[1], u.q[2], u.q[3], u.q[4])
    A = rot(q)                                   # rows: body axes in inertial coordinates
    engine_axis_i = -SVector{3, Float64}(A[3, 1], A[3, 2], A[3, 3])   # thrust pushes along body -z
    return thrust * engine_axis_i, model.actuators.torque_nm[i]
end

function calcControlMassFlowRate(model::ApolloDescentControlModel, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)::Float64
    n = length(model.actuators.thrust_n)
    (1 <= i <= n) || return 0.0
    cfg = model.config
    dps = model.actuators.thrust_n[i] / (cfg.dps_isp_s * _DESCENT_G0)
    τ = model.actuators.torque_nm[i]
    rcs = (abs(τ[1]) + abs(τ[2]) + abs(τ[3])) / (cfg.rcs_moment_arm_m * cfg.rcs_isp_s * _DESCENT_G0)
    return -(dps + rcs)
end

"""
    touchdown_spec(effector)

Hook for the engine's event registry: `nothing` for effectors without a
landing, else a NamedTuple with the terrain model, the reference radius
(m), the touchdown height (m) and an `on_touchdown(t, r_p, v_p, i)` record
callback. The registry then replaces the impact callback with a touchdown
event on the terrain.
"""
touchdown_spec(effector) = nothing
function touchdown_spec(model::ApolloDescentControlModel)
    frame = model.state.site_frame
    radius = model.terrain isa DEMTerrainModel ? model.terrain.reference_radius_m : NaN
    function on_touchdown(t::Float64, r_p::SVector{3, Float64}, v_p::SVector{3, Float64}, i::Int)
        st = model.state
        st.touchdown_s[i] = t
        st.phase[i] = :landed
        st.phase_start_s[4, i] = t
        st.thrust_cmd_n[i] = 0.0; st.throttle[i] = 0.0
        model.actuators.thrust_n[i] = 0.0
        model.actuators.torque_nm[i] = SVector{3, Float64}(0.0, 0.0, 0.0)
        fill!(model.actuators.thruster_level[i], 0.0)
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
