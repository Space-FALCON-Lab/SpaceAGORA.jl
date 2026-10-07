# Standalone MuJoCo proximity-scene runner (Stage 1).
#
# Frame: the MuJoCo world origin is a virtual chief point that follows the scene's mass-weighted center of
# mass under SpaceAGORA's gravity plus the feedback acceleration `a_fb` (mass-averaged non-gravitational
# acceleration). World axes are parallel to the inertial (J2000/ECI) axes and `opt.gravity` is zero, so the
# frame translates but does not rotate and no Coriolis terms arise. Each body receives
#     F_i = m_i * dg_i + F_ext_i - m_i * a_fb            at its center of mass, world axes,
# where dg_i = g(R + rho_i) - g(R) is SpaceAGORA's Encke relative gravity and rho_i is the body COM position
# relative to the chief. The wrench is evaluated from the kinematics at the start of the step and held over
# one `dt`, which is MuJoCo's documented `mj_step1`/`mj_step2` control pattern (first-order in dt).

using LinearAlgebra
using StaticArrays
import SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const V3 = SVector{3, Float64}

# RK4 is not offered: the runner steps with mj_step1/mj_step2 so it can write wrenches from fresh
# kinematics, and MuJoCo's mj_step2 integrates RK4 models with Euler. Accepting :rk4 would silently
# run Euler, so it is refused until a full-step (mj_step) path supports it.
const _INTEGRATORS = (euler = Binding.INT_EULER, implicit = Binding.INT_IMPLICIT,
    implicitfast = Binding.INT_IMPLICITFAST)

"""
    SceneBodyState(body, r, v; q=(1,0,0,0), ω=(0,0,0))

Absolute initial state of one free-root body of a scene: inertial position `r` [m] and velocity `v` [m/s]
of the body center of mass, attitude quaternion `q` (MuJoCo order `(w, x, y, z)`, body to world) and
body-frame angular velocity `ω` [rad/s]. `body` is the MJCF body name.
"""
struct SceneBodyState
    body::String
    r::V3
    v::V3
    q::SVector{4, Float64}
    ω::V3
end
SceneBodyState(body::AbstractString, r, v; q=(1.0, 0.0, 0.0, 0.0), ω=(0.0, 0.0, 0.0)) =
    SceneBodyState(String(body), V3(r), V3(v), normalize(SVector{4, Float64}(q)), V3(ω))

"""
    body_state_from_initial_condition(body, ic::SpaceAGORA.CartesianInitialCondition)

Map a SpaceAGORA spacecraft initial condition to a [`SceneBodyState`](@ref). Only position and velocity are
mapped; the attitude convention between SpaceAGORA and MuJoCo is Stage 2 work, so a non-identity attitude
or a nonzero angular rate is rejected rather than silently dropped.
"""
function body_state_from_initial_condition(body::AbstractString, ic)
    isapprox(ic.q, SVector(0.0, 0.0, 0.0, 1.0); atol = 1e-12) && iszero(ic.ang_vel) || throw(ArgumentError(
        "body_state_from_initial_condition: attitude/angular-rate mapping is not implemented; use identity attitude and zero rate, or build the SceneBodyState directly."))
    return SceneBodyState(body, ic.pos, ic.vel)
end

"""
    SceneState

Complete restorable state of a scene: step counter `n`, chief position `R` and velocity `V`, and MuJoCo's
`mjSTATE_INTEGRATION` vector (time, qpos, qvel, act, warmstart, user inputs). Plain data, safe to keep or copy.
"""
struct SceneState
    n::Int
    R::V3
    V::V3
    mj::Vector{Float64}
end

"""
    ProximityScene(; mjcf_path | mjcf_xml, dt, planet, gravity_effectors, initial_states,
                     integrator=:implicitfast, planet_rotation=nothing, t0=0.0)

A MuJoCo scene (MJCF file or string) flown around a virtual chief point. The scene step `dt` has no default
and must be given. `gravity_effectors` is a tuple of SpaceAGORA position-only gravity effectors
(for example `(SpaceAGORA.InverseSquaredJ2GravityModel(),)`); `planet` is the SpaceAGORA planet they use.
`planet_rotation(t) -> SMatrix{3,3}` returns the inertial-to-planet-fixed rotation at scene time `t` and is
required when an effector needs the planet frame (J2, harmonics). `initial_states` holds one
[`SceneBodyState`](@ref) per free-root body.

Every scene owns its `mjModel` and `mjData`; nothing is shared between scenes and no global state is
written, so scenes can be stepped on separate tasks. Native memory is freed by finalizers.
"""
mutable struct ProximityScene{E <: Tuple, F}
    model::Binding.MjModel
    data::Binding.MjData
    dt::Float64
    integrator::Symbol
    planet::Any
    effectors::E
    planet_rotation::F
    t0::Float64
    initial_states::Vector{SceneBodyState}
    n::Int
    R::V3
    V::V3
    fresh::Bool                      # kinematics (xipos, cvel) current with qpos/qvel
    nb::Int                          # bodies excluding world; body i has MuJoCo id i
    names::Vector{String}
    roots::Vector{Int}               # free-root body ids
    mass::Vector{Float64}            # view of body_mass, indexed id + 1
    total_mass::Float64
    qpos::Vector{Float64}
    qvel::Vector{Float64}
    ctrl::Vector{Float64}
    xfrc::Matrix{Float64}            # 6 x nbody
    xipos::Matrix{Float64}           # 3 x nbody
    xquat::Matrix{Float64}           # 4 x nbody
    vel6::Vector{Float64}            # scratch for mj_objectVelocity
end

function ProximityScene(; mjcf_path=nothing, mjcf_xml=nothing, dt, planet, gravity_effectors::Tuple,
        initial_states::AbstractVector{SceneBodyState}, integrator::Symbol=:implicitfast,
        planet_rotation=nothing, t0::Real=0.0)
    (mjcf_path === nothing) ⊻ (mjcf_xml === nothing) || throw(ArgumentError("give exactly one of mjcf_path or mjcf_xml"))
    isfinite(dt) && dt > 0 || throw(ArgumentError("dt must be positive and finite, got $dt"))
    integrator === :rk4 && throw(ArgumentError("integrator :rk4 is not supported: the scene steps with mj_step1/mj_step2, and MuJoCo integrates RK4 models with Euler on that path. Use :implicitfast, :implicit or :euler."))
    haskey(_INTEGRATORS, integrator) || throw(ArgumentError("integrator must be one of $(keys(_INTEGRATORS)), got :$integrator"))
    _check_effectors(gravity_effectors, planet_rotation)
    model = mjcf_path === nothing ? Binding.load_xml_string(mjcf_xml) : Binding.load_xml_file(String(mjcf_path))
    return _assemble(model, Float64(dt), integrator, planet, gravity_effectors, planet_rotation, Float64(t0), collect(initial_states))
end

function _check_effectors(effectors::Tuple, planet_rotation)
    isempty(effectors) && throw(ArgumentError("gravity_effectors is empty"))
    all(e -> SM.gravity_backbone_structure(e) === :position_only_static_gravity, effectors) || throw(ArgumentError(
        "every gravity effector must declare gravity_backbone_structure == :position_only_static_gravity"))
    _needs_frame(effectors) && planet_rotation === nothing && throw(ArgumentError(
        "a gravity effector needs the planet frame (J2/harmonics): pass planet_rotation(t) -> inertial-to-planet-fixed rotation"))
    return nothing
end

function _assemble(model, dt, integrator, planet, effectors, planet_rotation, t0, initial_states)
    B = Binding
    B.set_timestep!(model, dt)
    B.set_integrator!(model, _INTEGRATORS[integrator])
    B.set_gravity!(model, (0.0, 0.0, 0.0))
    data = B.MjData(model)
    nbody = B.nbody(model)
    nbody >= 2 || throw(ArgumentError("the MJCF scene has no bodies"))
    mass = B.body_mass(model)
    nb = nbody - 1
    total = sum(@view mass[2:end])
    total > 0 || throw(ArgumentError("the scene has zero total mass"))
    rootid = B.body_rootid(model); dofnum = B.body_dofnum(model); jnum = B.body_jntnum(model)
    jadr = B.body_jntadr(model); jtype = B.jnt_type(model)
    roots = Int[]
    names = String[B.id2name(model, B.OBJ_BODY, i) for i in 1:nb]
    for i in 1:nb
        rootid[i + 1] == i || continue
        (jnum[i + 1] == 1 && jtype[jadr[i + 1] + 1] == B.JNT_FREE) || throw(ArgumentError(
            "root body '$(names[i])' must have exactly one free joint; bodies fixed to the world are not supported"))
        push!(roots, i)
    end
    scene = ProximityScene(model, data, dt, integrator, planet, effectors, planet_rotation, t0, SceneBodyState[],
        0, zero(V3), zero(V3), false, nb, names, roots, mass, total,
        B.qpos(model, data), B.qvel(model, data), B.ctrl(model, data), B.xfrc_applied(model, data),
        B.xipos(model, data), B.xquat(model, data), zeros(6))
    scene_reset!(scene, initial_states)
    return scene
end

"""
    scene_reset!(scene[, initial_states])

Return the scene to step zero. Without `initial_states` the states given at construction are reused, so
repeated resets reproduce the same trajectory bit for bit. Passing new states replaces them (the hook a
vectorized runner uses to randomize episodes without rebuilding models).
"""
function scene_reset!(scene::ProximityScene, states::AbstractVector{SceneBodyState}=scene.initial_states)
    B = Binding
    m, d = scene.model, scene.data
    by_body = Dict{Int, SceneBodyState}()
    for s in states
        id = B.name2id(m, B.OBJ_BODY, s.body)
        id in scene.roots || throw(ArgumentError("'$(s.body)' is not a free-root body of the scene (roots: $(scene.names[scene.roots]))"))
        haskey(by_body, id) && throw(ArgumentError("duplicate initial state for body '$(s.body)'"))
        by_body[id] = s
    end
    length(by_body) == length(scene.roots) || throw(ArgumentError(
        "initial states are required for every free-root body: $(scene.names[scene.roots])"))
    B.reset_data!(m, d)
    jadr = B.body_jntadr(m); qadr = B.jnt_qposadr(m); dadr = B.jnt_dofadr(m)
    first_state = by_body[first(scene.roots)]
    r_ref, v_ref = first_state.r, first_state.v
    # Stage the free joints around a nearby reference so the large inertial coordinates never enter MuJoCo.
    for id in scene.roots
        s = by_body[id]; j = jadr[id + 1] + 1; qa = qadr[j]; da = dadr[j]
        for k in 1:3
            scene.qpos[qa + k] = s.r[k] - r_ref[k]
            scene.qvel[da + k] = s.v[k] - v_ref[k]
            scene.qvel[da + 3 + k] = s.ω[k]
        end
        for k in 1:4
            scene.qpos[qa + 3 + k] = s.q[k]
        end
    end
    B.forward!(m, d)
    # Body origin versus center of mass: correct both exactly (the offset depends on attitude only).
    for id in scene.roots
        s = by_body[id]; j = jadr[id + 1] + 1; qa = qadr[j]; da = dadr[j]
        B.body_velocity!(scene.vel6, m, d, id)
        for k in 1:3
            scene.qpos[qa + k] += (s.r[k] - r_ref[k]) - scene.xipos[k, id + 1]
            scene.qvel[da + k] += (s.v[k] - v_ref[k]) - scene.vel6[3 + k]
        end
    end
    B.forward!(m, d)
    com, vcom = _com(scene)
    for id in scene.roots
        j = jadr[id + 1] + 1; qa = qadr[j]; da = dadr[j]
        for k in 1:3
            scene.qpos[qa + k] -= com[k]
            scene.qvel[da + k] -= vcom[k]
        end
    end
    B.forward!(m, d)
    scene.R = r_ref + com
    scene.V = v_ref + vcom
    scene.initial_states = collect(states)
    scene.n = 0
    scene.fresh = true
    return scene
end

# Mass-weighted center of mass and its velocity in the scene frame (kinematics must be fresh).
function _com(scene::ProximityScene)
    B = Binding
    c = zero(V3); v = zero(V3)
    for i in 1:scene.nb
        mi = scene.mass[i + 1]
        c += mi * V3(scene.xipos[1, i + 1], scene.xipos[2, i + 1], scene.xipos[3, i + 1])
        B.body_velocity!(scene.vel6, scene.model, scene.data, i)
        v += mi * V3(scene.vel6[4], scene.vel6[5], scene.vel6[6])
    end
    return c / scene.total_mass, v / scene.total_mass
end

"""Scene time `t0 + n*dt` [s], from the integer step counter (no accumulated rounding)."""
scene_time(scene::ProximityScene) = scene.t0 + scene.n * scene.dt

"""Chief position and velocity `(R, V)` in the inertial frame [m, m/s]."""
scene_chief(scene::ProximityScene) = (scene.R, scene.V)

"""MJCF names of the scene bodies; body `i` is `scene_body_names(scene)[i]` (the world body is excluded)."""
scene_body_names(scene::ProximityScene) = copy(scene.names)

"""View of MuJoCo's `ctrl` vector, for actuator commands."""
scene_ctrl(scene::ProximityScene) = scene.ctrl

function _refresh!(scene::ProximityScene)
    scene.fresh || (Binding.forward!(scene.model, scene.data); scene.fresh = true)
    return nothing
end

"""
    scene_body_state(scene, i) -> (r, v, q, ω)

Absolute inertial state of body `i`: COM position `R + xipos`, COM velocity `V + v`, attitude quaternion
`(w, x, y, z)` body to world, and angular velocity in world axes.
"""
function scene_body_state(scene::ProximityScene, i::Integer)
    1 <= i <= scene.nb || throw(BoundsError(scene.names, i))
    _refresh!(scene)
    Binding.body_velocity!(scene.vel6, scene.model, scene.data, i)
    r = scene.R + V3(scene.xipos[1, i + 1], scene.xipos[2, i + 1], scene.xipos[3, i + 1])
    v = scene.V + V3(scene.vel6[4], scene.vel6[5], scene.vel6[6])
    q = SVector{4, Float64}(scene.xquat[1, i + 1], scene.xquat[2, i + 1], scene.xquat[3, i + 1], scene.xquat[4, i + 1])
    return (r = r, v = v, q = q, ω = V3(scene.vel6[1], scene.vel6[2], scene.vel6[3]))
end
function scene_body_state(scene::ProximityScene, name::AbstractString)
    i = findfirst(==(name), scene.names)
    i === nothing && throw(ArgumentError("no body named '$name'"))
    return scene_body_state(scene, i)
end

# --- gravity through SpaceAGORA's hooks --------------------------------------------------------------

@inline _needs_frame(::Tuple{}) = false
@inline _needs_frame(e::Tuple) = SM.environment_requirements(first(e)).planet_frame || _needs_frame(Base.tail(e))

@inline _g_sum(::Tuple{}, x, env, t) = zero(V3)
@inline _g_sum(e::Tuple, x, env, t) = SM.gravity_backbone_acceleration_ii(first(e), x, env, t) + _g_sum(Base.tail(e), x, env, t)

@inline _dg_sum(::Tuple{}, xb, xf, eb, ef, ρ, t) = zero(V3)
@inline _dg_sum(e::Tuple, xb, xf, eb, ef, ρ, t) =
    SM.gravity_backbone_relative_acceleration_ii(first(e), xb, ρ, xf, eb, ef, t) + _dg_sum(Base.tail(e), xb, xf, eb, ef, ρ, t)

@inline _lpi(scene::ProximityScene, t) = _needs_frame(scene.effectors) ? scene.planet_rotation(t) : nothing

@inline function _env(scene::ProximityScene, r::V3, l_pi)
    frame = l_pi === nothing ? nothing : SE.sample_planet_frame_with_lpi((pos_ii = r, vel_ii = zero(V3)), scene.planet, l_pi)
    return SM.EnvironmentSample(scene.planet; planet_frame = frame)
end

@inline _sample(r::V3) = SM.StateSample(r, zero(V3), 1.0)

@inline function _gravity(scene::ProximityScene, r::V3, t::Float64, l_pi)
    return _g_sum(scene.effectors, _sample(r), _env(scene, r, l_pi), t)
end

# Fixed-step RK4 of the chief under gravity + a_fb (held over the step).
function _advance_chief!(scene::ProximityScene, afb::V3)
    h = scene.dt; t = scene_time(scene)
    R = scene.R; V = scene.V
    l0 = _lpi(scene, t); lh = _lpi(scene, t + h / 2); l1 = _lpi(scene, t + h)
    k1v = _gravity(scene, R, t, l0) + afb
    k2r = V + (h / 2) * k1v
    k2v = _gravity(scene, R + (h / 2) * V, t + h / 2, lh) + afb
    k3r = V + (h / 2) * k2v
    k3v = _gravity(scene, R + (h / 2) * k2r, t + h / 2, lh) + afb
    k4r = V + h * k3v
    k4v = _gravity(scene, R + h * k3r, t + h, l1) + afb
    scene.R = R + (h / 6) * (V + 2k2r + 2k3r + k4r)
    scene.V = V + (h / 6) * (k1v + 2k2v + 2k3v + k4v)
    return nothing
end

const _UNSTABLE = (Binding.WARN_BADQACC, Binding.WARN_BADQPOS, Binding.WARN_BADQVEL)
const _UNSTABLE_NAMES = ("mjWARN_BADQACC", "mjWARN_BADQPOS", "mjWARN_BADQVEL")

"""
    scene_step!(scene; external_forces=nothing)

Advance the scene one `dt`. `external_forces` is an optional `3 x nbodies` matrix of non-gravitational
forces [N] in inertial axes, applied at each body's COM; their mass average is the feedback acceleration
`a_fb` given to the chief (and removed from every body), as in mjorbit.

Order: `mj_step1` (fresh kinematics) -> wrenches from those kinematics and the chief position `R_n` ->
`xfrc_applied` -> `mj_step2` -> chief RK4 from `R_n` to `R_{n+1}` -> `n += 1`.
"""
function scene_step!(scene::ProximityScene; external_forces::Union{Nothing, AbstractMatrix{<:Real}}=nothing)
    B = Binding
    nb = scene.nb
    if external_forces !== nothing
        size(external_forces) == (3, nb) || throw(DimensionMismatch("external_forces must be 3 x $nb"))
    end
    before = map(w -> B.warning_count(scene.data, w), _UNSTABLE)
    B.step1!(scene.model, scene.data)
    afb = zero(V3)
    if external_forces !== nothing
        for i in 1:nb
            afb += V3(external_forces[1, i], external_forces[2, i], external_forces[3, i])
        end
        afb /= scene.total_mass
    end
    t = scene_time(scene)
    R = scene.R
    l_pi = _lpi(scene, t)
    xb = _sample(R); eb = _env(scene, R, l_pi)
    for i in 1:nb
        mi = scene.mass[i + 1]
        ρ = V3(scene.xipos[1, i + 1], scene.xipos[2, i + 1], scene.xipos[3, i + 1])
        rf = R + ρ
        dg = _dg_sum(scene.effectors, xb, _sample(rf), eb, _env(scene, rf, l_pi), ρ, t)
        f = mi * (dg - afb)
        if external_forces !== nothing
            f += V3(external_forces[1, i], external_forces[2, i], external_forces[3, i])
        end
        scene.xfrc[1, i + 1] = f[1]; scene.xfrc[2, i + 1] = f[2]; scene.xfrc[3, i + 1] = f[3]
        scene.xfrc[4, i + 1] = 0.0; scene.xfrc[5, i + 1] = 0.0; scene.xfrc[6, i + 1] = 0.0
    end
    B.step2!(scene.model, scene.data)
    for (k, w) in enumerate(_UNSTABLE)
        # MuJoCo has already reset its state on a bad acceleration; continuing would run on reset data.
        B.warning_count(scene.data, w) == before[k] || error(
            "MuJoCo reported an unstable simulation ($(_UNSTABLE_NAMES[k])) at scene time $(t) s, step $(scene.n + 1); ",
            "the scene state was reset by MuJoCo and is no longer valid")
    end
    _advance_chief!(scene, afb)
    scene.n += 1
    scene.fresh = false
    return scene
end

"""Snapshot the restorable state of the scene as a [`SceneState`](@ref)."""
function scene_state(scene::ProximityScene)
    buf = Vector{Float64}(undef, Binding.state_size(scene.model))
    Binding.get_state!(buf, scene.model, scene.data)
    return SceneState(scene.n, scene.R, scene.V, buf)
end

"""Restore a [`SceneState`](@ref) taken from a scene of the same MJCF model."""
function scene_set_state!(scene::ProximityScene, st::SceneState)
    length(st.mj) == Binding.state_size(scene.model) || throw(DimensionMismatch(
        "state has $(length(st.mj)) MuJoCo entries, the scene needs $(Binding.state_size(scene.model))"))
    Binding.set_state!(scene.model, scene.data, copy(st.mj))
    scene.n = st.n; scene.R = st.R; scene.V = st.V
    scene.fresh = false
    return scene
end

"""An independent scene with its own copied `mjModel` and fresh `mjData`, at the same state."""
function Base.copy(scene::ProximityScene)
    other = _assemble(Binding.copy_model(scene.model), scene.dt, scene.integrator, scene.planet, scene.effectors,
        scene.planet_rotation, scene.t0, scene.initial_states)
    scene_set_state!(other, scene_state(scene))
    return other
end
