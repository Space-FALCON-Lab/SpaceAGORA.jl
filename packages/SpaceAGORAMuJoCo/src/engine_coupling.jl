# Stage 2: a ProximityScene as the external propagator of spacecraft inside `run_simulation`.
#
# The scene owns the listed spacecraft for the whole run. The engine keeps them in `u.sc` as shadow entries
# (translation follows the chief acceleration between syncs; attitude advances by plain kinematics), steps the
# scene by whole `dt` after every accepted solver step through `external_sync!`, and overwrites the shadow
# entries from the absolute scene state. Wrenches beyond differential gravity (thrusters, wheels, joint
# motors, environment loads) are not applied yet; a control effector on an owned spacecraft is refused.

const EP = SM.ExternalPropagation

"""
    ProximitySceneDynamics(scene, spacecraft_index => "body", ...; mass_rtol=nothing)

Hand the listed spacecraft (1-based indices into `dynamics_model.spacecraft`) to `scene`: each is paired with
the MJCF free-root body that carries it, and every free-root body of the scene must be paired. Pass it in
`SimulationConfiguration.external_propagators`, for example
`SimConfig._with_configuration(config; external_propagators = (ProximitySceneDynamics(scene, 1 => "chaser", 2 => "target"),))`.

`scene` is a template: every run copies it (its own `mjModel` and `mjData`) and starts the copy from the
initial state of the paired spacecraft, so the template is never stepped, can be shared by concurrent runs,
and `deepcopy` of a configuration copies the native model rather than its pointers. The scene's own
`initial_states` are therefore replaced by the spacecraft initial conditions. The scene's gravity effectors and
planet drive the owned bodies and should match the run's dynamics; the engine's dynamic effectors do not act on
owned spacecraft. `mass_rtol`, when given, refuses a run whose spacecraft mass differs from the MJCF body mass
by more than that relative tolerance (no check otherwise).
"""
struct ProximitySceneDynamics{S <: ProximityScene} <: EP.AbstractExternalPropagator
    scene::S
    spacecraft::Vector{Int}
    bodies::Vector{String}
    mass_rtol::Union{Nothing, Float64}
end

function ProximitySceneDynamics(scene::ProximityScene, pairs::Pair{<:Integer, <:AbstractString}...; mass_rtol::Union{Nothing, Real}=nothing)
    isempty(pairs) && throw(ArgumentError("ProximitySceneDynamics needs at least one spacecraft => body pair"))
    sc = Int[first(p) for p in pairs]
    bodies = String[last(p) for p in pairs]
    allunique(sc) || throw(ArgumentError("duplicate spacecraft index in $sc"))
    allunique(bodies) || throw(ArgumentError("a body can carry only one spacecraft; got $bodies"))
    all(>=(1), sc) || throw(ArgumentError("spacecraft indices are 1-based; got $sc"))
    for b in bodies
        id = findfirst(==(b), scene.names)
        id === nothing && throw(ArgumentError("the scene has no body named '$b' (bodies: $(scene.names))"))
        id in scene.roots || throw(ArgumentError("body '$b' is not a free-root body; only free-root bodies can carry a spacecraft (articulated links stay in the scene)"))
    end
    unpaired = setdiff(scene.names[scene.roots], bodies)
    isempty(unpaired) || throw(ArgumentError("free-root bodies without a spacecraft: $unpaired; the scene needs an initial state for every free-root body"))
    mass_rtol === nothing || (isfinite(mass_rtol) && mass_rtol >= 0) || throw(ArgumentError("mass_rtol must be nonnegative and finite"))
    return ProximitySceneDynamics(scene, sc, bodies, mass_rtol === nothing ? nothing : Float64(mass_rtol))
end

"""Run-owned state of a [`ProximitySceneDynamics`](@ref): the scene copy and the chief acceleration."""
mutable struct ProximitySceneRuntime
    scene::ProximityScene
    body_ids::Vector{Int}
    names::Vector{String}
    l_pi::Any                      # inertial-to-planet rotation at the scene's current time (nothing when unused)
end

EP.external_spacecraft(ep::ProximitySceneDynamics) = copy(ep.spacecraft)
EP.external_step(ep::ProximitySceneDynamics) = ep.scene.dt

function EP.external_preflight(ep::ProximitySceneDynamics, args)
    ep.mass_rtol === nothing && return nothing
    for (i, b) in zip(ep.spacecraft, ep.bodies)
        sc = args.dynamics_model.spacecraft[i]
        m_sc = sc.dry_mass + sc.prop_mass
        m_mj = ep.scene.mass[findfirst(==(b), ep.scene.names) + 1]
        abs(m_sc - m_mj) <= ep.mass_rtol * m_mj || throw(ArgumentError(
            "spacecraft $i has mass $m_sc kg but MJCF body '$b' has $m_mj kg (mass_rtol = $(ep.mass_rtol))"))
    end
    return nothing
end

# A scene copy with its own native model, ready to be reset to the run's initial state.
_fresh_copy(s::ProximityScene) = _assemble(Binding.copy_model(s.model), s.dt, s.integrator, s.planet, s.effectors,
    s.planet_rotation, s.t0, s.initial_states)

function EP.external_prepare(ep::ProximitySceneDynamics, args, u0)
    scene = _fresh_copy(ep.scene)
    states = SceneBodyState[]
    for (i, b) in zip(ep.spacecraft, ep.bodies)
        sc = u0.sc[i]
        q = hasproperty(sc, :q) ? sa_to_mujoco_quaternion(sc.q) : SVector(1.0, 0.0, 0.0, 0.0)
        ω = hasproperty(sc, :ω) ? V3(sc.ω) : zero(V3)
        push!(states, SceneBodyState(b, V3(sc.pos), V3(sc.vel); q, ω))
    end
    scene_reset!(scene, states)
    return ProximitySceneRuntime(scene, Int[findfirst(==(b), scene.names) for b in ep.bodies], copy(ep.bodies), _lpi(scene, scene_time(scene)))
end

function EP.external_sync!(rt::ProximitySceneRuntime, t::Float64)
    scene = rt.scene
    target = floor(Int, t / scene.dt + 1.0e-6)      # whole scene steps within t, tolerant of tstop rounding
    n = target - scene.n
    n > 0 || return 0
    for _ in 1:n
        scene_step!(scene)
    end
    rt.l_pi = _lpi(scene, scene_time(scene))
    return n
end

EP.external_time(rt::ProximitySceneRuntime) = scene_time(rt.scene)

# Hamilton product of scalar-last quaternions (x, y, z, w).
@inline function _quat_mul(a::SVector{4, Float64}, b::SVector{4, Float64})
    av = V3(a[1], a[2], a[3]); bv = V3(b[1], b[2], b[3])
    v = a[4] * bv + b[4] * av + cross(av, bv)
    return SVector{4, Float64}(v[1], v[2], v[3], a[4] * b[4] - dot(av, bv))
end

# SpaceAGORA's gravity at the chief, displaced to engine time t = (scene time) + τ by the chief's own velocity.
# The planet-frame rotation is the one at the scene's last step: it turns 7e-5 rad/s, so a stale value over a
# scene step perturbs the acceleration by far less than the chief displacement it replaces.
@inline function _chief_gravity(rt::ProximitySceneRuntime, τ::Float64, lag::Float64)
    sc = rt.scene
    return _gravity(sc, sc.R + sc.V * lag, scene_time(sc) + τ, rt.l_pi)
end

# The scene lags engine time `t` by τ in [0, dt). The shadow entry is written at `t`, so the scene state is
# advanced by τ with the shadow's own model: the chief's acceleration, taken at the midpoint of the interval, for
# translation, and constant body rate for attitude (q (x) exp(ω τ), the kinematics the engine integrates). At a
# scene-step boundary τ = 0 and the scene state is written unchanged (a remainder below 1e-9 dt is rounding of a
# tstop, not a lag).
function EP.external_state(rt::ProximitySceneRuntime, k::Int, t::Float64)
    s = scene_body_state(rt.scene, rt.body_ids[k])
    τ = t - rt.scene.n * rt.scene.dt
    q = mujoco_to_sa_quaternion(s.q)
    τ > 1.0e-9 * rt.scene.dt || return (pos = s.r, vel = s.v, q = q, ω = s.ω)
    a = _chief_gravity(rt, τ / 2, τ / 2)
    φ = norm(s.ω) * τ
    dq = φ > 0 ? SVector{4, Float64}((sin(φ / 2) / norm(s.ω)) * s.ω..., cos(φ / 2)) : SVector(0.0, 0.0, 0.0, 1.0)
    return (pos = s.r + s.v * τ + a * (τ^2 / 2), vel = s.v + a * τ, q = _quat_mul(q, dq), ω = s.ω)
end

# The chief's gravity at engine time t (no feedback acceleration yet), evaluated once per right-hand-side call.
# ponytail: one gravity evaluation per owned spacecraft per RHS call; share one per owner if a harmonics field makes it matter.
function EP.external_acceleration(rt::ProximitySceneRuntime, ::Int, t::Float64)
    τ = max(t - rt.scene.n * rt.scene.dt, 0.0)
    return _chief_gravity(rt, τ, τ)
end

"""
    scene_body_pose_save_field(name) -> SaveField

A save field `scene_pose_<name>` with the absolute pose of the MJCF body `name` (for example an arm link that
no spacecraft carries): inertial COM position (3) and the SpaceAGORA attitude quaternion (4, scalar-last, body
to inertial), written as columns `scene_pose_<name>_1..7`. The value is the pose at the last completed scene
step at or before the save time.
"""
function scene_body_pose_save_field(name::AbstractString)
    nm = String(name)
    getter = function (u, t, integrator)
        for rt in integrator.p.shared_buffers.external_runtimes
            rt isa ProximitySceneRuntime || continue
            i = findfirst(==(nm), rt.scene.names)
            i === nothing && continue
            s = scene_body_state(rt.scene, i)
            return Float64[s.r[1], s.r[2], s.r[3], mujoco_to_sa_quaternion(s.q)...]
        end
        throw(ArgumentError("no scene body named '$nm' in the run's external propagators"))
    end
    return SM.SaveField(Symbol("scene_pose_", nm), getter; per_satellite = false, column_prefix = "scene_pose_" * nm)
end
