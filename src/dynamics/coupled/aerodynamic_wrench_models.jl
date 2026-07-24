using SpecialFunctions

const sqrt_π = sqrt(π)
const inv_sqrt_π = 1 / sqrt(π)

@inline _parse_bool_env(name::String, default::Bool)::Bool = ParallelPolicy.parse_bool_env(name, default)

@inline function _multibody_parallel_mode()::Symbol
    return ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_MULTIBODY_PARALLEL")
end

@inline function _multibody_thread_threshold()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_MULTIBODY_THREAD_THRESHOLD", 4)
end

@inline function _multibody_max_threads()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_MULTIBODY_MAX_THREADS", 4)
end

@inline function _threadid_capacity()::Int
    # Compatibility shim for legacy tests/callers after migrating scratch indexing
    # from thread ids to stable worker ids.
    return max(Threads.maxthreadid(), _multibody_max_threads())
end

@inline function _multibody_outer_parallel_hint()::Bool
    # Explicit hint for disabling intra-satellite threading under outer parallel execution.
    return ParallelPolicy.outer_parallel_active()
end

@inline function _multibody_use_threads(num_items::Int; heavy_work::Bool=true)::Bool
    return _multibody_thread_decision(num_items; heavy_work=heavy_work).use_threads
end

@inline function _multibody_thread_decision(
    num_items::Int;
    heavy_work::Bool=true,
    env::Union{Nothing,ParallelPolicy.PolicyDecisionEnvConfig}=nothing,
)
    if env !== nothing && env.inner_thread_budget <= 1
        return (
            use_threads=false,
            allotment=1,
            mode=:off,
            policy_applied=false,
        )
    end
    mode = _multibody_parallel_mode()
    threshold = _multibody_thread_threshold()
    outer_active = env === nothing ? _multibody_outer_parallel_hint() : env.outer_parallel_active
    allow_with_outer = _parse_bool_env("SPACEAGORA_MULTIBODY_PARALLEL_ALLOW_WITH_OUTER", false)
    heavy_only = _parse_bool_env("SPACEAGORA_MULTIBODY_PARALLEL_HEAVY_ONLY", true)
    policy = ParallelPolicy.thread_policy_decision(
        num_items;
        mode=mode,
        threshold=threshold,
        heavy_work=heavy_work,
        heavy_only=heavy_only,
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=:multibody,
        env=env,
    )
    allotment = min(policy.allotment, _multibody_max_threads())
    use_threads = policy.use_threads && allotment > 1
    return (
        use_threads=use_threads,
        allotment=use_threads ? allotment : 1,
        mode=mode,
        policy_applied=true,
    )
end

# Aerodynamic models
@kwdef struct AerodynamicCoefficientConstant <: AbstractForceTorqueModel

end

@kwdef struct AerodynamicCoefficientfM <: AbstractForceTorqueModel
    # Per-model opt-in for per-link atmosphere sampling (preferred over the
    # process-wide set_per_link_atmosphere!, which leaks across simulations).
    per_link_atmosphere::Bool = false
end

@kwdef struct AerodynamicCoefficientNoBallisticFlight <: AbstractForceTorqueModel
    # Per-model opt-in for per-link atmosphere sampling (preferred over the
    # process-wide set_per_link_atmosphere!, which leaks across simulations).
    per_link_atmosphere::Bool = false
end

@inline function _make_aero_scratch_workspace(n_threads::Int)::AeroScratchWorkspace
    n_threads >= 1 || throw(ArgumentError("Aerodynamic scratch workspace thread count must be >= 1, got $n_threads"))
    thread_force = [MVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n_threads]
    thread_cl = zeros(Float64, n_threads)
    thread_cd = zeros(Float64, n_threads)
    thread_area = zeros(Float64, n_threads)
    return AeroScratchWorkspace(thread_force, thread_cl, thread_cd, thread_area)
end

@inline function _ensure_aero_workspace_capacity!(
    workspace::AeroScratchWorkspace,
    n_threads::Int
)::AeroScratchWorkspace
    n_threads >= 1 || throw(ArgumentError("Aerodynamic scratch workspace thread count must be >= 1, got $n_threads"))
    if length(workspace.thread_force) < n_threads
        old_len = length(workspace.thread_force)
        resize!(workspace.thread_force, n_threads)
        @inbounds for tid in (old_len + 1):n_threads
            workspace.thread_force[tid] = MVector{3, Float64}(0.0, 0.0, 0.0)
        end
    end
    if length(workspace.thread_cl) < n_threads
        resize!(workspace.thread_cl, n_threads)
        resize!(workspace.thread_cd, n_threads)
        resize!(workspace.thread_area, n_threads)
    end
    return workspace
end

@inline function _aero_workspace_for_sat!(
    param::ODEParams,
    sat_idx::Int,
    n_threads::Int
)::AeroScratchWorkspace
    workspaces = param.shared_buffers.aero_workspaces
    if sat_idx > length(workspaces)
        return _ensure_aero_workspace_capacity!(_make_aero_scratch_workspace(n_threads), n_threads)
    end

    sat_entry = @inbounds workspaces[sat_idx]
    workspace = if sat_entry === nothing
        created = _make_aero_scratch_workspace(n_threads)
        @inbounds workspaces[sat_idx] = created
        created
    else
        sat_entry::AeroScratchWorkspace
    end
    _ensure_aero_workspace_capacity!(workspace, n_threads)
    return workspace
end

@inline function _simulation_model_module_for_aero()
    mod = @__MODULE__
    while true
        if isdefined(mod, :planet_frame_lpi)
            return mod
        end
        parent = parentmodule(mod)
        parent === mod && break
        mod = parent
    end
    error("SimulationModel.planet_frame_lpi not found in module ancestry for aerodynamic_wrench_models.jl")
end

@inline function _store_vector_cache!(
    cache::Vector{SVector{3, Float64}},
    sat_idx::Int,
    value::SVector{3, Float64}
)::Nothing
    if length(cache) < sat_idx
        resize!(cache, sat_idx)
        @inbounds for idx in eachindex(cache)
            if !isassigned(cache, idx)
                cache[idx] = SVector{3, Float64}(0.0, 0.0, 0.0)
            end
        end
    end
    @inbounds cache[sat_idx] = value
    return nothing
end

@inline function _store_aero_caches!(
    param::ODEParams,
    sat_idx::Int,
    drag_ii::SVector{3, Float64},
    lift_ii::SVector{3, Float64},
    cross_ii::SVector{3, Float64}
)::Nothing
    _store_vector_cache!(param.save_cache.drag_cache, sat_idx, drag_ii)
    _store_vector_cache!(param.save_cache.lift_cache, sat_idx, lift_ii)
    _store_vector_cache!(param.save_cache.cross_cache, sat_idx, cross_ii)
    return nothing
end

function collect_and_reset_link_wrenches!(bodies)
    # Collect into fresh vectors to avoid aliasing when there is only one link.
    force_acc = MVector{3, Float64}(0.0, 0.0, 0.0)
    torque_acc = MVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for b in bodies
        force_acc .+= SVector{3, Float64}(b.net_force)
        torque_acc .+= SVector{3, Float64}(b.net_torque)
        b.net_force .= SVector{3, Float64}(0.0, 0.0, 0.0)
        b.net_torque .= SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    return SVector{3, Float64}(force_acc), SVector{3, Float64}(torque_acc)
end

@inline environment_requirements(::AerodynamicCoefficientConstant) = EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
@inline environment_requirements(::AerodynamicCoefficientfM) = EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
@inline environment_requirements(::AerodynamicCoefficientNoBallisticFlight) = EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
@inline solver_partition(::AerodynamicCoefficientConstant) = :implicit
@inline solver_partition(::AerodynamicCoefficientfM) = :implicit
@inline solver_partition(::AerodynamicCoefficientNoBallisticFlight) = :implicit

@inline function _constant_drag_coefficient(alpha_rad::Float64)::Float64
    return 2 * (2.2 - 0.8) / pi * alpha_rad + 0.8
end

@inline function _aero_link_angles(
    spacecraft,
    body,
    root_index::Int,
    orientation_sim::Bool,
    vel_pi::SVector{3, Float64},
    θ_body::Float64,
)::Tuple{Float64, Float64, Union{Nothing, SMatrix{3, 3, Float64, 9}}}
    if orientation_sim
        R = rotate_to_inertial(spacecraft, body, root_index)
        body_frame_velocity = R' * vel_pi
        α_body = atan(body_frame_velocity[1], body_frame_velocity[3])
        β_body = atan(body_frame_velocity[2], hypot(body_frame_velocity[1], body_frame_velocity[3]))
        return α_body, β_body, SMatrix{3, 3, Float64, 9}(R)
    end
    α_body = body.root ? (pi / 2) : atan((rot(body.q) * SVector{3, Float64}(1.0, 0.0, 0.0))[1], (rot(body.q) * SVector{3, Float64}(1.0, 0.0, 0.0))[3])
    return α_body, 0.0, nothing
end

const _AERO_ZERO3 = SVector{3, Float64}(0.0, 0.0, 0.0)
const _AERO_ZERO5 = (_AERO_ZERO3, _AERO_ZERO3, _AERO_ZERO3, _AERO_ZERO3, _AERO_ZERO3)

# Opt-in switch for per-link atmosphere sampling in the ODEParams-aware `wrench`
# methods below. Default `false` preserves the long-standing single-sample
# behavior (one atmosphere sample for the whole spacecraft), so default
# simulation physics is identical to the pre-feature code path; callers that
# want the per-link treatment (e.g. the CYGNSS attitude/torque reconstruction
# study) enable it explicitly via `set_per_link_atmosphere!(true)`.
const PER_LINK_ATMOSPHERE_ENABLED = Ref(false)
# NOTE: process-wide switch, kept for compatibility; it affects every
# aerodynamic effector in the Julia process, including concurrent or later
# simulations. Prefer the per-model field
# `AerodynamicCoefficientfM(per_link_atmosphere=true)` (idem
# NoBallisticFlight), which scopes the behavior to one effector instance.
set_per_link_atmosphere!(flag::Bool) = (PER_LINK_ATMOSPHERE_ENABLED[] = flag; nothing)

@inline _per_link_enabled(model)::Bool =
    (hasfield(typeof(model), :per_link_atmosphere) && model.per_link_atmosphere) ||
    PER_LINK_ATMOSPHERE_ENABLED[]


# Returns (force_ii, torque_body, drag_ii, lift_ii, cross_ii) all in the inertial frame.
#
# `link_atmosphere_fn`, when provided, is called as `link_atmosphere_fn(pos_pp_link)`
# for every non-root link and must return `(rho_kg_m3, temperature_k, wind_pp)` for
# that link's own planet-fixed position, rather than reusing the single spacecraft-level
# sample for every link. Panel offsets are normally tiny compared to atmospheric density
# scale heights, so this only matters for larger multi-body spacecraft or coarse density
# fields, but it is the physically correct per-body treatment. Pass `nothing` (the
# default) to keep the original single-sample-for-the-whole-spacecraft behavior, which
# is all that's available to callers without access to the live density model (i.e. the
# plain `wrench` methods below, which don't receive `p`).
function _aero_pure_wrench(
    coefficient_mode::Symbol,
    x::StateSample,
    env::EnvironmentSample,
    link_atmosphere_fn=nothing,
)::NTuple{5, SVector{3, Float64}}
    x.spacecraft === nothing && throw(ArgumentError("Aerodynamic wrench evaluation requires StateSample.spacecraft."))
    planet_frame = env.planet_frame
    atmosphere = env.atmosphere
    planet_frame === nothing && throw(ArgumentError("Aerodynamic wrench evaluation requires env.planet_frame."))
    atmosphere === nothing && throw(ArgumentError("Aerodynamic wrench evaluation requires env.atmosphere."))

    rho = atmosphere.rho_kg_m3
    T = atmosphere.temperature_k
    wind = atmosphere.wind_pp
    if !isfinite(rho) || rho <= eps(Float64) || !isfinite(T) || T <= 0.0
        return _AERO_ZERO5
    end

    spacecraft = x.spacecraft
    bodies = spacecraft.links
    root_index = 1
    orientation_sim = x.q_ib !== nothing
    planet = env.planet

    vel_pp = planet_frame.vel_pp
    h_pp = cross(planet_frame.pos_pp, vel_pp)
    h_pp_mag = norm(h_pp)
    if !isfinite(h_pp_mag) || h_pp_mag <= eps(Float64)
        return _AERO_ZERO5
    end

    bank_angle = deg2rad(0.0)
    uD, uN, uE = latlongtoNED((planet_frame.alt_m, planet_frame.lat_rad, planet_frame.lon_rad))
    wE, wN, wU = wind
    wind_pp = wN * uN + wE * uE - wU * uD
    # Airspeed is spacecraft velocity minus the atmosphere's own velocity.
    vel_pp_rw = vel_pp - wind_pp
    vel_pp_rw_mag = norm(vel_pp_rw)
    if vel_pp_rw_mag <= eps(Float64)
        return _AERO_ZERO5
    end
    # Free-molecular coefficients use the same wind-relative flow as the force
    # direction and dynamic pressure (the legacy calcForceTorque path already
    # does); mach_body/S_body are computed per link inside the loop below.

    vel_pp_rw_hat = vel_pp_rw / vel_pp_rw_mag
    h_pp_hat = h_pp / h_pp_mag
    lift_pp_hat = normalize(cross(h_pp_hat, vel_pp_rw_hat))
    drag_pp_hat = -vel_pp_rw_hat
    cross_pp_hat = cross(drag_pp_hat, lift_pp_hat)
    l_pi_t = planet_frame.l_pi'
    vel_pi = orientation_sim ? l_pi_t * vel_pp_rw : SVector{3, Float64}(0.0, 0.0, 0.0)
    θ_body = acos(clamp(vel_pp_rw[1] / vel_pp_rw_mag, -1.0, 1.0))

    force_ii  = MVector{3, Float64}(0.0, 0.0, 0.0)
    torque_body = MVector{3, Float64}(0.0, 0.0, 0.0)
    drag_ii   = MVector{3, Float64}(0.0, 0.0, 0.0)
    lift_ii   = MVector{3, Float64}(0.0, 0.0, 0.0)
    cross_ii  = MVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for body in bodies
        α_body, β_body, R_body_to_inertial = _aero_link_angles(spacecraft, body, root_index, orientation_sim, vel_pi, θ_body)

        rho_body, T_body = rho, T
        if link_atmosphere_fn !== nothing && orientation_sim && !body.root && R_body_to_inertial !== nothing
            pos_ii_body = x.pos_ii + R_body_to_inertial * SVector{3, Float64}(body.r)
            pos_pp_body = planet_frame.l_pi * pos_ii_body
            rho_link, T_link, wind_link = link_atmosphere_fn(pos_pp_body)
            if isfinite(rho_link) && rho_link > eps(Float64) && isfinite(T_link) && T_link > 0.0
                rho_body, T_body = rho_link, T_link
                # Per-link sampling applies density/temperature only; a
                # link-local WIND sample is discarded here (Mach and dynamic
                # pressure keep the spacecraft-level wind-relative velocity).
                # Warn exactly when nonzero wind data is actually thrown away,
                # so the maxlog budget cannot be consumed by runs that never
                # sample a link (Codex review on PR #62).
                if wind_link !== nothing && any(!iszero, wind_link)
                    @warn "per-link atmosphere discards the link-local WIND sample: Mach/dynamic pressure use the spacecraft-level wind-relative velocity (per-link sampling covers density/temperature only)." maxlog = 1
                end
            end
        end
        sound_velocity_body = sqrt(planet.γ * planet.R * T_body)
        mach_body = vel_pp_rw_mag / sound_velocity_body
        S_body = sqrt(planet.γ * 0.5) * mach_body
        q_body = 0.5 * rho_body * vel_pp_rw_mag^2
        lift_scale_body = q_body * cos(bank_angle)

        CL_body, CD_body, CS_body = if coefficient_mode == :fm
            coeffs = aerodynamic_coefficient_fM(body, T_body, S_body, α_body, β_body, θ_body)
            coeffs[1], coeffs[2], coeffs[3]
        else
            0.0, _constant_drag_coefficient(α_body), 0.0
        end

        drag_pp_body = q_body * CD_body * body.ref_area * drag_pp_hat
        lift_pp_body = lift_scale_body * CL_body * body.ref_area * lift_pp_hat
        cross_pp_body = orientation_sim ? (q_body * CS_body * body.ref_area * cross_pp_hat) : SVector{3, Float64}(0.0, 0.0, 0.0)
        drag_ii_body  = l_pi_t * drag_pp_body
        lift_ii_body  = l_pi_t * lift_pp_body
        cross_ii_body = l_pi_t * cross_pp_body
        force_ii  .+= drag_ii_body + lift_ii_body + cross_ii_body
        drag_ii   .+= drag_ii_body
        lift_ii   .+= lift_ii_body
        cross_ii  .+= cross_ii_body
        if orientation_sim && !(R_body_to_inertial === nothing)
            force_body = SVector{3, Float64}(drag_ii_body + lift_ii_body + cross_ii_body)
            torque_body .+= cross(SVector{3, Float64}(body.r), R_body_to_inertial' * force_body)
        end
    end

    return (
        SVector{3, Float64}(force_ii),
        SVector{3, Float64}(torque_body),
        SVector{3, Float64}(drag_ii),
        SVector{3, Float64}(lift_ii),
        SVector{3, Float64}(cross_ii),
    )
end

# Direct (uncached) per-link density query used by `link_atmosphere_fn`. Bypasses the
# satellite-level GRAM tracking/extrapolation cache in
# SimulationCallbacks._density_state_from_kinematics! deliberately: that cache is keyed
# and warmed per satellite trajectory point, and reusing it for a handful of link-offset
# positions per RHS call would corrupt its trajectory-continuity assumptions. This costs
# one extra raw density-model call per non-root link per wrench evaluation.
@inline function _aero_link_atmosphere_query(p, sat_idx::Int, t::Float64, pos_pp_link::SVector{3, Float64}, planet)
    alt, lat, lon = rtolatlong(pos_pp_link, planet)
    density_model = SimulationModel.SimulationCallbacks._density_model_for_sat(p, sat_idx)
    return SimulationModel.getDensity(density_model, alt, lat, lon, t, true, p)
end

@inline function wrench(
    model::AerodynamicCoefficientConstant,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    force, torque, _, _, _ = _aero_pure_wrench(:constant, x, env)
    return force, torque
end

@inline function wrench(
    model::AerodynamicCoefficientfM,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    force, torque, _, _, _ = _aero_pure_wrench(:fm, x, env)
    return force, torque
end

@inline function wrench(
    model::AerodynamicCoefficientNoBallisticFlight,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    force, torque, _, _, _ = _aero_pure_wrench(:constant, x, env)
    return force, torque
end

@inline function wrench_caching!(
    model::AerodynamicCoefficientConstant,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
    p::ODEParams,
    sat_idx::Int,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    link_atmosphere_fn = _per_link_enabled(model) ?
        (pos_pp_body -> _aero_link_atmosphere_query(p, sat_idx, t, pos_pp_body, env.planet)) : nothing
    force, torque, drag_ii, lift_ii, cross_ii = _aero_pure_wrench(:constant, x, env, link_atmosphere_fn)
    _store_aero_caches!(p, sat_idx, drag_ii, lift_ii, cross_ii)
    return force, torque
end

@inline function wrench_caching!(
    model::AerodynamicCoefficientfM,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
    p::ODEParams,
    sat_idx::Int,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    link_atmosphere_fn = _per_link_enabled(model) ?
        (pos_pp_body -> _aero_link_atmosphere_query(p, sat_idx, t, pos_pp_body, env.planet)) : nothing
    force, torque, drag_ii, lift_ii, cross_ii = _aero_pure_wrench(:fm, x, env, link_atmosphere_fn)
    _store_aero_caches!(p, sat_idx, drag_ii, lift_ii, cross_ii)
    return force, torque
end

@inline function wrench_caching!(
    model::AerodynamicCoefficientNoBallisticFlight,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
    p::ODEParams,
    sat_idx::Int,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    link_atmosphere_fn = _per_link_enabled(model) ?
        (pos_pp_body -> _aero_link_atmosphere_query(p, sat_idx, t, pos_pp_body, env.planet)) : nothing
    force, torque, drag_ii, lift_ii, cross_ii = _aero_pure_wrench(:constant, x, env, link_atmosphere_fn)
    _store_aero_caches!(p, sat_idx, drag_ii, lift_ii, cross_ii)
    return force, torque
end

# Calculate force/torque functions
function calcForceTorque(model::AerodynamicCoefficientConstant, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    m = param.m
    cnf = param.cnf
    orientation_sim = param.orientation_sim

    bodies, root_index = traverse_bodies(m.body, m.body.roots[1]) # Get all bodies in the simulation

    pos_ii = SVector{3, Float64}(x[1], x[2], x[3])
    vel_ii = SVector{3, Float64}(x[4], x[5], x[6])

    h_ii = cross(pos_ii, vel_ii)    # Inertial angular momentum vector [m ^ 2 / s]

    h_ii_mag = norm(h_ii)           # Magnitude of the inertial angular momentum [m ^ 2 / s]

    # Inertial to planet relative transformation
    pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, m.planet, cnf.et) # Position vector planet / planet[m] # Velocity vector planet / planet[m / s]
    pos_pp_mag = norm(pos_pp) # Magnitude of the planet relative position

    vel_pp_mag = norm(vel_pp)

    h_pp = cross(pos_pp, vel_pp)
    
    h_pp_mag = norm(h_pp)
    h_pp_hat = normalize(h_pp) # Unit vector of the planet relative angular momentum
    
    bank_angle = deg2rad(0.0)

    lift_pp_hat = normalize(cross(h_pp_hat, vel_pp_rw_hat))
    # lift_pp_hat /= norm(lift_pp_hat) # Normalize the lift vector in planet relative frame
    drag_pp_hat = -vel_pp_rw_hat # Planet relative drag force direction
    cross_pp_hat = cross(drag_pp_hat, lift_pp_hat) # Cross product of the drag and lift vectors in planet relative frame

    if orientation_sim
        Rot = [MMatrix{3,3,Float64}(zeros(3, 3)) for i in eachindex(bodies)] # Rotation matrix from the root body to the spacecraft link
        @inbounds for (i, b) in enumerate(bodies)
            Rot[i] .= rotate_to_inertial(m.body, b, root_index) # Rotation matrix from the spacecraft link to the inertial frame
        end
    end
    
    CL, CD = 0.0, 0.0 # Initialize aerodynamic coefficients
    total_area = 0.0 # Initialize total area

    α = MVector{length(bodies), Float64}(zeros(length(bodies))) # Initialize angle of attack vector
    β = MVector{length(bodies), Float64}(zeros(length(bodies))) # Initialize sideslip angle vector
    R = MMatrix{3, 3, Float64}(zeros(3, 3)) # Rotation matrix from the root body to the spacecraft link
    # Determine angle of attack (α) and sideslip angle (β)
    # Vehicle Aerodynamic Forces
    # CL and CD
    @inbounds for (i, b) in enumerate(bodies)
        if orientation_sim
            R .= Rot[i] # Rotation matrix from the spacecraft link to the inertial frame
            body_frame_velocity = R' * m.planet.L_PI' * vel_pp_rw # Velocity of the spacecraft link in inertial frame
            
            α_body = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack in radians
            β_body = atan(body_frame_velocity[2], norm([body_frame_velocity[1], body_frame_velocity[3]])) # Sideslip angle in radians
            α[i] = α_body # Angle of attack for the spacecraft link
            β[i] = β_body # Sideslip angle for the spacecraft link
            b.α = α_body
            b.β = β_body
            b.θ = acos(clamp(vel_pp_rw[1]/norm(vel_pp_rw), -1.0, 1.0)) # Elevation angle for the spacecraft link
        else
            # TODO: Change this so that it just uses above code even with orientation_sim = false
            if b.root
                α[i] = pi/2
                b.α = pi/2 # Angle of attack for the root body
            else
                body_frame_velocity = rot(b.q) * SVector{3, Float64}(1.0, 0.0, 0.0) # Velocity of the spacecraft link in inertial frame
                α[i] = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack for the spacecraft link
                b.α = α[i] # Angle of attack for the spacecraft link
            end
        end

        CL_body = 0.0
        CD_body = 2 * (2.2 - 0.8)/pi * args.α + 0.8

        drag_pp_body = q * CD_body * b.ref_area * drag_pp_hat                       # Planet relative drag force vector
        lift_pp_body = q * CL_body * b.ref_area * lift_pp_hat * cos(bank_angle)     # Planet relative lift force vector

        if orientation_sim
            cross_pp_body = q * CS_body * b.ref_area * cross_pp_hat # Planet relative cross force vector
            cross_body = m.planet.L_PI' * cross_pp_body # Inertial cross force vector
        else
            cross_pp_body = SVector{3, Float64}(0.0, 0.0, 0.0) # Planet relative cross force vector
            cross_body = SVector{3, Float64}(0.0, 0.0, 0.0) # Inertial cross force vector
        end

        drag_body = m.planet.L_PI' * drag_pp_body   # Inertial drag force vector
        lift_body = m.planet.L_PI' * lift_pp_body   # Inertial lift force vector

        # Update the force on the spacecraft link
        b.net_force .+= drag_body + lift_body + cross_body # Update the force on the spacecraft link, inertial frame
        b.net_torque .+= cross(b.r, rot_body_to_inertial' * (drag_body + lift_body + cross_body)) # Update the torque on the spacecraft link, body frame
        # Update the total CL/CD
        CL += CL_body * b.ref_area
        CD += CD_body * b.ref_area
        total_area += b.ref_area # Update the total area
        drag_ii += drag_body # Update the total drag force
        lift_ii += lift_body # Update the total lift force
        drag_pp += drag_pp_body # Update the total drag force in planet relative frame
        lift_pp += lift_pp_body # Update the total lift force in planet relative frame
    end
    
    # Normalize the aerodynamic coefficients
    CL = CL / total_area
    CD = CD / total_area

    force_ii, torque_ii = collect_and_reset_link_wrenches!(bodies)

    return force_ii, torque_ii
end

function calcForceTorque(model::AerodynamicCoefficientfM, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    planet = param.args.environment_model.planet
    ephemerides_model = param.args.environment_model.ephemerides_model
    orientation_sim = param.args.mission_configuration.orientation_sim
    spacecraft = param.args.dynamics_model.spacecraft[i]
    bodies = spacecraft.links # Include the root body of the spacecraft
    root_index = 1
    env_state = SimulationModel.SimulationCallbacks._stage_environment_state(
        x,
        param,
        i,
        param.shared_buffers.current_time[];
        write_buffers=true,
    )
    rho = env_state.rho
    T = env_state.T
    wind = env_state.wind

    # Skip expensive aerodynamic geometry when the flow is effectively vacuum.
    if !isfinite(rho) || rho <= eps(Float64)
        _store_aero_caches!(param, i, SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    sound_velocity = sqrt(planet.γ * planet.R * T) # Speed of sound at the current altitude based on the temperature from the shared buffer updated by the callback function
    pos_ii = env_state.pos_ii
    vel_ii = env_state.vel_ii

    h_ii = cross(pos_ii, vel_ii)    # Inertial angular momentum vector [m ^ 2 / s]

    h_ii_mag = norm(h_ii)           # Magnitude of the inertial angular momentum [m ^ 2 / s]

    pos_pp = env_state.pos_pp
    vel_pp = env_state.vel_pp
    pos_pp_mag = norm(pos_pp) # Magnitude of the planet relative position

    vel_pp_mag = norm(vel_pp)
    mach = vel_pp_mag / sound_velocity # Mach number of the flow at the current altitude
    S = sqrt(planet.γ * 0.5) * mach # Molecular speed ratio
    h_pp = cross(pos_pp, vel_pp)
    
    h_pp_mag = norm(h_pp)
    h_pp_hat = normalize(h_pp) # Unit vector of the planet relative angular momentum
    
    bank_angle = deg2rad(0.0)
    alt,lat,lon = env_state.alt, env_state.lat, env_state.lon

    uD, uN, uE = latlongtoNED((alt, lat, lon))
    wE, wN, wU = wind # positive to the east , m / s
    wind_pp = wN * uN + wE * uE - wU * uD         # wind velocity in pp frame, m / s
    vel_pp_rw = vel_pp - wind_pp                  # airspeed: spacecraft velocity minus atmosphere velocity, m / s
    vel_pp_rw_mag = norm(vel_pp_rw)
    if vel_pp_rw_mag <= eps(Float64)
        _store_aero_caches!(param, i, SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    mach = vel_pp_rw_mag / sound_velocity # Mach number of the flow at the current altitude
    S = sqrt(planet.γ * 0.5) * mach # Molecular speed ratio
    vel_pp_rw_hat = vel_pp_rw / vel_pp_rw_mag     # relative wind unit vector
    vel_pp_rw_inv_mag = inv(vel_pp_rw_mag)
    lift_pp_hat = normalize(cross(h_pp_hat, vel_pp_rw_hat))
    # lift_pp_hat /= norm(lift_pp_hat) # Normalize the lift vector in planet relative frame
    drag_pp_hat = -vel_pp_rw_hat # Planet relative drag force direction
    cross_pp_hat = cross(drag_pp_hat, lift_pp_hat) # Cross product of the drag and lift vectors in planet relative frame
    q = 0.5 * rho * vel_pp_rw_mag^2 # Dynamic pressure in planet relative frame, using the stage-consistent density value
    L_PI_t = env_state.l_pi'
    vel_pi = orientation_sim ? L_PI_t * vel_pp_rw : SVector{3, Float64}(0.0, 0.0, 0.0)
    lift_scale = q * cos(bank_angle)

    CL, CD = 0.0, 0.0 # Initialize aerodynamic coefficients
    total_area = 0.0 # Initialize total area
    θ_body = acos(clamp(vel_pp_rw[1] * vel_pp_rw_inv_mag, -1.0, 1.0))
    decision = _multibody_thread_decision(
        length(bodies);
        heavy_work=true,
        env=param.shared_buffers.policy_env_config[],
    )
    use_threads = decision.use_threads
    n_workers = use_threads ? ParallelPolicy.thread_worker_count(length(bodies), decision.allotment) : 1

    zero_vec = SVector{3, Float64}(0.0, 0.0, 0.0)

    function compute_link_wrench!(idx::Int)
        b = bodies[idx]
        if orientation_sim
            R = rotate_to_inertial(spacecraft, b, root_index) # Rotation matrix from the spacecraft link to the inertial frame
            body_frame_velocity = R' * vel_pi # Velocity of the spacecraft link in inertial frame

            α_body = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack in radians
            β_body = atan(body_frame_velocity[2], hypot(body_frame_velocity[1], body_frame_velocity[3])) # Sideslip angle in radians
            b.α = α_body
            b.β = β_body
            b.θ = θ_body
        else
            # TODO: Change this so that it just uses above code even with orientation_sim = false
            if b.root
                b.α = pi / 2 # Angle of attack for the root body
            else
                body_frame_velocity = rot(b.q) * SVector{3, Float64}(1.0, 0.0, 0.0) # Velocity of the spacecraft link in inertial frame
                b.α = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack for the spacecraft link
            end
        end

        CL_body, CD_body, CS_body, _, _, _ = aerodynamic_coefficient_fM(b, T, S)

        drag_pp_body = q * CD_body * b.ref_area * drag_pp_hat                       # Planet relative drag force vector
        lift_pp_body = lift_scale * CL_body * b.ref_area * lift_pp_hat     # Planet relative lift force vector

        if orientation_sim
            cross_pp_body = q * CS_body * b.ref_area * cross_pp_hat # Planet relative cross force vector
            cross_body = L_PI_t * cross_pp_body # Inertial cross force vector
        else
            cross_pp_body = zero_vec # Planet relative cross force vector
            cross_body = zero_vec # Inertial cross force vector
        end

        drag_body = L_PI_t * drag_pp_body   # Inertial drag force vector
        lift_body = L_PI_t * lift_pp_body   # Inertial lift force vector
        force_body = drag_body + lift_body + cross_body

        return force_body, drag_body, lift_body, cross_body, CL_body * b.ref_area, CD_body * b.ref_area, b.ref_area
    end

    force_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    drag_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    lift_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    cross_ii = MVector{3, Float64}(0.0, 0.0, 0.0)

    started_ns = decision.policy_applied ? time_ns() : UInt64(0)
    if use_threads
        workspace = _aero_workspace_for_sat!(param, i, n_workers)
        thread_force = workspace.thread_force
        thread_drag = [MVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n_workers]
        thread_lift = [MVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n_workers]
        thread_cross = [MVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n_workers]
        thread_cl = workspace.thread_cl
        thread_cd = workspace.thread_cd
        thread_area = workspace.thread_area
        @inbounds for worker_id in 1:n_workers
            thread_force[worker_id] .= 0.0
            thread_drag[worker_id] .= 0.0
            thread_lift[worker_id] .= 0.0
            thread_cross[worker_id] .= 0.0
            thread_cl[worker_id] = 0.0
            thread_cd[worker_id] = 0.0
            thread_area[worker_id] = 0.0
        end

        ParallelPolicy.threaded_foreach_worker_persistent(:rhs_aero, length(bodies), decision.allotment) do worker_id, idx
            force_body, drag_body, lift_body, cross_body, cl_area, cd_area, area = compute_link_wrench!(idx)
            thread_force[worker_id] .+= force_body
            thread_drag[worker_id] .+= drag_body
            thread_lift[worker_id] .+= lift_body
            thread_cross[worker_id] .+= cross_body
            thread_cl[worker_id] += cl_area
            thread_cd[worker_id] += cd_area
            thread_area[worker_id] += area
        end

        @inbounds for worker_id in 1:n_workers
            force_ii .+= thread_force[worker_id]
            drag_ii .+= thread_drag[worker_id]
            lift_ii .+= thread_lift[worker_id]
            cross_ii .+= thread_cross[worker_id]
            CL += thread_cl[worker_id]
            CD += thread_cd[worker_id]
            total_area += thread_area[worker_id]
        end
    else
        @inbounds for idx in eachindex(bodies)
            force_body, drag_body, lift_body, cross_body, cl_area, cd_area, area = compute_link_wrench!(idx)
            force_ii .+= force_body
            drag_ii .+= drag_body
            lift_ii .+= lift_body
            cross_ii .+= cross_body
            CL += cl_area
            CD += cd_area
            total_area += area
        end
    end
    if decision.policy_applied
        ParallelPolicy.record_policy_observation!(
            :multibody;
            mode=decision.mode,
            num_items=length(bodies),
            use_threads=use_threads,
            elapsed_ns=(time_ns() - started_ns)
        )
    end

    if total_area > 0.0
        CL = CL / total_area
        CD = CD / total_area
    end

    _store_aero_caches!(param, i, SVector{3, Float64}(drag_ii), SVector{3, Float64}(lift_ii), SVector{3, Float64}(cross_ii))
    return SVector{3, Float64}(force_ii), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function calcForceTorque(model::AerodynamicCoefficientNoBallisticFlight, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    planet = param.args.environment_model.planet
    ephemerides_model = param.args.environment_model.ephemerides_model
    orientation_sim = param.args.mission_configuration.orientation_sim
    spacecraft = param.args.dynamics_model.spacecraft[i]
    bodies = spacecraft.links # Include the root body of the spacecraft
    root = spacecraft.root
    root_index = 1
    pos_ii = SVector{3, Float64}(x[1], x[2], x[3])
    vel_ii = SVector{3, Float64}(x[4], x[5], x[6])

    h_ii = cross(pos_ii, vel_ii)    # Inertial angular momentum vector [m ^ 2 / s]

    h_ii_mag = norm(h_ii)           # Magnitude of the inertial angular momentum [m ^ 2 / s]

    # Inertial to planet relative transformation
    et = param.shared_buffers.et_start[] + param.shared_buffers.current_time[]
    println("ephemerides_model: ", ephemerides_model)
    pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, planet, et, ephemerides_model) # Position vector planet / planet[m] # Velocity vector planet / planet[m / s]
    pos_pp_mag = norm(pos_pp) # Magnitude of the planet relative position

    vel_pp_mag = norm(vel_pp)

    h_pp = cross(pos_pp, vel_pp)
    
    h_pp_mag = norm(h_pp)
    h_pp_hat = normalize(h_pp) # Unit vector of the planet relative angular momentum
    
    bank_angle = deg2rad(0.0)
        
    lift_pp_hat = normalize(cross(h_pp_hat, vel_pp_rw_hat))
    # lift_pp_hat /= norm(lift_pp_hat) # Normalize the lift vector in planet relative frame
    drag_pp_hat = -vel_pp_rw_hat # Planet relative drag force direction
    cross_pp_hat = cross(drag_pp_hat, lift_pp_hat) # Cross product of the drag and lift vectors in planet relative frame

    if orientation_sim
        Rot = [MMatrix{3,3,Float64}(zeros(3, 3)) for i in eachindex(bodies)] # Rotation matrix from the root body to the spacecraft link
        @inbounds for (i, b) in enumerate(bodies)
            Rot[i] .= rotate_to_inertial(spacecraft, b, root_index) # Rotation matrix from the spacecraft link to the inertial frame
        end
    end
    
    CL, CD = 0.0, 0.0 # Initialize aerodynamic coefficients
    total_area = 0.0 # Initialize total area

    α = MVector{length(bodies), Float64}(zeros(length(bodies))) # Initialize angle of attack vector
    β = MVector{length(bodies), Float64}(zeros(length(bodies))) # Initialize sideslip angle vector
    R = MMatrix{3, 3, Float64}(zeros(3, 3)) # Rotation matrix from the root body to the spacecraft link
    # Determine angle of attack (α) and sideslip angle (β)
    # Vehicle Aerodynamic Forces
    # CL and CD
    @inbounds for (i, b) in enumerate(bodies)
        if orientation_sim
            R .= Rot[i] # Rotation matrix from the spacecraft link to the inertial frame
            body_frame_velocity = R' * planet.L_PI' * vel_pp_rw # Velocity of the spacecraft link in inertial frame
            
            α_body = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack in radians
            β_body = atan(body_frame_velocity[2], norm([body_frame_velocity[1], body_frame_velocity[3]])) # Sideslip angle in radians
            α[i] = α_body # Angle of attack for the spacecraft link
            β[i] = β_body # Sideslip angle for the spacecraft link
            b.α = α_body
            b.β = β_body
            b.θ = acos(clamp(vel_pp_rw[1]/norm(vel_pp_rw), -1.0, 1.0)) # Elevation angle for the spacecraft link
        else
            # TODO: Change this so that it just uses above code even with orientation_sim = false
            if b.root
                α[i] = pi/2
                b.α = pi/2 # Angle of attack for the root body
            else
                body_frame_velocity = rot(b.q) * SVector{3, Float64}(1.0, 0.0, 0.0) # Velocity of the spacecraft link in inertial frame
                α[i] = atan(body_frame_velocity[1], body_frame_velocity[3]) # Angle of attack for the spacecraft link
                b.α = α[i] # Angle of attack for the spacecraft link
            end
        end

        CL_body = 0.0
        CD_body = 2 * (2.2 - 0.8)/pi * args.α + 0.8
        
        drag_pp_body = q * CD_body * b.ref_area * drag_pp_hat                       # Planet relative drag force vector
        lift_pp_body = q * CL_body * b.ref_area * lift_pp_hat * cos(bank_angle)     # Planet relative lift force vector

        if orientation_sim
            cross_pp_body = q * CS_body * b.ref_area * cross_pp_hat # Planet relative cross force vector
            cross_body = planet.L_PI' * cross_pp_body # Inertial cross force vector
        else
            cross_pp_body = SVector{3, Float64}(0.0, 0.0, 0.0) # Planet relative cross force vector
            cross_body = SVector{3, Float64}(0.0, 0.0, 0.0) # Inertial cross force vector
        end

        drag_body = planet.L_PI' * drag_pp_body   # Inertial drag force vector
        lift_body = planet.L_PI' * lift_pp_body   # Inertial lift force vector

        # Update the force on the spacecraft link
        b.net_force .+= drag_body + lift_body + cross_body # Update the force on the spacecraft link, inertial frame
        b.net_torque .+= cross(b.r, rot_body_to_inertial' * (drag_body + lift_body + cross_body)) # Update the torque on the spacecraft link, body frame
        # Update the total CL/CD
        CL += CL_body * b.ref_area
        CD += CD_body * b.ref_area
        total_area += b.ref_area # Update the total area
        drag_ii += drag_body # Update the total drag force
        lift_ii += lift_body # Update the total lift force
        drag_pp += drag_pp_body # Update the total drag force in planet relative frame
        lift_pp += lift_pp_body # Update the total lift force in planet relative frame
    end
    
    force_ii, torque_ii = collect_and_reset_link_wrenches!(bodies)

    return force_ii, torque_ii
end


function aerodynamic_coefficient_constant(α, body, T, S, args, montecarlo=false)
    """

    """

    CL_body = 0.0
    CD_body = 2 * (2.2 - 0.8)/pi * args.α + 0.8

    if montecarlo == true
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)
    end

    return CL_body, CD_body
end

function aerodynamic_coefficient_fM(body, T::Float64, S::Float64)
    return aerodynamic_coefficient_fM(body, T, S, body.α, body.β, body.θ)
end

function aerodynamic_coefficient_fM(
    body,
    T::Float64,
    S::Float64,
    α_input::Float64,
    β_input::Float64,
    θ_input::Float64,
)
    """
    Calculate the aerodynamic coefficients for a blunt body in ballistic flight using the F.M. model.
        Angle of attack and side slip are calculated as angles between CA (normal to flat plate), and wind-relative velocity vector.
    # Arguments
    - `body`: Body object containing dimensions and properties
    - `T`: Temperature, K
    - `S`: Molecular speed ratio
    - `args`: Dictionary containing additional parameters like reflection coefficient
    - `montecarlo`: Boolean flag to apply Monte Carlo perturbations
    # Returns
    - `CL`: Lift coefficient
    - `CD`: Drag coefficient
    - `CS`: Sideslip coefficient
    # Notes
    - Equations are from Hart et al. (2017, https://doi.org/10.2514/1.A33606), for a rectangular prism. 
    """

    α = α_input
    β = β_input
    θ = θ_input
    # Adjust angles to match model assumptions
    α -= π/2
    σ = body.reflection_coefficient
    Tw = T
    lx = body.dims[1]
    ly = body.dims[2]
    lz = body.dims[3]
    lx_over_ly = lx / ly
    lx_over_lz = lx / lz
    σN = σ
    σT = σ
    cosα = cos(α)
    cosβ = cos(β)
    sinβ = sin(β)
    sinα = sin(α)
    # The Hart et al. coefficients are derived in body axes, then rotated into wind axes below.
    CA = ((2-σN)/(S*sqrt_π)*cosα*cosβ+sign(cosα*cosβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*cosα^2*cosβ^2) +
            (2-σN)*(cosα^2*cosβ^2+1/(2*S^2)) * (sign(cosα*cosβ)+erf(S*cosα*cosβ)) + 
            (σN/(2*S)*cosα*cosβ*√(π*Tw/T)) * (1+sign(cosα*cosβ)*erf(S*cosα*cosβ)) +
            σT*cosα*cosβ*(lx_over_ly*(1/(S*sqrt_π)*exp(-S^2*sinβ^2)+sinβ*(sign(sinβ)+erf(S*sinβ))) +
            lx_over_lz*(1/(S*sqrt_π)*exp(-S^2*sinα^2*cosβ^2)+sinα*cosβ*(sign(sinα*cosβ)+erf(S*sinα*cosβ))))
    CS = lx_over_ly*(((2-σN)/(S*sqrt_π)*sinβ+sign(sinβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*sinβ^2) +
                (2-σN)*(sinβ^2+1/(2*S^2)) * (sign(sinβ)+erf(S*sinβ)) + 
                (σN/(2*S)*sinβ*√(π*Tw/T)) * (1+sign(sinβ)*erf(S*sinβ))) +
                σT*sinβ*(1/(S*sqrt_π)*exp(-S^2*cosα^2*cosβ^2) + cosα*cosβ*(sign(cosα*cosβ)+erf(S*cosα*cosβ)) + 
                lx_over_lz*(1/(S*sqrt_π)*exp(-S^2*sinα^2*cosβ^2) + sinα*cosβ*(sign(sinα*cosβ)+erf(S*sinα*cosβ))))
    CN = lx_over_lz*((((2-σN)/(S*sqrt_π)*sinα*cosβ+sign(sinα*cosβ)*σN/(2*S^2)*√(Tw/T))*exp(-S^2*sinα^2*cosβ^2) +
                (2-σN)*(sinα^2*cosβ^2+1/(2*S^2)) * (sign(sinα*cosβ)+erf(S*sinα*cosβ)) + 
                (σN/(2*S)*sinα*cosβ*√(π*Tw/T)) * (1+sign(sinα*cosβ)*erf(S*sinα*cosβ)))) +
                σT*sinα*cosβ*(lx_over_ly*(1/(S*sqrt_π)*exp(-S^2*sinβ^2) + sinβ*(sign(sinβ)+erf(S*sinβ))) + 
                (1/(S*sqrt_π)*exp(-S^2*cosα^2*cosβ^2) + cosα*cosβ*(sign(cosα*cosβ)+erf(S*cosα*cosβ))))

    CL = -sinα*CA + cosα*CN
    CD = cosα*cosβ*CA + sinβ*CS + sinα*cosβ*CN
    CS = -cosα*sinβ*CA + cosβ*CS - sinα*sinβ*CN

    return CL, CD, CS, 0.0, 0.0, 0.0
end

function aerodynamic_coefficient_no_ballistic_flight(α, body, args, T=0, S=0, a=0, montecarlo=false)
    """

    """

    # Newtonian flow
    NoseRadius = body.nose_radius
    BaseRadius = body.base_radius

    k = NoseRadius / BaseRadius

    Cp_max = 2
    δ = body.δ

    CA_body = (1 - sin(δ)^4)*k^2 + (2*sin(δ)^2*cos(α)^2 + cos(δ)^2*sin(α)^2) * (1 - (k*cos(δ))^2)
    CN_body = (1 - (k*cos(δ))^2) * cos(δ^2) *sin(2*α)
    
    CD_body = CA_body*cos(α) + CN_body*sin(α) - 0.15
    CL_body = CN_body*cos(α) - CA_body*sin(α)

    if montecarlo == true
        CL_body, CD_body = monte_carlo_aerodynamics(CL_body, CD_body, args)
    end

    return CL_body, CD_body
end

"""
    calculate_CA(t1, t2, t3, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf, θ)

Calculates the expression C_A from the provided image.
This version incorporates the clarification that θ is a scalar variable (an angle).

Arguments:
- t1, t2, t3: t₁, t₂, t₃
- cos_α, sin_α, cos_β, sin_β: trig functions of angle of attack and sideslip angles (angles)
- θ: theta (flow angle)
- σ_T, σ_N: sigma_T, sigma_N
- s: s
- Tw: T_w (Wall Temperature)
- T_inf: T_∞ (Temperature at infinity)
"""
function calculate_CA(t1, t2, t3, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf, θ)
    # --- Pre-calculate common values ---
    s_sq = s^2
    sqrt_Tw_Tinf = sqrt(Tw / T_inf)

    # --- Group 1: Terms 1 and 2 ---
    s_sin_β = s * sin_β
    erf_s_sin_β = erf(s_sin_β)
    exp_s_sin_β_sq = exp(-s_sin_β^2)

    paren_term_minus_1 = sqrt_π * s_sin_β * (erf_s_sin_β - 1) + exp_s_sin_β_sq
    paren_term_plus_1 = sqrt_π * s_sin_β * (erf_s_sin_β + 1) + exp_s_sin_β_sq
    
    # Denominator simplifies because sqrt(sec(β)²) * abs(cos(β)) = 1
    den1_2 = sqrt_π * s * t2
    common_factor_1_2 = t1 * cos_α * cos_β * σ_T * θ
    
    term1 = (common_factor_1_2 * (-sin_β) * paren_term_minus_1) / den1_2
    term2_frac = (common_factor_1_2 * sin_β * paren_term_plus_1) / den1_2

    # --- Group 2: Terms 3 and 4 (no denominators) ---
    s_ca_cb = s * cos_α * cos_β
    erf_s_cacb = erf(s_ca_cb)
    exp_s_cacb_sq = exp(-s_ca_cb^2)
    s_sq_ca_sq_cb_sq = s_sq * cos_α^2 * cos_β^2
    
    # Term 3
    exp_bracket3 = (s_ca_cb * (σ_N - 2)) / sqrt_π + 0.5 * σ_N * sqrt_Tw_Tinf
    main_bracket3 = (1 - erf_s_cacb) * ((2 - σ_N) * (s_sq_ca_sq_cb_sq + 0.5) - 
                    0.5 * sqrt_π * s_ca_cb * σ_N * sqrt_Tw_Tinf) + 
                    exp_s_cacb_sq * exp_bracket3
    term3 = 1/s_sq * θ * (-cos_α) * cos_β * main_bracket3
    
    # Term 4
    exp_bracket4 = 0.5 * σ_N * sqrt_Tw_Tinf - (s_ca_cb * (σ_N - 2)) / sqrt_π
    main_bracket4 = (erf_s_cacb + 1) * ((2 - σ_N) * (s_sq_ca_sq_cb_sq + 0.5) +
                    0.5 * sqrt_π * s_ca_cb * σ_N * sqrt_Tw_Tinf) +
                    exp_s_cacb_sq * exp_bracket4
    term4 = 1/s_sq * θ * cos_α * cos_β * main_bracket4

    # --- Group 3: Terms 5 and 6 ---
    s_sa_cb = s * sin_α * cos_β
    erf_s_sacb = erf(s_sa_cb)
    exp_s_sacb_sq = exp(-s_sa_cb^2)

    paren_term5 = sqrt_π * s_sa_cb * (erf_s_sacb - 1) + exp_s_sacb_sq
    paren_term6 = sqrt_π * s_sa_cb * (erf_s_sacb + 1) + exp_s_sacb_sq

    den5_6 = sqrt_π * s * t3
    common_factor_5_6 = t1 * cos_α * cos_β * σ_T * θ
    
    term5 = (common_factor_5_6 * (-cos_β * sin_α) * paren_term5) / den5_6
    term6 = (common_factor_5_6 * (cos_β * sin_α) * paren_term6) / den5_6

    # --- Final Combination ---
    # Based on the signs shown in the image:
    # T1 + T2_frac - 1/s² - T3 + T4 + T5 + T6
    Ca_result = term1 + term2_frac - term3 + term4 + term5 + term6
    
    return Ca_result
end

"""
    calculate_CS(t1, t2, t3, θ, α, β, σ_T, σ_N, s, Tw, T_inf)

Calculates the expression C_S from the provided image.

Arguments:
- t1, t2, t3: t₁, t₂, t₃
- θ, α, β: theta, alpha, beta (angles)
- σ_T, σ_N: sigma_T, sigma_N
- s: s
- Tw: T_w (Wall Temperature)
- T_inf: T_∞ (Temperature at infinity)
"""
function calculate_CS(t1, t2, t3, θ, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf)
    # --- Pre-calculate common values ---
    s_sq = s^2
    sqrt_Tw_Tinf = sqrt(Tw / T_inf)

    # --- Group 1: Terms 1 and 2 ---
    s_sin_β = s * sin_β
    erf_s_sin_β = erf(s_sin_β)
    exp_s_sin_β_sq = exp(-s_sin_β^2)
    s_sq_sin_β_sq = s_sq * sin_β^2

    # Term 1
    exp_bracket1 = (s_sin_β * (σ_N - 2)) / sqrt_π + 0.5 * σ_N * sqrt_Tw_Tinf
    main_bracket1 = (1 - erf_s_sin_β) * ((2 - σ_N) * s_sq_sin_β_sq + 0.5) -
                    0.5 * sqrt_π * s_sin_β * σ_N * sqrt_Tw_Tinf +
                    exp_s_sin_β_sq * exp_bracket1
    num1 = t1 * θ * (-sin_β) * main_bracket1
    den1_2 = s_sq * t2
    term1 = num1 / den1_2

    # Term 2
    exp_bracket2 = 0.5 * σ_N * sqrt_Tw_Tinf - (s_sin_β * (σ_N - 2)) / sqrt_π
    main_bracket2 = (erf_s_sin_β + 1) * ((2 - σ_N) * s_sq_sin_β_sq + 0.5) +
                    0.5 * sqrt_π * s_sin_β * σ_N * sqrt_Tw_Tinf +
                    exp_s_sin_β_sq * exp_bracket2
    num2 = t1 * θ * sin_β * main_bracket2
    term2 = num2 / den1_2

    # --- Group 2: Terms 3 and 4 ---
    s_ca_cb = s * cos_α * cos_β
    erf_s_cacb = erf(s_ca_cb)
    exp_s_cacb_sq = exp(-s_ca_cb^2)
    
    # Term 3
    paren_term3 = sqrt_π * s_ca_cb * (erf_s_cacb - 1) + exp_s_cacb_sq
    num3 = sin_β * σ_T * θ * (-cos_α * cos_β) * paren_term3
    den3_4 = sqrt_π * s
    term3 = num3 / den3_4

    # Term 4
    paren_term4 = sqrt_π * s_ca_cb * (erf_s_cacb + 1) + exp_s_cacb_sq
    num4 = sin_β * σ_T * θ * (cos_α * cos_β) * paren_term4
    term4 = num4 / den3_4

    # --- Group 3: Terms 5 and 6 ---
    s_sa_cb = s * sin_α * cos_β
    erf_s_sacb = erf(s_sa_cb)
    exp_s_sacb_sq = exp(-s_sa_cb^2)

    paren_term5 = sqrt_π * s_sa_cb * (erf_s_sacb - 1) + exp_s_sacb_sq
    paren_term6 = sqrt_π * s_sa_cb * (erf_s_sacb + 1) + exp_s_sacb_sq

    den5_6 = sqrt_π * s * t3
    common_factor_5_6_num = t1 * sin_β * σ_T * θ
    
    num5 = common_factor_5_6_num * (-cos_β * sin_α) * paren_term5
    term5 = num5 / den5_6

    num6 = common_factor_5_6_num * (cos_β * sin_α) * paren_term6
    term6 = num6 / den5_6
    
    # --- Final Combination ---
    # Based on the signs shown at the start and end of each line in the image
    Cs_result = -term1 + term2 + term3 + term4 + term5 + term6
    
    return Cs_result
end

"""
    calculate_CN(t1, t2, t3, θ, α, β, σ_T, σ_N, s, Tw, T_inf)

Calculates the expression C_N from the provided image.

Arguments:
- t1, t2, t3: t₁, t₂, t₃
- θ, α, β: theta, alpha, beta (angles)
- σ_T, σ_N: sigma_T, sigma_N
- s: s
- Tw: T_w (Wall Temperature)
- T_inf: T_∞ (Temperature at infinity)
"""
function calculate_CN(t1, t2, t3, θ, sin_α, cos_α, sin_β, cos_β, σ_T, σ_N, s, Tw, T_inf)
    # --- Pre-calculate common values ---
    s_sq = s^2
    sqrt_Tw_Tinf = sqrt(Tw / T_inf)

    # --- Group 1: Terms 1 and 2 ---
    s_sin_β = s * sin_β
    erf_s_sin_β = erf(s_sin_β)
    exp_s_sin_β_sq = exp(-s_sin_β^2)

    paren_term_minus_1 = sqrt_π * s_sin_β * (erf_s_sin_β - 1) + exp_s_sin_β_sq
    paren_term_plus_1 = sqrt_π * s_sin_β * (erf_s_sin_β + 1) + exp_s_sin_β_sq
    
    # Denominator simplifies because sqrt(sec(β)²) * abs(cos(β)) = 1
    den1_2 = sqrt_π * s * t2
    common_factor_1_2_num = t1 * sin_α * cos_β * σ_T * θ
    
    num1 = common_factor_1_2_num * (-sin_β) * paren_term_minus_1
    term1 = num1 / den1_2

    num2_frac = common_factor_1_2_num * sin_β * paren_term_plus_1
    term2_frac = num2_frac / den1_2
    term2_sub = 1 / (s_sq * t3)

    # --- Group 2: Terms 3 and 4 ---
    s_sa_cb = s * sin_α * cos_β
    erf_s_sacb = erf(s_sa_cb)
    exp_s_sacb_sq = exp(-s_sa_cb^2)
    s_sq_sa_sq_cb_sq = s_sq * sin_α^2 * cos_β^2

    # Term 3
    exp_bracket3 = (s_sa_cb * (σ_N - 2)) / sqrt_π + 0.5 * σ_N * sqrt_Tw_Tinf + 1 / (s_sq * t3)
    main_bracket3 = (1 - erf_s_sacb) * ((2 - σ_N) * s_sq_sa_sq_cb_sq - 0.5) -
                    0.5 * sqrt_π * s_sa_cb * σ_N * sqrt_Tw_Tinf +
                    exp_s_sacb_sq * exp_bracket3
    term3 = t1 * θ * (-cos_β * sin_α) * main_bracket3
    
    # Term 4
    exp_bracket4 = 0.5 * σ_N * sqrt_Tw_Tinf - (s_sa_cb * (σ_N - 2)) / sqrt_π
    main_bracket4 = (erf_s_sacb + 1) * ((2 - σ_N) * s_sq_sa_sq_cb_sq + 0.5) +
                    0.5 * sqrt_π * s_sa_cb * σ_N * sqrt_Tw_Tinf +
                    exp_s_sacb_sq * exp_bracket4
    term4 = t1 * θ * (cos_β * sin_α) * main_bracket4

    # --- Group 3: Terms 5 and 6 ---
    s_ca_cb = s * cos_α * cos_β
    erf_s_cacb = erf(s_ca_cb)
    exp_s_cacb_sq = exp(-s_ca_cb^2)

    paren_term5 = sqrt_π * s_ca_cb * (erf_s_cacb - 1) + exp_s_cacb_sq
    paren_term6 = sqrt_π * s_ca_cb * (erf_s_cacb + 1) + exp_s_cacb_sq
    
    den5_6 = sqrt_π * s
    common_factor_5_6_num = sin_α * cos_β * σ_T * θ

    num5 = common_factor_5_6_num * (-cos_α * cos_β) * paren_term5
    term5 = num5 / den5_6
    
    num6 = common_factor_5_6_num * (cos_α * cos_β) * paren_term6
    term6 = num6 / den5_6
    
    # --- Final Combination ---
    # Based on the signs shown in the image: 
    # T1 + T2_frac - 1/(s²t₃) - T3 + T4 + T5 + T6
    Cn_result = term1 + term2_frac - term2_sub - term3 + term4 + term5 + term6
    
    return Cn_result
end

"""
    calculate_Cl(σ_T, s, α, β, θ)

Calculates the expression C_L from the provided image.

Arguments:
- σ_T: sigma_T
- s: s
- α, β, θ: alpha, beta, theta (angles)
"""
function calculate_Cl(σ_T, s, sin_α, cos_α, sin_β, cos_β, θ)
    # --- Pre-calculate common values ---

    # --- Calculation for the "Top Part" (first three lines) ---
    # This part is structured as: ( sin(β) * [ ... ] ) - 1/(sqrt(π)*|cos(β)|)
    
    s_sa_cb = s * sin_α * cos_β
    erf_s_sacb = erf(s_sa_cb)
    exp_s_sacb_sq = exp(-s_sa_cb^2)

    # The two θ terms can be combined and simplified algebraically:
    # θ(-x)A + θ(x)B = -θx A + θx B = θx(B - A)
    paren_A = s_sa_cb * erf_s_sacb + exp_s_sacb_sq / sqrt_π - s_sa_cb
    paren_B = -s_sa_cb * erf_s_sacb - exp_s_sacb_sq / sqrt_π - s_sa_cb
    
    term_in_paren_top = paren_B - paren_A # This simplifies to -2*s_sa_cb*erf - 2*exp/sqrt_π
    
    # Note: θ(-cos(β)sin(α)) becomes -θ*cos(β)sin(α) etc.
    sub_part_1 = sin_β * (θ * cos_β * sin_α * term_in_paren_top)
    
    sub_part_2 = 1 / (sqrt_π * abs(cos_β))
    
    top_part = sub_part_1 #- sub_part_2

    # --- Calculation for the "Bottom Part" (last three lines) ---
    # This is structured as: [prefactors] * ( θ(-sinβ)[...] - θ(sinβ)[...] )
    s_sin_β = s * sin_β
    exp_s_sin_β_sq_neg = exp(-s_sin_β^2)
    
    # Note: sqrt(sec(β)²) simplifies to 1/abs(cos(β))
    prefactor_bottom = sin_α * cos_β^3 * (1 / abs(cos_β)) * exp_s_sin_β_sq_neg

    # The bracketed terms can also be simplified:
    exp_s_sin_β_sq_pos = exp(s_sin_β^2)
    erf_s_sin_β = erf(s_sin_β)
    common_in_brackets = sqrt_π * s_sin_β * exp_s_sin_β_sq_pos * erf_s_sin_β
    
    bracket_A = common_in_brackets + sqrt_π * s_sin_β * (-exp_s_sin_β_sq_neg + 1)
    bracket_B = common_in_brackets + sqrt_π * s_sin_β * (exp_s_sin_β_sq_neg + 1)
    
    # θ(-sinβ)A - θ(sinβ)B = -θ*sinβ*A - θ*sinβ*B = -θ*sinβ*(A + B)
    sum_of_brackets = bracket_A + bracket_B # Simplifies to 2*common + 2*sqrt(π)*s_sin_β

    main_paren_bottom = -θ * sin_β * sum_of_brackets

    bottom_part = prefactor_bottom * main_paren_bottom
    
    # --- Final Combination ---
    total_inside_braces = top_part - sub_part_2 * bottom_part
    
    Cl_result = (1 / (2 * s)) * σ_T * total_inside_braces
    
    return Cl_result
end

"""
    calculate_Cm(σ_T, s, α, β, θ)

Calculates the expression C_m from the provided image.

Arguments:
- σ_T: sigma_T
- s: s
- α, β, θ: alpha, beta, theta (angles)
"""
function calculate_Cm(σ_T, s, sin_α, cos_α, sin_β, cos_β, θ)
    # --- Pre-calculate common values ---

    # --- Define arguments for the two main patterns ---
    arg1 = s * sin_α * cos_β
    arg2 = s * cos_α * cos_β

    # --- Term 1 ---
    erf_arg1 = erf(arg1)
    exp_arg1_sq = exp(-arg1^2)
    
    # NOTE: The structure of this parenthesis is unusual compared to the others.
    # It is translated literally as it appears in the image.
    paren1 = sin_α * cos_β - sin_α * cos_β * erf_arg1 - exp_arg1_sq / (sqrt_π * s)
    # term1 = cos_α * (θ * (-cos_β * sin_α)) * paren1

    # --- Term 2 ---
    paren2 = arg1 * erf_arg1 + exp_arg1_sq / sqrt_π + arg1
    # term2 = (1/s) * cos_α * (θ * (cos_β * sin_α)) * paren2

    # --- Term 3 ---
    erf_arg2 = erf(arg2)
    exp_arg2_sq = exp(-arg2^2)
    
    paren3 = arg2 * erf_arg2 + exp_arg2_sq / sqrt_π - arg2
    # term3 = sin_α * (θ * (-cos_α * cos_β)) * paren3

    # --- Term 4 ---
    # NOTE: This term is not prefixed by sin(α), unlike its counterpart, Term 3.
    # This asymmetry is preserved from the image.
    paren4 = -arg2 * erf_arg2 - exp_arg2_sq / sqrt_π - arg2
    term4 = (θ * (cos_α * cos_β)) * paren4

    # --- Final Combination ---
    total_sum = cos_α*θ*(-cos_β*sin_α)*paren1 + 1/s*(cos_α*θ*cos_β*cos_α*paren2 + sin_α*(θ * (-cos_α * cos_β) * paren3 + term4))
    
    Cm_result = 0.5 * cos_β * σ_T * total_sum
    
    return Cm_result
end

"""
    calculate_Cn(σ_T, s, α, β, θ)

Calculates the expression C_n from the provided image.

Arguments:
- σ_T: sigma_T
- s: s
- α, β, θ: alpha, beta, theta (angles)
"""
function calculate_Cn(σ_T, s, sin_α, cos_α, sin_β, cos_β, θ)
    # --- Pre-calculate common values ---
    s_sq = s^2
    abs_cos_β = abs(cos_β)

    # --- Define common exponential terms ---
    exp_arg_cos_term = s_sq * cos_α^2 * cos_β^2
    exp_arg_sin_term = s_sq * sin_β^2
    
    exp_cos_pos = exp(exp_arg_cos_term)
    exp_cos_neg = exp(-exp_arg_cos_term)
    exp_sin_pos = exp(exp_arg_sin_term)
    
    exp_main_pos = exp(exp_arg_cos_term + exp_arg_sin_term)
    exp_main_neg = exp(-(exp_arg_cos_term + exp_arg_sin_term))

    # --- Outer Factor ---
    outer_factor = (1 / (2 * sqrt_π * s * abs_cos_β)) * σ_T * exp_main_neg

    # --- Part 1: Terms factored by (-sin(β)|cos(β)|) ---
    s_ca_cb = s * cos_α * cos_β
    
    # The two θ terms inside the first main parenthesis can be combined and simplified.
    # The structure is prop. to θ*cosαcosβ * (Paren_B - Paren_A)
    # (Paren_B - Paren_A) simplifies to: 2 * sqrt(π) * s_ca_cb * exp(-s²cos²αcos²β)
    diff_paren_AB = 2 * sqrt_π * s_ca_cb * exp_cos_neg
    combo_AB = θ * cos_α * cos_β * diff_paren_AB
    
    part1 = (-sin_β * abs_cos_β) * combo_AB
    
    # --- Part 2: Terms factored by (cos³(β)sqrt(sec²(β))) ---
    s_sin_β = s * sin_β
    erf_s_sinb = erf(s_sin_β)

    # The two θ terms inside the second main parenthesis can also be simplified.
    # The structure is prop. to -θ*sinβ * (Paren_C + Paren_D)
    # (Paren_C + Paren_D) simplifies to: 
    # 2 * (sqrt(π)*s*sinβ*erf(s*sinβ)*exp_main_pos + exp(s²cos²αcos²β))
    sum_paren_CD = 2 * (sqrt_π * s_sin_β * erf_s_sinb * exp_main_pos + exp_cos_pos)
    combo_CD = θ * (-sin_β) * sum_paren_CD
    
    # Prefactor for Part 2 simplifies: cos³(β)sqrt(sec²(β)) = cos³(β)/|cos(β)|
    prefactor2 = cos_β^3 / abs_cos_β
    part2 = prefactor2 * combo_CD

    # --- Final Combination ---
    total_sum_in_brackets = part1 + part2
    
    Cn_result = outer_factor * total_sum_in_brackets
    
    return Cn_result
end
