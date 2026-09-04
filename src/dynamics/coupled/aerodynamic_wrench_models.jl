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

@inline function _multibody_thread_decision(num_items::Int; heavy_work::Bool=true)
    mode = _multibody_parallel_mode()
    threshold = _multibody_thread_threshold()
    outer_active = _multibody_outer_parallel_hint()
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
        source=:multibody
    )
    allotment = min(policy.allotment, _multibody_max_threads())
    use_threads = policy.use_threads && allotment > 1
    return (use_threads=use_threads, allotment=use_threads ? allotment : 1, mode=mode)
end

# Aerodynamic models
@kwdef struct AerodynamicCoefficientConstant <: AbstractForceTorqueModel

end

"""
Free-molecular aerodynamics from the Hart et al. closed forms (rectangular
prism, doi 10.2514/1.A33606), evaluated per link and summed.

`fixed_attitude_incidence` selects how link incidence is treated when
`orientation_sim=false` (it has no effect when attitude is simulated):

- `:max_drag` (default, historical behavior): every link is treated as
  flow-normal and charged its full `ref_area` — the spacecraft permanently
  flies its maximum-drag attitude. For multi-link vehicles this overstates
  drag by roughly (sum of face areas)/(true projected area). Retained as
  the default so previously calibrated scenarios stay bit-identical.
- `:attitude`: each link's angle of attack is derived from its configured
  quaternion, so the Hart closed forms' incidence dependence is honored
  for a user-specified fixed attitude. The quaternion is interpreted
  relative to a flow-aligned frame (velocity along reference +x), giving a
  constant incidence around the orbit — an attitude held relative to the
  velocity direction, not an inertially fixed one. Child quaternions are
  root-relative and compose with the root attitude, so a rigidly mounted
  identity-quaternion child shares the root's incidence. Coincides with
  `:max_drag` when all link quaternions are identity. Sideslip remains
  zero in this mode, so pitch-plane attitudes are represented exactly
  while a yaw-only attitude (flow into body ±y) degenerates to
  `:max_drag`; full alpha/beta requires `orientation_sim`.
- `:tumbling_average`: uncontrolled-spacecraft mode — normal-incidence
  coefficients charged on each link's mean projected area per Cauchy's
  theorem (total surface area / 4, exact for convex bodies; for a thin
  panel this reduces to one-sided area / 2 plus edge terms).

Convention notes (documented, defaults unchanged): the spacecraft
builder's `bus_ram_face=:legacy` reference area is inconsistent with the
Hart normalization face for ram-aligned buses (an aspect-ratio-class
effect); and `planet.R`/`planet.γ` carry sea-level air values, which
overstates the molecular speed ratio at exospheric altitudes by a few
percent (under-drag direction).
"""
@kwdef struct AerodynamicCoefficientfM <: AbstractForceTorqueModel
    # Per-model opt-in for per-link atmosphere sampling (preferred over the
    # process-wide set_per_link_atmosphere!, which leaks across simulations).
    per_link_atmosphere::Bool = false
    # Fixed-attitude incidence mode (see docstring): :max_drag | :attitude |
    # :tumbling_average. Only consulted when orientation_sim=false.
    fixed_attitude_incidence::Symbol = :max_drag
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

# The constant model's linear CD(alpha) is physical only on [0, pi/2]
# (CD 0.8 -> 2.2). The signed atan2 incidence a propagated attitude produces
# spans (-pi, pi]; fed in raw it goes as low as CD = -0.6 at alpha = -pi/2 —
# drag turned into thrust, pumping orbital energy into a tumbling spacecraft
# (Codex, PR #86). Fold by the box's x/z symmetry: reversed flow hits the
# same face. The fM model keeps the signed angles it needs.
@inline function _fold_constant_incidence(alpha_rad::Float64)::Float64
    a = abs(alpha_rad)
    return min(a, pi - a)
end

# `rot(q)` maps reference -> body, so this vector is the assumed flow
# direction (reference +x) expressed in body coordinates — NOT the body x-axis
# in the reference frame (the transpose reading flips the sign of lift).
# Uncomposed on purpose: this is the historical :max_drag child formula, which
# ignores that child quaternions are root-relative. Must stay bit-identical.
@inline function _quaternion_link_alpha(body)::Float64
    flow_body = rot(body.q) * SVector{3, Float64}(1.0, 0.0, 0.0)
    return atan(flow_body[1], flow_body[3])
end

# :attitude-mode incidence. Child quaternions are stored relative to the root
# frame (see vehicle/kinematics), so children compose with the root attitude —
# a rigidly mounted identity-quaternion child shares the root's incidence,
# matching the orientation_sim=true composition rot(child.q)*rot(root.q)*v.
@inline function _attitude_link_alpha(body, root)::Float64
    e1 = SVector{3, Float64}(1.0, 0.0, 0.0)
    flow_body = body.root ? rot(body.q) * e1 : rot(body.q) * (rot(root.q) * e1)
    return atan(flow_body[1], flow_body[3])
end

@inline function _validate_fm_incidence(incidence::Symbol)
    (incidence === :max_drag || incidence === :attitude || incidence === :tumbling_average) ||
        throw(ArgumentError(
            "AerodynamicCoefficientfM fixed_attitude_incidence must be :max_drag, :attitude, or :tumbling_average, got :$(incidence)."
        ))
    return nothing
end

# Effective per-link area for the fixed-attitude fM force: the configured
# ref_area, except in :tumbling_average mode where Cauchy's theorem gives the
# mean projected area of a convex body as total surface area / 4.
@inline function _aero_link_area(body, incidence::Symbol)::Float64
    if incidence === :tumbling_average
        d = body.dims
        return 0.5 * (d[1] * d[2] + d[1] * d[3] + d[2] * d[3])
    end
    return body.ref_area
end

@inline function _aero_link_angles(
    spacecraft,
    body,
    root_index::Int,
    orientation_sim::Bool,
    vel_pi::SVector{3, Float64},
    θ_body::Float64,
    fixed_attitude_incidence::Symbol=:max_drag,
    q_root_ib::Union{Nothing, SVector{4, Float64}}=nothing,
)::Tuple{Float64, Float64, Union{Nothing, SMatrix{3, 3, Float64, 9}}}
    if orientation_sim
        # The root attitude must be the PROPAGATED state (StateSample.q_ib),
        # not the static Link.q stored on the spacecraft model — the stored
        # quaternion never tracks the attitude ODE, so using it froze the aero
        # geometry at the initial attitude for the whole simulation. Child
        # attitudes stay the stored root-relative Link.q (articulated panels).
        q_root_ib === nothing && throw(ArgumentError(
            "orientation_sim aero evaluation requires the propagated root attitude q_root_ib."))
        R_root = rot(q_root_ib)'
        R = body.root ? SMatrix{3, 3, Float64, 9}(R_root) :
            SMatrix{3, 3, Float64, 9}(R_root * rot(SVector{4, Float64}(body.q...))')
        body_frame_velocity = R' * vel_pi
        α_body = atan(body_frame_velocity[1], body_frame_velocity[3])
        β_body = atan(body_frame_velocity[2], hypot(body_frame_velocity[1], body_frame_velocity[3]))
        return α_body, β_body, R
    end
    α_body = if fixed_attitude_incidence === :attitude
        _attitude_link_alpha(body, spacecraft.root)
    elseif fixed_attitude_incidence === :tumbling_average
        pi / 2
    else # :max_drag — historical behavior
        body.root ? (pi / 2) : _quaternion_link_alpha(body)
    end
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


# Returns (force_ii, torque_body, drag_ii, lift_ii, cross_ii): the forces in the
# inertial frame, torque_body in the ROOT BODY frame (the wrench torque contract,
# matching every other torque-returning effector).
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
    fixed_attitude_incidence::Symbol=:max_drag,
)::NTuple{5, SVector{3, Float64}}
    _validate_fm_incidence(fixed_attitude_incidence)
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
        α_body, β_body, R_body_to_inertial = _aero_link_angles(spacecraft, body, root_index, orientation_sim, vel_pi, θ_body, fixed_attitude_incidence, x.q_ib)
        link_area = orientation_sim ? body.ref_area : _aero_link_area(body, fixed_attitude_incidence)

        rho_body, T_body = rho, T
        if link_atmosphere_fn !== nothing && orientation_sim && !body.root && R_body_to_inertial !== nothing
            # body.r is a ROOT-frame offset; rotate it with the propagated root
            # attitude (not the link frame, which differs for canted panels).
            pos_ii_body = x.pos_ii + rot(x.q_ib)' * SVector{3, Float64}(body.r)
            pos_pp_body = planet_frame.l_pi * pos_ii_body
            rho_link, T_link, wind_link = link_atmosphere_fn(pos_pp_body)
            if isfinite(rho_link) && rho_link > eps(Float64) && isfinite(T_link) && T_link > 0.0
                rho_body, T_body = rho_link, T_link
                # Per-link sampling applies density/temperature only; the
                # link-local WIND sample is discarded here (Mach and dynamic
                # pressure keep the spacecraft-level wind-relative velocity).
                # Warn exactly when the discard changes the answer — the link
                # wind DISAGREES with the spacecraft-level sample. A nonzero
                # test alone misses a calm link inside a nonzero spacecraft
                # wind (the link then wrongly inherits the spacecraft wind)
                # and false-alarms when the two agree, where the physics is
                # correct (Codex reviews on PRs #62/#63).
                # The triples are local E/N/U components at each sample point,
                # but link offsets are spacecraft-scale (body.r, meters), so
                # the two bases differ by <= |body.r|/R_planet ~ 1e-6 rad —
                # comparing raw components is exact to ~ppm, and a
                # frame-resolved comparison would put a per-link geodetic
                # conversion in the aero hot loop for a maxlog=1 diagnostic.
                if wind_link !== nothing && !all(wind_link .== wind)
                    @warn "per-link atmosphere discards the link-local WIND sample where it differs from the spacecraft-level wind: Mach/dynamic pressure use the spacecraft-level wind-relative velocity (per-link sampling covers density/temperature only)." maxlog = 1
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
            0.0, _constant_drag_coefficient(_fold_constant_incidence(α_body)), 0.0
        end

        drag_pp_body = q_body * CD_body * link_area * drag_pp_hat
        lift_pp_body = lift_scale_body * CL_body * link_area * lift_pp_hat
        cross_pp_body = orientation_sim ? (q_body * CS_body * link_area * cross_pp_hat) : SVector{3, Float64}(0.0, 0.0, 0.0)
        drag_ii_body  = l_pi_t * drag_pp_body
        lift_ii_body  = l_pi_t * lift_pp_body
        cross_ii_body = l_pi_t * cross_pp_body
        force_ii  .+= drag_ii_body + lift_ii_body + cross_ii_body
        drag_ii   .+= drag_ii_body
        lift_ii   .+= lift_ii_body
        cross_ii  .+= cross_ii_body
        if orientation_sim && !(R_body_to_inertial === nothing)
            # Torque contract is the ROOT BODY frame: rotate the inertial link
            # force there and build the lever arm there. body.r is already a
            # root-frame offset (zero for the root link); the center-of-pressure
            # offset is a link-frame vector and composes through the child
            # attitude. The previous code crossed the root-frame body.r with a
            # LINK-frame force — frame-inconsistent for canted panels — and had
            # no CoP term, so root-only spacecraft could never carry aero torque.
            force_ii_link = SVector{3, Float64}(drag_ii_body + lift_ii_body + cross_ii_body)
            force_root = rot(x.q_ib) * force_ii_link
            cop = SVector{3, Float64}(body.cop_offset_b)
            lever_root = body.root ? cop :
                SVector{3, Float64}(body.r) + rot(SVector{4, Float64}(body.q...))' * cop
            torque_body .+= cross(lever_root, force_root)
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
    force, torque, _, _, _ = _aero_pure_wrench(:fm, x, env, nothing, model.fixed_attitude_incidence)
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
    force, torque, drag_ii, lift_ii, cross_ii = _aero_pure_wrench(:fm, x, env, link_atmosphere_fn, model.fixed_attitude_incidence)
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
    incidence = model.fixed_attitude_incidence
    _validate_fm_incidence(incidence)

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
    decision = _multibody_thread_decision(length(bodies); heavy_work=true)
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
            # Fixed-attitude incidence per the model's fixed_attitude_incidence
            # mode (see the AerodynamicCoefficientfM docstring). :max_drag is
            # the historical behavior: root pinned flow-normal.
            if incidence === :attitude
                b.α = _attitude_link_alpha(b, spacecraft.root)
            elseif incidence === :tumbling_average
                b.α = pi / 2
            else # :max_drag
                if b.root
                    b.α = pi / 2 # Angle of attack for the root body
                else
                    b.α = _quaternion_link_alpha(b)
                end
            end
        end

        CL_body, CD_body, CS_body, _, _, _ = aerodynamic_coefficient_fM(b, T, S)
        link_area = orientation_sim ? b.ref_area : _aero_link_area(b, incidence)

        drag_pp_body = q * CD_body * link_area * drag_pp_hat                       # Planet relative drag force vector
        lift_pp_body = lift_scale * CL_body * link_area * lift_pp_hat     # Planet relative lift force vector

        if orientation_sim
            cross_pp_body = q * CS_body * link_area * cross_pp_hat # Planet relative cross force vector
            cross_body = L_PI_t * cross_pp_body # Inertial cross force vector
        else
            cross_pp_body = zero_vec # Planet relative cross force vector
            cross_body = zero_vec # Inertial cross force vector
        end

        drag_body = L_PI_t * drag_pp_body   # Inertial drag force vector
        lift_body = L_PI_t * lift_pp_body   # Inertial lift force vector
        force_body = drag_body + lift_body + cross_body

        return force_body, drag_body, lift_body, cross_body, CL_body * link_area, CD_body * link_area, link_area
    end

    force_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    drag_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    lift_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    cross_ii = MVector{3, Float64}(0.0, 0.0, 0.0)

    started_ns = time_ns()
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
    ParallelPolicy.record_policy_observation!(
        :multibody;
        mode=decision.mode,
        num_items=length(bodies),
        use_threads=use_threads,
        elapsed_ns=(time_ns() - started_ns)
    )

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
