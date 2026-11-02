using SpecialFunctions

const sqrt_π = sqrt(π)
const inv_sqrt_π = 1 / sqrt(π)

# Historical bindings remain aliases to the neutral parallel-policy owner.
const _parse_bool_env = ParallelPolicy._parse_bool_env
const _MULTIBODY_MODE_CACHE = ParallelPolicy._MULTIBODY_MODE_CACHE
const refresh_multibody_parallel_mode! = ParallelPolicy.refresh_multibody_parallel_mode!
const _multibody_parallel_mode = ParallelPolicy._multibody_parallel_mode
const _multibody_thread_threshold = ParallelPolicy._multibody_thread_threshold
const _multibody_max_threads = ParallelPolicy._multibody_max_threads
const _threadid_capacity = ParallelPolicy._threadid_capacity
const _multibody_outer_parallel_hint = ParallelPolicy._multibody_outer_parallel_hint
const _multibody_use_threads = ParallelPolicy._multibody_use_threads
const _multibody_thread_decision = ParallelPolicy._multibody_thread_decision

# Aerodynamic models
@kwdef struct AerodynamicCoefficientConstant <: AbstractForceTorqueModel

end

"""
Aerodynamic drag model driven by articulated exposed area with a constant drag
coefficient. Uncontrolled links contribute their full reference area;
controlled links contribute `ref_area*sin(alpha)`. This matches the exposed
area and force law used by the KS aerobraking MPC.
"""
Base.@kwdef struct AerodynamicCommandedAreaDragModel <: AbstractForceTorqueModel
    drag_coefficient::Float64 = 2.2
    controlled_link_indices::Tuple{Vararg{Int}} = (2, 3)
end

"""
Free-molecular aerodynamics from the Hart et al. closed forms (rectangular
prism, doi 10.2514/1.A33606), evaluated per link and summed.

`fixed_attitude_incidence` selects the link incidence used by both aerodynamics
and heating when `orientation_sim=false`. When attitude is propagated, both
use the current spacecraft and link geometry instead:

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

@inline function _make_aero_scratch_workspace(n_links::Int)::AeroScratchWorkspace
    n_links >= 1 || throw(ArgumentError("Aerodynamic scratch workspace link count must be >= 1, got $n_links"))
    zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
    return AeroScratchWorkspace(
        fill(zero3, n_links), fill(zero3, n_links), fill(zero3, n_links), fill(zero3, n_links),
        zeros(Float64, n_links), zeros(Float64, n_links), zeros(Float64, n_links)
    )
end

@inline function _ensure_aero_workspace_capacity!(
    workspace::AeroScratchWorkspace,
    n_links::Int
)::AeroScratchWorkspace
    n_links >= 1 || throw(ArgumentError("Aerodynamic scratch workspace link count must be >= 1, got $n_links"))
    if length(workspace.link_force) < n_links
        zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
        for v in (workspace.link_force, workspace.link_drag, workspace.link_lift, workspace.link_cross)
            old_len = length(v)
            resize!(v, n_links)
            @inbounds for idx in (old_len + 1):n_links
                v[idx] = zero3
            end
        end
        for v in (workspace.link_cl_area, workspace.link_cd_area, workspace.link_area)
            old_len = length(v)
            resize!(v, n_links)
            @inbounds for idx in (old_len + 1):n_links
                v[idx] = 0.0
            end
        end
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
        if isdefined(mod, :ephemerides_frame_lpi)
            return mod
        end
        parent = parentmodule(mod)
        parent === mod && break
        mod = parent
    end
    error("SimulationModel.ephemerides_frame_lpi not found in module ancestry for aerodynamic_wrench_models.jl")
end

# Guards the lazy-growth path below, and nothing else.
#
# Not on the hot path: `_initialize_save_cache_buffers!` (engine/setup.jl) sizes
# these caches to the satellite count before any threaded dispatch can begin, so
# in a solve the growth branch is never taken and this lock is never acquired.
const _SAVE_CACHE_GROW_LOCK = ReentrantLock()

@inline function _store_vector_cache!(
    cache::Vector{SVector{3, Float64}},
    sat_idx::Int,
    value::SVector{3, Float64}
)::Nothing
    if length(cache) < sat_idx
        # Growing a SHARED vector from inside a parallel region.
        #
        # AerodynamicCoefficientfM's wrench is dispatched across persistent
        # worker threads by the flat-constellation-effector-queue route, so two
        # satellites on different workers can both see the cache too short and
        # both call `resize!` on the same Vector. Julia 1.12 catches that as a
        # ConcurrencyViolationError; at larger satellite counts, where the race
        # window is tighter, it has instead corrupted memory badly enough to
        # segfault.
        #
        # `_initialize_save_cache_buffers!` pre-sizes the caches and is called by
        # `run_simulation`, which is why solves do not hit this. But safety then
        # rests on every construction path remembering one call, and the paths
        # that build `ODEParams` directly to drive the RHS -- the cost-model
        # validation script and the contention probes -- did not. A missed call
        # should cost a lock acquisition on a cold path, not a segfault at N=1024.
        #
        # Double-checked: the length is re-tested under the lock, so concurrent
        # arrivals grow once and the losers fall through to the indexed write.
        #
        # This serialises growth against growth. It does NOT serialise growth
        # against a concurrent indexed write lower down the same vector, which
        # `resize!` can invalidate by reallocating -- guarding every write would
        # put a lock on the hot path to fix a case that should not arise. So
        # this is a survivability net, not a licence to skip initialisation:
        # pre-sizing remains the only correct configuration, and callers that
        # build `ODEParams` themselves must still call
        # `_initialize_save_cache_buffers!`.
        lock(_SAVE_CACHE_GROW_LOCK) do
            if length(cache) < sat_idx
                old_len = length(cache)
                resize!(cache, sat_idx)
                @inbounds for idx in (old_len + 1):sat_idx
                    cache[idx] = SVector{3, Float64}(0.0, 0.0, 0.0)
                end
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
@inline environment_requirements(::AerodynamicCommandedAreaDragModel) = EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
@inline environment_requirements(::AerodynamicCoefficientfM) = EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
@inline environment_requirements(::AerodynamicCoefficientNoBallisticFlight) = EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
@inline solver_partition(::AerodynamicCoefficientConstant) = :implicit
@inline solver_partition(::AerodynamicCommandedAreaDragModel) = :implicit
@inline solver_partition(::AerodynamicCoefficientfM) = :implicit
@inline solver_partition(::AerodynamicCoefficientNoBallisticFlight) = :implicit

# Internal bridge to thermal sampling. Built-in aero models own geometric
# incidence. With no recognized owner, heating retains stored-angle inputs.
@inline _thermal_incidence_mode(::Any) = nothing
@inline _thermal_incidence_mode(::AerodynamicCoefficientConstant) = :max_drag
@inline _thermal_incidence_mode(model::AerodynamicCoefficientfM) = model.fixed_attitude_incidence
@inline _thermal_incidence_mode(::AerodynamicCoefficientNoBallisticFlight) = :max_drag

# Peel the heterogeneous tuple so each effector keeps its concrete type; a
# runtime loop boxes the models on every thermal/RHS sample.
@inline _thermal_incidence_mode(effectors::Tuple, orientation_sim::Bool) =
    _thermal_incidence_mode(effectors, orientation_sim, nothing)
@inline _thermal_incidence_mode(::Tuple{}, ::Bool, mode) = mode

@inline function _thermal_incidence_mode(effectors::Tuple, orientation_sim::Bool, mode)
    candidate = _thermal_incidence_mode(first(effectors))
    if candidate !== nothing
        # Fixed-attitude policies are irrelevant when attitude is propagated.
        orientation_sim && return :attitude
        _validate_fm_incidence(candidate)
        if mode !== nothing && mode !== candidate
            throw(ArgumentError("Thermal incidence is ambiguous: aerodynamic effectors specify different fixed-attitude incidence modes."))
        end
        mode = candidate
    end
    return _thermal_incidence_mode(Base.tail(effectors), orientation_sim, mode)
end

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
    # direction and dynamic pressure; mach_body/S_body are computed per link
    # inside the loop below.

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

@inline function _commanded_area_drag_wrench(
    model::AerodynamicCommandedAreaDragModel,
    x::StateSample,
    env::EnvironmentSample,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    x.spacecraft === nothing && throw(ArgumentError(
        "Commanded-area drag requires StateSample.spacecraft."))
    planet_frame = env.planet_frame
    atmosphere = env.atmosphere
    planet_frame === nothing && throw(ArgumentError(
        "Commanded-area drag requires env.planet_frame."))
    atmosphere === nothing && throw(ArgumentError(
        "Commanded-area drag requires env.atmosphere."))
    rho = atmosphere.rho_kg_m3
    (!isfinite(rho) || rho <= eps(Float64)) && return _AERO_ZERO3, _AERO_ZERO3

    uD, uN, uE = latlongtoNED((
        planet_frame.alt_m, planet_frame.lat_rad, planet_frame.lon_rad))
    wE, wN, wU = atmosphere.wind_pp
    wind_pp = wN * uN + wE * uE - wU * uD
    relative_velocity_pp = planet_frame.vel_pp - wind_pp
    speed = norm(relative_velocity_pp)
    speed <= eps(Float64) && return _AERO_ZERO3, _AERO_ZERO3

    controlled = model.controlled_link_indices
    area_m2 = 0.0
    @inbounds for (index, link) in pairs(x.spacecraft.links)
        if index in controlled
            area_m2 += link.ref_area * sin(clamp(Float64(link.α), 0.0, pi / 2))
        else
            area_m2 += link.ref_area
        end
    end
    drag_pp = -0.5 * rho * speed * model.drag_coefficient * area_m2 *
        relative_velocity_pp
    drag_ii = planet_frame.l_pi' * drag_pp
    return SVector{3, Float64}(drag_ii), _AERO_ZERO3
end

# Direct (uncached) per-link density query used by `link_atmosphere_fn`. Bypasses the
# satellite-level GRAM tracking/extrapolation cache in
# SimulationCallbacks._density_state_from_kinematics! deliberately: that cache is keyed
# and warmed per satellite trajectory point, and reusing it for a handful of link-offset
# positions per RHS call would corrupt its trajectory-continuity assumptions. This costs
# one extra raw density-model call per non-root link per wrench evaluation.
@inline function _aero_link_atmosphere_query(p, sat_idx::Int, t::Float64, pos_pp_link::SVector{3, Float64}, planet)
    # A history-dependent model (GRAM perturbed winds) under freeze-per-step:
    # a per-link native query from inside a (possibly threaded) RHS would
    # advance the shared random walk in scheduling order. Use the satellite's
    # once-per-step sample instead, as every other RHS-side atmosphere read
    # does in that mode.
    cb_env = SimulationModel.SimulationCallbacks._callback_env_config(p)
    if cb_env.density_history_dependent && cb_env.density_freeze_per_step
        sb = p.shared_buffers
        if sat_idx <= length(sb.density_sample_t) && isfinite(sb.density_sample_t[sat_idx])
            return sb.densities[sat_idx], sb.temperatures[sat_idx], sb.winds[sat_idx]
        end
    end
    alt, lat, lon = rtolatlong(pos_pp_link, planet)
    density_model = SimulationModel.SimulationCallbacks._density_model_for_sat(p, sat_idx)
    EM = SimulationModel.EnvironmentModels
    wind_requested = EM._environment_wind_enabled(p)
    rho, T, wind_vec = SimulationModel.getDensity(density_model, alt, lat, lon, t, wind_requested, p)
    # Same opt-in GRAM perturbation factor as the satellite-level path; returns
    # (rho, T, wind) unchanged when no mode is installed.
    rho, T, wind_vec = SimulationModel.SimulationCallbacks._apply_gram_density_perturbation(p, sat_idx, t, alt, rho, T, wind_vec)
    return rho, T, EM._environment_wind(wind_requested, wind_vec)
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
    model::AerodynamicCommandedAreaDragModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return _commanded_area_drag_wrench(model, x, env)
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
    model::AerodynamicCommandedAreaDragModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
    p::ODEParams,
    sat_idx::Int,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    force, torque = _commanded_area_drag_wrench(model, x, env)
    _store_aero_caches!(p, sat_idx, force, _AERO_ZERO3, _AERO_ZERO3)
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

# ---------------------------------------------------------------------------
# Per-link kernel at a live link pose (articulated spacecraft, opt-in; see `link_wrench`)
# ---------------------------------------------------------------------------
#
# The orientation_sim branch of `_aero_pure_wrench` for ONE link, evaluated at the link's live COM
# position, velocity and attitude. The relative wind therefore includes `ω × r` of the carrying body
# (aerodynamic damping). Link `r`/`q` are never read: the pose is `xl`. `env.planet_frame` is the
# planet-frame sample at the link COM; `env.atmosphere` is the spacecraft-level sample unless the
# model's `per_link_atmosphere` opt-in is set (then a non-root link queries density and temperature at
# its own position, as `_aero_pure_wrench` does). Force at the link COM and torque about it, inertial.
function _aero_link_pure_wrench(
    coefficient_mode::Symbol,
    link,
    xl::LinkStateSample,
    env::EnvironmentSample,
    per_link_density::Bool,
    t::Float64,
    p,
    sat_idx::Int,
)::NTuple{5, SVector{3, Float64}}
    planet_frame = env.planet_frame
    atmosphere = env.atmosphere
    planet_frame === nothing && throw(ArgumentError("Aerodynamic link wrench evaluation requires env.planet_frame."))
    atmosphere === nothing && throw(ArgumentError("Aerodynamic link wrench evaluation requires env.atmosphere."))
    planet = env.planet
    rho = atmosphere.rho_kg_m3
    T = atmosphere.temperature_k
    wind = atmosphere.wind_pp
    if per_link_density && p !== nothing && !link.root
        rho_l, T_l, _ = _aero_link_atmosphere_query(p, sat_idx, t, planet_frame.pos_pp, planet)
        if isfinite(rho_l) && rho_l > eps(Float64) && isfinite(T_l) && T_l > 0.0
            rho, T = rho_l, T_l
        end
    end
    if !isfinite(rho) || rho <= eps(Float64) || !isfinite(T) || T <= 0.0
        return _AERO_ZERO5
    end
    vel_pp = planet_frame.vel_pp
    h_pp = cross(planet_frame.pos_pp, vel_pp)
    h_pp_mag = norm(h_pp)
    if !isfinite(h_pp_mag) || h_pp_mag <= eps(Float64)
        return _AERO_ZERO5
    end
    uD, uN, uE = latlongtoNED((planet_frame.alt_m, planet_frame.lat_rad, planet_frame.lon_rad))
    wE, wN, wU = wind
    wind_pp = wN * uN + wE * uE - wU * uD
    vel_pp_rw = vel_pp - wind_pp
    vel_pp_rw_mag = norm(vel_pp_rw)
    if vel_pp_rw_mag <= eps(Float64)
        return _AERO_ZERO5
    end
    vel_pp_rw_hat = vel_pp_rw / vel_pp_rw_mag
    h_pp_hat = h_pp / h_pp_mag
    lift_pp_hat = normalize(cross(h_pp_hat, vel_pp_rw_hat))
    drag_pp_hat = -vel_pp_rw_hat
    cross_pp_hat = cross(drag_pp_hat, lift_pp_hat)
    l_pi_t = planet_frame.l_pi'
    vel_pi = l_pi_t * vel_pp_rw
    θ_body = acos(clamp(vel_pp_rw[1] / vel_pp_rw_mag, -1.0, 1.0))
    R_link = SMatrix{3, 3, Float64, 9}(rot(xl.q_ib)')       # link frame -> inertial
    body_frame_velocity = R_link' * vel_pi
    α_body = atan(body_frame_velocity[1], body_frame_velocity[3])
    β_body = atan(body_frame_velocity[2], hypot(body_frame_velocity[1], body_frame_velocity[3]))
    link_area = link.ref_area
    mach_body = vel_pp_rw_mag / sqrt(planet.γ * planet.R * T)
    S_body = sqrt(planet.γ * 0.5) * mach_body
    q_body = 0.5 * rho * vel_pp_rw_mag^2
    CL_body, CD_body, CS_body = if coefficient_mode == :fm
        coeffs = aerodynamic_coefficient_fM(link, T, S_body, α_body, β_body, θ_body)
        coeffs[1], coeffs[2], coeffs[3]
    else
        0.0, _constant_drag_coefficient(_fold_constant_incidence(α_body)), 0.0
    end
    drag_ii = l_pi_t * (q_body * CD_body * link_area * drag_pp_hat)
    lift_ii = l_pi_t * (q_body * CL_body * link_area * lift_pp_hat)
    cross_ii = l_pi_t * (q_body * CS_body * link_area * cross_pp_hat)
    force_ii = drag_ii + lift_ii + cross_ii
    cop = SVector{3, Float64}(link.cop_offset_b)
    torque_ii = cross(R_link * cop, force_ii)
    return force_ii, torque_ii, drag_ii, lift_ii, cross_ii
end

@inline link_wrench_capable(::AerodynamicCoefficientConstant) = true
@inline link_wrench_capable(::AerodynamicCoefficientfM) = true
@inline link_wrench_capable(::AerodynamicCoefficientNoBallisticFlight) = true

@inline link_wrench(model::AerodynamicCoefficientConstant, link, xl::LinkStateSample, env::EnvironmentSample, t::Float64, p, sat_idx::Int) =
    _aero_link_pure_wrench(:constant, link, xl, env, _per_link_enabled(model), t, p, sat_idx)
@inline link_wrench(model::AerodynamicCoefficientfM, link, xl::LinkStateSample, env::EnvironmentSample, t::Float64, p, sat_idx::Int) =
    _aero_link_pure_wrench(:fm, link, xl, env, _per_link_enabled(model), t, p, sat_idx)
@inline link_wrench(model::AerodynamicCoefficientNoBallisticFlight, link, xl::LinkStateSample, env::EnvironmentSample, t::Float64, p, sat_idx::Int) =
    _aero_link_pure_wrench(:constant, link, xl, env, _per_link_enabled(model), t, p, sat_idx)

@inline link_wrench_store!(::Union{AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight},
    p::ODEParams, sat_idx::Int, drag_ii, lift_ii, cross_ii) = _store_aero_caches!(p, sat_idx, drag_ii, lift_ii, cross_ii)


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
