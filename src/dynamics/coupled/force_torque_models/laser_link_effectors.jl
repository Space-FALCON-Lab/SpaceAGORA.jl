module LaserLinkEffectors

using LinearAlgebra
using StaticArrays
using DiffEqBase

using ...AbstractTypes: AbstractForceTorqueModel
import ..DynamicEffectors: calcForceTorque, solver_partition

export OpenCavityLaserLinkModel, laser_link_scheduler_callback
export laser_link_force_magnitude, laser_link_pair_force, laser_link_active_pair
export update_laser_link_schedule!, accumulate_laser_link_forces!
export LaserImpulseTracker, laser_impulse_callback, tracked_dv_at

const SPEED_OF_LIGHT_MPS = 299_792_458.0
const SUPPORTED_LASER_SCHEDULES = (
    :naive_next_entering,
    :positive_along_track,
    :gve_sma,   # maximise semi-major axis rate
    :gve_ecc,   # maximise eccentricity rate
    :gve_inc,   # maximise inclination rate
    :gve_raan,  # maximise RAAN rate
    :gve_argp,  # maximise argument-of-periapsis rate
    :external_schedule,  # explicit precomputed (interval × helper) on/off table; supports multiple simultaneous helpers
)

"""
    OpenCavityLaserLinkModel(; kwargs...)

Open-cavity laser interlink force model for one target satellite and a helper
population. The model applies at most one active target-helper link at a time,
with equal-and-opposite inertial-frame forces along the instantaneous
inter-satellite line of sight.
"""
mutable struct OpenCavityLaserLinkModel <: AbstractForceTorqueModel
    target_idx::Int
    helper_indices::Vector{Int}
    range_m::Float64
    power_w::Float64
    magnification::Float64
    beta::Float64
    eta::Float64
    schedule::Symbol
    active_helper_idx::Int
    active_helper_indices::Vector{Int}   # used by :external_schedule; may hold >1 simultaneously-active helper
    schedule_starts::Vector{Float64}     # :external_schedule interval start times [s], sorted ascending
    schedule_active::BitMatrix           # :external_schedule (n_intervals × n_helpers) on/off table
    previous_in_range::Vector{Bool}
    link_activation_count::Int
    active_link_step_count::Int
end

function OpenCavityLaserLinkModel(;
    target_idx::Integer=1,
    helper_indices::AbstractVector{<:Integer}=Int[],
    range_m::Real=200e3,
    power_w::Real=10_000.0,
    magnification::Real=100.0,
    beta::Real=1.0,
    eta::Real=1.0,
    schedule::Symbol=:naive_next_entering,
    active_helper_idx::Integer=0,
    active_helper_indices::AbstractVector{<:Integer}=Int[],
    schedule_starts::AbstractVector{<:Real}=Float64[],
    schedule_active::AbstractMatrix{Bool}=falses(0, 0),
    previous_in_range::AbstractVector{Bool}=Bool[],
    link_activation_count::Integer=0,
    active_link_step_count::Integer=0,
)
    model = OpenCavityLaserLinkModel(
        Int(target_idx),
        collect(Int, helper_indices),
        Float64(range_m),
        Float64(power_w),
        Float64(magnification),
        Float64(beta),
        Float64(eta),
        schedule,
        Int(active_helper_idx),
        collect(Int, active_helper_indices),
        collect(Float64, schedule_starts),
        BitMatrix(schedule_active),
        collect(Bool, previous_in_range),
        Int(link_activation_count),
        Int(active_link_step_count),
    )
    _validate_laser_link_model!(model)
    _ensure_laser_link_state!(model)
    return model
end

# Validates all model fields at construction and before scheduling; throws ArgumentError on any violation.
# Input: model instance.
# Output: nothing (throws on invalid state).
function _validate_laser_link_model!(model::OpenCavityLaserLinkModel)::Nothing
    model.target_idx > 0 || throw(ArgumentError("OpenCavityLaserLinkModel.target_idx must be positive."))
    isempty(model.helper_indices) && throw(ArgumentError("OpenCavityLaserLinkModel requires at least one helper index."))
    any(==(model.target_idx), model.helper_indices) &&
        throw(ArgumentError("OpenCavityLaserLinkModel helper_indices cannot include target_idx."))
    length(unique(model.helper_indices)) == length(model.helper_indices) ||
        throw(ArgumentError("OpenCavityLaserLinkModel helper_indices must be unique."))
    model.range_m >= 0.0 || throw(ArgumentError("OpenCavityLaserLinkModel.range_m must be >= 0.0."))
    model.power_w >= 0.0 || throw(ArgumentError("OpenCavityLaserLinkModel.power_w must be >= 0.0."))
    model.magnification >= 0.0 || throw(ArgumentError("OpenCavityLaserLinkModel.magnification must be >= 0.0."))
    model.beta >= 0.0 || throw(ArgumentError("OpenCavityLaserLinkModel.beta must be >= 0.0."))
    model.eta >= 0.0 || throw(ArgumentError("OpenCavityLaserLinkModel.eta must be >= 0.0."))
    model.schedule in SUPPORTED_LASER_SCHEDULES ||
        throw(ArgumentError("Unsupported OpenCavityLaserLinkModel.schedule=$(repr(model.schedule)). Supported schedules: $(SUPPORTED_LASER_SCHEDULES)."))
    if model.schedule === :external_schedule
        isempty(model.schedule_starts) &&
            throw(ArgumentError("OpenCavityLaserLinkModel schedule=:external_schedule requires non-empty schedule_starts."))
        issorted(model.schedule_starts) ||
            throw(ArgumentError("OpenCavityLaserLinkModel.schedule_starts must be sorted ascending."))
        size(model.schedule_active) == (length(model.schedule_starts), length(model.helper_indices)) ||
            throw(ArgumentError("OpenCavityLaserLinkModel.schedule_active must be sized (n_intervals, n_helpers) = ($(length(model.schedule_starts)), $(length(model.helper_indices)))."))
    end
    return nothing
end

# Ensures previous_in_range has the correct length, resizing and zeroing if the helper count changed.
# Input: model instance.
# Output: nothing (mutates model.previous_in_range in-place).
@inline function _ensure_laser_link_state!(model::OpenCavityLaserLinkModel)::Nothing
    if length(model.previous_in_range) != length(model.helper_indices)
        resize!(model.previous_in_range, length(model.helper_indices))
        fill!(model.previous_in_range, false)
    end
    return nothing
end

# Declares this effector as explicit for IMEX splitting — laser thrust is non-stiff.
# Input: OpenCavityLaserLinkModel (type dispatch only).
# Output: Symbol :explicit.
@inline solver_partition(::OpenCavityLaserLinkModel) = :explicit

# Required interface stub — laser force is applied via callback, not through the ODE RHS.
# Input: spacecraft state view x, parameter struct p, satellite index i.
# Output: zero force and zero torque (3-vectors each); callback handles the real kick.
function calcForceTorque(
    ::OpenCavityLaserLinkModel,
    x,
    p,
    i::Int64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end


# Extracts position and velocity arrays for all spacecraft from the ODE state.
# Input: ODE state u with u.sc array of per-spacecraft state views.
# Output: two Vector{SVector{3,Float64}} — positions and velocities for all N spacecraft.
function _state_vectors(u)::Tuple{Vector{SVector{3, Float64}}, Vector{SVector{3, Float64}}}
    sc_state = u.sc
    pos = Vector{SVector{3, Float64}}(undef, length(sc_state))
    vel = Vector{SVector{3, Float64}}(undef, length(sc_state))
    @inbounds for idx in eachindex(sc_state)
        sc = sc_state[idx]
        pos[idx] = SVector{3, Float64}(sc[1], sc[2], sc[3])
        vel[idx] = SVector{3, Float64}(sc[4], sc[5], sc[6])
    end
    return pos, vel
end

# Maps a global spacecraft index to its position in model.helper_indices.
# Input: model, helper_idx (global spacecraft index).
# Output: slot integer (1-based position in helper_indices), or nothing if not found.
@inline _helper_slot(model::OpenCavityLaserLinkModel, helper_idx::Int)::Union{Nothing, Int} =
    findfirst(==(helper_idx), model.helper_indices)

# Determines which helpers are currently within laser range of the target.
# Input: reusable flags buffer, model (for target index and range_m), pos array for all spacecraft.
# Output: flags Vector{Bool} mutated in-place — true means that helper slot is in range.
function _in_range_flags!(
    flags::Vector{Bool},
    model::OpenCavityLaserLinkModel,
    pos::AbstractVector{SVector{3, Float64}},
)::Vector{Bool}
    _ensure_laser_link_state!(model)
    if length(flags) != length(model.helper_indices)
        resize!(flags, length(model.helper_indices))
    end
    target_pos = pos[model.target_idx]
    @inbounds for slot in eachindex(model.helper_indices)
        helper_idx = model.helper_indices[slot]
        flags[slot] = norm(pos[helper_idx] - target_pos) <= model.range_m
    end
    return flags
end

# Computes the Radial-Transverse-Normal (RTN) unit vector triad for a spacecraft.
# Input: target position pos_t and velocity vel_t in ECI frame (SVector{3}).
# Output: rhat (radial, away from Earth), that (along-track), nhat (orbit-normal) — all unit SVector{3}.
@inline function _rtn_basis(pos_t::SVector{3, Float64}, vel_t::SVector{3, Float64})
    rnorm = norm(pos_t)
    rhat = rnorm > 0.0 ? pos_t / rnorm : SVector{3, Float64}(1.0, 0.0, 0.0)
    h = cross(pos_t, vel_t)
    hnorm = norm(h)
    nhat = hnorm > 0.0 ? h / hnorm : SVector{3, Float64}(0.0, 0.0, 1.0)
    that = cross(nhat, rhat)
    return rhat, that, nhat
end

# Scores how well a given helper's laser direction aligns with the target's along-track axis.
# Input: model, helper_idx, pos and vel arrays for all spacecraft.
# Output: scalar projection — positive means the laser pushes the target in the +T (prograde) direction.
@inline function _along_track_projection(
    model::OpenCavityLaserLinkModel,
    helper_idx::Int,
    pos::AbstractVector{SVector{3, Float64}},
    vel::AbstractVector{SVector{3, Float64}},
)::Float64
    target_pos = pos[model.target_idx]
    rel = target_pos - pos[helper_idx]   # force direction: helper → target
    rho = norm(rel)
    rho > 0.0 || return 0.0
    _, that, _ = _rtn_basis(target_pos, vel[model.target_idx])
    return dot(rel / rho, that)
end

const _MU_EARTH_GVE = 3.986004418e14  # Earth gravitational parameter [m³/s²]

# Scores the instantaneous rate of change of orbital element `elem` if the laser fires to `helper_idx`.
# Input: elem symbol (:gve_sma/:gve_ecc/:gve_inc/:gve_raan/:gve_argp), model, helper_idx, pos/vel arrays.
# Output: scalar GVE score — higher means this helper maximises the chosen element's rate of change.
function _gve_score(
    elem::Symbol,
    model::OpenCavityLaserLinkModel,
    helper_idx::Int,
    pos::AbstractVector{SVector{3, Float64}},
    vel::AbstractVector{SVector{3, Float64}},
)::Float64
    tgt_pos = pos[model.target_idx]
    tgt_vel = vel[model.target_idx]

    # Unit force direction on target (helper → target)
    rel = tgt_pos - pos[helper_idx]
    rho = norm(rel)
    rho > 0.0 || return 0.0
    f̂ = rel / rho

    # RTN decomposition of force direction
    rhat, that, nhat = _rtn_basis(tgt_pos, tgt_vel)
    aR = dot(f̂, rhat)
    aT = dot(f̂, that)
    aN = dot(f̂, nhat)

    # --- Orbital state of target ---
    r  = norm(tgt_pos)
    r > 0.0 || return 0.0
    v2 = dot(tgt_vel, tgt_vel)
    a  = -_MU_EARTH_GVE / (v2 - 2.0 * _MU_EARTH_GVE / r)   # vis-viva
    a > 0.0 || return 0.0                                     # skip hyperbolic

    h_vec  = cross(tgt_pos, tgt_vel)
    h_sq   = dot(h_vec, h_vec)
    h_norm = sqrt(h_sq)
    h_norm > 0.0 || return 0.0
    p_slr  = h_sq / _MU_EARTH_GVE                            # semi-latus rectum

    e_vec = cross(tgt_vel, h_vec) / _MU_EARTH_GVE - tgt_pos / r
    e     = norm(e_vec)
    e_sq  = clamp(e * e, 0.0, 1.0 - 1e-12)
    sqrt_1me2 = sqrt(1.0 - e_sq)
    n_mean    = sqrt(_MU_EARTH_GVE / (a * a * a))

    # True anomaly ν
    ν = acos(clamp(dot(e_vec / max(e, 1e-12), tgt_pos / r), -1.0, 1.0))
    dot(tgt_pos, tgt_vel) < 0.0 && (ν = 2π - ν)

    # --- GVE scoring per element ---
    if elem === :gve_sma
        # ȧ = 2/(n√(1-e²)) * (e sinν · aR  +  p/r · aT)
        return (2.0 / (n_mean * sqrt_1me2)) * (e * sin(ν) * aR + (p_slr / r) * aT)

    elseif elem === :gve_ecc
        # ė = √(1-e²)/(na) * [sinν · aR  +  (cosν + (e+cosν)/(1+e cosν)) · aT]
        coeff_T = cos(ν) + (e + cos(ν)) / (1.0 + e * cos(ν))
        return (sqrt_1me2 / (n_mean * a)) * (sin(ν) * aR + coeff_T * aT)

    else
        # :gve_inc, :gve_raan, :gve_argp all need inclination and argument of latitude u
        i_rad  = acos(clamp(h_vec[3] / h_norm, -1.0, 1.0))
        sin_i  = sin(i_rad)
        cos_i  = cos(i_rad)

        # Ascending node vector (n_asc = k̂ × ĥ)
        n_asc  = cross(SVector(0.0, 0.0, 1.0), h_vec)
        n_mag  = norm(n_asc)

        # Argument of latitude u = ν + ω
        u = if n_mag > 1e-12 && e > 1e-12
            # General case
            ω = acos(clamp(dot(n_asc / n_mag, e_vec / e), -1.0, 1.0))
            e_vec[3] < 0.0 && (ω = 2π - ω)
            ν + ω
        elseif n_mag > 1e-12
            # Circular orbit: use angle from node to position
            u_tmp = acos(clamp(dot(n_asc / n_mag, tgt_pos / r), -1.0, 1.0))
            tgt_pos[3] < 0.0 ? 2π - u_tmp : u_tmp
        else
            # Equatorial: use true longitude from x-axis
            atan(tgt_pos[2], tgt_pos[1])
        end

        denom = n_mean * a * a * sqrt_1me2

        if elem === :gve_inc
            # i̇ = r cos(u) / (na²√(1-e²)) · aN
            return (r * cos(u) / denom) * aN

        elseif elem === :gve_raan
            # Ω̇ = r sin(u) / (na²√(1-e²) sin i) · aN   [singular at i = 0]
            abs(sin_i) < 1e-6 && return 0.0
            return (r * sin(u) / (denom * sin_i)) * aN

        else  # :gve_argp
            # ω̇ = √(1-e²)/(nae) [-cosν · aR + (1+r/p) sinν · aT]
            #       - r sin(u) cos(i) / (na²√(1-e²) sin i) · aN
            abs(e) < 1e-6    && return 0.0   # circular: ω undefined
            abs(sin_i) < 1e-6 && return 0.0   # equatorial: ω undefined
            term_RT = (sqrt_1me2 / (n_mean * a * e)) *
                      (-cos(ν) * aR + (1.0 + r / p_slr) * sin(ν) * aT)
            term_N  = -(r * sin(u) * cos_i / (denom * sin_i)) * aN
            return term_RT + term_N
        end
    end
end

# Returns the closest in-range helper; when entering_only=true, skips helpers already in range last step.
# Input: model, pos array for all spacecraft, in_range Bool flags, entering_only flag.
# Output: global spacecraft index of the best matching helper, or 0 if none.
function _closest_helper(
    model::OpenCavityLaserLinkModel,
    pos::AbstractVector{SVector{3, Float64}},
    in_range::AbstractVector{Bool};
    entering_only::Bool=false,
)::Int
    target_pos = pos[model.target_idx]
    best_idx = 0
    best_range = Inf
    @inbounds for slot in eachindex(model.helper_indices)
        in_range[slot] || continue
        entering_only && model.previous_in_range[slot] && continue
        helper_idx = model.helper_indices[slot]
        rho = norm(pos[helper_idx] - target_pos)
        if rho < best_range
            best_idx = helper_idx
            best_range = rho
        end
    end
    return best_idx
end

# Sets the active helper, incrementing link_activation_count if the selection changes.
# Input: model, helper_idx (0 means deactivate).
# Output: nothing (mutates model.active_helper_idx and model.link_activation_count).
function _activate_helper!(model::OpenCavityLaserLinkModel, helper_idx::Int)::Nothing
    if helper_idx > 0 && model.active_helper_idx != helper_idx
        model.link_activation_count += 1
    end
    model.active_helper_idx = helper_idx
    return nothing
end

# Runs the scheduling policy to select which helper fires this step, mutating model.active_helper_idx.
# Input: model, pos and vel arrays for all spacecraft (pre-extracted from ODE state).
# Output: active_helper_idx after update (0 if none selected).
function update_laser_link_schedule!(
    model::OpenCavityLaserLinkModel,
    pos::AbstractVector{SVector{3, Float64}},
    vel::AbstractVector{SVector{3, Float64}};
    t::Float64=NaN,
)::Int
    _validate_laser_link_model!(model)
    in_range = _in_range_flags!(Bool[], model, pos)

    if model.schedule === :external_schedule
        # Precomputed (e.g. MPC) plan is the sole authority — no in-range/LOS re-gating here.
        empty!(model.active_helper_indices)
        if !isnan(t) && !isempty(model.schedule_starts)
            k = clamp(searchsortedlast(model.schedule_starts, t), 1, length(model.schedule_starts))
            @inbounds for slot in eachindex(model.helper_indices)
                model.schedule_active[k, slot] || continue
                push!(model.active_helper_indices, model.helper_indices[slot])
            end
        end
        model.active_helper_idx = isempty(model.active_helper_indices) ? 0 : model.active_helper_indices[1]
    elseif model.schedule === :naive_next_entering
        if model.active_helper_idx > 0
            slot = _helper_slot(model, model.active_helper_idx)
            if slot === nothing || !in_range[slot]
                model.active_helper_idx = 0
                helper_idx = _closest_helper(model, pos, in_range; entering_only=true)
                helper_idx > 0 && _activate_helper!(model, helper_idx)
            end
        else
            helper_idx = any(model.previous_in_range) ?
                _closest_helper(model, pos, in_range; entering_only=true) :
                _closest_helper(model, pos, in_range)
            helper_idx > 0 && _activate_helper!(model, helper_idx)
        end
    elseif model.schedule === :positive_along_track
        if model.active_helper_idx > 0
            slot = _helper_slot(model, model.active_helper_idx)
            valid = slot !== nothing &&
                in_range[slot] &&
                _along_track_projection(model, model.active_helper_idx, pos, vel) > 0.0
            valid || (model.active_helper_idx = 0)
        end
        if model.active_helper_idx == 0
            best_idx = 0
            best_projection = 0.0
            @inbounds for slot in eachindex(model.helper_indices)
                in_range[slot] || continue
                helper_idx = model.helper_indices[slot]
                projection = _along_track_projection(model, helper_idx, pos, vel)
                if projection > best_projection
                    best_projection = projection
                    best_idx = helper_idx
                end
            end
            best_idx > 0 && _activate_helper!(model, best_idx)
        end
    elseif model.schedule === :maximize_sma || model.schedule in (:gve_sma, :gve_ecc, :gve_inc, :gve_raan, :gve_argp)
        # GVE-optimal scheduling: always fire the helper that maximises the chosen element rate.
        # No sticky retention — globally re-evaluated every scheduler step.
        elem = model.schedule === :maximize_sma ? :gve_sma : model.schedule
        best_idx   = 0
        best_score = 0.0
        @inbounds for slot in eachindex(model.helper_indices)
            in_range[slot] || continue
            helper_idx = model.helper_indices[slot]
            score = _gve_score(elem, model, helper_idx, pos, vel)
            if score > best_score
                best_score = score
                best_idx   = helper_idx
            end
        end
        _activate_helper!(model, best_idx)   # switches to best, or deactivates if none qualify
    end

    model.previous_in_range .= in_range
    has_active = model.schedule === :external_schedule ?
        !isempty(model.active_helper_indices) : model.active_helper_idx > 0
    has_active && (model.active_link_step_count += 1)
    return model.active_helper_idx
end

# Convenience overload that unpacks the ODE state u before calling the pos/vel overload.
# Input: model, ODE state u (integrator.u), optional simulation time t (required by :external_schedule).
# Output: active_helper_idx after update (0 if none selected).
function update_laser_link_schedule!(model::OpenCavityLaserLinkModel, u, t::Real=NaN)::Int
    pos, vel = _state_vectors(u)
    return update_laser_link_schedule!(model, pos, vel; t=Float64(t))
end

# Returns the (target_idx, active_helper_idx) pair identifying the current active link.
# Input: model.
# Output: Tuple{Int,Int} — (target index, active helper index); second element is 0 if inactive.
@inline laser_link_active_pair(model::OpenCavityLaserLinkModel)::Tuple{Int, Int} =
    (model.target_idx, model.active_helper_idx)


# Computes the scalar laser force magnitude F = η·β·M·P/c.
# Input: model (reads eta, beta, magnification, power_w).
# Output: force magnitude in Newtons.
@inline laser_link_force_magnitude(model::OpenCavityLaserLinkModel)::Float64 =
    model.eta * model.beta * model.magnification * model.power_w / SPEED_OF_LIGHT_MPS

# Computes the laser force vector on the target from a given helper position.
# Input: model, target_pos and helper_pos as SVector{3,Float64} in ECI frame.
# Output: force SVector{3,Float64} in direction helper→target, zero if out of range.
function laser_link_pair_force(
    model::OpenCavityLaserLinkModel,
    target_pos::SVector{3, Float64},
    helper_pos::SVector{3, Float64},
)::SVector{3, Float64}
    rel = target_pos - helper_pos
    rho = norm(rel)
    (rho > 0.0 && rho <= model.range_m) || return SVector{3, Float64}(0.0, 0.0, 0.0)
    return laser_link_force_magnitude(model) * rel / rho
end

# Adds equal-and-opposite laser forces for the active link into the shared force totals matrix.
# Input: totals matrix (6×N), model, pos array, active_flags (Bool per satellite).
# Output: nothing (mutates totals in-place; rows 1-3 are force components per satellite column).
function accumulate_laser_link_forces!(
    totals::AbstractMatrix{Float64},
    model::OpenCavityLaserLinkModel,
    pos::AbstractVector{SVector{3, Float64}},
    active_flags,
)::Nothing
    target_idx = model.target_idx
    target_idx <= length(pos) || return nothing
    active_flags[target_idx] || return nothing

    helper_ids = model.schedule === :external_schedule ? model.active_helper_indices :
        (model.active_helper_idx > 0 ? (model.active_helper_idx,) : ())
    @inbounds for helper_idx in helper_ids
        helper_idx <= length(pos) || continue
        active_flags[helper_idx] || continue
        force = laser_link_pair_force(model, pos[target_idx], pos[helper_idx])
        totals[1, target_idx] += force[1]
        totals[2, target_idx] += force[2]
        totals[3, target_idx] += force[3]
        totals[1, helper_idx] -= force[1]
        totals[2, helper_idx] -= force[2]
        totals[3, helper_idx] -= force[3]
    end
    return nothing
end

# Finds laser models matching the template in the integrator's effector list and runs their scheduler.
# Input: template model (used to match by identity or target/helper indices), integrator.
# Output: nothing (mutates matching models' active_helper_idx via update_laser_link_schedule!).
function _update_matching_laser_models!(template::OpenCavityLaserLinkModel, integrator)::Nothing
    for effector in integrator.p.args.dynamics_model.dynamic_effectors
        effector isa OpenCavityLaserLinkModel || continue
        if effector === template ||
           (effector.target_idx == template.target_idx && effector.helper_indices == template.helper_indices)
            update_laser_link_schedule!(effector, integrator.u, integrator.t)
        end
    end
    return nothing
end

# Builds a DiscreteCallback that runs the link scheduler at every accepted ODE step.
# Input: model (captured by closure; used to find and update matching effectors in the integrator).
# Output: DiffEqBase.DiscreteCallback that mutates model.active_helper_idx each step.
function laser_link_scheduler_callback(model::OpenCavityLaserLinkModel)
    condition(u, t, integrator) = true
    affect!(integrator) = _update_matching_laser_models!(model, integrator)
    initialize = (cb, u, t, integrator) -> _update_matching_laser_models!(model, integrator)
    return DiffEqBase.DiscreteCallback(condition, affect!; initialize=initialize)
end

# Running accumulator for laser ΔV in RTN components; also stores the full time-series history.
Base.@kwdef mutable struct LaserImpulseTracker
    t_prev::Float64            = 0.0
    dv_R::Float64              = 0.0
    dv_T::Float64              = 0.0
    dv_N::Float64              = 0.0
    t_hist::Vector{Float64}    = Float64[]
    dv_R_hist::Vector{Float64} = Float64[]
    dv_T_hist::Vector{Float64} = Float64[]
    dv_N_hist::Vector{Float64} = Float64[]
end

# Builds a DiscreteCallback that applies a discrete velocity kick and accumulates RTN ΔV each step.
# Input: model (for active helper and force params), tracker (accumulates ΔV history), mass_kg of target.
# Output: DiffEqBase.DiscreteCallback that mutates integrator.u velocities and tracker fields.
function laser_impulse_callback(
    model::OpenCavityLaserLinkModel,
    tracker::LaserImpulseTracker,
    mass_kg::Float64,
)
    function affect!(integrator)
        dt = integrator.t - tracker.t_prev
        if dt > 0.0
            helper_ids = model.schedule === :external_schedule ? model.active_helper_indices :
                (model.active_helper_idx > 0 ? (model.active_helper_idx,) : ())
            for helper_idx in helper_ids
                sc      = integrator.u.sc
                tgt_pos = SVector{3, Float64}(sc[model.target_idx].pos)
                tgt_vel = SVector{3, Float64}(sc[model.target_idx].vel)
                hlp_pos = SVector{3, Float64}(sc[helper_idx].pos)
                rel     = tgt_pos - hlp_pos
                rho     = norm(rel)
                if rho > 0.0 && rho <= model.range_m
                    F_mag = model.eta * model.beta * model.magnification *
                            model.power_w / SPEED_OF_LIGHT_MPS
                    force = F_mag * rel / rho
                    rhat, that, nhat = _rtn_basis(tgt_pos, tgt_vel)
                    accel = force / mass_kg
                    tracker.dv_R += dot(accel, rhat) * dt
                    tracker.dv_T += dot(accel, that) * dt
                    tracker.dv_N += dot(accel, nhat) * dt
                    # apply equal-and-opposite velocity kick to ODE state
                    dv = accel * dt
                    integrator.u.sc[model.target_idx].vel .+= dv
                    integrator.u.sc[helper_idx].vel .-= dv
                    DiffEqBase.u_modified!(integrator, true)
                end
            end
        end
        tracker.t_prev = integrator.t
        push!(tracker.t_hist,    integrator.t)
        push!(tracker.dv_R_hist, tracker.dv_R)
        push!(tracker.dv_T_hist, tracker.dv_T)
        push!(tracker.dv_N_hist, tracker.dv_N)
    end
    return DiffEqBase.DiscreteCallback(
        (u, t, integrator) -> true,
        affect!;
        save_positions=(false, false),
    )
end

# Looks up the cumulative RTN ΔV at time t by binary search in the tracker history.
# Input: tracker, query time t in seconds.
# Output: (dv_R, dv_T, dv_N) tuple in m/s at the last recorded step ≤ t.
function tracked_dv_at(tracker::LaserImpulseTracker, t::Float64)
    isempty(tracker.t_hist) && return (0.0, 0.0, 0.0)
    k = clamp(searchsortedlast(tracker.t_hist, t), 1, length(tracker.t_hist))
    return (tracker.dv_R_hist[k], tracker.dv_T_hist[k], tracker.dv_N_hist[k])
end

end
