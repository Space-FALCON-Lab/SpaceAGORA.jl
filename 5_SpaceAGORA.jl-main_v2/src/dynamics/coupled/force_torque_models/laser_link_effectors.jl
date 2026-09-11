module LaserLinkEffectors

using LinearAlgebra
using StaticArrays
using DiffEqBase

using ...AbstractTypes: AbstractForceTorqueModel
using ...EffectorSampling: StateSample, EnvironmentSample
import ..DynamicEffectors: wrench

include(joinpath(@__DIR__, "..", "..", "..", "core", "interfaces", "reference_system.jl"))

export OpenCavityLaserLinkModel_struct, build_OpenCavityLaserLinkModel, laser_link_scheduler_callback
export choose_active_helper!
export LaserImpulseTracker, laser_impulse_callback

const SPEED_OF_LIGHT_MPS = 299_792_458.0
const _MU_EARTH_GVE = 3.986004418e14  # Earth gravitational parameter [m³/s²]

# ── Model construction ──────────────────────────────────────────────────────
mutable struct OpenCavityLaserLinkModel_struct <: AbstractForceTorqueModel
    target_idx::Int
    helper_indices::Vector{Int}
    range_m::Float64
    power_w::Float64
    magnification::Float64
    beta::Float64
    eta::Float64
    schedule::Symbol
    active_helper_idx::Int
    previous_in_range::Vector{Bool}
    active_link_step_count::Int
end

function build_OpenCavityLaserLinkModel(;
    target_idx::Integer=1,
    helper_indices::AbstractVector{<:Integer}=Int[],
    range_m::Real=200e3,
    power_w::Real=10_000.0,
    magnification::Real=100.0,
    beta::Real=1.0,
    eta::Real=2.0,
    schedule::Symbol=:naive_next_entering,
    active_helper_idx::Integer=0,
    previous_in_range::AbstractVector{Bool}=Bool[],
    active_link_step_count::Integer=0,
)
    # Step 1: build the model, copying the two array fields so the model owns its own storage.
    model = OpenCavityLaserLinkModel_struct(
        target_idx,
        collect(Int, helper_indices),
        range_m,
        power_w,
        magnification,
        beta,
        eta,
        schedule,
        active_helper_idx,
        collect(Bool, previous_in_range),
        active_link_step_count,
    )
    # Step 2: initialize the model's helper-state buffer, resizing and zeroing if the helper count changed.
    if length(model.previous_in_range) != length(model.helper_indices)
        resize!(model.previous_in_range, length(model.helper_indices))
        fill!(model.previous_in_range, false)
    end
    # Step 3: return the ready-to-use model.
    return model
end

# Required interface stub — laser force is applied via callback, not through the ODE RHS.
# Input: spacecraft state/environment samples (type dispatch only).
# Output: zero force and zero torque (3-vectors each); callback handles the real kick.
@inline function wrench(
    ::OpenCavityLaserLinkModel_struct,
    ::StateSample,
    ::EnvironmentSample,
    ::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    # Step 1: return zero force and torque because the callback applies the kick.
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end


# ── Helper-scoring (one scoring function per scheduling policy) ────────────────
# Returns the closest in-range helper; when entering_only=true, skips helpers already in range last step.
# Input: model, pos array for all spacecraft, in_range Bool flags, entering_only flag.
# Output: global spacecraft index of the best matching helper, or 0 if none.
function _closest_helper(
    model::OpenCavityLaserLinkModel_struct,
    pos::AbstractVector{SVector{3, Float64}},
    in_range::AbstractVector{Bool};
    entering_only::Bool=false,
)::Int
    # Step 1: scan in-range helpers and track the closest eligible candidate.
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
    # Step 2: return the best helper index, or zero when none qualifies.
    return best_idx
end

# Scores how well a given helper's laser direction aligns with the target's along-track axis.
# Input: model, helper_idx, pos and vel arrays for all spacecraft.
# Output: scalar projection — positive means the laser pushes the target in the +T (prograde) direction.
@inline function _along_track_projection(
    model::OpenCavityLaserLinkModel_struct,
    helper_idx::Int,
    pos::AbstractVector{SVector{3, Float64}},
    vel::AbstractVector{SVector{3, Float64}},
)::Float64
    # Step 1: compute the helper-to-target line-of-sight vector.
    target_pos = pos[model.target_idx]
    rel = target_pos - pos[helper_idx]   # force direction: helper → target
    rho = norm(rel)
    rho > 0.0 || return 0.0
    # Step 2: project that direction onto the target's along-track axis.
    that = rtn_dcm_from_inertial(target_pos, vel[model.target_idx])[:, 2]
    return dot(rel / rho, that)
end

# Scores the instantaneous rate of change of orbital element `elem` if the laser fires to `helper_idx`.
# Input: elem symbol (:gve_sma/:gve_ecc/:gve_inc/:gve_raan/:gve_argp), model, helper_idx, pos/vel arrays.
# Output: scalar GVE score — higher means this helper maximises the chosen element's rate of change.
function _gve_score(
    elem::Symbol,
    model::OpenCavityLaserLinkModel_struct,
    helper_idx::Int,
    pos::AbstractVector{SVector{3, Float64}},
    vel::AbstractVector{SVector{3, Float64}},
)::Float64
    # Step 1: compute the unit laser direction and its RTN components.
    tgt_pos = pos[model.target_idx]
    tgt_vel = vel[model.target_idx]

    # Unit force direction on target (helper → target)
    rel = tgt_pos - pos[helper_idx]
    rho = norm(rel)
    rho > 0.0 || return 0.0
    f̂ = rel / rho

    # RTN decomposition of force direction
    C = rtn_dcm_from_inertial(tgt_pos, tgt_vel)
    rhat, that, nhat = C[:, 1], C[:, 2], C[:, 3]
    aR = dot(f̂, rhat)
    aT = dot(f̂, that)
    aN = dot(f̂, nhat)

    # Step 2: compute the target's orbital elements and anomaly.
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

    # Step 3: calculate the selected Gauss variational-equation score.
    # Shared inclination/argument-of-latitude terms needed by gve_inc, gve_raan, gve_argp.
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

    if elem === :gve_sma
        # ȧ = 2/(n√(1-e²)) * (e sinν · aR  +  p/r · aT)
        return (2.0 / (n_mean * sqrt_1me2)) * (e * sin(ν) * aR + (p_slr / r) * aT)

    elseif elem === :gve_ecc
        # ė = √(1-e²)/(na) * [sinν · aR  +  (cosν + (e+cosν)/(1+e cosν)) · aT]
        coeff_T = cos(ν) + (e + cos(ν)) / (1.0 + e * cos(ν))
        return (sqrt_1me2 / (n_mean * a)) * (sin(ν) * aR + coeff_T * aT)

    elseif elem === :gve_inc
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

# ── Scheduling: choose which helper fires each step, and wire it into the ODE solver ──
# Determines which helpers are currently within laser range of the target.
# Input: reusable flags buffer, model (for target index and range_m), pos array for all spacecraft.
# Output: flags Vector{Bool} mutated in-place — true means that helper slot is in range.
function _in_range_flags!(
    flags::Vector{Bool},
    model::OpenCavityLaserLinkModel_struct,
    pos::AbstractVector{SVector{3, Float64}},
)::Vector{Bool}
    # Step 1: resize the reusable flags buffer to match the helpers.
    if length(flags) != length(model.helper_indices)                            # buffer length changed since last call?
        resize!(flags, length(model.helper_indices))                            # grow/shrink it to match the helper count
    end
    # Step 2: mark each helper that is within the configured range.
    target_pos = pos[model.target_idx]                                          # the target spacecraft's position
    @inbounds for slot in eachindex(model.helper_indices)                       # loop over every configured helper
        helper_idx = model.helper_indices[slot]                                 # this helper's global spacecraft index
        flags[slot] = norm(pos[helper_idx] - target_pos) <= model.range_m       # true if within laser range
    end
    # Step 3: return the updated flags.
    return flags
end

# Extracts position and velocity arrays for all spacecraft from the ODE state.
# Input: ODE state u with u.sc array of per-spacecraft state views.
# Output: two Vector{SVector{3,Float64}} — positions and velocities for all N spacecraft.
function _state_vectors(u)::Tuple{Vector{SVector{3, Float64}}, Vector{SVector{3, Float64}}}
    # Step 1: allocate position and velocity vectors for all spacecraft.
    sc_state = u.sc
    pos = Vector{SVector{3, Float64}}(undef, length(sc_state))
    vel = Vector{SVector{3, Float64}}(undef, length(sc_state))
    # Step 2: copy each spacecraft's position and velocity from the ODE state.
    @inbounds for idx in eachindex(sc_state)
        sc = sc_state[idx]
        pos[idx] = SVector{3, Float64}(sc[1], sc[2], sc[3])
        vel[idx] = SVector{3, Float64}(sc[4], sc[5], sc[6])
    end
    # Step 3: return the extracted state vectors.
    return pos, vel
end

# Finds laser models matching the template in the integrator's effector list and runs the scheduling policy on each.
# Input: template model (used to match by identity or target/helper indices), integrator.
# Output: nothing (mutates matching models' active_helper_idx).
function choose_active_helper!(template::OpenCavityLaserLinkModel_struct, integrator)::Nothing
    # Step 1: extract positions and velocities once, shared by all matching effectors.
    pos, vel = _state_vectors(integrator.u)
    # Step 2: inspect every dynamic effector in the integrator.
    for model in integrator.p.args.dynamics_model.dynamic_effectors                               # loop over all effectors
        model isa OpenCavityLaserLinkModel_struct || continue                                     # skip non-laser effectors
        (model === template ||                                                                    # same instance as the template
         (model.target_idx == template.target_idx && model.helper_indices == template.helper_indices)) || continue  # or same target/helpers

        # Step 3: calculate current in-range helpers.
        in_range = _in_range_flags!(Bool[], model, pos)                                      # true/false per helper: in laser range or not

        # Step 4: apply the configured scheduling policy.
        if model.schedule === :naive_next_entering
            if model.active_helper_idx > 0                                                   # a helper is currently firing
                slot = findfirst(==(model.active_helper_idx), model.helper_indices)          # find that helper's slot in the list
                if slot === nothing || !in_range[slot]                                       # helper removed, or now out of range
                    model.active_helper_idx = 0                                              # deactivate the current link
                    helper_idx = _closest_helper(model, pos, in_range; entering_only=true)   # look for a newly-entering helper
                    helper_idx > 0 && (model.active_helper_idx = helper_idx)                 # activate it if one was found
                end
            else
                helper_idx = any(model.previous_in_range) ?                                  # were any helpers already in range last step?
                    _closest_helper(model, pos, in_range; entering_only=true) :              # yes: pick the closest newly-entering helper
                    _closest_helper(model, pos, in_range)                                    # no: pick the closest in-range helper
                helper_idx > 0 && (model.active_helper_idx = helper_idx)                     # activate the chosen helper, if any
            end
        elseif model.schedule === :positive_along_track
            if model.active_helper_idx > 0                                                   # a helper is currently firing
                slot = findfirst(==(model.active_helper_idx), model.helper_indices)          # find that helper's slot in the list
                valid = slot !== nothing &&                                                  # still a configured helper
                    in_range[slot] &&                                                        # still within laser range
                    _along_track_projection(model, model.active_helper_idx, pos, vel) > 0.0  # still pushing prograde
                valid || (model.active_helper_idx = 0)                                       # deactivate if any check failed
            end
            if model.active_helper_idx == 0                                                  # no helper currently active
                best_idx = 0                                                                 # best candidate found so far
                best_projection = 0.0                                                        # its along-track projection score
                @inbounds for slot in eachindex(model.helper_indices)                        # scan every configured helper
                    in_range[slot] || continue                                               # skip helpers out of range
                    helper_idx = model.helper_indices[slot]                                  # this helper's global spacecraft index
                    projection = _along_track_projection(model, helper_idx, pos, vel)        # how well it pushes prograde
                    if projection > best_projection                                          # better than the current best?
                        best_projection = projection                                         # remember the new best score
                        best_idx = helper_idx                                                # remember the new best helper
                    end
                end
                best_idx > 0 && (model.active_helper_idx = best_idx)                         # activate the best helper, if any
            end
        elseif model.schedule in (:gve_sma, :gve_ecc, :gve_inc, :gve_raan, :gve_argp)        # GVE-optimal scheduling: always fire the helper that maximises the chosen element rate.
            elem = model.schedule                                                            # which orbital element to optimize
            best_idx   = 0                                                                   # best candidate found so far
            best_score = 0.0                                                                 # its GVE score
            @inbounds for slot in eachindex(model.helper_indices)                            # scan every configured helper
                in_range[slot] || continue                                                   # skip helpers out of range
                helper_idx = model.helper_indices[slot]                                      # this helper's global spacecraft index
                score = _gve_score(elem, model, helper_idx, pos, vel)                        # how much it improves the chosen element
                if score > best_score                                                        # better than the current best?
                    best_score = score                                                       # remember the new best score
                    best_idx   = helper_idx                                                  # remember the new best helper
                end
            end
            model.active_helper_idx = best_idx                                               # switches to best, or deactivates if none qualify
        end

        # Step 5: save range history and count active-link steps.
        model.previous_in_range .= in_range                                                  # remember this step's in-range flags
        model.active_helper_idx > 0 && (model.active_link_step_count += 1)                   # count this step if a link is active
    end
    # Step 6: report completion after updating matching laser models.
    return nothing
end

# Builds a DiscreteCallback that runs the link scheduler at every accepted ODE step.
# Input: model (captured by closure; used to find and update matching effectors in the integrator).
# Output: DiffEqBase.DiscreteCallback that mutates model.active_helper_idx each step.
function laser_link_scheduler_callback(model::OpenCavityLaserLinkModel_struct)
    # Step 1: define a condition that runs the scheduler at every accepted step.
    condition(u, t, integrator) = true                                                             # always trigger
    # Step 2: update matching laser models when the callback fires or initializes.
    affect!(integrator) = choose_active_helper!(model, integrator)                                  # run the scheduler
    initialize = (cb, u, t, integrator) -> choose_active_helper!(model, integrator)                # also run it at setup
    # Step 3: return the configured discrete callback.
    return DiffEqBase.DiscreteCallback(condition, affect!; initialize=initialize)
end

# ── Impulse application: apply the velocity kick and track cumulative ΔV ──────────

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
    model::OpenCavityLaserLinkModel_struct,
    tracker::LaserImpulseTracker,
    mass_kg::Float64,
)
    # Step 1: define the callback that processes each accepted integration step.
    function affect!(integrator)
        # Step 2: calculate elapsed time since the previous callback.
        dt = integrator.t - tracker.t_prev
        if dt > 0.0                                                              # time actually elapsed?
            helper_idx = model.active_helper_idx                                 # currently active helper, if any
            if helper_idx > 0                                                    # a link is active
                sc      = integrator.u.sc                                        # all spacecraft states
                tgt_pos = SVector{3, Float64}(sc[model.target_idx].pos)          # target position
                tgt_vel = SVector{3, Float64}(sc[model.target_idx].vel)          # target velocity
                hlp_pos = SVector{3, Float64}(sc[helper_idx].pos)                # helper position
                rel     = tgt_pos - hlp_pos                                      # helper-to-target vector
                rho     = norm(rel)                                              # distance between them
                if rho > 0.0 && rho <= model.range_m                             # still within laser range
                    # Step 3: calculate the active laser force and RTN acceleration.
                    force = (model.eta * model.beta * model.magnification * model.power_w / SPEED_OF_LIGHT_MPS) * rel / rho
                    C = rtn_dcm_from_inertial(tgt_pos, tgt_vel)                  # target's RTN basis
                    rhat, that, nhat = C[:, 1], C[:, 2], C[:, 3]                 # radial/along-track/normal axes
                    accel = force / mass_kg                                      # force to acceleration
                    # Step 4: RTN delta-V computing
                    tracker.dv_R += dot(accel, rhat) * dt                        # add radial ΔV
                    tracker.dv_T += dot(accel, that) * dt                        # add along-track ΔV
                    tracker.dv_N += dot(accel, nhat) * dt                        # add normal ΔV
                    # Step 5: apply kick in integrator for target and the helper
                    dv = accel * dt                                              # velocity change this step
                    integrator.u.sc[model.target_idx].vel .+= dv                 # push the target
                    integrator.u.sc[helper_idx].vel .-= dv                       # recoil on the helper
                    DiffEqBase.u_modified!(integrator, true)                     # tell the solver state changed
                end
            end
        end
        # Step 5: record callback time and cumulative delta-V history.
        tracker.t_prev = integrator.t                                            # remember this callback's time
        push!(tracker.t_hist,    integrator.t)                                   # log the time
        push!(tracker.dv_R_hist, tracker.dv_R)                                   # log cumulative radial ΔV
        push!(tracker.dv_T_hist, tracker.dv_T)                                   # log cumulative along-track ΔV
        push!(tracker.dv_N_hist, tracker.dv_N)                                   # log cumulative normal ΔV
    end
    # Step 6: return the configured impulse callback.
    return DiffEqBase.DiscreteCallback(
        (u, t, integrator) -> true,
        affect!;
        save_positions=(false, false),
    )
end

end
