module LaserLinkEffectors

using LinearAlgebra
using StaticArrays
using DiffEqBase

using ...AbstractTypes: AbstractForceTorqueModel
using ...EffectorSampling: StateSample, EnvironmentSample
using ...Constellations: constellation_struct
import ..DynamicEffectors: wrench

include(joinpath(@__DIR__, "..", "..", "..", "core", "interfaces", "reference_system.jl"))

export LaserThrusterParams, LaserCommunicationParams, LaserPowerTransferParams
export LaserLinkModel, build_LaserLinkModel, laser_link_scheduler_callback
export choose_active_links!
export LaserImpulseTracker, laser_impulse_callback

const SPEED_OF_LIGHT_MPS = 299_792_458.0
const _MU_EARTH_GVE = 3.986004418e14  # Earth gravitational parameter [m³/s²]

# ── Model construction ──────────────────────────────────────────────────────
# Params type picked by laser_type (e.g. LaserThrusterParams for :thruster).
abstract type AbstractLaserLinkParams end

# Thruster physics, used when laser_type === :thruster.
struct LaserThrusterParams <: AbstractLaserLinkParams
    power_w::Float64
    magnification::Float64
    beta::Float64
    eta::Float64
end

# Comms physics, used when laser_type === :communication.
struct LaserCommunicationParams <: AbstractLaserLinkParams
    power_w::Float64
end

# Power-transfer physics, used when laser_type === :power_transfer.
struct LaserPowerTransferParams <: AbstractLaserLinkParams
    power_w::Float64
end

# One laser link between two satellites. Shared scheduling state lives on constellation_struct, not here.
mutable struct LaserLinkModel <: AbstractForceTorqueModel
    sat_i::Int
    sat_j::Int
    range_m::Float64
    schedule::Symbol
    laser_type::Symbol                              # :communication, :power_transfer, :thruster
    params::Union{Nothing, AbstractLaserLinkParams} # picked by laser_type
end

function build_LaserLinkModel(;
    sat_i::Integer,
    sat_j::Integer,
    range_m::Real=200e3,
    schedule::Symbol=:naive_next_entering,
    laser_type::Symbol=:thruster,
    power_w::Real=10_000.0,
    magnification::Real=100.0,
    beta::Real=1.0,
    eta::Real=2.0,
)
    # pick params by laser_type, then build the link.
    params = if laser_type === :thruster
        LaserThrusterParams(power_w, magnification, beta, eta)
    elseif laser_type === :communication
        LaserCommunicationParams(power_w)
    elseif laser_type === :power_transfer
        LaserPowerTransferParams(power_w)
    else
        nothing
    end
    return LaserLinkModel(Int(sat_i), Int(sat_j), Float64(range_m), schedule, laser_type, params)
end

# Laser force is applied via callback, not the ODE RHS, so wrench is always zero.
@inline function wrench(
    ::LaserLinkModel,
    ::StateSample,
    ::EnvironmentSample,
    ::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end


# ── Helper-scoring (one scoring function per scheduling policy) ────────────────
# Projection of the emitter->receiver laser direction onto the receiver's along-track axis;
# positive means the laser pushes the receiver prograde.
@inline function _along_track_projection(
    receiver::Int,
    emitter::Int,
    pos::AbstractVector{SVector{3, Float64}},
    vel::AbstractVector{SVector{3, Float64}},
)::Float64
    receiver_pos = pos[receiver]
    rel = receiver_pos - pos[emitter]       # emitter → receiver
    rho = norm(rel)
    that = rtn_dcm_from_inertial(receiver_pos, vel[receiver])[:, 2]
    return dot(rel / rho, that)
end

# Scores the instantaneous rate of change of orbital element `elem` if the laser fires emitter -> receiver.
function _gve_score(
    elem::Symbol,
    receiver::Int,
    emitter::Int,
    pos::AbstractVector{SVector{3, Float64}},
    vel::AbstractVector{SVector{3, Float64}},
)::Float64
    tgt_pos = pos[receiver]
    tgt_vel = vel[receiver]

    # Unit force direction on receiver (emitter → receiver)
    rel = tgt_pos - pos[emitter]
    rho = norm(rel)
    f̂ = rel / rho

    # RTN decomposition of force direction
    C = rtn_dcm_from_inertial(tgt_pos, tgt_vel)
    rhat, that, nhat = C[:, 1], C[:, 2], C[:, 3]
    aR = dot(f̂, rhat)
    aT = dot(f̂, that)
    aN = dot(f̂, nhat)

    # Step 2: compute the target's orbital elements and anomaly.
    r  = norm(tgt_pos)
    v2 = dot(tgt_vel, tgt_vel)
    a  = -_MU_EARTH_GVE / (v2 - 2.0 * _MU_EARTH_GVE / r)   # vis-viva

    h_vec  = cross(tgt_pos, tgt_vel)
    h_sq   = dot(h_vec, h_vec)
    h_norm = sqrt(h_sq)
    p_slr  = h_sq / _MU_EARTH_GVE                            # semi-latus rectum

    e_vec = cross(tgt_vel, h_vec) / _MU_EARTH_GVE - tgt_pos / r
    e     = norm(e_vec)
    e_sq  = clamp(e * e, 0.0, 1.0 - 1e-12)
    sqrt_1me2 = sqrt(1.0 - e_sq)
    n_mean    = sqrt(_MU_EARTH_GVE / (a * a * a))

    # True anomaly ν
    ν = acos(clamp(dot(e_vec / max(e, 1e-12), tgt_pos / r), -1.0, 1.0))
    dot(tgt_pos, tgt_vel) < 0.0 && (ν = 2π - ν)

    # Inclination/argument-of-latitude terms needed by gve_inc, gve_raan, gve_argp.
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
        # Ω̇ = r sin(u) / (na²√(1-e²) sin i) · aN   [singular at i = 0 — kept: our default scenario is equatorial]
        abs(sin_i) < 1e-6 && return 0.0
        return (r * sin(u) / (denom * sin_i)) * aN

    else  # :gve_argp
        # ω̇ = √(1-e²)/(nae) [-cosν · aR + (1+r/p) sinν · aT]
        #       - r sin(u) cos(i) / (na²√(1-e²) sin i) · aN   [singular for e=0 or i=0 — kept: our defaults]
        abs(e) < 1e-6    && return 0.0
        abs(sin_i) < 1e-6 && return 0.0
        term_RT = (sqrt_1me2 / (n_mean * a * e)) *
                  (-cos(ν) * aR + (1.0 + r / p_slr) * sin(ν) * aT)
        term_N  = -(r * sin(u) * cos_i / (denom * sin_i)) * aN
        return term_RT + term_N
    end
end

# ── Scheduling: choose which links fire each step, and wire it into the ODE solver ──
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

# Greedily assigns satellite-disjoint links from `candidates` (highest priority first) on top of `kept`.
# Input: kept (already-selected links, endpoints pre-claimed), candidates (priority-ordered pool), used (claimed satellites).
# Output: Vector{Tuple{Int,Int}} — kept plus every candidate whose endpoints were still free when its turn came.
function _assign_disjoint_links!(
    kept::Vector{Tuple{Int, Int}},
    candidates::AbstractVector{Tuple{Int, Int}},
    used::Set{Int},
)::Vector{Tuple{Int, Int}}
    # Step 1: walk the priority-ordered candidates, claiming any whose endpoints are still free.
    for (r, e) in candidates
        (r in used || e in used) && continue                                    # one or both satellites already committed this tick
        push!(kept, (r, e))                                                     # claim the link
        push!(used, r)
        push!(used, e)                                                          # mark both endpoints as committed
    end
    # Step 2: return the combined selection.
    return kept
end

# Runs the constellation's scheduling policy and updates its shared link state.
# Each satellite may appear in at most one entry of constellation.active_links at a time (enforced by the
# matching below). Assumes every registered LaserLinkModel shares the same schedule.
# Input: constellation (owns possible_links/active_links/previous_in_range_links), integrator.
# Output: nothing (mutates constellation.active_links and constellation.previous_in_range_links).
function choose_active_links!(constellation::constellation_struct, integrator)::Nothing
    # Step 1: extract positions and velocities once, shared by every registered link.
    pos, vel = _state_vectors(integrator.u)                                       # positions/velocities for every spacecraft
    models = filter(m -> m isa LaserLinkModel, integrator.p.args.dynamics_model.dynamic_effectors)  # all registered links
    schedule = models[1].schedule                                                # shared by every registered link

    # Step 2: calculate the currently in-range subset of every registered link.
    in_range = [
        (m.sat_i, m.sat_j) for m in models
        if norm(pos[m.sat_j] - pos[m.sat_i]) <= m.range_m                        # true if within this link's own range
    ]
    in_range_set = Set(in_range)                                                 # fast membership test for the matching below

    used = Set{Int}()                                                            # satellites already claimed this tick
    selected = Tuple{Int, Int}[]                                                 # this tick's constellation-wide active-link selection

    # Step 3: apply the configured scheduling policy as a satellite-disjoint matching over in_range.
    if schedule === :naive_next_entering
        # Sticky: keep any currently active link that is still in range.
        kept = [link for link in constellation.active_links if link in in_range_set]  # links carried over from last tick
        for (r, e) in kept
            push!(used, r); push!(used, e)                                       # claim both endpoints of every kept link
        end
        # Rank unclaimed in-range links: newly-entering ones outrank already-in-range ones; ties broken by distance.
        entering  = Tuple{Int, Int}[]                                            # in-range links not seen last tick
        remaining = Tuple{Int, Int}[]                                            # in-range links already seen last tick
        for link in in_range
            link in kept && continue                                             # already carried over above
            if link in constellation.previous_in_range_links
                push!(remaining, link)                                           # was already in range last tick
            else
                push!(entering, link)                                            # newly entered range this tick
            end
        end
        sort!(entering;  by = link -> norm(pos[link[2]] - pos[link[1]]))         # closest entering link first
        sort!(remaining; by = link -> norm(pos[link[2]] - pos[link[1]]))         # closest already-in-range link first
        _assign_disjoint_links!(kept, entering, used)                            # claim entering links first
        selected = _assign_disjoint_links!(kept, remaining, used)                # then fill remaining free satellites
    elseif schedule === :positive_along_track
        # Sticky: keep any currently active link that is still in range and still pushing prograde.
        kept = [
            link for link in constellation.active_links
            if link in in_range_set && _along_track_projection(link[1], link[2], pos, vel) > 0.0
        ]                                                                        # links carried over from last tick
        for (r, e) in kept
            push!(used, r); push!(used, e)                                       # claim both endpoints of every kept link
        end
        scored = [
            (link, _along_track_projection(link[1], link[2], pos, vel))
            for link in in_range if !(link in kept)
        ]                                                                                 # (link, projection) for unclaimed in-range links
        filter!(x -> x[2] > 0.0, scored)                                                  # only prograde-pushing candidates qualify
        sort!(scored; by = x -> x[2], rev = true)                                         # best projection first
        selected = _assign_disjoint_links!(kept, first.(scored), used)                    # claim links in descending projection order
    elseif schedule in (:gve_sma, :gve_ecc, :gve_inc, :gve_raan, :gve_argp)               # GVE-optimal scheduling: not sticky, fully re-matched every tick.
        scored = [
            (link, _gve_score(schedule, link[1], link[2], pos, vel))
            for link in in_range
        ]                                                                                 # (link, GVE score) for every in-range link
        filter!(x -> x[2] > 0.0, scored)                                                  # only links that actually help qualify
        sort!(scored; by = x -> x[2], rev = true)                                         # best score first
        selected = _assign_disjoint_links!(Tuple{Int, Int}[], first.(scored), used)       # claim links in descending score order
    end

    # Step 4: commit the constellation-wide selection and range history for next tick's entering-detection.
    constellation.active_links = selected                                        # publish this tick's active links
    constellation.previous_in_range_links = in_range                             # remember this tick's in-range links
    return nothing
end

# Builds a DiscreteCallback that runs the link scheduler at every accepted ODE step.
# Input: constellation (captured by closure; mutated in place by choose_active_links!).
# Output: DiffEqBase.DiscreteCallback that mutates constellation.active_links each step.
function laser_link_scheduler_callback(constellation::constellation_struct)
    # Step 1: define a condition that runs the scheduler at every accepted step.
    condition(u, t, integrator) = true                                                             # always trigger
    # Step 2: update the constellation's link state when the callback fires or initializes.
    affect!(integrator) = choose_active_links!(constellation, integrator)                          # run the scheduler
    initialize = (cb, u, t, integrator) -> choose_active_links!(constellation, integrator)          # also run it at setup
    # Step 3: return the configured discrete callback.
    return DiffEqBase.DiscreteCallback(condition, affect!; initialize=initialize)
end

# ── Impulse application: apply the velocity kick and track cumulative ΔV ──────────

# Running accumulator for laser ΔV in RTN components; also stores the full time-series history.
Base.@kwdef mutable struct LaserImpulseTracker
    t_prev::Float64            = 0.0                                             # time of the previous callback invocation
    dv_R::Dict{Tuple{Int, Int}, Float64}              = Dict{Tuple{Int, Int}, Float64}()             # cumulative radial ΔV per link
    dv_T::Dict{Tuple{Int, Int}, Float64}              = Dict{Tuple{Int, Int}, Float64}()             # cumulative along-track ΔV per link
    dv_N::Dict{Tuple{Int, Int}, Float64}              = Dict{Tuple{Int, Int}, Float64}()             # cumulative normal ΔV per link
    active_link_steps::Int     = 0                                               # count of steps with at least one active link
    t_hist::Vector{Float64}    = Float64[]                                       # time at each recorded step
    dv_R_hist::Dict{Tuple{Int, Int}, Vector{Float64}} = Dict{Tuple{Int, Int}, Vector{Float64}}()     # radial ΔV time series per link
    dv_T_hist::Dict{Tuple{Int, Int}, Vector{Float64}} = Dict{Tuple{Int, Int}, Vector{Float64}}()     # along-track ΔV time series per link
    dv_N_hist::Dict{Tuple{Int, Int}, Vector{Float64}} = Dict{Tuple{Int, Int}, Vector{Float64}}()     # normal ΔV time series per link
end

# Records `value` into link's history vector, backfilling with zeros if this is the link's first appearance.
function _push_link_hist!(
    hist::Dict{Tuple{Int, Int}, Vector{Float64}},
    link::Tuple{Int, Int},
    value::Float64,
    n_steps::Int,
)::Nothing
    v = get!(() -> zeros(n_steps - 1), hist, link)                               # existing history, or zero-padded if new
    push!(v, value)                                                              # append this step's cumulative value
    return nothing
end

# Finds the registered thruster-type LaserLinkModel matching `link`, if any.
# Input: dynamic_effectors (from the integrator), link (a (sat_i, sat_j) pair).
# Output: the matching LaserLinkModel, or nothing (e.g. link is :communication / :power_transfer, not yet modeled).
function _find_thruster_model(dynamic_effectors, link::Tuple{Int, Int})::Union{Nothing, LaserLinkModel}
    for model in dynamic_effectors
        model isa LaserLinkModel || continue                                     # skip non-laser effectors
        model.laser_type === :thruster || continue                               # skip comms/power-transfer links
        (model.sat_i, model.sat_j) == link && return model                       # matching link found
    end
    return nothing                                                               # no thruster-type model registered for this link
end

# Builds a DiscreteCallback that applies a discrete velocity kick and accumulates RTN ΔV each step.
# Input: constellation (for active_links), tracker (accumulates ΔV history, keyed per link), mass_kg of target.
# Output: DiffEqBase.DiscreteCallback that mutates integrator.u velocities and tracker fields.
function laser_impulse_callback(
    constellation::constellation_struct,
    tracker::LaserImpulseTracker,
    mass_kg::Float64,
)
    # Step 1: define the callback that processes each accepted integration step.
    function affect!(integrator)
        # Step 2: calculate elapsed time since the previous callback.
        dt = integrator.t - tracker.t_prev
        if dt > 0.0                                                              # time actually elapsed?
            isempty(constellation.active_links) || (tracker.active_link_steps += 1)  # count this step if any link is active
            for link in constellation.active_links                              # apply the kick for every currently active link
                model = _find_thruster_model(integrator.p.args.dynamics_model.dynamic_effectors, link)
                model === nothing && continue                                    # not a thruster link — no mechanical kick (yet)
                (receiver, emitter) = link
                sc      = integrator.u.sc                                        # all spacecraft states
                tgt_pos = SVector{3, Float64}(sc[receiver].pos)                  # receiver position
                tgt_vel = SVector{3, Float64}(sc[receiver].vel)                  # receiver velocity
                hlp_pos = SVector{3, Float64}(sc[emitter].pos)                   # emitter position
                rel     = tgt_pos - hlp_pos                                      # emitter-to-receiver vector
                rho     = norm(rel)                                              # distance between them
                if rho <= model.range_m                                          # still within laser range
                    # Step 3: calculate the active laser force and RTN acceleration.
                    tp = model.params::LaserThrusterParams
                    force = (tp.eta * tp.beta * tp.magnification * tp.power_w / SPEED_OF_LIGHT_MPS) * rel / rho
                    C = rtn_dcm_from_inertial(tgt_pos, tgt_vel)                  # receiver's RTN basis
                    rhat, that, nhat = C[:, 1], C[:, 2], C[:, 3]                 # radial/along-track/normal axes
                    accel = force / mass_kg                                      # force to acceleration
                    # Step 4: RTN delta-V computing, accumulated per link.
                    tracker.dv_R[link] = get(tracker.dv_R, link, 0.0) + dot(accel, rhat) * dt   # add radial ΔV
                    tracker.dv_T[link] = get(tracker.dv_T, link, 0.0) + dot(accel, that) * dt   # add along-track ΔV
                    tracker.dv_N[link] = get(tracker.dv_N, link, 0.0) + dot(accel, nhat) * dt   # add normal ΔV
                    # Step 5: apply kick in integrator for receiver and the emitter
                    dv = accel * dt                                              # velocity change this step
                    integrator.u.sc[receiver].vel .+= dv                         # push the receiver
                    integrator.u.sc[emitter].vel .-= dv                          # recoil on the emitter
                    DiffEqBase.u_modified!(integrator, true)                     # tell the solver state changed
                end
            end
        end
        # Step 5: record callback time and cumulative delta-V history for every link seen so far.
        tracker.t_prev = integrator.t                                            # remember this callback's time
        push!(tracker.t_hist, integrator.t)                                      # log the time
        n_steps = length(tracker.t_hist)
        for link in keys(tracker.dv_R)                                          # every link that has ever fired
            _push_link_hist!(tracker.dv_R_hist, link, tracker.dv_R[link], n_steps)   # log cumulative radial ΔV
            _push_link_hist!(tracker.dv_T_hist, link, tracker.dv_T[link], n_steps)   # log cumulative along-track ΔV
            _push_link_hist!(tracker.dv_N_hist, link, tracker.dv_N[link], n_steps)   # log cumulative normal ΔV
        end
    end
    # Step 6: return the configured impulse callback.
    return DiffEqBase.DiscreteCallback(
        (u, t, integrator) -> true,
        affect!;
        save_positions=(false, false),
    )
end

end
