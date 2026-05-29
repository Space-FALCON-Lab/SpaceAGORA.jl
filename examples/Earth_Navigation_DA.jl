include(joinpath(@__DIR__, "common.jl"))
using .SimulationModel
using Dates
using SPICE
using StaticArrays
using LinearAlgebra
using Random
using CSV
using DataFrames
const RuntimeServices = SpaceAGORA.RuntimeServices
const EARTH_HARMONICS_FILE = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
const MISSION_TIME_SEC = 500.0
const THRESHOLD_DISTANCE_KM = 500.0 
const COMMUNICATION_RANGE_KM = 300.0
const SIGMA_THETA_RAD = 0.0001 
const NAVIGATION_RATE_SEC = 5.0
const NAVIGATION_DT_TOL_SEC = 0.25
const EKF_CONSENSUS_ITERS = 5
const EKF_PROCESS_Q_DIAG = SVector{6, Float64}(5e-1, 5e-1, 5e-1, 5e-3, 5e-3, 5e-3)
const ORPHAN_ATTACH_MAX_DTHETA_RAD = 0.05   
const IOD_NEIGHBOR_MISS_DISTANCE_MAX_M = 200.0
const CONSENSUS_MATCH_MAHAL_IOD_MAX_D2 = 2.0
# const CONSENSUS_MATCH_MAHAL_FILTER_MAX_D2 = 100.0
const CONSENSUS_MATCH_MAHAL_FILTER_MAX_D2 = 27.86
const CONSENSUS_MATCH_COV_REG_EPS = 1e-6
# const MEAS_ASSOC_MAHAL_MAX_D2 = 300.0
const MEAS_ASSOC_MAHAL_MAX_D2 = 21.11
const MEAS_ASSOC_DISAMBIG_RATIO_MIN = 1.5
const MEAS_ASSOC_A_INCLUDE_SELF_SCORE = true
const CONSENSUS_GROUP_FALLBACK_ENABLE = false
const DEKF_INIT_ALLOW_SINGLETON_IF_NO_GROUP = true
const DEKF_CONSENSUS_UPDATE_MODE = :kcf # :information, :icf, or :kcf
const DEKF_KCF_EPSILON = 0.05
const TRACK_CLOSE_AFTER_MISSED_MEASUREMENTS = 5
const ENABLE_TRACK_CLOSE_AFTER_SOLO_TARGET_MEASURE_STREAK = true
const SOLO_TARGET_MEASURE_STREAK_CLOSE_STEPS = 10
const ENABLE_OBSERVER_OD_PERTURBATION = true
const OBSERVER_OD_POS_STD_M = 5.0 # usually for LEO 1-10 cm with GNSS, 1-5m on-board KF 
const OBSERVER_OD_VEL_STD_MPS = 0.05 # usually 0.1-1 mm/s with GNSS,  2-5 cm/s on board
# const OBSERVER_OD_RNG_SEED = 12345
const NAN_LOS = SVector{3, Float64}(NaN, NaN, NaN)
const NAN_STATE6 = SVector{6, Float64}(ntuple(_ -> NaN, 6))
@inline _nan_cov6() = fill(NaN, 6, 6)
(DEKF_CONSENSUS_UPDATE_MODE == :information || DEKF_CONSENSUS_UPDATE_MODE == :icf || DEKF_CONSENSUS_UPDATE_MODE == :kcf) ||
    error("Invalid DEKF_CONSENSUS_UPDATE_MODE=$(DEKF_CONSENSUS_UPDATE_MODE). Use :information, :icf, or :kcf.")

Base.@kwdef mutable struct LOSMeasurement
    t::Float64
    observer::Int
    target::Int
    range_m::Float64
    los_unit::SVector{3, Float64}
end

Base.@kwdef mutable struct OpticalLOSSensorModel
    observer_idxs::Vector{Int}
    target_idxs::Vector{Int}
    detection_range_m::Float64
    sigma_theta_rad::Float64
    counts::Matrix{Int}
    latest::Dict{Tuple{Int, Int}, LOSMeasurement}
    previous::Dict{Tuple{Int, Int}, LOSMeasurement}
    measurements_now::Vector{Vector{LOSMeasurement}}
end

function OpticalLOSSensorModel(
    observer_idxs::Vector{Int},
    target_idxs::Vector{Int},
    detection_range_m::Float64,
    sigma_theta_rad::Float64,
    num_sats::Int
)
    return OpticalLOSSensorModel(
        observer_idxs=observer_idxs,
        target_idxs=target_idxs,
        detection_range_m=detection_range_m,
        sigma_theta_rad=sigma_theta_rad,
        counts=zeros(Int, num_sats, num_sats),
        latest=Dict{Tuple{Int, Int}, LOSMeasurement}(),
        previous=Dict{Tuple{Int, Int}, LOSMeasurement}(),
        measurements_now=[LOSMeasurement[] for _ in 1:num_sats]
    )
end

@inline function _perturb_los(
    los_true::SVector{3, Float64},
    sigma_theta_rad::Float64,
    rng::AbstractRNG=Random.default_rng()
)::SVector{3, Float64}
    sigma_theta_rad <= 0.0 && return los_true
    theta = randn(rng) * sigma_theta_rad

    v = randn(rng, 3)
    v = v - dot(v, los_true) * los_true
    v = v / norm(v) 

    K = [ 0.0 -v[3] v[2];
        v[3] 0.0 -v[1];
       -v[2] v[1] 0.0 ]

    R = I + sin(theta) * K + (1 - cos(theta)) * (v * v')

    return R * los_true / norm(R * los_true)
end

function SimulationModel.calcNavigationEffect!(
    model::OpticalLOSSensorModel,
    u,
    p,
    t::Float64,
    sat_idx::Int
)
    sat_idx in model.observer_idxs || return nothing # Only compute measurements for observer satellites.
    empty!(model.measurements_now[sat_idx]) # Clear measurements from previous nav tick for this observer.
    observer_pos = SVector{3, Float64}(u.sc[sat_idx].pos)

    @inbounds for target_idx in model.target_idxs # Iterate over all target satellites, Skip if observer and target are the same satellite.
        target_idx == sat_idx && continue
        rel_pos = SVector{3, Float64}(u.sc[target_idx].pos) - observer_pos
        dist = norm(rel_pos)
        if dist <= model.detection_range_m
            los_true = rel_pos / dist
            los_meas = _perturb_los(los_true, model.sigma_theta_rad)
            model.counts[sat_idx, target_idx] += 1
            key = (sat_idx, target_idx) 
            # Before updating the latest measurement, move the current latest to previous if it exists. 
            if haskey(model.latest, key) 
                model.previous[key] = model.latest[key] 
            end
            meas = LOSMeasurement(
                t=t,
                observer=sat_idx,
                target=target_idx,
                range_m=dist,
                los_unit=los_meas
            )
            model.latest[key] = meas
            push!(model.measurements_now[sat_idx], meas)
        end
    end

    return nothing
end

Base.@kwdef mutable struct InterAgentCommunicationModel
    participant_idxs::Vector{Int}
    comm_range_m::Float64
    adjacency::BitMatrix
    neighbors::Vector{Vector{Int}}
end

function InterAgentCommunicationModel(
    participant_idxs::Vector{Int},
    comm_range_m::Float64,
    num_sats::Int
)
    return InterAgentCommunicationModel(
        participant_idxs=participant_idxs,
        comm_range_m=comm_range_m,
        adjacency=falses(num_sats, num_sats),
        neighbors=[Int[] for _ in 1:num_sats]
    )
end

@inline function _is_consecutive_measure_pair(t_last::Float64, t_prev::Float64)::Bool
    dt = t_last - t_prev
    return (dt > 0.0) && (abs(dt - NAVIGATION_RATE_SEC) <= NAVIGATION_DT_TOL_SEC)
end

# This function updates the adjacency matrix and neighbor lists based on the current positions of the agents. 
function _update_neighbors!(comms::InterAgentCommunicationModel, u, p)
    fill!(comms.adjacency, false) # Reset adjacency matrix
    @inbounds for idx in eachindex(comms.neighbors)
        empty!(comms.neighbors[idx])
    end
    participants = comms.participant_idxs
    n = length(participants)
    @inbounds for a in 1:n
        i = participants[a]
        p.is_active[i] || continue
        pos_i = SVector{3, Float64}(u.sc[i].pos)
        for b in (a + 1):n
            j = participants[b]
            p.is_active[j] || continue
            pos_j = SVector{3, Float64}(u.sc[j].pos)
            if norm(pos_i - pos_j) <= comms.comm_range_m
                comms.adjacency[i, j] = true
                comms.adjacency[j, i] = true
                push!(comms.neighbors[i], j)
                push!(comms.neighbors[j], i)
            end
        end
    end
    return nothing
end

@inline function _safe_unit(v::SVector{3, Float64})::SVector{3, Float64}
    n = norm(v)
    return n > 1e-12 ? (v / n) : SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline _is_finite_state(x::SVector{6, Float64}) = all(isfinite, x)
@inline _is_finite_cov(P::Matrix{Float64}) = all(isfinite, P)

@inline function _skew(v::SVector{3, Float64})::SMatrix{3, 3, Float64, 9}
    return SMatrix{3, 3, Float64, 9}(
         0.0,   -v[3],  v[2],
         v[3],   0.0,  -v[1],
        -v[2],   v[1],  0.0
    )
end

function _build_agent_equations(
    r::SVector{3, Float64},
    v::SVector{3, Float64},
    los::SVector{3, Float64},
    los_rate::SVector{3, Float64}
)
    r_norm = norm(r)
    if r_norm <= 0.0
        return nothing
    end
    R_unit = -r / r_norm
    n1 = cross(los, R_unit)
    n2 = cross(n1, los)
    dR_unit = -v / r_norm + r * (dot(r, v) / r_norm^3)
    dn1 = cross(los_rate, R_unit) + cross(los, dR_unit)
    dn2 = cross(dn1, los) + cross(n1, los_rate)
    A = [n1'; n2']
    b = [dot(n1, r); dot(n2, r)]
    dA = [dn1'; dn2']
    db = [dot(dn1, r) + dot(n1, v); dot(dn2, r) + dot(n2, v)]
    return (A=A, b=b, dA=dA, db=db, n1=n1, R_unit=R_unit, dR_unit=dR_unit, dn1=dn1)
end

function _compute_state_covariance(
    nodes,
    x_est::SVector{3, Float64},
    x_dot_est::SVector{3, Float64},
    sigma_theta::Float64,
    dt::Float64
)::Union{Nothing, Matrix{Float64}}
    m = length(nodes)
    m == 0 && return nothing
    H = Matrix{Float64}(undef, 4 * m, 6)
    Cov_e = zeros(4 * m, 4 * m)
    for (k, nd) in enumerate(nodes)
        r = nd.r
        v = nd.v
        l = nd.los
        l_dot = nd.los_rate
        eq = _build_agent_equations(r, v, l, l_dot)
        eq === nothing && return nothing
        row = 4 * (k - 1)
        H[row + 1:row + 2, 1:3] .= eq.A
        H[row + 1:row + 2, 4:6] .= 0.0
        H[row + 3:row + 4, 1:3] .= eq.dA
        H[row + 3:row + 4, 4:6] .= eq.A

        rho = r - x_est
        eta = v - x_dot_est
        R = eq.R_unit
        Rdot = eq.dR_unit
        n1 = eq.n1
        n1dot = eq.dn1
        SR = _skew(R)
        SRdot = _skew(Rdot)
        Sl = _skew(l)
        Sldot = _skew(l_dot)
        Sn1 = _skew(n1)
        Sn1dot = _skew(n1dot)

        row1_l = -(rho' * SR)
        row2_l = rho' * (-Sl * SR + Sn1)
        row3_l = -(rho' * SRdot + eta' * SR)
        row3_ldot = -(rho' * SR)
        row4_l = rho' * (Sl * SRdot + Sn1dot + Sldot * SR) + eta' * (-Sl * SR + Sn1)
        row4_ldot = rho' * (Sl * SR + Sn1)
        J = zeros(4, 6)
        J[1, 1:3] .= vec(row1_l)
        J[2, 1:3] .= vec(row2_l)
        J[3, 1:3] .= vec(row3_l)
        J[3, 4:6] .= vec(row3_ldot)
        J[4, 1:3] .= vec(row4_l)
        J[4, 4:6] .= vec(row4_ldot)

        P = Matrix(I, 3, 3) - l * l'
        cov_l = sigma_theta^2 * P
        cov_l_dot = (2.0 * sigma_theta^2 / dt^2) * P
        cov_cross = (sigma_theta^2 / dt) * P
        Sigma_i = [cov_l cov_cross; cov_cross cov_l_dot]
        Cov_i = J * Sigma_i * J'
        Cov_e[row + 1:row + 4, row + 1:row + 4] .= Cov_i
    end
    HtH = H' * H
    HtH_inv = pinv(HtH)
    return HtH_inv * H' * Cov_e * H * HtH_inv
end

# Local track memory for each observer and slot, storing the latest two LOS measures and associated context for seeding and fusion.
Base.@kwdef mutable struct LocalTrack
    slot::Int
    last_meas::Union{Nothing, LOSMeasurement}
    prev_meas::Union{Nothing, LOSMeasurement}
    plot_target_id::Int
    status::Symbol
    seed_ready::Bool
    has_measure_now::Bool
    consecutive_missed::Int
    last_update_t::Float64
    state_estimate_now::SVector{6, Float64}
    covariance_estimate_now::Matrix{Float64}
    observer_pos_now::SVector{3, Float64}
    observer_pos_prev::SVector{3, Float64}
end


Base.@kwdef mutable struct ObserverNavigationModel
    sensor::OpticalLOSSensorModel
    comms::InterAgentCommunicationModel
    observer_idxs::Vector{Int}
    track_slots::Vector{Int} 
    local_tracks::Vector{Dict{Int, LocalTrack}} 
    od_perturb_enabled::Bool
    od_pos_std_m::Float64
    od_vel_std_mps::Float64
    od_rng::MersenneTwister
    known_observer_state::Vector{SVector{6, Float64}}
    known_observer_state_t::Vector{Float64}
    conflict_count_total::Vector{Int}
    skipped_assoc_total::Vector{Int}
    disambiguation_calls_total::Vector{Int}
    disambiguation_success_total::Vector{Int}
    disambiguation_ratio_success_correct_total::Vector{Int}
    disambiguation_ratio_success_wrong_total::Vector{Int}
    disambiguation_ratio_fail_total::Vector{Int}
    meas_commit_unique_total::Vector{Int}
    meas_commit_unique_correct_total::Vector{Int}
    meas_commit_unique_wrong_total::Vector{Int}
    meas_commit_ambig_total::Vector{Int}
    meas_commit_ambig_correct_total::Vector{Int}
    meas_commit_ambig_wrong_total::Vector{Int}
    meas_commit_ambig_dropped_total::Vector{Int}
    b_conflict_events_total::Vector{Int}
    skipped_collision_B_total::Vector{Int}
    next_local_slot::Vector{Int}       # counter for assigning new local slot IDs to orphan measurements.
end

function ObserverNavigationModel(
    sensor::OpticalLOSSensorModel,
    comms::InterAgentCommunicationModel,
    observer_idxs::Vector{Int},
    track_slots::Vector{Int},
    num_sats::Int;
    od_perturb_enabled::Bool=ENABLE_OBSERVER_OD_PERTURBATION,
    od_pos_std_m::Float64=OBSERVER_OD_POS_STD_M,
    od_vel_std_mps::Float64=OBSERVER_OD_VEL_STD_MPS
)
    local_slot_seed = max(maximum(track_slots), num_sats) + 1
    return ObserverNavigationModel(
        sensor=sensor,
        comms=comms,
        observer_idxs=observer_idxs,
        track_slots=track_slots,
        local_tracks=[Dict{Int, LocalTrack}() for _ in 1:num_sats],
        od_perturb_enabled=od_perturb_enabled,
        od_pos_std_m=od_pos_std_m,
        od_vel_std_mps=od_vel_std_mps,
        od_rng=MersenneTwister(),
        known_observer_state=fill(NAN_STATE6, num_sats),
        known_observer_state_t=fill(NaN, num_sats),
        conflict_count_total=zeros(Int, num_sats),
        skipped_assoc_total=zeros(Int, num_sats),
        disambiguation_calls_total=zeros(Int, num_sats),
        disambiguation_success_total=zeros(Int, num_sats),
        disambiguation_ratio_success_correct_total=zeros(Int, num_sats),
        disambiguation_ratio_success_wrong_total=zeros(Int, num_sats),
        disambiguation_ratio_fail_total=zeros(Int, num_sats),
        meas_commit_unique_total=zeros(Int, num_sats),
        meas_commit_unique_correct_total=zeros(Int, num_sats),
        meas_commit_unique_wrong_total=zeros(Int, num_sats),
        meas_commit_ambig_total=zeros(Int, num_sats),
        meas_commit_ambig_correct_total=zeros(Int, num_sats),
        meas_commit_ambig_wrong_total=zeros(Int, num_sats),
        meas_commit_ambig_dropped_total=zeros(Int, num_sats),
        b_conflict_events_total=zeros(Int, num_sats),
        skipped_collision_B_total=zeros(Int, num_sats),
        next_local_slot=fill(local_slot_seed, num_sats)
    )
end

@inline function _true_observer_state(u, observer_id::Int)::SVector{6, Float64}
    pos = SVector{3, Float64}(u.sc[observer_id].pos)
    vel = SVector{3, Float64}(u.sc[observer_id].vel)
    return SVector{6, Float64}(pos[1], pos[2], pos[3], vel[1], vel[2], vel[3])
end

@inline function _known_observer_state!(
    model::ObserverNavigationModel,
    observer_id::Int,
    u,
    t::Float64
)::SVector{6, Float64}
    true_state = _true_observer_state(u, observer_id)
    model.od_perturb_enabled || return true_state

    last_t = model.known_observer_state_t[observer_id]
    if !isfinite(last_t) || !isapprox(last_t, t; atol=1e-9, rtol=0.0)
        pos = SVector{3, Float64}(true_state[1], true_state[2], true_state[3])
        vel = SVector{3, Float64}(true_state[4], true_state[5], true_state[6])
        pos_noise = model.od_pos_std_m .* SVector{3, Float64}(
            randn(model.od_rng), randn(model.od_rng), randn(model.od_rng)
        )
        vel_noise = model.od_vel_std_mps .* SVector{3, Float64}(
            randn(model.od_rng), randn(model.od_rng), randn(model.od_rng)
        )
        pos_k = pos + pos_noise
        vel_k = vel + vel_noise
        model.known_observer_state[observer_id] = SVector{6, Float64}(
            pos_k[1], pos_k[2], pos_k[3], vel_k[1], vel_k[2], vel_k[3]
        )
        model.known_observer_state_t[observer_id] = t
    end

    return model.known_observer_state[observer_id]
end

@inline function _known_observer_pos!(
    model::ObserverNavigationModel,
    observer_id::Int,
    u,
    t::Float64
)::SVector{3, Float64}
    x = _known_observer_state!(model, observer_id, u, t)
    return SVector{3, Float64}(x[1], x[2], x[3])
end

@inline function _known_observer_vel!(
    model::ObserverNavigationModel,
    observer_id::Int,
    u,
    t::Float64
)::SVector{3, Float64}
    x = _known_observer_state!(model, observer_id, u, t)
    return SVector{3, Float64}(x[4], x[5], x[6])
end


# Create (if needed) and return the local track object for this observer/slot pair.
@inline function _ensure_local_track!(
    model::ObserverNavigationModel,
    observer_id::Int,
    slot_id::Int,
    t::Float64
)::LocalTrack
    tracks = model.local_tracks[observer_id]
    if !haskey(tracks, slot_id)
        tracks[slot_id] = LocalTrack(
            slot=slot_id,
            last_meas=nothing,
            prev_meas=nothing,
            plot_target_id=0,
            # Every new local track starts as tentative.
            # Status upgrades happen later through local association / IOD / filter init.
            status=:tentative,
            seed_ready=false,
            has_measure_now=false,
            consecutive_missed=0,
            last_update_t=t,
            state_estimate_now=NAN_STATE6,
            covariance_estimate_now=_nan_cov6(),
            observer_pos_now=NAN_LOS,
            observer_pos_prev=NAN_LOS
        )
        println("Local track init | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s")
    end
    return tracks[slot_id]
end

# Assigns a new local slot ID for an orphan measurement for this observer, and increments the counter for the next assignment. 
@inline function _next_local_slot!(model::ObserverNavigationModel, observer_id::Int)::Int
    slot_id = model.next_local_slot[observer_id]
    model.next_local_slot[observer_id] += 1
    return slot_id
end

@inline function _los_angle_rad(a::SVector{3, Float64}, b::SVector{3, Float64})::Float64
    au = _safe_unit(a)
    bu = _safe_unit(b)
    return acos(clamp(dot(au, bu), -1.0, 1.0))
end

@inline function _passes_orphan_attach_gate(
    track::LocalTrack,
    measurement::LOSMeasurement
)::Bool
    track.last_meas === nothing && return false
    dt = measurement.t - track.last_meas.t
    (dt > 0.0) || return false
    _is_consecutive_measure_pair(measurement.t, track.last_meas.t) || return false
    dtheta = _los_angle_rad(track.last_meas.los_unit, measurement.los_unit)
    return dtheta <= ORPHAN_ATTACH_MAX_DTHETA_RAD
end

function _associate_orphan_measurement_to_local_slot!(
    model::ObserverNavigationModel,
    observer_id::Int,
    measurement::LOSMeasurement
)::Int
    tracks = model.local_tracks[observer_id]
    best_slot = 0
    best_score = Inf

    # Try to attach orphan to an existing tentative track first.
    for (slot_id, track) in tracks
        (track.status in (:tentative, :seed_ready, :iod_initialized)) || continue
        _passes_orphan_attach_gate(track, measurement) || continue
        score = _los_angle_rad(track.last_meas.los_unit, measurement.los_unit)
        if score < best_score
            best_score = score
            best_slot = slot_id
        end
    end

    if best_slot != 0
        return best_slot
    end

    # No compatible tentative track found: create a new local slot.
    return _next_local_slot!(model, observer_id)
end

@inline function _refresh_track_tick_context!(
    track::LocalTrack,
    observer_pos::SVector{3, Float64}
)
    track.has_measure_now = false
    track.observer_pos_prev = track.observer_pos_now
    track.observer_pos_now = observer_pos
end

# True when the track has two finite LOS samples separated by one nav tick.
@inline function _track_has_two_measures(
    model::ObserverNavigationModel,
    observer_id::Int,
    slot_id::Int
)::Bool
    tracks = model.local_tracks[observer_id]
    haskey(tracks, slot_id) || return false
    tr = tracks[slot_id]
    tr.last_meas === nothing && return false
    tr.prev_meas === nothing && return false
    last = tr.last_meas
    prev = tr.prev_meas
    return _is_consecutive_measure_pair(last.t, prev.t)
end

# Returns LOS rate computed from the latest two measures for this track if they are consecutive and valid, otherwise returns nothing.
@inline function _track_los_rate(
    model::ObserverNavigationModel,
    observer_id::Int,
    slot_id::Int
)::Union{Nothing, SVector{3, Float64}}
    tracks = model.local_tracks[observer_id]
    haskey(tracks, slot_id) || return nothing
    tr = tracks[slot_id]
    tr.last_meas === nothing && return nothing
    tr.prev_meas === nothing && return nothing
    last = tr.last_meas
    prev = tr.prev_meas
    _is_consecutive_measure_pair(last.t, prev.t) || return nothing
    dt = last.t - prev.t
    return (last.los_unit - prev.los_unit) / dt
end

# Returns latest LOS sample for (observer, slot) only if fresh at current nav epoch.
@inline function _latest_track_measure(
    model::ObserverNavigationModel,
    observer_id::Int,
    slot_id::Int,
    t::Float64
)::Union{Nothing, LOSMeasurement}
    tracks = model.local_tracks[observer_id]
    haskey(tracks, slot_id) || return nothing
    tr = tracks[slot_id]
    tr.last_meas === nothing && return nothing
    # Use only measurements committed at the current nav tick.
    tr.has_measure_now || return nothing
    meas = tr.last_meas
    abs(t - meas.t) <= NAVIGATION_DT_TOL_SEC || return nothing
    return meas
end

# Selects neighbor observers and measurement pairs usable for IOD.
# Match each neighbor with the local track that is most coherent with the
# seed track LOS by ray miss-distance at current tick.
@inline function _line_param_interval_in_sphere(
    origin::SVector{3, Float64},
    dir_u::SVector{3, Float64},
    sphere_center::SVector{3, Float64},
    sphere_radius_m::Float64
)::Union{Nothing, Tuple{Float64, Float64}}
    sphere_radius_m > 0.0 || return nothing
    Δ = origin - sphere_center
    b = dot(dir_u, Δ)
    c = dot(Δ, Δ) - sphere_radius_m^2
    disc = b^2 - c
    disc < 0.0 && return nothing
    sdisc = sqrt(max(disc, 0.0))
    t1 = -b - sdisc
    t2 = -b + sdisc
    return (min(t1, t2), max(t1, t2))
end

@inline function _line_param_interval_in_visibility_intersection(
    origin::SVector{3, Float64},
    dir_u::SVector{3, Float64},
    own_visibility_radius_m::Float64,
    other_origin::SVector{3, Float64},
    other_visibility_radius_m::Float64
)::Union{Nothing, Tuple{Float64, Float64}}
    own_visibility_radius_m > 0.0 || return nothing
    other_visibility_radius_m > 0.0 || return nothing

    own_interval = (0.0, own_visibility_radius_m)
    other_interval = _line_param_interval_in_sphere(
        origin,
        dir_u,
        other_origin,
        other_visibility_radius_m
    )
    other_interval === nothing && return nothing

    lo = max(own_interval[1], other_interval[1])
    hi = min(own_interval[2], other_interval[2])
    return hi >= lo ? (lo, hi) : nothing
end

@inline function _los_ray_miss_distance(
    r1::SVector{3, Float64},
    l1::SVector{3, Float64},
    visibility_r1_m::Float64,
    r2::SVector{3, Float64},
    l2::SVector{3, Float64},
    visibility_r2_m::Float64
)::Float64
    l1_norm = norm(l1)
    l2_norm = norm(l2)
    (l1_norm > 1e-12 && l2_norm > 1e-12) || return Inf
    l1u = l1 / l1_norm
    l2u = l2 / l2_norm

    i1 = _line_param_interval_in_visibility_intersection(r1, l1u, visibility_r1_m, r2, visibility_r2_m)
    i2 = _line_param_interval_in_visibility_intersection(r2, l2u, visibility_r2_m, r1, visibility_r1_m)
    (i1 === nothing || i2 === nothing) && return Inf

    s1_lo, s1_hi = i1
    s2_lo, s2_hi = i2

    d12 = dot(l1u, l2u)
    denom = 1.0 - d12^2

    min_d = Inf

    # Candidate from interior stationary point (when lines are not parallel).
    if denom > 1e-12
        w0 = r1 - r2
        d = dot(l1u, w0)
        e = dot(l2u, w0)
        s1_star = (d12 * e - d) / denom
        s2_star = (e - d12 * d) / denom
        if s1_lo <= s1_star <= s1_hi && s2_lo <= s2_star <= s2_hi
            p1 = r1 + s1_star * l1u
            p2 = r2 + s2_star * l2u
            min_d = min(min_d, norm(p1 - p2))
        end
    end

    # Edge candidates on the rectangle [s1_lo,s1_hi] x [s2_lo,s2_hi].
    for s1 in (s1_lo, s1_hi)
        p1 = r1 + s1 * l1u
        s2 = clamp(dot(p1 - r2, l2u), s2_lo, s2_hi)
        p2 = r2 + s2 * l2u
        min_d = min(min_d, norm(p1 - p2))
    end
    for s2 in (s2_lo, s2_hi)
        p2 = r2 + s2 * l2u
        s1 = clamp(dot(p2 - r1, l1u), s1_lo, s1_hi)
        p1 = r1 + s1 * l1u
        min_d = min(min_d, norm(p1 - p2))
    end

    return min_d
end

function _select_iod_neighbors(
    nav_model::ObserverNavigationModel,
    observer_id::Int,
    slot_id::Int,
    neighbor_ids::Vector{Int},
    t::Float64,
    u
)::Tuple{Vector{Tuple{Int, Int}}, Int, Int}
    # extract track latest meas and pos for the seed track, return empty if not available or not fresh enough for IOD at current tick.
    seed_meas = _latest_track_measure(nav_model, observer_id, slot_id, t)
    seed_meas === nothing && return (Tuple{Int, Int}[], 0, 0)
    seed_los = _safe_unit(seed_meas.los_unit)
    norm(seed_los) > 1e-12 || return (Tuple{Int, Int}[], 0, 0)
    # use perturbed observer pos for IOD neighbor selection.
    seed_pos = _known_observer_pos!(nav_model, observer_id, u, t)
    seed_visibility = nav_model.sensor.detection_range_m

    # For each neighbor, find the local track with the best LOS coherence with the seed track, and if the miss distance is within the threshold, consider it a valid neighbor for IOD.
    valid_neighbors = Tuple{Int, Int}[]
    multi_pass_events = 0
    multi_pass_extra = 0
    for neighbor_id in neighbor_ids
        best_slot = 0
        best_miss = Inf
        pass_count = 0
        for neighbor_slot in keys(nav_model.local_tracks[neighbor_id])
            # check if the neighbor track has the potential to be an IOD seed based on its own local measures.
            _track_has_two_measures(nav_model, neighbor_id, neighbor_slot) || continue
            # extract neighbor track latest measure and pos, skip if not available or not fresh enough for IOD at current tick.
            neighbor_meas = _latest_track_measure(nav_model, neighbor_id, neighbor_slot, t)
            neighbor_meas === nothing && continue
            neighbor_los = _safe_unit(neighbor_meas.los_unit)
            norm(neighbor_los) > 1e-12 || continue
            neighbor_pos = _known_observer_pos!(nav_model, neighbor_id, u, t)
            # calculate miss distance between seed and neighbor track LOS rays, and keep track of the best (lowest) miss distance for this neighbor.
            miss_d = _los_ray_miss_distance(
                seed_pos,
                seed_los,
                seed_visibility,
                neighbor_pos,
                neighbor_los,
                nav_model.sensor.detection_range_m
            )
            miss_d <= IOD_NEIGHBOR_MISS_DISTANCE_MAX_M && (pass_count += 1)
            if miss_d < best_miss
                best_miss = miss_d
                best_slot = neighbor_slot
            end
        end
        if pass_count > 1
            multi_pass_events += 1
            multi_pass_extra += (pass_count - 1)
        end

        # check if the min miss distance satisfies the treshold. If so this neighbor and its best slot are valid for IOD.
        if best_slot != 0 && best_miss <= IOD_NEIGHBOR_MISS_DISTANCE_MAX_M
            push!(valid_neighbors, (neighbor_id, best_slot))
        end
    end
    return valid_neighbors, multi_pass_events, multi_pass_extra
end

# Resolve which estimator mode is available for this local track.
@inline function _track_match_mode(
    model,
    observer_id::Int,
    slot_id::Int
)::Union{Nothing, Symbol}
    cap = _slot_capacity(model)
    slot_id <= cap || return nothing

    if model.filter_initialized[observer_id, slot_id]
        x_try = model.state_pred[observer_id, slot_id]
        P_try = model.covariance_pred[observer_id, slot_id]
        if !(_is_finite_state(x_try) && _is_finite_cov(P_try))
            x_try = model.state[observer_id, slot_id]
            P_try = model.covariance[observer_id, slot_id]
        end
        (_is_finite_state(x_try) && _is_finite_cov(P_try)) || return nothing
        return :filter
    end

    if model.iod_initialized[observer_id, slot_id] &&
       _is_finite_state(model.iod_estimate_state[observer_id, slot_id]) &&
       _is_finite_cov(model.iod_estimate_covariance[observer_id, slot_id])
        return :iod
    end

    return nothing
end

@inline function _track_state_cov_for_matching(
    model,
    observer_id::Int,
    slot_id::Int,
    mode::Symbol
)::Union{Nothing, Tuple{SVector{6, Float64}, Matrix{Float64}}}
    cap = _slot_capacity(model)
    slot_id <= cap || return nothing

    if mode === :filter
        model.filter_initialized[observer_id, slot_id] || return nothing
        x_try = model.state_pred[observer_id, slot_id]
        P_try = model.covariance_pred[observer_id, slot_id]
        if !(_is_finite_state(x_try) && _is_finite_cov(P_try))
            x_try = model.state[observer_id, slot_id]
            P_try = model.covariance[observer_id, slot_id]
        end
        (_is_finite_state(x_try) && _is_finite_cov(P_try)) || return nothing
        return (x_try, P_try)
    end

    if mode === :iod
        # In IOD bootstrap mode we still use IOD state/cov, even if this local
        # track was already promoted to filter_initialized earlier in the same
        # nav tick. This avoids order-dependent split init across observers.
        model.iod_initialized[observer_id, slot_id] || return nothing
        x_try = model.iod_estimate_state[observer_id, slot_id]
        P_try = model.iod_estimate_covariance[observer_id, slot_id]
        (_is_finite_state(x_try) && _is_finite_cov(P_try)) || return nothing
        return (x_try, P_try)
    end

    return nothing
end

@inline function _mahalanobis_distance_sq(
    dx::SVector{6, Float64},
    P_sum::Matrix{Float64}
)::Float64
    (size(P_sum, 1) == 6 && size(P_sum, 2) == 6) || return Inf
    all(isfinite, P_sum) || return Inf

    P_sym = 0.5 * (P_sum + P_sum')
    all(isfinite, P_sym) || return Inf

    reg_levels = (
        CONSENSUS_MATCH_COV_REG_EPS,
        1e-4,
        1e-2,
        1.0
    )
    dx_vec = Vector{Float64}(dx)
    for eps in reg_levels
        P_reg = P_sym + eps * Matrix(I, 6, 6)
        chol = cholesky(Hermitian(P_reg); check=false)
        isposdef(chol) || continue
        y = chol \ dx_vec
        return dot(y, y)
    end
    return Inf
end

function _best_slot_match_by_mahalanobis(
    model,
    source_observer_id::Int,
    source_slot_id::Int,
    candidate_observer_id::Int,
    mode::Symbol
)::Tuple{Int, Float64, Int}
    source_state_cov = _track_state_cov_for_matching(
        model,
        source_observer_id,
        source_slot_id,
        mode
    )
    source_state_cov === nothing && return (0, Inf, 0)
    x_source, P_source = source_state_cov

    tracks = model.nav.local_tracks[candidate_observer_id]
    best_slot = 0
    best_d2 = Inf
    pass_count = 0
    for candidate_slot_id in keys(tracks)
        candidate_state_cov = _track_state_cov_for_matching(
            model,
            candidate_observer_id,
            candidate_slot_id,
            mode
        )
        candidate_state_cov === nothing && continue
        x_candidate, P_candidate = candidate_state_cov
        d2 = _mahalanobis_distance_sq(x_source - x_candidate, P_source + P_candidate)
        d2 <= (mode === :iod ? CONSENSUS_MATCH_MAHAL_IOD_MAX_D2 : CONSENSUS_MATCH_MAHAL_FILTER_MAX_D2) && (pass_count += 1)
        if d2 < best_d2
            best_d2 = d2
            best_slot = candidate_slot_id
        end
    end
    return (best_slot, best_d2, pass_count)
end

@inline function _neighbor_has_same_target_candidate(
    model,
    candidate_observer_id::Int,
    source_gid::Int,
    mode::Symbol
)::Bool
    source_gid > 0 || return false
    tracks = model.nav.local_tracks[candidate_observer_id]
    for candidate_slot_id in keys(tracks)
        tr = tracks[candidate_slot_id]
        (tr.last_meas !== nothing && tr.last_meas.target == source_gid) || continue
        candidate_state_cov = _track_state_cov_for_matching(
            model,
            candidate_observer_id,
            candidate_slot_id,
            mode
        )
        candidate_state_cov === nothing && continue
        return true
    end
    return false
end

# Cross-track one-shot association between two observers.
# Rows = source observer tracks, Cols = neighbor observer tracks.
# Gate on local Mahalanobis distance, then degree-1 reductions, then Hungarian on residual.
function _cross_track_match_for_consensus_pairwise(
    fusion_model,
    source_observer_id::Int,
    source_slot_ids::Vector{Int},
    candidate_observer_id::Int,
    candidate_slot_ids::Vector{Int}
)::Tuple{Dict{Int, Int}, Int, Int}
    matched_slot_by_source_slot = Dict{Int, Int}()
    n_rows = length(source_slot_ids)
    n_cols = length(candidate_slot_ids)
    n_rows == 0 && return (matched_slot_by_source_slot, 0, 0)

    row_modes = fill(:none, n_rows)
    row_source_gid = zeros(Int, n_rows)
    row_attempted = falses(n_rows)
    pass_count_by_row = zeros(Int, n_rows)
    local_d2 = fill(Inf, n_rows, n_cols)

    # populate the local Mahalanobis distance matrix for all pairs.
    for r in 1:n_rows
        source_slot_id = source_slot_ids[r]
        source_track = get(fusion_model.nav.local_tracks[source_observer_id], source_slot_id, nothing)
        row_source_gid[r] = (source_track === nothing || source_track.last_meas === nothing) ? 0 : source_track.last_meas.target

        mode = _track_match_mode(fusion_model, source_observer_id, source_slot_id)
        mode === nothing && continue
        row_modes[r] = mode
        row_attempted[r] = true
        fusion_model.tt_attempt_total += 1
        d2_gate = mode === :iod ? CONSENSUS_MATCH_MAHAL_IOD_MAX_D2 : CONSENSUS_MATCH_MAHAL_FILTER_MAX_D2

        source_state_cov = _track_state_cov_for_matching(
            fusion_model,
            source_observer_id,
            source_slot_id,
            mode
        )
        source_state_cov === nothing && continue
        x_source, P_source = source_state_cov

        for c in 1:n_cols
            candidate_slot_id = candidate_slot_ids[c]
            candidate_state_cov = _track_state_cov_for_matching(
                fusion_model,
                candidate_observer_id,
                candidate_slot_id,
                mode
            )
            candidate_state_cov === nothing && continue
            x_candidate, P_candidate = candidate_state_cov
            d2 = _mahalanobis_distance_sq(x_source - x_candidate, P_source + P_candidate)
            if isfinite(d2) && d2 <= d2_gate
                local_d2[r, c] = d2
                pass_count_by_row[r] += 1
            end
        end
    end

    # Count multi-pass events, where a single source track passes the gate with multiple candidate tracks. 
    multi_pass_events = 0
    multi_pass_extra = 0
    for r in 1:n_rows
        row_attempted[r] || continue
        if pass_count_by_row[r] > 1
            multi_pass_events += 1
            multi_pass_extra += (pass_count_by_row[r] - 1)
        end
    end

    assigned_pairs = Tuple{Int, Int}[]
    if n_cols > 0
        active_rows = collect(row_attempted)
        active_cols = collect(trues(n_cols))
        dummy_orphans = collect(falses(n_cols))
        # First resolve all degree-1 reductions iteratively, then run Hungarian on the residual.
        while _assign_degree1_and_remove!(active_rows, active_cols, local_d2, assigned_pairs, dummy_orphans)
        end

        residual_rows = findall(active_rows)
        residual_cols = findall(active_cols)
        if !isempty(residual_rows) && !isempty(residual_cols)
            residual_cost = local_d2[residual_rows, residual_cols]
            valid_row_mask = [any(isfinite, residual_cost[ri, :]) for ri in 1:size(residual_cost, 1)]
            valid_col_mask = [any(isfinite, residual_cost[:, ci]) for ci in 1:size(residual_cost, 2)]

            residual_rows = [residual_rows[ri] for ri in 1:length(residual_rows) if valid_row_mask[ri]]
            residual_cols = [residual_cols[ci] for ci in 1:length(residual_cols) if valid_col_mask[ci]]

            if !isempty(residual_rows) && !isempty(residual_cols)
                residual_cost = residual_cost[valid_row_mask, valid_col_mask]
                hungarian_assignment = _hungarian_assign_min_cost(residual_cost) # row -> col (0 if unmatched)
                for (ri, ci) in enumerate(hungarian_assignment)
                    ci == 0 && continue
                    push!(assigned_pairs, (residual_rows[ri], residual_cols[ci]))
                end
            end
        end
    end

    matched_rows = Set{Int}()
    for (r, c) in assigned_pairs
        source_slot_id = source_slot_ids[r]
        candidate_slot_id = candidate_slot_ids[c]
        matched_slot_by_source_slot[source_slot_id] = candidate_slot_id
        push!(matched_rows, r)

        fusion_model.tt_commit_total += 1
        source_gid = row_source_gid[r]
        candidate_track = get(fusion_model.nav.local_tracks[candidate_observer_id], candidate_slot_id, nothing)
        candidate_gid = (candidate_track === nothing || candidate_track.last_meas === nothing) ? 0 : candidate_track.last_meas.target
        if source_gid > 0 && candidate_gid > 0
            if source_gid == candidate_gid
                fusion_model.tt_commit_correct_total += 1
            else
                fusion_model.tt_commit_wrong_total += 1
            end
        else
            fusion_model.tt_commit_unknown_total += 1
        end
    end

    for r in 1:n_rows
        row_attempted[r] || continue
        r in matched_rows && continue

        fusion_model.tt_skipped_total += 1
        if pass_count_by_row[r] == 0
            fusion_model.tt_no_candidate_total += 1
        else
            fusion_model.tt_gate_fail_total += 1
        end

        source_gid = row_source_gid[r]
        mode = row_modes[r]
        has_same_target = mode == :none ? false : _neighbor_has_same_target_candidate(
            fusion_model,
            candidate_observer_id,
            source_gid,
            mode
        )
        if has_same_target
            fusion_model.tt_skipped_same_target_present_total += 1
        elseif source_gid > 0
            fusion_model.tt_skipped_no_same_target_total += 1
        else
            fusion_model.tt_skipped_unknown_source_target_total += 1
        end
    end

    return matched_slot_by_source_slot, multi_pass_events, multi_pass_extra
end

Base.@kwdef struct MatchGroup
    tracks::Vector{Tuple{Int, Int}} = Tuple{Int, Int}[]
    selected_tracks::Vector{Tuple{Int, Int}} = Tuple{Int, Int}[]
end

@inline function _first_track_per_observer(
    keys::Vector{Tuple{Int, Int}}
)::Vector{Tuple{Int, Int}}
    selected = Tuple{Int, Int}[]
    seen = Set{Int}()
    for key in keys
        observer_id = key[1]
        observer_id in seen && continue
        push!(seen, observer_id)
        push!(selected, key)
    end
    return selected
end

function _build_match_groups(
    model,
    active_tracks::Vector{Tuple{Int, Int}},
    neighbor_map::Dict{Int, Vector{Int}}
)::Tuple{Vector{MatchGroup}, Int, Int, Int}
    active_unique = sort(unique(active_tracks))
    isempty(active_unique) && return (MatchGroup[], 0, 0, 0)

    # every active track creates a node in the graph, and we add undirected edges between nodes that are matched through cross-track consensus. 
    #Then each connected component of the graph is a match group.
    active_set = Set(active_unique)
    graph = Dict{Tuple{Int, Int}, Set{Tuple{Int, Int}}}(key => Set{Tuple{Int, Int}}() for key in active_unique)

    groups = MatchGroup[]
    fallback_count = 0
    mahal_multi_pass_events = 0
    mahal_multi_pass_extra = 0

    function _record_group_target_consistency!(selected_tracks::Vector{Tuple{Int, Int}})
        model.consensus_group_total += 1
        gids = Int[]
        has_unknown = false

        for (obs, slot) in selected_tracks
            
            tr = get(model.nav.local_tracks[obs], slot, nothing)
            # get the global ID from the latest measure of this track, if available.
            gid = (tr === nothing || tr.last_meas === nothing) ? 0 : tr.last_meas.target
            if gid > 0
                push!(gids, gid)
            else
                has_unknown = true
            end
        end
        if isempty(gids) || has_unknown
            model.consensus_group_unknown_total += 1
            return nothing
        end
        if length(unique(gids)) == 1
            model.consensus_group_same_target_total += 1
        else
            model.consensus_group_mixed_target_total += 1
        end
        return nothing
    end

    active_slots_by_observer = Dict{Int, Vector{Int}}()
    for (observer_id, slot_id) in active_unique
        push!(get!(active_slots_by_observer, observer_id, Int[]), slot_id)
    end
    for slot_ids in values(active_slots_by_observer)
        sort!(slot_ids)
    end

    # For each directed observer-neighbor pair, solve TT association in one shot.
    for observer_id in sort(collect(keys(active_slots_by_observer)))
        source_slots = active_slots_by_observer[observer_id]
        neighbor_ids = sort(unique(get(neighbor_map, observer_id, Int[])))
        for candidate_observer in neighbor_ids
            candidate_observer == observer_id && continue
            haskey(active_slots_by_observer, candidate_observer) || continue
            candidate_slots = active_slots_by_observer[candidate_observer]

            matched_slot_by_source_slot, pass_events, pass_extra = _cross_track_match_for_consensus_pairwise(
                model,
                observer_id,
                source_slots,
                candidate_observer,
                candidate_slots
            )
            mahal_multi_pass_events += pass_events
            mahal_multi_pass_extra += pass_extra

            for (source_slot, matched_slot) in matched_slot_by_source_slot
                key = (observer_id, source_slot)
                matched_key = (candidate_observer, matched_slot)
                key in active_set || continue
                matched_key in active_set || continue
                push!(graph[key], matched_key)
                push!(graph[matched_key], key)
            end
        end
    end

    # build connected components of the graph.
    # this is built in a centralized way here for simplicity, but it can be done in a distributed way as well if needed.
    visited = Set{Tuple{Int, Int}}()
    for start_key in active_unique
        start_key in visited && continue

        stack = [start_key]
        component_tracks = Tuple{Int, Int}[]
        while !isempty(stack)
            key = pop!(stack)
            key in visited && continue
            push!(visited, key)
            push!(component_tracks, key)
            for neighbor_key in graph[key]
                neighbor_key in visited && continue
                push!(stack, neighbor_key)
            end
        end
        sort!(component_tracks)

        # Baseline behavior: connected-component grouping + deterministic
        # one-track-per-observer selection (no degree heuristic).
        selected_tracks = _first_track_per_observer(component_tracks)

        # Optional global-id fallback: if enabled and this connected component
        # mixes multiple true IDs, split by ID and log a fallback event.
        if CONSENSUS_GROUP_FALLBACK_ENABLE
            gids = Dict{Tuple{Int, Int}, Int}()
            for key in component_tracks
                obs, slot = key
                meas = model.nav.local_tracks[obs][slot].last_meas
                gids[key] = (meas === nothing) ? 0 : meas.target
            end
            nonzero_ids = sort(unique([gid for gid in values(gids) if gid > 0]))
            if length(nonzero_ids) > 1
                fallback_count += 1
                for gid in nonzero_ids
                    gid_tracks = sort([key for key in component_tracks if gids[key] == gid])
                    isempty(gid_tracks) && continue
                    gid_selected = _first_track_per_observer(gid_tracks)
                    _record_group_target_consistency!(gid_selected)
                    push!(groups, MatchGroup(tracks=gid_tracks, selected_tracks=gid_selected))
                end
                unknown_tracks = sort([key for key in component_tracks if gids[key] == 0])
                if !isempty(unknown_tracks)
                    unknown_selected = _first_track_per_observer(unknown_tracks)
                    _record_group_target_consistency!(unknown_selected)
                    push!(groups, MatchGroup(tracks=unknown_tracks, selected_tracks=unknown_selected))
                end
                continue
            end
        end

        _record_group_target_consistency!(selected_tracks)
        push!(groups, MatchGroup(tracks=component_tracks, selected_tracks=selected_tracks))
    end
    return groups, fallback_count, mahal_multi_pass_events, mahal_multi_pass_extra
end

Base.@kwdef struct MeasAssocCandidate
    slot_id::Int
    d2::Float64
    id_match::Bool
end

@inline function _measurement_mahalanobis_d2(
    model::ObserverNavigationModel,
    track::LocalTrack,
    measurement::LOSMeasurement
)::Float64
    return _measurement_mahalanobis_d2_from_state_cov(
        model,
        track.state_estimate_now,
        track.covariance_estimate_now,
        track.observer_pos_now,
        measurement
    )
end

@inline function _measurement_mahalanobis_d2_from_state_cov(
    model::ObserverNavigationModel,
    x::SVector{6, Float64},
    P::Matrix{Float64},
    observer_pos::SVector{3, Float64},
    measurement::LOSMeasurement
)::Float64
    (_is_finite_state(x) && _is_finite_cov(P)) || return Inf
    (size(P, 1) == 6 && size(P, 2) == 6) || return Inf

    all(isfinite, observer_pos) || return Inf

    z = measurement.los_unit
    all(isfinite, z) || return Inf
    z_unit = _safe_unit(z)
    norm(z_unit) > 1e-12 || return Inf
    h = _safe_unit(x[1:3] - observer_pos)
    norm(h) > 1e-12 || return Inf

    H = _measurement_jacobian(x, observer_pos)
    σ2 = model.sensor.sigma_theta_rad^2
    I3 = Matrix(I, 3, 3)
    ν = Vector(z_unit - h)
    h_vec = Vector(h)
    R = σ2 * (I3 - h_vec * h_vec')
    S = H * P * H' + R
    S_sym = 0.5 * (S + S')
    all(isfinite, S_sym) || return Inf

    S_pinv = try
        pinv(S_sym)
    catch
        return Inf
    end
    all(isfinite, S_pinv) || return Inf

    y = S_pinv * ν
    d2 = dot(ν, y)
    return isfinite(d2) ? d2 : Inf
end

@inline function _collaborative_score_for_candidate(
    model::ObserverNavigationModel,
    observer_id::Int,
    observer_pos::SVector{3, Float64},
    measurement::LOSMeasurement,
    candidate::MeasAssocCandidate
)::Union{Nothing, Float64}
    neighbor_ids = [
        nid for nid in model.comms.neighbors[observer_id]
        if (nid != observer_id) && (nid in model.observer_idxs)
    ]

    tracks = model.local_tracks[observer_id]
    haskey(tracks, candidate.slot_id) || return nothing
    track = tracks[candidate.slot_id]
    x_src = track.state_estimate_now
    P_src = track.covariance_estimate_now
    (_is_finite_state(x_src) && _is_finite_cov(P_src)) || return nothing

    support_scores = Float64[]
    MEAS_ASSOC_A_INCLUDE_SELF_SCORE && push!(support_scores, candidate.d2)
    for neighbor_id in neighbor_ids
        best_neighbor_slot = 0
        best_neighbor_d2 = Inf
        for (neighbor_slot, neighbor_track) in model.local_tracks[neighbor_id]
            (neighbor_track.status == :filter_initialized) || continue
            x_nb = neighbor_track.state_estimate_now
            P_nb = neighbor_track.covariance_estimate_now
            (_is_finite_state(x_nb) && _is_finite_cov(P_nb)) || continue
            d2_track = _mahalanobis_distance_sq(x_src - x_nb, P_src + P_nb)
            if d2_track < best_neighbor_d2
                best_neighbor_d2 = d2_track
                best_neighbor_slot = neighbor_slot
            end
        end
        (best_neighbor_slot != 0 && best_neighbor_d2 <= CONSENSUS_MATCH_MAHAL_FILTER_MAX_D2) || continue
        neighbor_track = model.local_tracks[neighbor_id][best_neighbor_slot]
        d2_meas = _measurement_mahalanobis_d2_from_state_cov(
            model,
            neighbor_track.state_estimate_now,
            neighbor_track.covariance_estimate_now,
            observer_pos,
            measurement
        )
        isfinite(d2_meas) && push!(support_scores, d2_meas)
    end
    isempty(support_scores) && return nothing
    return sum(support_scores) / length(support_scores)
end

# Local measurement-to-filter-initialized-track candidates by Mahalanobis gate.
@inline function _local_association_measure_to_track_candidates(
    model::ObserverNavigationModel,
    observer_id::Int,
    measurement::LOSMeasurement
)::Vector{MeasAssocCandidate}
    tracks = model.local_tracks[observer_id]
    candidates = MeasAssocCandidate[]

    for (slot_id, track) in tracks
        # Keep filter tracks eligible even if they were marked :assigned in a
        # previous observer pass before fusion status sync.
        (track.status == :filter_initialized) || continue
        d2 = _measurement_mahalanobis_d2(model, track, measurement)
        (isfinite(d2) && d2 <= MEAS_ASSOC_MAHAL_MAX_D2) || continue
        id_match = (track.last_meas !== nothing) && (track.last_meas.target == measurement.target)
        push!(candidates, MeasAssocCandidate(slot_id=slot_id, d2=d2, id_match=id_match))
    end
    return candidates
end

# Deterministic local picker: prioritize global-ID match, then best d2.
@inline function _pick_candidate_for_measurement(
    candidates::Vector{MeasAssocCandidate}
)::MeasAssocCandidate
    best = candidates[1]
    for c in candidates[2:end]
        if (c.id_match != best.id_match && c.id_match) || (c.id_match == best.id_match && c.d2 < best.d2)
            best = c
        end
    end
    return best
end

function _hungarian_rows_leq_cols(cost::Matrix{Float64})::Vector{Int}
    n_rows, n_cols = size(cost)
    assignment = zeros(Int, n_rows)
    n_rows == 0 && return assignment

    u = zeros(Float64, n_rows + 1)
    v = zeros(Float64, n_cols + 1)
    p = zeros(Int, n_cols + 1)
    way = zeros(Int, n_cols + 1)

    for i in 1:n_rows
        p[1] = i
        j0 = 1
        minv = fill(Inf, n_cols + 1)
        used = falses(n_cols + 1)
        while true
            used[j0] = true
            i0 = p[j0]
            delta = Inf
            j1 = 1
            for j in 2:(n_cols + 1)
                used[j] && continue
                cur = cost[i0, j - 1] - u[i0 + 1] - v[j]
                if cur < minv[j]
                    minv[j] = cur
                    way[j] = j0
                end
                if minv[j] < delta
                    delta = minv[j]
                    j1 = j
                end
            end
            for j in 1:(n_cols + 1)
                if used[j]
                    u[p[j] + 1] += delta
                    v[j] -= delta
                else
                    minv[j] -= delta
                end
            end
            j0 = j1
            p[j0] == 0 && break
        end

        while true
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            j0 == 1 && break
        end
    end

    for j in 2:(n_cols + 1)
        row = p[j]
        row > 0 || continue
        assignment[row] = j - 1
    end
    return assignment
end

function _hungarian_assign_min_cost(cost::Matrix{Float64})::Vector{Int}
    n_rows, n_cols = size(cost)
    n_rows == 0 && return Int[]
    n_cols == 0 && return zeros(Int, n_rows)

    finite_vals = [x for x in cost if isfinite(x)]
    if isempty(finite_vals)
        return zeros(Int, n_rows)
    end

    max_finite = maximum(finite_vals)
    big_cost = max(1e12, max_finite * 1e6 + 1.0)
    safe_cost = [isfinite(x) ? x : big_cost for x in cost]

    if n_rows <= n_cols
        raw = _hungarian_rows_leq_cols(safe_cost)
        for r in eachindex(raw)
            c = raw[r]
            if c == 0 || !isfinite(cost[r, c]) || cost[r, c] >= big_cost / 10.0
                raw[r] = 0
            end
        end
        return raw
    end

    raw_t = _hungarian_rows_leq_cols(permutedims(safe_cost))
    assignment = zeros(Int, n_rows)
    for col in 1:n_cols
        row = raw_t[col]
        row == 0 && continue
        isfinite(cost[row, col]) || continue
        assignment[row] = col
    end
    return assignment
end

function _assign_degree1_and_remove!(
    active_rows::AbstractVector{Bool},
    active_cols::AbstractVector{Bool},
    local_d2::Matrix{Float64},
    assigned_pairs::Vector{Tuple{Int, Int}},
    orphan_cols::AbstractVector{Bool}
)::Bool
    changed = false
    n_rows, n_cols = size(local_d2)

    for c in 1:n_cols
        active_cols[c] || continue
        deg = 0
        row_idx = 0
        for r in 1:n_rows
            active_rows[r] || continue
            isfinite(local_d2[r, c]) || continue
            deg += 1
            row_idx = r
            deg > 1 && break
        end
        if deg == 0
            active_cols[c] = false
            orphan_cols[c] = true
            changed = true
        elseif deg == 1 && active_rows[row_idx] && active_cols[c]
            push!(assigned_pairs, (row_idx, c))
            active_rows[row_idx] = false
            active_cols[c] = false
            changed = true
        end
    end

    for r in 1:n_rows
        active_rows[r] || continue
        deg = 0
        col_idx = 0
        for c in 1:n_cols
            active_cols[c] || continue
            isfinite(local_d2[r, c]) || continue
            deg += 1
            col_idx = c
            deg > 1 && break
        end
        if deg == 0
            active_rows[r] = false
            changed = true
        elseif deg == 1 && active_rows[r] && active_cols[col_idx]
            push!(assigned_pairs, (r, col_idx))
            active_rows[r] = false
            active_cols[col_idx] = false
            changed = true
        end
    end

    return changed
end

function _split_assigned_and_orphan_measurements(
    model::ObserverNavigationModel,
    observer_id::Int,
    raw_measurements::Vector{LOSMeasurement},
    t::Float64,
    observer_pos::SVector{3, Float64}
)::Tuple{Dict{Int, LOSMeasurement}, Vector{LOSMeasurement}, Int, Int}
    assigned_by_slot = Dict{Int, LOSMeasurement}()
    orphan_measurements = LOSMeasurement[]
    conflict_count = 0
    skipped_count = 0

    # Count skipped associations only after filter bootstrap starts for this observer.
    count_skips_this_tick = any(
        track -> track.status == :filter_initialized,
        values(model.local_tracks[observer_id])
    )

    # only consider filter-initialized tracks as candidates for local association
    tracks = model.local_tracks[observer_id]
    slot_ids = sort([
        slot_id for (slot_id, track) in tracks
        if track.status == :filter_initialized
    ])

    # build association matrix. Rows = filter-initialized track slots, Cols = raw measurements. Cell = Mahalanobis d2 or Inf if gated out.
    n_rows = length(slot_ids)
    n_cols = length(raw_measurements)
    if n_cols == 0
        return assigned_by_slot, orphan_measurements, conflict_count, skipped_count
    end

    if n_rows == 0
        append!(orphan_measurements, raw_measurements)
        return assigned_by_slot, orphan_measurements, conflict_count, skipped_count
    end

    local_d2 = fill(Inf, n_rows, n_cols)
    local_id_match = falses(n_rows, n_cols)
    for r in 1:n_rows
        slot_id = slot_ids[r]
        track = tracks[slot_id]
        for c in 1:n_cols
            measurement = raw_measurements[c]
            d2 = _measurement_mahalanobis_d2(model, track, measurement)
            if isfinite(d2) && d2 <= MEAS_ASSOC_MAHAL_MAX_D2
                local_d2[r, c] = d2
                local_id_match[r, c] = (track.last_meas !== nothing) && (track.last_meas.target == measurement.target)
            end
        end
    end

    # First assign 1-1 matches and compute residual matrix.
    active_rows = collect(trues(n_rows))
    active_cols = collect(trues(n_cols))
    orphan_cols = collect(falses(n_cols))
    assigned_pairs = Tuple{Int, Int}[]

    while _assign_degree1_and_remove!(active_rows, active_cols, local_d2, assigned_pairs, orphan_cols)
    end

    # compute collaborative scores for residual matrix.
    residual_rows = findall(active_rows)
    residual_cols = findall(active_cols)
    if !isempty(residual_rows) && !isempty(residual_cols)
        residual_scores = fill(Inf, length(residual_rows), length(residual_cols))
        for (ri, r) in enumerate(residual_rows)
            slot_id = slot_ids[r]
            for (ci, c) in enumerate(residual_cols)
                d2_local = local_d2[r, c]
                isfinite(d2_local) || continue
                measurement = raw_measurements[c]
                candidate = MeasAssocCandidate(slot_id=slot_id, d2=d2_local, id_match=local_id_match[r, c])
                score = _collaborative_score_for_candidate(model, observer_id, observer_pos, measurement, candidate)
                score === nothing && continue
                residual_scores[ri, ci] = score
            end
        end

        # Hungarian assignment on residual score matrix.
        valid_row_mask = [any(isfinite, residual_scores[ri, :]) for ri in 1:size(residual_scores, 1)]
        valid_col_mask = [any(isfinite, residual_scores[:, ci]) for ci in 1:size(residual_scores, 2)]
        residual_rows = [residual_rows[ri] for ri in 1:length(residual_rows) if valid_row_mask[ri]]
        residual_cols = [residual_cols[ci] for ci in 1:length(residual_cols) if valid_col_mask[ci]]
        if !isempty(residual_rows) && !isempty(residual_cols)
            residual_scores = residual_scores[valid_row_mask, valid_col_mask]
            hungarian_assignment = _hungarian_assign_min_cost(residual_scores) # row -> col (0 if unmatched)
            for (ri, ci) in enumerate(hungarian_assignment)
                ci == 0 && continue
                push!(assigned_pairs, (residual_rows[ri], residual_cols[ci]))
            end

            # Columns still active but unmatched become orphan measurements.
            matched_cols = Set([residual_cols[ci] for ci in hungarian_assignment if ci > 0])
            for c in residual_cols
                c in matched_cols && continue
                orphan_cols[c] = true
            end
        else
            for c in residual_cols
                orphan_cols[c] = true
            end
        end
    end

    # check ratio gate for ambigous assignments and update counters.
    for (r, c) in assigned_pairs
        measurement = raw_measurements[c]
        slot_id = slot_ids[r]

        candidate_rows = Int[]
        for rr in 1:n_rows
            isfinite(local_d2[rr, c]) && push!(candidate_rows, rr)
        end

        chosen_candidate = MeasAssocCandidate(slot_id=slot_id, d2=local_d2[r, c], id_match=local_id_match[r, c])
        chosen = chosen_candidate
        ratio_failed = false
        ratio_used = false
        is_ambiguous = (length(candidate_rows) > 1)

        if is_ambiguous
            conflict_count += 1
            model.disambiguation_calls_total[observer_id] += 1

            score_pairs = Tuple{MeasAssocCandidate, Float64}[]
            for rr in candidate_rows
                cand = MeasAssocCandidate(
                    slot_id=slot_ids[rr],
                    d2=local_d2[rr, c],
                    id_match=local_id_match[rr, c]
                )
                sc = _collaborative_score_for_candidate(model, observer_id, observer_pos, measurement, cand)
                sc === nothing && continue
                push!(score_pairs, (cand, sc))
            end

            if length(score_pairs) >= 2
                chosen_score = Inf
                alt_best = Inf
                for (cand, sc) in score_pairs
                    if cand.slot_id == chosen_candidate.slot_id
                        chosen_score = sc
                    else
                        alt_best = min(alt_best, sc)
                    end
                end
                if isfinite(chosen_score) && isfinite(alt_best)
                    ratio = chosen_score <= 1e-12 ? (alt_best > 1e-12 ? Inf : 1.0) : (alt_best / chosen_score)
                    ratio_used = true
                    if ratio < MEAS_ASSOC_DISAMBIG_RATIO_MIN
                        ratio_failed = true
                        chosen = nothing
                    end
                end
            end

            if ratio_failed
                model.disambiguation_ratio_fail_total[observer_id] += 1
            elseif ratio_used && chosen !== nothing
                model.disambiguation_success_total[observer_id] += 1
                if chosen.id_match
                    model.disambiguation_ratio_success_correct_total[observer_id] += 1
                else
                    model.disambiguation_ratio_success_wrong_total[observer_id] += 1
                end
            end
        end

        if chosen === nothing
            if is_ambiguous
                model.meas_commit_ambig_dropped_total[observer_id] += 1
            end
            if count_skips_this_tick
                skipped_count += 1
            end
            orphan_cols[c] = true
            continue
        end

        # One-shot output is already 1-to-1; this is just a sanity guard.
        if haskey(assigned_by_slot, chosen.slot_id)
            conflict_count += 1
            model.b_conflict_events_total[observer_id] += 1
            if is_ambiguous
                model.meas_commit_ambig_dropped_total[observer_id] += 1
            end
            if count_skips_this_tick
                skipped_count += 1
                model.skipped_collision_B_total[observer_id] += 1
            end
            orphan_cols[c] = true
            continue
        end

        assigned_by_slot[chosen.slot_id] = measurement
        if is_ambiguous
            model.meas_commit_ambig_total[observer_id] += 1
            if chosen.id_match
                model.meas_commit_ambig_correct_total[observer_id] += 1
            else
                model.meas_commit_ambig_wrong_total[observer_id] += 1
            end
        else
            model.meas_commit_unique_total[observer_id] += 1
            if chosen.id_match
                model.meas_commit_unique_correct_total[observer_id] += 1
            else
                model.meas_commit_unique_wrong_total[observer_id] += 1
            end
        end
    end

    for c in 1:n_cols
        if orphan_cols[c]
            push!(orphan_measurements, raw_measurements[c])
        end
    end

    return assigned_by_slot, orphan_measurements, conflict_count, skipped_count
end

# Commit the selected measurements to their associated local tracks, creating/updating tracks as needed, and update the assigned count for logging/analysis.
function _commit_assigned_measurements_to_local_tracks!(
    model::ObserverNavigationModel,
    observer_id::Int,
    selected_by_slot::Dict{Int, LOSMeasurement},
    t::Float64,
    observer_pos::SVector{3, Float64}
)
    tracks = model.local_tracks[observer_id]
    for (slot_id, measurement) in selected_by_slot
        # By construction selected slots come from existing filter tracks.
        track = tracks[slot_id] 
        track.prev_meas = track.last_meas
        track.last_meas = measurement
        track.has_measure_now = true
        track.consecutive_missed = 0
        track.last_update_t = t
        track.status = :filter_initialized
        track.observer_pos_prev = track.observer_pos_now
        track.observer_pos_now = observer_pos
    end
end

# For measurements that were not selected for assignment to any existing track, create or update tentative local tracks as needed.
function _create_or_update_tentative_tracks_from_orphans!(
    model::ObserverNavigationModel,
    observer_id::Int,
    orphan_measurements::Vector{LOSMeasurement},
    t::Float64,
    observer_pos::SVector{3, Float64}
)::Int
    orphan_count = 0
    for measurement in orphan_measurements
        orphan_count += 1
        slot_id = _associate_orphan_measurement_to_local_slot!(model, observer_id, measurement)
        track = _ensure_local_track!(model, observer_id, slot_id, t)
        track.prev_meas = track.last_meas
        track.last_meas = measurement
        track.has_measure_now = true
        track.consecutive_missed = 0
        track.last_update_t = t
        track.status = :tentative
        track.observer_pos_prev = track.observer_pos_now
        track.observer_pos_now = observer_pos
        track.seed_ready = _track_has_two_measures(model, observer_id, slot_id)
        if track.seed_ready
            track.status = :seed_ready
        end
    end
    return orphan_count
end


function SimulationModel.calcNavigationEffect!(
    model::ObserverNavigationModel,
    u,
    p,
    t::Float64,
    sat_idx::Int
)   
    # Only for observers.
    sat_idx in model.observer_idxs || return nothing
    observer_id = sat_idx
    # perturb observer position with additive gaussian noise to replicate self-OD inaccuracies.
    observer_pos = _known_observer_pos!(model, observer_id, u, t)

    # Refresh per-tick local-track context for this observer.
    observer_tracks = model.local_tracks[observer_id]
    for track in values(observer_tracks)
        _refresh_track_tick_context!(track, observer_pos)
    end

    # Raw detections at this navigation tick.
    raw_measurements = model.sensor.measurements_now[observer_id]

    # Local track update pipeline:
    # 1) Split raw measurements into assigned (max 1 per slot) and orphan.
    # 2) Commit assigned measurements to local tracks.
    # 3) Attach orphan measurements to tentative tracks if possible, otherwise create new tentative tracks.

    # Step 1: split assigned vs orphan measurements via modified Hungarian algorithm.
    assigned_by_slot, orphan_measurements, conflict_count, skipped_count = _split_assigned_and_orphan_measurements(
        model,
        observer_id,
        raw_measurements,
        t,
        observer_pos
    )
    model.conflict_count_total[observer_id] += conflict_count
    model.skipped_assoc_total[observer_id] += skipped_count

    # Commit assigned measurements to local tracks.
    _commit_assigned_measurements_to_local_tracks!(
        model,
        observer_id,
        assigned_by_slot,
        t,
        observer_pos
    )

    # Create or update tentative tracks from orphan measurements. 
    # Here we use the angle constraint to check if we need to update existing track or init. a new one.
    _create_or_update_tentative_tracks_from_orphans!(
        model,
        observer_id,
        orphan_measurements,
        t,
        observer_pos
    )

    return nothing
end

Base.@kwdef mutable struct DistributedFusionModel
    nav::ObserverNavigationModel
    observer_idxs::Vector{Int}
    track_slots::Vector{Int}
    min_neighbor_count::Int
    sigma_theta_rad::Float64
    num_consensus_iter::Int
    process_q_diag::SVector{6, Float64}
    μ::Float64
    J2::Float64
    R_ref_m::Float64
    iod_estimate_state::Matrix{SVector{6, Float64}}
    iod_estimate_covariance::Matrix{Matrix{Float64}}
    iod_estimate_time_s::Matrix{Float64}
    iod_used_neighbors::Matrix{Int}
    iod_triangulation_ready::BitMatrix
    iod_initialized::BitMatrix
    state::Matrix{SVector{6, Float64}}
    covariance::Matrix{Matrix{Float64}}
    state_pred::Matrix{SVector{6, Float64}}
    covariance_pred::Matrix{Matrix{Float64}}
    filter_initialized::BitMatrix
    solo_target_measure_streak::Matrix{Int}
    last_update_t::Matrix{Float64}
    last_no_measure_warning_t::Vector{Float64}
    grouping_fallback_total::Int
    miss_multi_pass_total::Int
    miss_multi_pass_extra_total::Int
    mahal_multi_pass_total::Int
    mahal_multi_pass_extra_total::Int
    tt_attempt_total::Int
    tt_commit_total::Int
    tt_commit_correct_total::Int
    tt_commit_wrong_total::Int
    tt_commit_unknown_total::Int
    tt_skipped_total::Int
    tt_no_candidate_total::Int
    tt_gate_fail_total::Int
    tt_skipped_same_target_present_total::Int
    tt_skipped_no_same_target_total::Int
    tt_skipped_unknown_source_target_total::Int
    consensus_group_total::Int
    consensus_group_same_target_total::Int
    consensus_group_mixed_target_total::Int
    consensus_group_unknown_total::Int
end

function DistributedFusionModel(
    nav::ObserverNavigationModel,
    observer_idxs::Vector{Int},
    track_slots::Vector{Int},
    min_neighbor_count::Int,
    sigma_theta_rad::Float64,
    num_consensus_iter::Int,
    process_q_diag::SVector{6, Float64},
    μ::Float64,
    J2::Float64,
    R_ref_m::Float64,
    num_sats::Int
)
    nan_cov = () -> fill(NaN, 6, 6)
    return DistributedFusionModel(
        nav=nav,
        observer_idxs=observer_idxs,
        track_slots=track_slots,
        min_neighbor_count=min_neighbor_count,
        sigma_theta_rad=sigma_theta_rad,
        num_consensus_iter=num_consensus_iter,
        process_q_diag=process_q_diag,
        μ=μ,
        J2=J2,
        R_ref_m=R_ref_m,
        iod_estimate_state=[NAN_STATE6 for _ in 1:num_sats, _ in 1:num_sats],
        iod_estimate_covariance=[nan_cov() for _ in 1:num_sats, _ in 1:num_sats],
        iod_estimate_time_s=fill(NaN, num_sats, num_sats),
        iod_used_neighbors=zeros(Int, num_sats, num_sats),
        iod_triangulation_ready=falses(num_sats, num_sats),
        iod_initialized=falses(num_sats, num_sats),
        state=[NAN_STATE6 for _ in 1:num_sats, _ in 1:num_sats],
        covariance=[nan_cov() for _ in 1:num_sats, _ in 1:num_sats], # n x n matrix of 6x6 covariances (one per sat x slot)
        state_pred=[NAN_STATE6 for _ in 1:num_sats, _ in 1:num_sats],
        covariance_pred=[nan_cov() for _ in 1:num_sats, _ in 1:num_sats],
        filter_initialized=falses(num_sats, num_sats),
        solo_target_measure_streak=zeros(Int, num_sats, num_sats),
        last_update_t=fill(NaN, num_sats, num_sats),
        last_no_measure_warning_t=fill(NaN, num_sats),
        grouping_fallback_total=0,
        miss_multi_pass_total=0,
        miss_multi_pass_extra_total=0,
        mahal_multi_pass_total=0,
        mahal_multi_pass_extra_total=0,
        tt_attempt_total=0,
        tt_commit_total=0,
        tt_commit_correct_total=0,
        tt_commit_wrong_total=0,
        tt_commit_unknown_total=0,
        tt_skipped_total=0,
        tt_no_candidate_total=0,
        tt_gate_fail_total=0,
        tt_skipped_same_target_present_total=0,
        tt_skipped_no_same_target_total=0,
        tt_skipped_unknown_source_target_total=0,
        consensus_group_total=0,
        consensus_group_same_target_total=0,
        consensus_group_mixed_target_total=0,
        consensus_group_unknown_total=0
    )
end

@inline function _slot_capacity(model::DistributedFusionModel)::Int
    return size(model.filter_initialized, 2)
end

function _ensure_fusion_slot_capacity!(
    model::DistributedFusionModel,
    slot_id::Int
)
    slot_id <= _slot_capacity(model) && return nothing

    n_rows = size(model.filter_initialized, 1)
    old_cols = _slot_capacity(model)
    new_cols = slot_id

    new_iod_state = fill(NAN_STATE6, n_rows, new_cols)
    new_iod_state[:, 1:old_cols] = model.iod_estimate_state
    model.iod_estimate_state = new_iod_state

    new_iod_cov = [_nan_cov6() for _ in 1:n_rows, _ in 1:new_cols]
    for i in 1:n_rows
        for j in 1:old_cols
            new_iod_cov[i, j] = copy(model.iod_estimate_covariance[i, j])
        end
    end
    model.iod_estimate_covariance = new_iod_cov

    new_iod_time = fill(NaN, n_rows, new_cols)
    new_iod_time[:, 1:old_cols] = model.iod_estimate_time_s
    model.iod_estimate_time_s = new_iod_time

    new_iod_neighbors = zeros(Int, n_rows, new_cols)
    new_iod_neighbors[:, 1:old_cols] = model.iod_used_neighbors
    model.iod_used_neighbors = new_iod_neighbors

    new_iod_ready = falses(n_rows, new_cols)
    new_iod_ready[:, 1:old_cols] = model.iod_triangulation_ready
    model.iod_triangulation_ready = new_iod_ready

    new_iod_init = falses(n_rows, new_cols)
    new_iod_init[:, 1:old_cols] = model.iod_initialized
    model.iod_initialized = new_iod_init

    new_state = fill(NAN_STATE6, n_rows, new_cols)
    new_state[:, 1:old_cols] = model.state
    model.state = new_state

    new_cov = [_nan_cov6() for _ in 1:n_rows, _ in 1:new_cols]
    for i in 1:n_rows
        for j in 1:old_cols
            new_cov[i, j] = copy(model.covariance[i, j])
        end
    end
    model.covariance = new_cov

    new_state_pred = fill(NAN_STATE6, n_rows, new_cols)
    new_state_pred[:, 1:old_cols] = model.state_pred
    model.state_pred = new_state_pred

    new_cov_pred = [_nan_cov6() for _ in 1:n_rows, _ in 1:new_cols]
    for i in 1:n_rows
        for j in 1:old_cols
            new_cov_pred[i, j] = copy(model.covariance_pred[i, j])
        end
    end
    model.covariance_pred = new_cov_pred

    new_filter_init = falses(n_rows, new_cols)
    new_filter_init[:, 1:old_cols] = model.filter_initialized
    model.filter_initialized = new_filter_init

    new_solo_target_streak = zeros(Int, n_rows, new_cols)
    new_solo_target_streak[:, 1:old_cols] = model.solo_target_measure_streak
    model.solo_target_measure_streak = new_solo_target_streak

    new_last_update = fill(NaN, n_rows, new_cols)
    new_last_update[:, 1:old_cols] = model.last_update_t
    model.last_update_t = new_last_update

    old_warn_len = length(model.last_no_measure_warning_t)
    resize!(model.last_no_measure_warning_t, new_cols)
    for idx in (old_warn_len + 1):new_cols
        model.last_no_measure_warning_t[idx] = NaN
    end

    return nothing
end

@inline function _measurement_jacobian(x::SVector{6, Float64}, r_agent::SVector{3, Float64})
    ρ = x[1:3] - r_agent
    ρn = norm(ρ)
    if ρn <= 0.0
        return zeros(3, 6)
    end
    H_pos = Matrix(I, 3, 3) / ρn - (ρ * ρ') / ρn^3
    return hcat(H_pos, zeros(3, 3))
end

# Forward declarations avoid Julia 1.12 world-age depwarns when reloading with Revise.
function _gravity_accel_j2 end
function _gravity_pos_jacobian end

@inline function _kepler_dynamics(
    x::AbstractVector{<:Real},
    μ::Float64,
    J2::Float64,
    R_ref_m::Float64
)::SVector{6, Float64}
    r = SVector{3, Float64}(x[1], x[2], x[3])
    v = SVector{3, Float64}(x[4], x[5], x[6])
    rn = norm(r)
    if !(isfinite(rn) && rn > 1.0)
        return SVector{6, Float64}(v..., 0.0, 0.0, 0.0)
    end
    a = _gravity_accel_j2(r, μ, J2, R_ref_m)
    return SVector{6, Float64}(v..., a...)
end

@inline function _propagate_keplerian(
    x::AbstractVector{<:Real},
    μ::Float64,
    dt::Float64,
    J2::Float64,
    R_ref_m::Float64
)::SVector{6, Float64}
    x0 = SVector{6, Float64}(x)
    k1 = _kepler_dynamics(x0, μ, J2, R_ref_m)
    k2 = _kepler_dynamics(x0 + 0.5 * dt * k1, μ, J2, R_ref_m)
    k3 = _kepler_dynamics(x0 + 0.5 * dt * k2, μ, J2, R_ref_m)
    k4 = _kepler_dynamics(x0 + dt * k3, μ, J2, R_ref_m)
    return x0 + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
end

@inline function _gravity_accel_j2(
    r::AbstractVector{<:Real},
    μ::Float64,
    J2::Float64,
    R_ref_m::Float64
)::SVector{3, Float64}
    rv = SVector{3, Float64}(r[1], r[2], r[3])
    rn = norm(rv)
    if !(isfinite(rn) && rn > 1.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    a_central = -μ * rv / rn^3
    if abs(J2) <= 0.0 || R_ref_m <= 0.0
        return a_central
    end

    r2 = rn^2
    z2 = rv[3]^2
    j2_scale = 1.5 * J2 * μ * R_ref_m^2 / rn^5
    xy_factor = 5.0 * (z2 / r2) - 1.0
    z_factor = 5.0 * (z2 / r2) - 3.0
    a_j2 = SVector{3, Float64}(
        j2_scale * rv[1] * xy_factor,
        j2_scale * rv[2] * xy_factor,
        j2_scale * rv[3] * z_factor
    )
    return a_central + a_j2
end

@inline function _gravity_pos_jacobian(
    r::AbstractVector{<:Real},
    μ::Float64,
    J2::Float64,
    R_ref_m::Float64
)::Matrix{Float64}
    rv = SVector{3, Float64}(r[1], r[2], r[3])
    rn = norm(rv)
    if !(isfinite(rn) && rn > 1.0)
        return zeros(3, 3)
    end

    I3 = Matrix(I, 3, 3)
    Ω_central = -μ * (I3 / rn^3 - 3.0 * (rv * rv') / rn^5)
    if abs(J2) <= 0.0 || R_ref_m <= 0.0
        return Ω_central
    end

    x = rv[1]
    y = rv[2]
    z = rv[3]
    r2 = rn^2
    r5 = rn^5
    r7 = rn^7
    r9 = rn^9
    z2 = z^2
    z3 = z^3
    z4 = z2^2
    x2 = x^2
    y2 = y^2

    k = -1.5 * J2 * μ * R_ref_m^2

    dax_dx = k * (1.0 / r5 - 5.0 * (x2 + z2) / r7 + 35.0 * x2 * z2 / r9)
    dax_dy = k * (-5.0 * x * y / r7 + 35.0 * x * y * z2 / r9)
    dax_dz = k * (-15.0 * x * z / r7 + 35.0 * x * z3 / r9)

    day_dy = k * (1.0 / r5 - 5.0 * (y2 + z2) / r7 + 35.0 * y2 * z2 / r9)
    day_dz = k * (-15.0 * y * z / r7 + 35.0 * y * z3 / r9)

    daz_dz = k * (3.0 / r5 - 30.0 * z2 / r7 + 35.0 * z4 / r9)

    Ω_j2 = zeros(3, 3)
    Ω_j2[1, 1] = dax_dx
    Ω_j2[1, 2] = dax_dy
    Ω_j2[1, 3] = dax_dz
    Ω_j2[2, 1] = dax_dy
    Ω_j2[2, 2] = day_dy
    Ω_j2[2, 3] = day_dz
    Ω_j2[3, 1] = dax_dz
    Ω_j2[3, 2] = day_dz
    Ω_j2[3, 3] = daz_dz

    return Ω_central + Ω_j2
end

@inline function _compute_process_jacobian(
    x::AbstractVector{<:Real},
    μ::Float64,
    dt::Float64,
    J2::Float64,
    R_ref_m::Float64
)
    r = SVector{3, Float64}(x[1], x[2], x[3])
    if norm(r) <= 0.0
        return Matrix(I, 6, 6)
    end
    I3 = Matrix(I, 3, 3)
    Ω = _gravity_pos_jacobian(r, μ, J2, R_ref_m)
    F = zeros(6, 6)
    F[1:3, 1:3] .= I3
    F[1:3, 4:6] .= I3 * dt
    F[4:6, 1:3] .= Ω * dt
    F[4:6, 4:6] .= I3
    return F
end

function _low_pass_consensus(
    local_values::Dict{Int, Any},
    observer_idxs::Vector{Int},
    neighbor_map::Dict{Int, Vector{Int}},
    num_iter::Int
)
    consensus_values = Dict{Int, Any}(observer_id => deepcopy(local_values[observer_id]) for observer_id in observer_idxs)
    for _ in 1:num_iter
        new_values = Dict{Int, Any}()
        for observer_id in observer_idxs
            neighbors = get(neighbor_map, observer_id, Int[])
            degree_i = length(neighbors)
            weight_sum = 0.0
            sum_val = zero(consensus_values[observer_id])
            for nid in neighbors
                degree_j = length(get(neighbor_map, nid, Int[]))
                w_ij = 1.0 / (max(degree_i, degree_j) + 1.0)
                weight_sum += w_ij
                sum_val += w_ij * consensus_values[nid]
            end
            w_ii = 1.0 - weight_sum
            sum_val += w_ii * consensus_values[observer_id]
            new_values[observer_id] = sum_val
        end
        consensus_values = new_values
    end
    return consensus_values
end

@inline function _symmetrize_psd(P::Matrix{Float64})::Matrix{Float64}
    n = size(P, 1)
    fallback = Matrix(Diagonal(fill(1e6, n)))
    (size(P, 1) == size(P, 2) && all(isfinite, P)) || return fallback

    Ps = 0.5 * (P + P')
    all(isfinite, Ps) || return fallback

    eig = try
        eigen(Symmetric(Ps))
    catch
        return fallback
    end
    λ = eig.values
    V = eig.vectors
    (all(isfinite, λ) && all(isfinite, V)) || return fallback

    λ_clamped = max.(λ, 1e-12)
    P_psd = Matrix(V * Diagonal(λ_clamped) * V')
    all(isfinite, P_psd) || return fallback
    return 0.5 * (P_psd + P_psd')
end

function _initialize_group_filter_from_iod_consensus!(
    model::DistributedFusionModel,
    slot_match_by_observer::Dict{Int, Int},
    neighbor_map::Dict{Int, Vector{Int}},
    t::Float64,
    u
)
    active = Int[]
    for (observer_id, slot_id) in slot_match_by_observer
        if model.iod_initialized[observer_id, slot_id] &&
           isfinite(model.iod_estimate_time_s[observer_id, slot_id]) &&
           (abs(t - model.iod_estimate_time_s[observer_id, slot_id]) <= NAVIGATION_DT_TOL_SEC) &&
           _is_finite_state(model.iod_estimate_state[observer_id, slot_id]) &&
           _is_finite_cov(model.iod_estimate_covariance[observer_id, slot_id])
            push!(active, observer_id)
        end
    end
    required_group_size = max(2, model.min_neighbor_count + 1)
    has_full_group = length(active) >= required_group_size
    allow_singleton_bootstrap = DEKF_INIT_ALLOW_SINGLETON_IF_NO_GROUP && (length(active) == 1)
    (has_full_group || allow_singleton_bootstrap) || return nothing

    # Distributed IOD-consensus init in information form:
    # 1) each node builds Y_i = P_i^{-1}, y_i = Y_i*x_i
    # 2) low-pass consensus gives averages of Y_i and y_i
    # 3) each node reconstructs sum-terms via N*avg and recovers x,P locally
    active_set = Set(active)
    active_neighbor_map = Dict{Int, Vector{Int}}(
        observer_id => [n for n in get(neighbor_map, observer_id, Int[]) if n in active_set] for observer_id in active
    )
    Y_local = Dict{Int, Any}()
    y_local = Dict{Int, Any}()
    info_active = Int[]
    for observer_id in active
        slot_id = slot_match_by_observer[observer_id]
        x_i = model.iod_estimate_state[observer_id, slot_id]
        P_i = model.iod_estimate_covariance[observer_id, slot_id]
        iod_t_i = model.iod_estimate_time_s[observer_id, slot_id]
        (isfinite(iod_t_i) && abs(t - iod_t_i) <= NAVIGATION_DT_TOL_SEC) || continue
        (_is_finite_state(x_i) && _is_finite_cov(P_i)) || continue
        Y_i = pinv(P_i)
        all(isfinite, Y_i) || continue
        y_i = Y_i * Vector{Float64}(x_i)
        all(isfinite, y_i) || continue
        Y_local[observer_id] = Y_i
        y_local[observer_id] = y_i
        push!(info_active, observer_id)
    end
    isempty(info_active) && return nothing

    required_info_group_size = max(2, model.min_neighbor_count + 1)
    has_info_group = length(info_active) >= required_info_group_size
    allow_info_singleton = DEKF_INIT_ALLOW_SINGLETON_IF_NO_GROUP && (length(info_active) == 1)
    (has_info_group || allow_info_singleton) || return nothing

    info_set = Set(info_active)
    info_neighbor_map = Dict{Int, Vector{Int}}(
        observer_id => [n for n in get(active_neighbor_map, observer_id, Int[]) if n in info_set] for observer_id in info_active
    )
    Y_cons = _low_pass_consensus(Y_local, info_active, info_neighbor_map, model.num_consensus_iter)
    y_cons = _low_pass_consensus(y_local, info_active, info_neighbor_map, model.num_consensus_iter)
    N_participants = length(info_active)

    for observer_id in info_active
        slot_id = slot_match_by_observer[observer_id]
        # If filter already initialized for this slot, skip.
        model.filter_initialized[observer_id, slot_id] && continue
        Y_sum_hat = N_participants * Matrix(Y_cons[observer_id])
        y_sum_hat = N_participants * Vector{Float64}(y_cons[observer_id])
        (all(isfinite, Y_sum_hat) && all(isfinite, y_sum_hat)) || continue
        P_init = _symmetrize_psd(pinv(Y_sum_hat))
        all(isfinite, P_init) || continue
        x_init_vec = P_init * y_sum_hat
        all(isfinite, x_init_vec) || continue
        x_init = SVector{6, Float64}(x_init_vec)

        model.state[observer_id, slot_id] = x_init
        model.covariance[observer_id, slot_id] = P_init
        model.state_pred[observer_id, slot_id] = model.state[observer_id, slot_id]
        model.covariance_pred[observer_id, slot_id] = model.covariance[observer_id, slot_id]
        model.filter_initialized[observer_id, slot_id] = true
        model.solo_target_measure_streak[observer_id, slot_id] = 0
        model.last_update_t[observer_id, slot_id] = t
        tracks = model.nav.local_tracks[observer_id]
        track = tracks[slot_id]
        track.seed_ready = true
        track.status = :filter_initialized
        if track.plot_target_id == 0 && (track.last_meas !== nothing)
            # Safety: set the plot label once if this slot had no prior IOD label.
            track.plot_target_id = track.last_meas.target
        end
        track.consecutive_missed = 0
        # Keep local track estimate aligned to x^- for next-tick local gating.
        track.state_estimate_now = model.state_pred[observer_id, slot_id]
        track.covariance_estimate_now = copy(model.covariance_pred[observer_id, slot_id])

        println("DEKF init (IOD+consensus) | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s")
    end
    return nothing
end

function _retire_stale_filter_tracks!(model::DistributedFusionModel, t::Float64)
    cap = _slot_capacity(model)
    close_on_missed = TRACK_CLOSE_AFTER_MISSED_MEASUREMENTS > 0
    close_on_solo = ENABLE_TRACK_CLOSE_AFTER_SOLO_TARGET_MEASURE_STREAK &&
                    SOLO_TARGET_MEASURE_STREAK_CLOSE_STEPS > 0
    !(close_on_missed || close_on_solo) && return nothing

    solo_target_observers = Dict{Int, Set{Int}}()
    if close_on_solo
        for observer_id in model.observer_idxs
            tracks = model.nav.local_tracks[observer_id]
            for (slot_id, track) in tracks
                slot_id <= cap || continue
                model.filter_initialized[observer_id, slot_id] || continue
                track.has_measure_now || continue
                track.last_meas === nothing && continue
                target_id = track.last_meas.target
                if !haskey(solo_target_observers, target_id)
                    solo_target_observers[target_id] = Set{Int}()
                end
                push!(solo_target_observers[target_id], observer_id)
            end
        end
    end

    for observer_id in model.observer_idxs
        tracks = model.nav.local_tracks[observer_id]
        for (slot_id, track) in tracks
            slot_id <= cap || continue
            model.filter_initialized[observer_id, slot_id] || continue

            close_reason = nothing

            if close_on_missed
                if track.has_measure_now
                    track.consecutive_missed = 0
                else
                    track.consecutive_missed += 1
                    if track.consecutive_missed >= TRACK_CLOSE_AFTER_MISSED_MEASUREMENTS
                        close_reason = :missed
                    end
                end
            end

            if close_reason === nothing && close_on_solo
                if track.has_measure_now && (track.last_meas !== nothing)
                    target_id = track.last_meas.target
                    observers_with_fresh = get(solo_target_observers, target_id, Set{Int}())
                    is_only_observer = (length(observers_with_fresh) == 1) && (observer_id in observers_with_fresh)
                    if is_only_observer
                        model.solo_target_measure_streak[observer_id, slot_id] += 1
                        streak = model.solo_target_measure_streak[observer_id, slot_id]
                        if streak == 1
                            println(
                                "Solo-target-measure streak | target=$(target_id) | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s"
                            )
                        end
                        if streak > SOLO_TARGET_MEASURE_STREAK_CLOSE_STEPS
                            close_reason = :solo_target_measure_streak
                        end
                    else
                        model.solo_target_measure_streak[observer_id, slot_id] = 0
                    end
                else
                    model.solo_target_measure_streak[observer_id, slot_id] = 0
                end
            end

            if close_reason !== nothing
                close_streak = model.solo_target_measure_streak[observer_id, slot_id]
                model.filter_initialized[observer_id, slot_id] = false
                model.solo_target_measure_streak[observer_id, slot_id] = 0
                model.iod_initialized[observer_id, slot_id] = false
                model.iod_triangulation_ready[observer_id, slot_id] = false
                model.iod_used_neighbors[observer_id, slot_id] = 0
                track.seed_ready = false
                track.status = :closed
                if close_reason === :missed
                    println("Local track closed | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s | missed=$(track.consecutive_missed) | reason=missed_measurements_threshold")
                else
                    target_id = (track.last_meas === nothing) ? 0 : track.last_meas.target
                    println("Local track closed | target=$(target_id) | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s | solo_streak=$(close_streak) | reason=solo_target_measure_streak_threshold")
                end
            end
        end
    end

    return nothing
end

@inline function _resolve_initialized_slot_for_plot_target(
    model::DistributedFusionModel,
    obs::Int,
    plot_target_id::Int
)::Int
    tracks = model.nav.local_tracks[obs]
    best_slot = 0
    best_t = -Inf
    cap = _slot_capacity(model)
    for (slot_id, track) in tracks
        slot_id <= cap || continue
        model.filter_initialized[obs, slot_id] || continue
        track.plot_target_id == plot_target_id || continue
        t_last = model.last_update_t[obs, slot_id]
        if isfinite(t_last) && t_last > best_t
            best_t = t_last
            best_slot = slot_id
        end
    end
    return best_slot
end

@inline function _fusion_slot_for_save(
    model::DistributedFusionModel,
    obs::Int,
    slot_or_target::Int,
    t::Float64
)::Int
    cap = _slot_capacity(model)
    slot_eff = slot_or_target

    # Save path: caller passes target id; routing uses track.plot_target_id (plot-only label).
    if !(slot_eff <= cap && model.filter_initialized[obs, slot_eff])
        slot_eff = _resolve_initialized_slot_for_plot_target(model, obs, slot_or_target)
        slot_eff == 0 && return 0
    end

    t_last = model.last_update_t[obs, slot_eff]
    if !isfinite(t_last) || abs(t - t_last) > NAVIGATION_DT_TOL_SEC
        return 0
    end
    return slot_eff
end

@inline function _fusion_state_for_save(
    model::DistributedFusionModel,
    obs::Int,
    slot::Int,
    t::Float64
)::SVector{6, Float64}
    slot_eff = _fusion_slot_for_save(model, obs, slot, t)
    slot_eff == 0 && return NAN_STATE6
    return model.state[obs, slot_eff]
end

@inline function _build_local_information_terms(
    model::DistributedFusionModel,
    x_prior::SVector{6, Float64},
    observer_pos::SVector{3, Float64},
    los_meas::SVector{3, Float64}
)::Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Bool}
    # Build local information-form measurement terms:
    #   u = H'R^-1 z
    #   uhat = H'R^-1 h(x)
    #   U = H'R^-1 H
    H = _measurement_jacobian(x_prior, observer_pos)
    all(isfinite, los_meas) || return (zeros(6), zeros(6), zeros(6, 6), false)
    los_unit = _safe_unit(los_meas)

    I3 = Matrix(I, 3, 3)
    R_i = (model.sigma_theta_rad^2) * (I3 - los_unit * los_unit')
    R_inv = pinv(R_i)
    all(isfinite, R_inv) || return (zeros(6), zeros(6), zeros(6, 6), false)
    u_local = H' * (R_inv * los_meas)
    all(isfinite, u_local) || return (zeros(6), zeros(6), zeros(6, 6), false)

    los_pred = _safe_unit(x_prior[1:3] - observer_pos)
    uhat_local = H' * (R_inv * los_pred)
    U_local = H' * (R_inv * H)
    if !(all(isfinite, uhat_local) && all(isfinite, U_local))
        return (zeros(6), zeros(6), zeros(6, 6), false)
    end
    return (u_local, uhat_local, U_local, true)
end

function SimulationModel.calcNavigationEffect!(
    model::DistributedFusionModel,
    u,
    p,
    t::Float64,
    sat_idx::Int
)
    # Execute one distributed navigation fusion cycle per nav tick.
    # Run on the last observer callback so sensor/local-track effectors
    # have already populated this tick's measurements for all observers.
    sat_idx == model.observer_idxs[end] || return nothing
    # update current comm. network based on current inter-agent distance.
    _update_neighbors!(model.nav.comms, u, p)
    # check if a track should be retired due to missed measurements or solo-target-measure streak.
    _retire_stale_filter_tracks!(model, t)

    observer_set = Set(model.observer_idxs) 
    neighbor_map = Dict{Int, Vector{Int}}() 
    for sid in model.observer_idxs
        neighbor_map[sid] = [n for n in model.nav.comms.neighbors[sid] if n in observer_set]
    end

    Q = Diagonal(model.process_q_diag)

    # Snapshot state at tick start so all observers process the same time layer.
    x_prior_snap = copy(model.state_pred) 
    P_prior_snap = deepcopy(model.covariance_pred)
    x_snap_cols = size(x_prior_snap, 2)
    P_snap_cols = size(P_prior_snap, 2)

    # Buffers for updated estimates and consensus; keyed by (observer_id, slot_id).
    x_upd_buf = Dict{Tuple{Int, Int}, SVector{6, Float64}}() 
    P_upd_buf = Dict{Tuple{Int, Int}, Matrix{Float64}}()     
    x_pred_buf = Dict{Tuple{Int, Int}, SVector{6, Float64}}()
    P_pred_buf = Dict{Tuple{Int, Int}, Matrix{Float64}}()    
    touched = Set{Tuple{Int, Int}}()
    all_active_tracks = Tuple{Int, Int}[] # collection of all the active tracks from all the obs.

    for observer_id in model.observer_idxs
        # Active track is defined as any track that is not closed, and either has an initialized filter or can seed IOD.
        for slot_id in keys(model.nav.local_tracks[observer_id])
            track = model.nav.local_tracks[observer_id][slot_id]
            track.status == :closed && continue
            can_seed_iod = _track_has_two_measures(model.nav, observer_id, slot_id)
            is_initialized = (slot_id <= _slot_capacity(model)) &&
                             model.filter_initialized[observer_id, slot_id] &&
                             (track.status == :filter_initialized)
            (is_initialized || can_seed_iod) || continue
            _ensure_fusion_slot_capacity!(model, slot_id)
            push!(all_active_tracks, (observer_id, slot_id))
        end
    end

    # STEP 1/4: IOD attempt pass (observer-centric, all observers first).
    for (observer_id, slot_id) in all_active_tracks
        # If already init, skip this track's IOD attempt
        model.iod_initialized[observer_id, slot_id] && continue

        tracks = model.nav.local_tracks[observer_id]
        track = tracks[slot_id]
        # make sure track as at least two consecutive measures.
        can_seed_iod = _track_has_two_measures(model.nav, observer_id, slot_id)

        if !can_seed_iod
            model.iod_triangulation_ready[observer_id, slot_id] = false
            model.iod_used_neighbors[observer_id, slot_id] = 0
            track.seed_ready = false
            track.status = :tentative
            continue
        end

        track.seed_ready = true
        track.status = :seed_ready

        valid_neighbors, miss_pass_events, miss_pass_extra = _select_iod_neighbors(
            model.nav,
            observer_id,
            slot_id,
            get(neighbor_map, observer_id, Int[]),
            t,
            u
        )
        model.miss_multi_pass_total += miss_pass_events
        model.miss_multi_pass_extra_total += miss_pass_extra
        model.iod_used_neighbors[observer_id, slot_id] = length(valid_neighbors)
        if length(valid_neighbors) < model.min_neighbor_count
            model.iod_triangulation_ready[observer_id, slot_id] = false
            continue
        end

        # use valid neighbors + self for building the IOD equations.
        nodes_idx = vcat([(observer_id, slot_id)], valid_neighbors)
        nodes = NamedTuple[]
        A_rows = Matrix{Float64}(undef, 0, 3)
        b_rows = Float64[]
        dA_rows = Matrix{Float64}(undef, 0, 3)
        db_rows = Float64[]

        for (idx, matched_slot) in nodes_idx
            meas = _latest_track_measure(model.nav, idx, matched_slot, t)
            meas === nothing && continue

            tracks_i = model.nav.local_tracks[idx]
            tr = tracks_i[matched_slot]
            tr.prev_meas === nothing && continue

            lrate = _track_los_rate(model.nav, idx, matched_slot)
            lrate === nothing && continue

            # use perturbed obs states for IOD.
            l_unit = _safe_unit(meas.los_unit)
            r = _known_observer_pos!(model.nav, idx, u, t)
            v = _known_observer_vel!(model.nav, idx, u, t)
            eq = _build_agent_equations(r, v, l_unit, lrate)
            eq === nothing && continue

            A_rows = vcat(A_rows, eq.A)
            append!(b_rows, eq.b)
            dA_rows = vcat(dA_rows, eq.dA)
            append!(db_rows, eq.db)
            push!(nodes, (r=r, v=v, los=l_unit, los_rate=lrate))
        end

        length(nodes) >= 2 || continue
        H = [A_rows zeros(size(A_rows, 1), 3); dA_rows A_rows]
        y = vcat(b_rows, db_rows)
        size(H, 1) >= 6 || continue

        z_est = pinv(H) * y
        all(isfinite, z_est) || continue

        x_est = SVector{3, Float64}(z_est[1], z_est[2], z_est[3])
        x_dot_est = SVector{3, Float64}(z_est[4], z_est[5], z_est[6])
        cov = _compute_state_covariance(nodes, x_est, x_dot_est, model.sigma_theta_rad, NAVIGATION_RATE_SEC)
        (cov !== nothing && all(isfinite, cov)) || continue


        # update tracks with IOD results.
        model.iod_estimate_state[observer_id, slot_id] = SVector{6, Float64}(z_est)
        model.iod_estimate_covariance[observer_id, slot_id] = copy(cov)
        model.iod_estimate_time_s[observer_id, slot_id] = t
        model.iod_triangulation_ready[observer_id, slot_id] = true
        model.iod_initialized[observer_id, slot_id] = true

        track.status = :iod_initialized
        if track.plot_target_id == 0 && (track.last_meas !== nothing)
            # Plot-only stable label: freeze at first successful IOD init.
            track.plot_target_id = track.last_meas.target
        end
        track.state_estimate_now = model.iod_estimate_state[observer_id, slot_id]
        track.covariance_estimate_now = copy(model.iod_estimate_covariance[observer_id, slot_id])

        println("IOD init | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s")
    end

    # STEP 2/4 + 3/4: build match groups and run one consensus pass per group.
    # match group is defined as a set of tracks across observers that are believed to correspond to the same target.
    match_groups, grouping_fallback_count, mahal_pass_events, mahal_pass_extra = _build_match_groups(model, all_active_tracks, neighbor_map)
    model.grouping_fallback_total += grouping_fallback_count
    model.mahal_multi_pass_total += mahal_pass_events
    model.mahal_multi_pass_extra_total += mahal_pass_extra
    
    for group in match_groups

        group_pairs = group.selected_tracks # Vector of (observer_id, slot_id) pairs that are in this consensus group
        isempty(group_pairs) && continue

        # Initialize filters for this whole match group from IOD consensus if not already initialized.
        group_slot_match = Dict{Int, Int}(observer_id => slot_id for (observer_id, slot_id) in group_pairs)
        _initialize_group_filter_from_iod_consensus!(model, group_slot_match, neighbor_map, t, u)

        u_local = Dict{Int, Any}()
        uhat_local = Dict{Int, Any}()
        U_local = Dict{Int, Any}()
        x_prior_by_observer = Dict{Int, SVector{6, Float64}}()
        P_prior_by_observer = Dict{Int, Matrix{Float64}}()

        # Build local information terms once per (observer, local-track) group member.
        for (observer_id, slot_id) in group_pairs
            u_local[observer_id] = zeros(6)
            uhat_local[observer_id] = zeros(6)
            U_local[observer_id] = zeros(6, 6)

            # skip if filter not init.
            if !model.filter_initialized[observer_id, slot_id]
                continue
            end

            # get prior from snapshot to ensure same time layer for all observers in this consensus group.
            if slot_id <= x_snap_cols
                x_prior = x_prior_snap[observer_id, slot_id]
            else
                x_prior = model.state_pred[observer_id, slot_id]
            end
            if slot_id <= P_snap_cols
                P_prior = P_prior_snap[observer_id, slot_id]
            else
                P_prior = model.covariance_pred[observer_id, slot_id]
            end
            if !(_is_finite_state(x_prior) && _is_finite_cov(P_prior))
                x_prior = model.state[observer_id, slot_id]
                P_prior = model.covariance[observer_id, slot_id]
            end
            if !(_is_finite_state(x_prior) && _is_finite_cov(P_prior))
                continue
            end
            x_prior_by_observer[observer_id] = x_prior
            P_prior_by_observer[observer_id] = P_prior

            # compute local information terms for this observer-track using the prior and the latest measurement.
            observer_pos = _known_observer_pos!(model.nav, observer_id, u, t)
            measurement = _latest_track_measure(model.nav, observer_id, slot_id, t)
            measurement === nothing && continue

            u_i, uhat_i, U_i, _ = _build_local_information_terms(
                model,
                x_prior,
                observer_pos,
                measurement.los_unit
            )
            u_local[observer_id] = u_i
            uhat_local[observer_id] = uhat_i
            U_local[observer_id] = U_i
        end

        # Identify consensus participants as all initialized tracks in the group.
        # Observers without a fresh measurement already carry null local info terms
        # (u=0, uhat=0, U=0), so they participate in the same average consensus.
        consensus_participants = [
            observer_id for (observer_id, slot_id) in group_pairs
            if model.filter_initialized[observer_id, slot_id]
        ]
        q_scale = (DEKF_CONSENSUS_UPDATE_MODE == :information) ? max(1, length(consensus_participants)) : 1
        participant_set = Set(consensus_participants)
        
        # Run consensus whenever we have at least one initialized participant.
        has_consensus = !isempty(consensus_participants)
        z_cons = Dict{Int, Any}()
        zhat_cons = Dict{Int, Any}()
        S_cons = Dict{Int, Any}()
        V_cons = Dict{Int, Any}()
        v_cons = Dict{Int, Any}()
        icf_participant_set = Set{Int}()
        icf_participant_count = 0
        kcf_sum_i = Dict{Int, Vector{Float64}}()
        kcf_sum_U = Dict{Int, Matrix{Float64}}()
        kcf_delta_x = Dict{Int, Vector{Float64}}()
        kcf_participant_set = Set{Int}()

        # Run consensus on the selected update mode.
        if has_consensus
            if DEKF_CONSENSUS_UPDATE_MODE == :information
                participant_neighbor_map = Dict{Int, Vector{Int}}(
                    observer_id => [n for n in get(neighbor_map, observer_id, Int[]) if n in participant_set]
                    for observer_id in consensus_participants
                )
                u_part = Dict{Int, Any}(observer_id => u_local[observer_id] for observer_id in consensus_participants)
                uhat_part = Dict{Int, Any}(observer_id => uhat_local[observer_id] for observer_id in consensus_participants)
                U_part = Dict{Int, Any}(observer_id => U_local[observer_id] for observer_id in consensus_participants)
                z_cons = _low_pass_consensus(u_part, consensus_participants, participant_neighbor_map, model.num_consensus_iter)
                zhat_cons = _low_pass_consensus(uhat_part, consensus_participants, participant_neighbor_map, model.num_consensus_iter)
                S_cons = _low_pass_consensus(U_part, consensus_participants, participant_neighbor_map, model.num_consensus_iter)
            elseif DEKF_CONSENSUS_UPDATE_MODE == :icf
                Y_prior_local = Dict{Int, Matrix{Float64}}()
                y_prior_local = Dict{Int, Vector{Float64}}()
                i_local = Dict{Int, Vector{Float64}}()

                icf_participants = Int[]
                for observer_id in consensus_participants
                    haskey(x_prior_by_observer, observer_id) || continue
                    haskey(P_prior_by_observer, observer_id) || continue
                    x_prior = x_prior_by_observer[observer_id]
                    P_prior = P_prior_by_observer[observer_id]
                    Y_prior = pinv(P_prior)
                    all(isfinite, Y_prior) || continue
                    y_prior = Vector{Float64}(Y_prior * x_prior)
                    i_term = Vector{Float64}((u_local[observer_id] - uhat_local[observer_id]) + U_local[observer_id] * x_prior)
                    if !(all(isfinite, y_prior) && all(isfinite, i_term))
                        continue
                    end
                    Y_prior_local[observer_id] = Y_prior
                    y_prior_local[observer_id] = y_prior
                    i_local[observer_id] = i_term
                    push!(icf_participants, observer_id)
                end

                icf_participant_count = length(icf_participants)
                if icf_participant_count > 0
                    icf_participant_set = Set(icf_participants)
                    icf_neighbor_map = Dict{Int, Vector{Int}}(
                        observer_id => [n for n in get(neighbor_map, observer_id, Int[]) if n in icf_participant_set]
                        for observer_id in icf_participants
                    )

                    V_part = Dict{Int, Any}()
                    v_part = Dict{Int, Any}()
                    for observer_id in icf_participants
                        V_part[observer_id] = (Y_prior_local[observer_id] / icf_participant_count) + U_local[observer_id]
                        v_part[observer_id] = (y_prior_local[observer_id] / icf_participant_count) + i_local[observer_id]
                    end

                    V_cons = _low_pass_consensus(V_part, icf_participants, icf_neighbor_map, model.num_consensus_iter)
                    v_cons = _low_pass_consensus(v_part, icf_participants, icf_neighbor_map, model.num_consensus_iter)
                end
            elseif DEKF_CONSENSUS_UPDATE_MODE == :kcf
                kcf_participants = [
                    observer_id for observer_id in consensus_participants
                    if haskey(x_prior_by_observer, observer_id) && haskey(P_prior_by_observer, observer_id)
                ]
                kcf_participant_set = Set(kcf_participants)
                for observer_id in kcf_participants
                    neighbors = [
                        n for n in get(neighbor_map, observer_id, Int[])
                        if n in kcf_participant_set
                    ]
                    support_nodes = vcat([observer_id], neighbors)

                    i_sum = zeros(6)
                    S_sum = zeros(6, 6)
                    for nid in support_nodes
                        x_nid = x_prior_by_observer[nid]
                        i_sum += Vector{Float64}((u_local[nid] - uhat_local[nid]) + U_local[nid] * x_nid)
                        S_sum += Matrix{Float64}(U_local[nid])
                    end

                    delta_x = zeros(6)
                    x_i = x_prior_by_observer[observer_id]
                    for nid in neighbors
                        delta_x += Vector{Float64}(x_prior_by_observer[nid] - x_i)
                    end

                    kcf_sum_i[observer_id] = i_sum
                    kcf_sum_U[observer_id] = S_sum
                    kcf_delta_x[observer_id] = delta_x
                end
            end

            for (_, slot_id) in group_pairs
                model.last_no_measure_warning_t[slot_id] = NaN
            end
        end

        # Apply correction/prediction for every initialized observer/local-track in this matched group.
        kcf_local_pairs = Tuple{Int, Int}[]
        kcf_x_upd_local = Dict{Tuple{Int, Int}, SVector{6, Float64}}()
        kcf_P_upd_local = Dict{Tuple{Int, Int}, Matrix{Float64}}()
        for (observer_id, slot_id) in group_pairs
            model.filter_initialized[observer_id, slot_id] || continue
            key = (observer_id, slot_id)
            if slot_id <= x_snap_cols
                x_prior = x_prior_snap[observer_id, slot_id]
            else
                x_prior = model.state_pred[observer_id, slot_id]
            end
            if slot_id <= P_snap_cols
                P_prior = P_prior_snap[observer_id, slot_id]
            else
                P_prior = model.covariance_pred[observer_id, slot_id]
            end
            if !(_is_finite_state(x_prior) && _is_finite_cov(P_prior))
                x_prior = model.state[observer_id, slot_id]
                P_prior = model.covariance[observer_id, slot_id]
            end
            if !(_is_finite_state(x_prior) && _is_finite_cov(P_prior))
                continue
            end

            x_upd = x_prior
            P_upd = P_prior

            if DEKF_CONSENSUS_UPDATE_MODE == :information
                # Information-domain update.
                if has_consensus && (observer_id in participant_set)
                    info_mat = pinv(P_prior) + S_cons[observer_id]
                    if all(isfinite, info_mat) && all(isfinite, z_cons[observer_id]) && all(isfinite, zhat_cons[observer_id])
                        P_upd_try = _symmetrize_psd(pinv(info_mat))
                        if all(isfinite, P_upd_try)
                            x_upd_try = x_prior + P_upd_try * (z_cons[observer_id] - zhat_cons[observer_id])
                            if all(isfinite, x_upd_try)
                                x_upd = SVector{6, Float64}(x_upd_try)
                                P_upd = P_upd_try
                            end
                        end
                    end
                end
            elseif DEKF_CONSENSUS_UPDATE_MODE == :icf
                # ICF update with consensus on V and v:
                # V_j = Y_j/N + I_j, v_j = y_j/N + i_j.
                # Recover Y and y as N*V_cons and N*v_cons before solving x = Y^-1 y.
                if has_consensus && (observer_id in icf_participant_set) && (icf_participant_count > 0)
                    Y_upd = icf_participant_count * V_cons[observer_id]
                    y_upd = icf_participant_count * v_cons[observer_id]
                    if all(isfinite, Y_upd) && all(isfinite, y_upd)
                        P_upd_try = _symmetrize_psd(pinv(Y_upd))
                        if all(isfinite, P_upd_try)
                            x_upd_try = P_upd_try * y_upd
                            if all(isfinite, x_upd_try)
                                x_upd = SVector{6, Float64}(x_upd_try)
                                P_upd = P_upd_try
                            end
                        end
                    end
                end
            elseif DEKF_CONSENSUS_UPDATE_MODE == :kcf
                # One-step KCF update:
                # i_i = sum_{j in N_i U {i}} i_j
                # S_i = sum_{j in N_i U {i}} U_j
                # x_i^+ = x_i^- + (Y_i^- + S_i)^-1 (i_i - S_i x_i^-)
                #         + gamma_i * sum_{j in N_i}(x_j^- - x_i^-)
                if has_consensus && (observer_id in kcf_participant_set) &&
                   haskey(kcf_sum_i, observer_id) && haskey(kcf_sum_U, observer_id) && haskey(kcf_delta_x, observer_id)
                    Y_prior = pinv(P_prior)
                    S_i = kcf_sum_U[observer_id]
                    i_i = kcf_sum_i[observer_id]
                    info_mat = Y_prior + S_i
                    if all(isfinite, info_mat) && all(isfinite, i_i)
                        P_upd_try = _symmetrize_psd(pinv(info_mat))
                        if all(isfinite, P_upd_try)
                            info_term = P_upd_try * (i_i - S_i * x_prior)
                            fro_sq = tr(P_prior' * P_prior)
                            gamma = DEKF_KCF_EPSILON / (sqrt(max(fro_sq, 0.0)) + 1.0)
                            cons_term = gamma * (P_prior * kcf_delta_x[observer_id])
                            x_upd_try = x_prior + info_term + cons_term
                            if all(isfinite, x_upd_try)
                                x_upd = SVector{6, Float64}(x_upd_try)
                                P_upd = P_upd_try
                            end
                        end
                    end
                end
                # run average consensus after the update step.
                # Store local post-update terms; consensus/prediction are done after
                # all group members are processed.
                kcf_x_upd_local[key] = x_upd
                kcf_P_upd_local[key] = P_upd
                push!(kcf_local_pairs, key)
                continue
            end

            x_pred_try = _propagate_keplerian(x_upd, model.μ, NAVIGATION_RATE_SEC, model.J2, model.R_ref_m)
            x_pred = _is_finite_state(x_pred_try) ? x_pred_try : x_upd
            F = _compute_process_jacobian(x_upd, model.μ, NAVIGATION_RATE_SEC, model.J2, model.R_ref_m)
            all(isfinite, F) || (F = Matrix(I, 6, 6))
            P_pred_raw = F * P_upd * F' + q_scale * Matrix(Q)
            all(isfinite, P_pred_raw) || (P_pred_raw = P_upd + q_scale * Matrix(Q))
            P_pred = _symmetrize_psd(P_pred_raw)

            x_upd_buf[key] = x_upd
            P_upd_buf[key] = P_upd
            x_pred_buf[key] = x_pred
            P_pred_buf[key] = P_pred
            push!(touched, key)
        end

        if DEKF_CONSENSUS_UPDATE_MODE == :kcf && !isempty(kcf_local_pairs)
            # Post-update average consensus on information states (Y, y) as in KCF.
            kcf_obs = Int[]
            obs_to_key = Dict{Int, Tuple{Int, Int}}()
            Y_post_local = Dict{Int, Any}()
            y_post_local = Dict{Int, Any}()

            for key in kcf_local_pairs
                observer_id = key[1]
                P_loc = kcf_P_upd_local[key]
                x_loc = kcf_x_upd_local[key]
                Y_loc = pinv(P_loc)
                y_loc = Vector{Float64}(Y_loc * x_loc)
                if all(isfinite, Y_loc) && all(isfinite, y_loc)
                    push!(kcf_obs, observer_id)
                    obs_to_key[observer_id] = key
                    Y_post_local[observer_id] = Y_loc
                    y_post_local[observer_id] = y_loc
                end
            end

            kcf_obs = sort(unique(kcf_obs))
            Y_post_cons = Dict{Int, Any}()
            y_post_cons = Dict{Int, Any}()
            if !isempty(kcf_obs)
                kcf_post_neighbor_map = Dict{Int, Vector{Int}}(
                    observer_id => [n for n in get(neighbor_map, observer_id, Int[]) if n in kcf_obs]
                    for observer_id in kcf_obs
                )
                Y_post_cons = _low_pass_consensus(Y_post_local, kcf_obs, kcf_post_neighbor_map, model.num_consensus_iter)
                y_post_cons = _low_pass_consensus(y_post_local, kcf_obs, kcf_post_neighbor_map, model.num_consensus_iter)
            end

            for key in kcf_local_pairs
                observer_id = key[1]
                x_upd = kcf_x_upd_local[key]
                P_upd = kcf_P_upd_local[key]

                if haskey(obs_to_key, observer_id) && haskey(Y_post_cons, observer_id) && haskey(y_post_cons, observer_id)
                    Yc = Y_post_cons[observer_id]
                    yc = y_post_cons[observer_id]
                    if all(isfinite, Yc) && all(isfinite, yc)
                        P_try = _symmetrize_psd(pinv(Yc))
                        if all(isfinite, P_try)
                            x_try = P_try * yc
                            if all(isfinite, x_try)
                                P_upd = P_try
                                x_upd = SVector{6, Float64}(x_try)
                            end
                        end
                    end
                end

                x_pred_try = _propagate_keplerian(x_upd, model.μ, NAVIGATION_RATE_SEC, model.J2, model.R_ref_m)
                x_pred = _is_finite_state(x_pred_try) ? x_pred_try : x_upd
                F = _compute_process_jacobian(x_upd, model.μ, NAVIGATION_RATE_SEC, model.J2, model.R_ref_m)
                all(isfinite, F) || (F = Matrix(I, 6, 6))
                P_pred_raw = F * P_upd * F' + q_scale * Matrix(Q)
                all(isfinite, P_pred_raw) || (P_pred_raw = P_upd + q_scale * Matrix(Q))
                P_pred = _symmetrize_psd(P_pred_raw)

                x_upd_buf[key] = x_upd
                P_upd_buf[key] = P_upd
                x_pred_buf[key] = x_pred
                P_pred_buf[key] = P_pred
                push!(touched, key)
            end
        end
    end

    # STEP 4/4: Global commit after full observer-centric pass.
    # This preserves synchronous update semantics across observers.
    for key in touched
        sid, slot = key
        model.state[sid, slot] = x_upd_buf[key]
        model.covariance[sid, slot] = P_upd_buf[key]
        model.state_pred[sid, slot] = x_pred_buf[key]
        model.covariance_pred[sid, slot] = P_pred_buf[key]
        model.last_update_t[sid, slot] = t
        tracks = model.nav.local_tracks[sid]
        track = tracks[slot]
        # Keep local track estimate aligned to x^- for next-tick local gating.
        track.state_estimate_now = model.state_pred[sid, slot]
        track.covariance_estimate_now = copy(model.covariance_pred[sid, slot])
        track.status = :filter_initialized
    end
    return nothing
end

@inline function mean_to_true_anomaly_rad(M_rad::Float64, e::Float64; tol::Float64=1e-13, max_iter::Int=30)
    M_norm = mod(M_rad + pi, 2pi) - pi
    E = e < 0.8 ? M_norm : pi
    for _ in 1:max_iter
        f = E - e * sin(E) - M_norm
        fp = 1.0 - e * cos(E)
        Δ = f / fp
        E -= Δ
        if abs(Δ) < tol
            break
        end
    end
    ν = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(E / 2.0))
    return mod(ν, 2pi)
end

@inline mean_to_true_anomaly_deg(M_deg::Float64, e::Float64) = rad2deg(mean_to_true_anomaly_rad(deg2rad(M_deg), e))

@inline function eclipse_factor_cylindrical(r_sc::SVector{3, Float64}, r_sun::SVector{3, Float64}, rp::Float64)::Float64
    if dot(r_sc, r_sun) >= 0.0
        return 1.0
    end
    d = norm(cross(r_sc, r_sun)) / norm(r_sun)
    return d <= rp ? 0.0 : 1.0
end

function moon_pos_analytic_j2000(et::Float64)::SVector{3, Float64}
    d = et / 86_400.0
    Ω = deg2rad(mod(125.1228 - 0.0529538083 * d, 360.0))
    i = deg2rad(5.1454)
    ω = deg2rad(mod(318.0634 + 0.1643573223 * d, 360.0))
    a = 384_400_000.0
    e = 0.0549
    M = deg2rad(mod(115.3654 + 13.0649929509 * d, 360.0))

    E = M
    for _ in 1:20
        Δ = (E - e * sin(E) - M) / (1.0 - e * cos(E))
        E -= Δ
        abs(Δ) < 1e-13 && break
    end

    ν = 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(E / 2.0))
    r = a * (1.0 - e * cos(E))
    u = ω + ν

    x_ecl = r * (cos(Ω) * cos(u) - sin(Ω) * sin(u) * cos(i))
    y_ecl = r * (sin(Ω) * cos(u) + cos(Ω) * sin(u) * cos(i))
    z_ecl = r * (sin(u) * sin(i))

    ϵ = deg2rad(23.439291)
    x_eq = x_ecl
    y_eq = y_ecl * cos(ϵ) - z_ecl * sin(ϵ)
    z_eq = y_ecl * sin(ϵ) + z_ecl * cos(ϵ)
    return SVector{3, Float64}(x_eq, y_eq, z_eq)
end

const _moon_fallback_warned = Ref(false)

function body_pos_wrt_earth_pci(body_name::String, et::Float64, planet::AbstractPlanet)::SVector{3, Float64}
    name = uppercase(body_name)
    if name == "MOON"
        try
            r_moon_j2000 = lock(RuntimeServices.SPICE_LOCK) do
                SVector{3, Float64}(spkpos("MOON", et, "J2000", "NONE", "EARTH")[1]) * 1e3
            end
            return r_moon_j2000
        catch
            if !_moon_fallback_warned[]
                _moon_fallback_warned[] = true
                @warn "MOON SPICE ephemeris is unavailable; using analytic lunar orbit fallback for third-body acceleration."
            end
            return moon_pos_analytic_j2000(et)
        end
    end

    r_body_j2000 = lock(RuntimeServices.SPICE_LOCK) do
        SVector{3, Float64}(spkpos(name, et, "J2000", "NONE", "EARTH")[1]) * 1e3
    end
    return r_body_j2000
end

struct SunMoonThirdBodyModel{P <: AbstractPlanet} <: AbstractForceTorqueModel
    planet::P
    body_names::NTuple{2, String}
    body_mus::NTuple{2, Float64}
end

SunMoonThirdBodyModel(planet::P) where {P <: AbstractPlanet} = SunMoonThirdBodyModel(
    planet,
    ("SUN", "MOON"),
    (1.3271244002331e20, 4.9028005821478e12)
)

function SimulationModel.calcForceTorque(
    model::SunMoonThirdBodyModel,
    x,
    p::SimulationModel.ODEParams,
    i::Int64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_sc = SVector{3, Float64}(x.pos)
    mass_sc = x.mass
    et = p.shared_buffers.et_start[] + p.shared_buffers.current_time[]
    a_tb = MVector{3, Float64}(0.0, 0.0, 0.0)

    @inbounds for k in eachindex(model.body_names)
        body_name = model.body_names[k]
        μ_k = model.body_mus[k]
        r_primary_k = body_pos_wrt_earth_pci(body_name, et, model.planet)
        r_sc_k = r_primary_k - pos_sc
        a_tb .+= μ_k * (r_sc_k / norm(r_sc_k)^3 - r_primary_k / norm(r_primary_k)^3)
    end

    return mass_sc * SVector{3, Float64}(a_tb), SVector{3, Float64}(0.0, 0.0, 0.0)
end

struct CannonballSRPModel{P <: AbstractPlanet} <: AbstractForceTorqueModel
    planet::P
    areas_m2::Vector{Float64}
    cr::Float64
    p_ref::Float64
    au_m::Float64
end

function SimulationModel.calcForceTorque(
    model::CannonballSRPModel,
    x,
    p::SimulationModel.ODEParams,
    i::Int64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_sc = SVector{3, Float64}(x.pos)
    et = p.shared_buffers.et_start[] + p.shared_buffers.current_time[]
    r_sun = body_pos_wrt_earth_pci("SUN", et, model.planet)
    r_sc_to_sun = r_sun - pos_sc
    d_sc_sun = norm(r_sc_to_sun)
    u_sun_to_sc = -r_sc_to_sun / d_sc_sun

    eclipse = eclipse_factor_cylindrical(pos_sc, r_sun, model.planet.Rp_e)
    p_srp = model.p_ref * (model.au_m / d_sc_sun)^2
    area = model.areas_m2[min(i, length(model.areas_m2))]
    force = p_srp * model.cr * area * eclipse * u_sun_to_sc
    return force, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function make_translational_spacecraft(
    id::Int64,
    a_m::Float64,
    e::Float64,
    i_deg::Float64,
    raan_deg::Float64,
    aop_deg::Float64,
    mean_anomaly_deg::Float64;
    mass_kg::Float64=220.0,
    area_m2::Float64=2.0
)::SpacecraftModel
    root = Link{0}(root=true, m=mass_kg, ref_area=area_m2, dims=MVector{3, Float64}(1.0, 1.0, 1.0))
    ν_deg = mean_to_true_anomaly_deg(mean_anomaly_deg, e)
    ic = InitialCondition(a_m, e, i_deg, aop_deg, raan_deg, ν_deg)
    return SpacecraftModel(
        Joint[],
        [root],
        root,
        true,
        root.m,
        0.0,
        root.inertia,
        0,
        0,
        ic,
        id
    )
end


function generate_observer_centered_debris_cluster(
    n_targets::Int,
    observer_a_m::Float64,
    observer_i_deg::Float64,
    observer_mean_anomalies_deg;
    observer_e::Float64=1e-4,
    observer_raan_deg::Float64=10.0,
    observer_aop_deg::Float64=14.0
)
    rng = MersenneTwister()
    # Center the target mean anomalies around the observer mean anomalies.
    m_center = sum(observer_mean_anomalies_deg) / length(observer_mean_anomalies_deg)
    target_defs = NamedTuple{(:a_m, :e, :i_deg, :raan_deg, :aop_deg, :M_deg), Tuple{Float64, Float64, Float64, Float64, Float64, Float64}}[]

    for _ in 1:n_targets
        # Semi-major axis
        a_offset_m = 200_000.0 * randn(rng)

        # Inclination
        i_offset_deg = 4.0 * randn(rng)

        # Eccentricity
        e_val = clamp(abs(observer_e + 0.0015 * randn(rng) + 0.0010 * abs(randn(rng))), 1e-5, 0.02)

        # RAAN
        raan_offset_deg = 8.0 * randn(rng)
        raan_deg = mod(observer_raan_deg + raan_offset_deg, 360.0)

        # Argument of perigee
        aop_offset_deg = 40.0 * randn(rng)
        aop_deg = mod(observer_aop_deg + aop_offset_deg, 360.0)

        # Mean anomaly
        m_offset =   1.5 * randn(rng)

        push!(target_defs, (
            a_m   = observer_a_m + a_offset_m,
            e = e_val,
            i_deg = observer_i_deg + i_offset_deg,
            raan_deg = raan_deg,
            aop_deg = aop_deg,
            M_deg = m_center + m_offset
        ))
    end

    return Tuple(target_defs)
end

@inline function _wrap_deg(x::Float64)::Float64
    y = mod(x, 360.0)
    return y < 0.0 ? y + 360.0 : y
end


planet = Earth("", SPICE_PATH)
isfile(EARTH_HARMONICS_FILE) || error("Earth harmonics file not found: $EARTH_HARMONICS_FILE")

e = 1e-4
raan_deg = 10.0
aop_deg = 14.0

observer_a_m = 6_963.0e3
target_a_m = 6_949.0e3

observer_i_deg = 85.0
target_primary_i_deg = 86.0
target_secondary_i_deg = 2.0 * observer_i_deg - target_primary_i_deg

observer_mean_anomalies_deg = (288.0, 290.0, 292.0)
target_primary_mean_anomaly_deg = 288.0
target_secondary_mean_anomaly_deg = 288.0


# # Challenging for meas association.
# observer_defs = (
#     (a_m=observer_a_m, i_deg=observer_i_deg, M_deg=observer_mean_anomalies_deg[1]),
#     (a_m=observer_a_m, i_deg=observer_i_deg, M_deg=observer_mean_anomalies_deg[2]),
#     (a_m=observer_a_m, i_deg=observer_i_deg, M_deg=observer_mean_anomalies_deg[3])
# )
# target_defs = (
#     (a_m=target_a_m, i_deg=target_primary_i_deg, M_deg=target_primary_mean_anomaly_deg),
#     # symmetric inclination with respect to observer orbit
#     (a_m=target_a_m, i_deg=target_secondary_i_deg, M_deg=target_secondary_mean_anomaly_deg),
#     (a_m=target_a_m, i_deg=85.5, M_deg=target_secondary_mean_anomaly_deg),
#     (a_m=target_a_m, i_deg=84.5, M_deg=target_secondary_mean_anomaly_deg),
#     # Late-entry target: starts outside observer range and typically enters
#     # visibility later (often first for one observer) with a shorter pass.
#     (a_m=target_a_m + 250.0e3, i_deg=72.0, M_deg=302.0),
#     # Companion challenging target with similar visibility window but
#     # genuinely different geometry (different altitude and inclination).
#     (a_m=target_a_m + 290.0e3, i_deg=68.0, M_deg=303.0)
# )

# # Create spacecraft models for all observers and targets.
# spacecraft = SpacecraftModel[]
# for (sat_id, obs) in enumerate(observer_defs)
#     push!(
#         spacecraft,
#         make_translational_spacecraft(
#             sat_id,
#             obs.a_m,
#             e,
#             obs.i_deg,
#             raan_deg,
#             aop_deg,
#             obs.M_deg
#         )
#     )
# end
# first_target_id = length(observer_defs) + 1
# for (k, tgt) in enumerate(target_defs)
#     sat_id = first_target_id + k - 1
#     push!(
#         spacecraft,
#         make_translational_spacecraft(
#             sat_id,
#             tgt.a_m,
#             e,
#             tgt.i_deg,
#             raan_deg,
#             aop_deg,
#             tgt.M_deg
#         )
#     )
# end

# observer_idxs = collect(1:length(observer_defs))
# target_idxs = collect(first_target_id:length(spacecraft))
# # end baseline scenario


# random debris cluster scenario.
# same observer configuration
observer_defs = (
    (a_m=observer_a_m, i_deg=observer_i_deg, M_deg=observer_mean_anomalies_deg[1]),
    (a_m=observer_a_m, i_deg=observer_i_deg, M_deg=observer_mean_anomalies_deg[2]),
    (a_m=observer_a_m, i_deg=observer_i_deg, M_deg=observer_mean_anomalies_deg[3])
)
N_CLUSTER_TARGETS = 100
target_defs = generate_observer_centered_debris_cluster(
    N_CLUSTER_TARGETS,
    observer_a_m,
    observer_i_deg,
    observer_mean_anomalies_deg;
    observer_e=e,
    observer_raan_deg=raan_deg,
    observer_aop_deg=aop_deg
)



spacecraft = SpacecraftModel[]
for (sat_id, obs) in enumerate(observer_defs)
    push!(
        spacecraft,
        make_translational_spacecraft(
            sat_id,
            obs.a_m,
            e,
            obs.i_deg,
            raan_deg,
            aop_deg,
            obs.M_deg
        )
    )
end
first_target_id = length(observer_defs) + 1
for (k, tgt) in enumerate(target_defs)
    sat_id = first_target_id + k - 1
    tgt_e = hasproperty(tgt, :e) ? tgt.e : e
    tgt_raan_deg = hasproperty(tgt, :raan_deg) ? tgt.raan_deg : raan_deg
    tgt_aop_deg = hasproperty(tgt, :aop_deg) ? tgt.aop_deg : aop_deg
    push!(
        spacecraft,
        make_translational_spacecraft(
            sat_id,
            tgt.a_m,
            tgt_e,
            tgt.i_deg,
            tgt_raan_deg,
            tgt_aop_deg,
            tgt.M_deg
        )
    )
end

observer_idxs = collect(1:length(observer_defs))
target_idxs = collect(first_target_id:length(spacecraft))


optical_sensor = OpticalLOSSensorModel(  
    observer_idxs,
    target_idxs,
    THRESHOLD_DISTANCE_KM * 1e3,
    SIGMA_THETA_RAD,
    length(spacecraft)
)
comms_model = InterAgentCommunicationModel(
    observer_idxs,
    COMMUNICATION_RANGE_KM * 1e3,
    length(spacecraft)
)
observer_nav = ObserverNavigationModel(
    optical_sensor,
    comms_model,
    observer_idxs,
    target_idxs,
    length(spacecraft);
    od_perturb_enabled=ENABLE_OBSERVER_OD_PERTURBATION,
    od_pos_std_m=OBSERVER_OD_POS_STD_M,
    od_vel_std_mps=OBSERVER_OD_VEL_STD_MPS
)
gh_model = GravitationalHarmonicsModel(4, 0, EARTH_HARMONICS_FILE, planet) # up to J4 (zonal)
fusion_ref_radius_m = gh_model.reference_radius_m
fusion_j2 = begin
    c20 = (size(gh_model.C, 1) >= 3) ? gh_model.C[3, 1] : NaN
    (isfinite(c20) && c20 != 0.0) ? (-sqrt(5.0) * c20) : planet.J2
end
fusion_model = DistributedFusionModel(
    observer_nav,
    observer_idxs,
    target_idxs,
    1,
    SIGMA_THETA_RAD,
    EKF_CONSENSUS_ITERS,
    EKF_PROCESS_Q_DIAG,
    planet.μ,
    fusion_j2,
    fusion_ref_radius_m,
    length(spacecraft)
)

navigation_effectors = Any[optical_sensor]
push!(navigation_effectors, observer_nav)
push!(navigation_effectors, fusion_model)
navigation_rates = fill(NAVIGATION_RATE_SEC, length(navigation_effectors))

areas = [sc.root.ref_area for sc in spacecraft]
dynamic_effectors = (
    InverseSquaredGravityModel(),
    gh_model,
    SunMoonThirdBodyModel(planet),
    CannonballSRPModel(planet, areas, 1.3, 4.56e-6, 149_597_870_700.0)
)

now_utc = Dates.now(Dates.UTC)
args = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=true,
        verbose=true,
        generate_plots=false,
        results_directory=joinpath(REPO_ROOT, "output_observer_noDA"),
        normalize=false
    ),
    mission_configuration=MissionConfiguration(
        mission_type=MissionTime,
        keplerian=true,
        number_of_orbits=1,
        mission_time=MISSION_TIME_SEC,
        orientation_sim=false,
        num_steps_to_save=2_000
    ),
    environment_model=EnvironmentModel(
        planet=planet,
        EI=120.0,
        density_model=NoAtmosphereModel(),
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        wind=false
    ),
    dynamics_model=DynamicsModel(spacecraft, dynamic_effectors),
    guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
    navigation_model=NavigationModel(
        navigation_effectors=tuple(navigation_effectors...),
        navigation_rates=navigation_rates
    ),
    control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
    initial_time=InitialTime(
        year=Dates.year(now_utc),
        month=Dates.month(now_utc),
        day=Dates.day(now_utc),
        hour=Dates.hour(now_utc),
        minute=Dates.minute(now_utc),
        second=Float32(Dates.second(now_utc))
    ),
    integration_tolerances=IntegrationTolerances(
        reltol_orbit=1e-9,
        abstol_orbit=1e-9,
        dt_max_orbit=20.0
    )
)

println("Mission start UTC: $(now_utc)")
println("Propagation time: $(MISSION_TIME_SEC / 3600.0) hours")
println("Observer idxs: $(observer_idxs)")
println("Target idxs: $(target_idxs)")

extra_save_fields = SaveField[]
for target_idx in target_idxs
    let tgt = target_idx
        push!(
            extra_save_fields,
            SaveField(
                Symbol("dekf_target$(tgt)_state"),
                (u, t, integrator) -> begin
                    return Dict("obs$(obs)" => _fusion_state_for_save(fusion_model, obs, tgt, t) for obs in observer_idxs)
                end;
                per_satellite=false,
                column_prefix="dekf_t$(tgt)_state"
            )
        )
        push!(
            extra_save_fields,
            SaveField(
                Symbol("dekf_target$(tgt)_slot"),
                (u, t, integrator) -> begin
                    return Dict("obs$(obs)" => Float64(_fusion_slot_for_save(fusion_model, obs, tgt, t)) for obs in observer_idxs)
                end;
                per_satellite=false,
                column_prefix="dekf_t$(tgt)_slot"
            )
        )
    end
end

save_fields = vcat(SimulationModel.default_save_fields(args), extra_save_fields)
t = @elapsed run_simulation(args; isolate_state=false, save_fields=save_fields)

total_conflicts = sum(observer_nav.conflict_count_total[obs] for obs in observer_idxs)
total_skipped = sum(observer_nav.skipped_assoc_total[obs] for obs in observer_idxs)
println("Association diagnostics (DA): total_conflicts=$(total_conflicts), total_skipped_associations=$(total_skipped)")
for obs in observer_idxs
    println("  obs$(obs): conflicts=$(observer_nav.conflict_count_total[obs]), skipped=$(observer_nav.skipped_assoc_total[obs])")
end
total_disambig_calls = sum(observer_nav.disambiguation_calls_total[obs] for obs in observer_idxs)
total_ratio_success = sum(observer_nav.disambiguation_success_total[obs] for obs in observer_idxs)
total_ratio_success_correct = sum(observer_nav.disambiguation_ratio_success_correct_total[obs] for obs in observer_idxs)
total_ratio_success_wrong = sum(observer_nav.disambiguation_ratio_success_wrong_total[obs] for obs in observer_idxs)
total_ratio_fail = sum(observer_nav.disambiguation_ratio_fail_total[obs] for obs in observer_idxs)
total_b_conflict_events = sum(observer_nav.b_conflict_events_total[obs] for obs in observer_idxs)
total_b_conflict_dropped = sum(observer_nav.skipped_collision_B_total[obs] for obs in observer_idxs)
println("Conflict summary (DA):")
println("  local A (meas->track): total=$(total_disambig_calls), include_self_score=$(MEAS_ASSOC_A_INCLUDE_SELF_SCORE)")
println("    ratio: fail=$(total_ratio_fail), success=$(total_ratio_success), success_correct=$(total_ratio_success_correct), success_wrong=$(total_ratio_success_wrong)")
println("  type B (many meas -> one slot): events=$(total_b_conflict_events), dropped_measurements=$(total_b_conflict_dropped)")
println("  miss-distance ambiguity (IOD neighbors): events=$(fusion_model.miss_multi_pass_total), extra_candidates=$(fusion_model.miss_multi_pass_extra_total)")
println("  local track-to-track ambiguity (consensus Mahalanobis): events=$(fusion_model.mahal_multi_pass_total), extra_candidates=$(fusion_model.mahal_multi_pass_extra_total)")
println("  consensus-group fallback (global-id split): enabled=$(CONSENSUS_GROUP_FALLBACK_ENABLE), count=$(fusion_model.grouping_fallback_total)")

meas_unique_total = sum(observer_nav.meas_commit_unique_total[obs] for obs in observer_idxs)
meas_unique_correct = sum(observer_nav.meas_commit_unique_correct_total[obs] for obs in observer_idxs)
meas_unique_wrong = sum(observer_nav.meas_commit_unique_wrong_total[obs] for obs in observer_idxs)
meas_ambig_total = sum(observer_nav.meas_commit_ambig_total[obs] for obs in observer_idxs)
meas_ambig_correct = sum(observer_nav.meas_commit_ambig_correct_total[obs] for obs in observer_idxs)
meas_ambig_wrong = sum(observer_nav.meas_commit_ambig_wrong_total[obs] for obs in observer_idxs)
meas_ambig_dropped = sum(observer_nav.meas_commit_ambig_dropped_total[obs] for obs in observer_idxs)
meas_total_committed = meas_unique_total + meas_ambig_total
meas_total_correct = meas_unique_correct + meas_ambig_correct
meas_total_wrong = meas_unique_wrong + meas_ambig_wrong
meas_commit_acc_pct = meas_total_committed > 0 ? 100.0 * meas_total_correct / meas_total_committed : NaN
meas_unique_acc_pct = meas_unique_total > 0 ? 100.0 * meas_unique_correct / meas_unique_total : NaN
meas_ambig_acc_pct = meas_ambig_total > 0 ? 100.0 * meas_ambig_correct / meas_ambig_total : NaN

tt_commit_total = fusion_model.tt_commit_total
tt_commit_correct = fusion_model.tt_commit_correct_total
tt_commit_wrong = fusion_model.tt_commit_wrong_total
tt_commit_unknown = fusion_model.tt_commit_unknown_total
tt_skipped = fusion_model.tt_skipped_total
tt_no_candidate = fusion_model.tt_no_candidate_total
tt_gate_fail = fusion_model.tt_gate_fail_total
tt_skipped_same_target_present = fusion_model.tt_skipped_same_target_present_total
tt_skipped_no_same_target = fusion_model.tt_skipped_no_same_target_total
tt_skipped_unknown_source_target = fusion_model.tt_skipped_unknown_source_target_total
tt_attempt_total = fusion_model.tt_attempt_total
tt_known_total = tt_commit_correct + tt_commit_wrong
tt_acc_pct = tt_known_total > 0 ? 100.0 * tt_commit_correct / tt_known_total : NaN

group_total = fusion_model.consensus_group_total
group_same = fusion_model.consensus_group_same_target_total
group_mixed = fusion_model.consensus_group_mixed_target_total
group_unknown = fusion_model.consensus_group_unknown_total
group_known = group_same + group_mixed
group_same_pct_known = group_known > 0 ? 100.0 * group_same / group_known : NaN
group_same_pct_all = group_total > 0 ? 100.0 * group_same / group_total : NaN

summary_rows = DataFrame(
    section=String[],
    metric=String[],
    value=Float64[]
)
push!(summary_rows, ("meas_assoc", "committed_total", Float64(meas_total_committed)))
push!(summary_rows, ("meas_assoc", "committed_correct", Float64(meas_total_correct)))
push!(summary_rows, ("meas_assoc", "committed_wrong", Float64(meas_total_wrong)))
push!(summary_rows, ("meas_assoc", "commit_accuracy_pct", meas_commit_acc_pct))
push!(summary_rows, ("meas_assoc", "unique_committed", Float64(meas_unique_total)))
push!(summary_rows, ("meas_assoc", "unique_correct", Float64(meas_unique_correct)))
push!(summary_rows, ("meas_assoc", "unique_wrong", Float64(meas_unique_wrong)))
push!(summary_rows, ("meas_assoc", "unique_accuracy_pct", meas_unique_acc_pct))
push!(summary_rows, ("meas_assoc", "ambig_committed", Float64(meas_ambig_total)))
push!(summary_rows, ("meas_assoc", "ambig_correct", Float64(meas_ambig_correct)))
push!(summary_rows, ("meas_assoc", "ambig_wrong", Float64(meas_ambig_wrong)))
push!(summary_rows, ("meas_assoc", "ambig_dropped", Float64(meas_ambig_dropped)))
push!(summary_rows, ("meas_assoc", "ambig_accuracy_pct", meas_ambig_acc_pct))

push!(summary_rows, ("track_assoc", "tt_attempt_total", Float64(tt_attempt_total)))
push!(summary_rows, ("track_assoc", "tt_committed_total", Float64(tt_commit_total)))
push!(summary_rows, ("track_assoc", "tt_committed_correct", Float64(tt_commit_correct)))
push!(summary_rows, ("track_assoc", "tt_committed_wrong", Float64(tt_commit_wrong)))
push!(summary_rows, ("track_assoc", "tt_committed_unknown", Float64(tt_commit_unknown)))
push!(summary_rows, ("track_assoc", "tt_skipped", Float64(tt_skipped)))
push!(summary_rows, ("track_assoc", "tt_skipped_no_candidate", Float64(tt_no_candidate)))
push!(summary_rows, ("track_assoc", "tt_skipped_gate_fail", Float64(tt_gate_fail)))
push!(summary_rows, ("track_assoc", "tt_skipped_same_target_present", Float64(tt_skipped_same_target_present)))
push!(summary_rows, ("track_assoc", "tt_skipped_no_same_target", Float64(tt_skipped_no_same_target)))
push!(summary_rows, ("track_assoc", "tt_skipped_unknown_source_target", Float64(tt_skipped_unknown_source_target)))
push!(summary_rows, ("track_assoc", "tt_accuracy_pct_known_only", tt_acc_pct))
push!(summary_rows, ("track_assoc", "consensus_group_total", Float64(group_total)))
push!(summary_rows, ("track_assoc", "consensus_group_same_target", Float64(group_same)))
push!(summary_rows, ("track_assoc", "consensus_group_mixed_target", Float64(group_mixed)))
push!(summary_rows, ("track_assoc", "consensus_group_unknown_target", Float64(group_unknown)))
push!(summary_rows, ("track_assoc", "consensus_group_same_target_pct_known_only", group_same_pct_known))
push!(summary_rows, ("track_assoc", "consensus_group_same_target_pct_all", group_same_pct_all))

println("Association quality table:")
show(summary_rows; allrows=true, allcols=true)
println()

results_dir = args.simulation_settings.results_directory
mkpath(results_dir)
CSV.write(joinpath(results_dir, "association_quality_table_noDA.csv"), summary_rows)

println("COMPUTATIONAL TIME = $(t) s")
