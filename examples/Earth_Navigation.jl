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
const MISSION_TIME_SEC = 4800
const THRESHOLD_DISTANCE_KM = 500.0 
const COMMUNICATION_RANGE_KM = 300.0
const SIGMA_THETA_RAD = 0.0001
const NAVIGATION_RATE_SEC = 5.0
const NAVIGATION_DT_TOL_SEC = 0.25
const EKF_CONSENSUS_ITERS = 30
const EKF_PROCESS_Q_DIAG = SVector{6, Float64}(5e-1, 5e-1, 5e-1, 5e-3, 5e-3, 5e-3)
const NAN_LOS = SVector{3, Float64}(NaN, NaN, NaN)
const NAN_STATE6 = SVector{6, Float64}(ntuple(_ -> NaN, 6))

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
    if !(isfinite(t_last) && isfinite(t_prev) && t_last > t_prev)
        return false
    end
    dt = t_last - t_prev
    return abs(dt - NAVIGATION_RATE_SEC) <= NAVIGATION_DT_TOL_SEC
end

@inline function _has_two_measures(sensor::OpticalLOSSensorModel, sat_idx::Int, target_idx::Int)::Bool
    key = (sat_idx, target_idx)
    haskey(sensor.latest, key) || return false 
    haskey(sensor.previous, key) || return false 
    last = sensor.latest[key]
    prev = sensor.previous[key]
    t_last = last.t
    t_prev = prev.t
    has_counts = sensor.counts[sat_idx, target_idx] >= 2
    has_los = all(isfinite, last.los_unit) && all(isfinite, prev.los_unit)
    return has_counts && has_los && _is_consecutive_measure_pair(t_last, t_prev)
end

@inline function _los_rate(sensor::OpticalLOSSensorModel, sat_idx::Int, target_idx::Int)::Union{Nothing, SVector{3, Float64}}
    key = (sat_idx, target_idx)
    haskey(sensor.latest, key) || return nothing
    haskey(sensor.previous, key) || return nothing
    last = sensor.latest[key]
    prev = sensor.previous[key]
    _is_consecutive_measure_pair(last.t, prev.t) || return nothing
    dt = last.t - prev.t
    return (last.los_unit - prev.los_unit) / dt
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

@inline function _safe_unit(v::SVector{3, Float64})::Union{Nothing, SVector{3, Float64}}
    n = norm(v)
    if !isfinite(n) || n <= 0.0
        return nothing
    end
    return v / n
end

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
    if !isfinite(r_norm) || r_norm <= 0.0
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

Base.@kwdef mutable struct LocalTrack
    slot::Int
    last_meas::Union{Nothing, LOSMeasurement}
    prev_meas::Union{Nothing, LOSMeasurement}
    has_measure_now::Bool
    last_update_t::Float64
end

Base.@kwdef mutable struct ObserverNavigationModel
    sensor::OpticalLOSSensorModel
    comms::InterAgentCommunicationModel
    observer_idxs::Vector{Int}
    track_slots::Vector{Int} 
    local_tracks::Vector{Dict{Int, LocalTrack}} 
    raw_count_now::Vector{Int}
    assigned_count_now::Vector{Int}
    unassigned_count_now::Vector{Int}
end

Base.@kwdef struct AssociationCandidate
    measurement_idx::Int
    measurement::LOSMeasurement
    slot_id::Int
    local_score::Float64
    passed_local_gate::Bool
end

Base.@kwdef struct AssociationDecision
    selected_by_slot::Dict{Int, LOSMeasurement}
    rejected_count::Int
end

function ObserverNavigationModel(
    sensor::OpticalLOSSensorModel,
    comms::InterAgentCommunicationModel,
    observer_idxs::Vector{Int},
    track_slots::Vector{Int},
    num_sats::Int
)
    return ObserverNavigationModel(
        sensor=sensor,
        comms=comms,
        observer_idxs=observer_idxs,
        track_slots=track_slots,
        local_tracks=[Dict{Int, LocalTrack}() for _ in 1:num_sats],
        raw_count_now=zeros(Int, num_sats),
        assigned_count_now=zeros(Int, num_sats),
        unassigned_count_now=zeros(Int, num_sats)
    )
end


# Create (if needed) and return the local track object for this observer/slot pair.
# Single allocation point for local track entries in noDA mode.
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
            has_measure_now=false,
            last_update_t=t
        )
    end
    return tracks[slot_id]
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
    all(isfinite, last.los_unit) || return false
    all(isfinite, prev.los_unit) || return false
    return _is_consecutive_measure_pair(last.t, prev.t)
end

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
    meas = tr.last_meas
    (t - meas.t) <= (NAVIGATION_RATE_SEC + NAVIGATION_DT_TOL_SEC) || return nothing
    (meas.t - t) <= NAVIGATION_DT_TOL_SEC || return nothing
    return meas
end

# Checks if two consecutive local measures can seed IOD for this local track.
# TODO(DA): replace slot-id rule with omega/angle-based same-object validation.
@inline function _can_seed_iod_from_local_consecutive_measures(
    nav_model::ObserverNavigationModel,
    observer_id::Int,
    slot_id::Int
)::Bool
    return _track_has_two_measures(nav_model, observer_id, slot_id)
end

# Selects neighbor observers and measurement pairs usable for IOD.
# TODO(DA): replace slot-id rule with miss-distance based collaboration test.
function _select_iod_neighbors(
    nav_model::ObserverNavigationModel,
    observer_id::Int,
    slot_id::Int,
    neighbor_ids::Vector{Int}
)::Vector{Int}
    valid_neighbors = Int[]
    for neighbor_id in neighbor_ids
        if _can_seed_iod_from_local_consecutive_measures(nav_model, neighbor_id, slot_id)
            push!(valid_neighbors, neighbor_id)
        end
    end
    return valid_neighbors
end

# Cross-track match for distributed consensus among neighbors.
# TODO(DA): replace slot-id match with Mahalanobis state-level matching.
function _cross_track_match_for_consensus(
    fusion_model,
    observer_id::Int,
    local_slot_id::Int,
    neighbor_ids::Vector{Int}
)::Dict{Int, Int}
    matched_slot_by_observer = Dict{Int, Int}()
    matched_slot_by_observer[observer_id] = local_slot_id
    for neighbor_id in neighbor_ids
        neighbor_tracks = fusion_model.nav.local_tracks[neighbor_id]
        if haskey(neighbor_tracks, local_slot_id)
            matched_slot_by_observer[neighbor_id] = local_slot_id
        end
    end
    return matched_slot_by_observer
end

function _build_association_candidates_mask(
    model::ObserverNavigationModel,
    observer_id::Int,
    raw_measurements::Vector{LOSMeasurement}
)::Vector{AssociationCandidate}
    candidates = AssociationCandidate[]
    for (measurement_idx, measurement) in pairs(raw_measurements)
        candidate_slots = _local_association_measure_to_track_candidates(model, observer_id, measurement)
        isempty(candidate_slots) && continue
        for slot_id in candidate_slots
            push!(candidates, AssociationCandidate(
                measurement_idx=measurement_idx,
                measurement=measurement,
                slot_id=slot_id,
                local_score=measurement.range_m,
                passed_local_gate=true
            ))
        end
    end
    return candidates
end

# Local measure-to-local-track association candidates.
# TODO(DA): replace truth-label mapping with local innovation/Mahalanobis gating.
@inline function _local_association_measure_to_track_candidates(
    model::ObserverNavigationModel,
    observer_id::Int,
    measurement::LOSMeasurement
)::Vector{Int}
    slot_id = measurement.target
    return (slot_id in model.track_slots) ? [slot_id] : Int[]
end

# Collaborative disambiguation for ambiguous local associations.
# TODO(DA): replace local score tie-break with distributed disambiguation metric.
@inline function _collaborative_disambiguation(
    model::ObserverNavigationModel,
    observer_id::Int,
    measurement::LOSMeasurement,
    candidates::Vector{AssociationCandidate}
)::Union{Nothing, AssociationCandidate}
    isempty(candidates) && return nothing
    best = candidates[1]
    for c in candidates[2:end]
        if c.local_score < best.local_score
            best = c
        end
    end
    return best
end

function _resolve_association_mask(
    model::ObserverNavigationModel,
    observer_id::Int,
    candidates::Vector{AssociationCandidate}
)::AssociationDecision
    candidates_by_measurement = Dict{Int, Vector{AssociationCandidate}}()
    for c in candidates
        bucket = get!(candidates_by_measurement, c.measurement_idx, AssociationCandidate[])
        push!(bucket, c)
    end

    selected_candidates = AssociationCandidate[]
    rejected_count = 0
    for bucket in values(candidates_by_measurement)
        valid_bucket = [c for c in bucket if c.passed_local_gate]
        if isempty(valid_bucket)
            rejected_count += 1
            continue
        elseif length(valid_bucket) == 1
            push!(selected_candidates, valid_bucket[1])
            continue
        end
        winner = _collaborative_disambiguation(model, observer_id, valid_bucket[1].measurement, valid_bucket)
        if winner === nothing
            rejected_count += 1
        else
            push!(selected_candidates, winner)
        end
    end

    selected_by_slot = Dict{Int, LOSMeasurement}()
    for candidate in selected_candidates
        slot_id = candidate.slot_id
        measurement = candidate.measurement
        if haskey(selected_by_slot, slot_id)
            # Keep one measurement per slot; tie-break by smaller range.
            if measurement.range_m < selected_by_slot[slot_id].range_m
                rejected_count += 1
                selected_by_slot[slot_id] = measurement
            else
                rejected_count += 1
            end
        else
            selected_by_slot[slot_id] = measurement
        end
    end

    return AssociationDecision(
        selected_by_slot=selected_by_slot,
        rejected_count=rejected_count
    )
end

function SimulationModel.calcNavigationEffect!(
    model::ObserverNavigationModel,
    u,
    p,
    t::Float64,
    sat_idx::Int
)   
    # only for obsservers
    sat_idx in model.observer_idxs || return nothing
    observer_id = sat_idx

    # Clear "measured this tick" flags for this observer local catalog.
    observer_tracks = model.local_tracks[observer_id]
    for track in values(observer_tracks)
        track.has_measure_now = false
    end

    # Raw detections at this navigation tick.
    raw_measurements = model.sensor.measurements_now[observer_id]
    model.raw_count_now[observer_id] = length(raw_measurements)
    model.assigned_count_now[observer_id] = 0
    model.unassigned_count_now[observer_id] = 0

    # DA bridge pipeline:
    # 1) Build candidates.
    # 2) Resolve one selected measurement per slot (if any) through local and collaborative disambiguation.
    # 3) keep track of orphan measurements that fail association for initialize new tracks.
    # 4) Commit selected assignments to local track memory.
    # Current behavior is intentionally equivalent to noDA label-based assignment.
    candidates = _build_association_candidates_mask(model, observer_id, raw_measurements)
    decision = _resolve_association_mask(model, observer_id, candidates)
    model.unassigned_count_now[observer_id] = decision.rejected_count

    # TODO(DA): For unassigned/orphan measurements, create tentative local tracks.
    # Keep iod_initialized=false and filter_initialized=false for these tracks.
    # Promote a tentative track only when IOD readiness conditions are met.

    # Update local track memory with selected measurements.
    for (slot_id, measurement) in decision.selected_by_slot
        track = _ensure_local_track!(model, observer_id, slot_id, t) # Create track if it doesn't exist.
        track.prev_meas = track.last_meas
        track.last_meas = measurement
        track.has_measure_now = true
        track.last_update_t = t
        model.assigned_count_now[observer_id] += 1
    end
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
    last_update_t::Matrix{Float64}
    last_no_measure_warning_t::Vector{Float64}
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
        iod_estimate_state=[NAN_STATE6 for _ in 1:num_sats, _ in 1:num_sats],
        iod_estimate_covariance=[nan_cov() for _ in 1:num_sats, _ in 1:num_sats],
        iod_estimate_time_s=fill(NaN, num_sats, num_sats),
        iod_used_neighbors=zeros(Int, num_sats, num_sats),
        iod_triangulation_ready=falses(num_sats, num_sats),
        iod_initialized=falses(num_sats, num_sats),
        state=[NAN_STATE6 for _ in 1:num_sats, _ in 1:num_sats],
        covariance=[nan_cov() for _ in 1:num_sats, _ in 1:num_sats],
        state_pred=[NAN_STATE6 for _ in 1:num_sats, _ in 1:num_sats],
        covariance_pred=[nan_cov() for _ in 1:num_sats, _ in 1:num_sats],
        filter_initialized=falses(num_sats, num_sats),
        last_update_t=fill(NaN, num_sats, num_sats),
        last_no_measure_warning_t=fill(NaN, num_sats)
    )
end

@inline function _measurement_jacobian(x::SVector{6, Float64}, r_agent::SVector{3, Float64})
    ρ = x[1:3] - r_agent
    ρn = norm(ρ)
    if !isfinite(ρn) || ρn <= 0.0
        return zeros(3, 6)
    end
    H_pos = Matrix(I, 3, 3) / ρn - (ρ * ρ') / ρn^3
    return hcat(H_pos, zeros(3, 3))
end

@inline function _kepler_dynamics(
    x::SVector{6, Float64},
    μ::Float64
)::SVector{6, Float64}
    r = x[1:3]
    v = x[4:6]
    rn = norm(r)
    a = -μ * r / rn^3
    return SVector{6, Float64}(v..., a...)
end

@inline function _propagate_keplerian(
    x::SVector{6, Float64},
    μ::Float64,
    dt::Float64
)::SVector{6, Float64}
    k1 = _kepler_dynamics(x, μ)
    k2 = _kepler_dynamics(x + 0.5 * dt * k1, μ)
    k3 = _kepler_dynamics(x + 0.5 * dt * k2, μ)
    k4 = _kepler_dynamics(x + dt * k3, μ)
    return x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
end

@inline function _compute_process_jacobian(
    x::SVector{6, Float64},
    μ::Float64,
    dt::Float64
)
    r = x[1:3]
    rn = norm(r)
    if !isfinite(rn) || rn <= 0.0
        return Matrix(I, 6, 6)
    end
    I3 = Matrix(I, 3, 3)
    Ω = -μ * (I3 / rn^3 - 3.0 * (r * r') / rn^5)
    F = zeros(6, 6)
    F[1:3, 1:3] .= I3
    F[1:3, 4:6] .= I3 * dt
    F[4:6, 1:3] .= Ω * dt
    F[4:6, 4:6] .= I3
    return F
end

@inline function _max_degree_in_neighborhood(observer_id::Int, neighbor_map::Dict{Int, Vector{Int}})
    d_self = length(get(neighbor_map, observer_id, Int[]))
    d_max = d_self
    for n in get(neighbor_map, observer_id, Int[])
        d_max = max(d_max, length(get(neighbor_map, n, Int[])))
    end
    return d_max
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
            max_deg = _max_degree_in_neighborhood(observer_id, neighbor_map)
            w_ii = 1.0 - degree_i / (max_deg + 1.0)
            sum_val = w_ii * consensus_values[observer_id]
            for nid in neighbors
                degree_j = length(get(neighbor_map, nid, Int[]))
                w_ij = 1.0 / (max(degree_i, degree_j) + 1.0)
                sum_val += w_ij * consensus_values[nid]
            end
            new_values[observer_id] = sum_val
        end
        consensus_values = new_values
    end
    return consensus_values
end

@inline function _symmetrize_psd(P::Matrix{Float64})::Matrix{Float64}
    Ps = 0.5 * (P + P')
    λ, V = eigen(Symmetric(Ps))
    λ_clamped = max.(λ, 1e-12)
    return Matrix(V * Diagonal(λ_clamped) * V')
end

@inline function _is_finite_state(x::SVector{6, Float64})::Bool
    return all(isfinite, x)
end

@inline function _is_finite_cov(P::Matrix{Float64})::Bool
    return size(P, 1) == 6 && size(P, 2) == 6 && all(isfinite, P)
end

function _initialize_slot_filter_from_iod_consensus!(
    model::DistributedFusionModel,
    slot::Int,
    neighbor_map::Dict{Int, Vector{Int}},
    t::Float64,
    u
)
    active = Int[]
    for observer_id in model.observer_idxs
        if model.iod_initialized[observer_id, slot] &&
           _is_finite_state(model.iod_estimate_state[observer_id, slot]) &&
           _is_finite_cov(model.iod_estimate_covariance[observer_id, slot])
            push!(active, observer_id)
        end
    end
    isempty(active) && return nothing

    active_set = Set(active)
    active_neighbor_map = Dict{Int, Vector{Int}}(
        observer_id => [n for n in get(neighbor_map, observer_id, Int[]) if n in active_set] for observer_id in active
    )
    x_local = Dict{Int, Any}(observer_id => model.iod_estimate_state[observer_id, slot] for observer_id in active)
    P_local = Dict{Int, Any}(observer_id => model.iod_estimate_covariance[observer_id, slot] for observer_id in active)

    x_cons = _low_pass_consensus(x_local, active, active_neighbor_map, model.num_consensus_iter)
    P_cons = _low_pass_consensus(P_local, active, active_neighbor_map, model.num_consensus_iter)

    target_true = SVector{6, Float64}(
        SVector{3, Float64}(u.sc[slot].pos)...,
        SVector{3, Float64}(u.sc[slot].vel)...
    )
    for observer_id in active
        model.filter_initialized[observer_id, slot] && continue
        model.state[observer_id, slot] = SVector{6, Float64}(x_cons[observer_id])
        model.covariance[observer_id, slot] = _symmetrize_psd(Matrix(P_cons[observer_id]))
        model.state_pred[observer_id, slot] = model.state[observer_id, slot]
        model.covariance_pred[observer_id, slot] = model.covariance[observer_id, slot]
        model.filter_initialized[observer_id, slot] = true
        model.last_update_t[observer_id, slot] = t

        err = model.state[observer_id, slot] - target_true
        println("DEKF init (IOD+consensus) | slot=$(slot) | obs=$(observer_id) | t=$(round(t; digits=3)) s")
        println("  state = ", model.state[observer_id, slot])
        println("  cov(6x6) = ")
        show(stdout, "text/plain", model.covariance[observer_id, slot])
        println("\n  |err_r|= $(round(norm(err[1:3]); digits=3)) m, |err_v|= $(round(norm(err[4:6]); digits=6)) m/s\n")
    end
    return nothing
end

@inline function _fusion_state_for_save(
    model::DistributedFusionModel,
    obs::Int,
    slot::Int,
    t::Float64
)::SVector{6, Float64}
    model.filter_initialized[obs, slot] || return NAN_STATE6
    t_last = model.last_update_t[obs, slot]
    if !isfinite(t_last) || abs(t - t_last) > NAVIGATION_DT_TOL_SEC
        return NAN_STATE6
    end
    return model.state[obs, slot]
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
    # Returns valid=false when measurement/prediction is not usable.
    H = _measurement_jacobian(x_prior, observer_pos)
    if !all(isfinite, los_meas)
        return (zeros(6), zeros(6), zeros(6, 6), false)
    end
    los_unit = _safe_unit(los_meas)
    if los_unit === nothing
        return (zeros(6), zeros(6), zeros(6, 6), false)
    end

    I3 = Matrix(I, 3, 3)
    R_i = (model.sigma_theta_rad^2) * (I3 - los_unit * los_unit')
    R_inv = pinv(R_i)
    u_local = H' * (R_inv * los_meas)

    los_pred = _safe_unit(x_prior[1:3] - observer_pos)
    if los_pred === nothing
        return (zeros(6), zeros(6), zeros(6, 6), false)
    end
    uhat_local = H' * (R_inv * los_pred)
    U_local = H' * (R_inv * H)
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
    _update_neighbors!(model.nav.comms, u, p)

    observer_set = Set(model.observer_idxs) 
    neighbor_map = Dict{Int, Vector{Int}}() 
    for sid in model.observer_idxs
        neighbor_map[sid] = [n for n in model.nav.comms.neighbors[sid] if n in observer_set]
    end

    Q = Diagonal(model.process_q_diag)
    num_agents = length(model.observer_idxs)

    # Snapshot state at tick start so all observers process the same time layer.
    x_prior_snap = copy(model.state_pred) 
    P_prior_snap = deepcopy(model.covariance_pred)

    # Buffers for updated estimates and consensus; keyed by (observer_id, slot_id).
    x_upd_buf = Dict{Tuple{Int, Int}, SVector{6, Float64}}() 
    P_upd_buf = Dict{Tuple{Int, Int}, Matrix{Float64}}()     
    x_pred_buf = Dict{Tuple{Int, Int}, SVector{6, Float64}}()
    P_pred_buf = Dict{Tuple{Int, Int}, Matrix{Float64}}()    
    touched = Set{Tuple{Int, Int}}()
    all_active_slots = Set{Int}()
    for observer_id in model.observer_idxs
        for slot_id in keys(model.nav.local_tracks[observer_id])
            slot_id in model.track_slots || continue
            push!(all_active_slots, slot_id)
        end
    end

    # STEP 1/4: IOD attempt pass (observer-centric, all observers first).
    for observer_id in model.observer_idxs
        for slot_id in keys(model.nav.local_tracks[observer_id])
            slot_id in model.track_slots || continue
            if !model.iod_initialized[observer_id, slot_id]
                if !_can_seed_iod_from_local_consecutive_measures(model.nav, observer_id, slot_id)
                    model.iod_triangulation_ready[observer_id, slot_id] = false
                    model.iod_used_neighbors[observer_id, slot_id] = 0
                else
                    valid_neighbors = _select_iod_neighbors(
                        model.nav,
                        observer_id,
                        slot_id,
                        get(neighbor_map, observer_id, Int[])
                    )
                    model.iod_used_neighbors[observer_id, slot_id] = length(valid_neighbors)
                    if length(valid_neighbors) < model.min_neighbor_count
                        model.iod_triangulation_ready[observer_id, slot_id] = false
                    else
                        nodes_idx = vcat([observer_id], valid_neighbors)
                        nodes = NamedTuple[]
                        A_rows = Matrix{Float64}(undef, 0, 3)
                        b_rows = Float64[]
                        dA_rows = Matrix{Float64}(undef, 0, 3)
                        db_rows = Float64[]
                        for idx in nodes_idx
                            meas = _latest_track_measure(model.nav, idx, slot_id, t)
                            meas === nothing && continue
                            tracks = model.nav.local_tracks[idx]
                            haskey(tracks, slot_id) || continue
                            tr = tracks[slot_id]
                            tr.prev_meas === nothing && continue
                            l = meas.los_unit
                            l_prev = tr.prev_meas.los_unit
                            lrate = _track_los_rate(model.nav, idx, slot_id)
                            lrate === nothing && continue
                            l_unit = _safe_unit(l)
                            l_prev_unit = _safe_unit(l_prev)
                            (l_unit === nothing || l_prev_unit === nothing) && continue
                            r = SVector{3, Float64}(u.sc[idx].pos)
                            v = SVector{3, Float64}(u.sc[idx].vel)
                            eq = _build_agent_equations(r, v, l_unit, lrate)
                            eq === nothing && continue
                            A_rows = vcat(A_rows, eq.A)
                            append!(b_rows, eq.b)
                            dA_rows = vcat(dA_rows, eq.dA)
                            append!(db_rows, eq.db)
                            push!(nodes, (r=r, v=v, los=l_unit, los_rate=lrate))
                        end
                        if length(nodes) >= 2
                            H = [A_rows zeros(size(A_rows, 1), 3); dA_rows A_rows]
                            y = vcat(b_rows, db_rows)
                            if size(H, 1) >= 6
                                z_est = pinv(H) * y
                                x_est = SVector{3, Float64}(z_est[1], z_est[2], z_est[3])
                                x_dot_est = SVector{3, Float64}(z_est[4], z_est[5], z_est[6])
                                cov = _compute_state_covariance(nodes, x_est, x_dot_est, model.sigma_theta_rad, NAVIGATION_RATE_SEC)
                                if cov !== nothing
                                    model.iod_estimate_state[observer_id, slot_id] = SVector{6, Float64}(z_est)
                                    model.iod_estimate_covariance[observer_id, slot_id] = copy(cov)
                                    model.iod_estimate_time_s[observer_id, slot_id] = t
                                    model.iod_triangulation_ready[observer_id, slot_id] = true
                                    model.iod_initialized[observer_id, slot_id] = true
                                    println("IOD init | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s")
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # STEP 2/4: Initialize filter from IOD consensus once per active slot.
    # Assumes same slot across observers corresponds to same object.
    # needs to be updated with DA-enabled cross-track matching based on miss distance/MD on IOD estimates
    for slot_id in all_active_slots
        _initialize_slot_filter_from_iod_consensus!(model, slot_id, neighbor_map, t, u)
    end

    # STEP 3/4: Slot-wise DEKF update/predict pass with a common consensus set.
    for slot_id in all_active_slots
        observers_with_slot = [obs for obs in model.observer_idxs if haskey(model.nav.local_tracks[obs], slot_id)]
        isempty(observers_with_slot) && continue

        u_local = Dict{Int, Any}()
        uhat_local = Dict{Int, Any}()
        U_local = Dict{Int, Any}()
        has_valid_local_measure = Dict{Int, Bool}()

        # Build local information terms once per (observer, slot).
        for observer_id in observers_with_slot
            if !model.filter_initialized[observer_id, slot_id]
                u_local[observer_id] = zeros(6)
                uhat_local[observer_id] = zeros(6)
                U_local[observer_id] = zeros(6, 6)
                has_valid_local_measure[observer_id] = false
                continue
            end

            x_prior = _is_finite_state(x_prior_snap[observer_id, slot_id]) ? x_prior_snap[observer_id, slot_id] : model.state_pred[observer_id, slot_id]
            observer_pos = SVector{3, Float64}(u.sc[observer_id].pos)
            measurement = _latest_track_measure(model.nav, observer_id, slot_id, t)
            if measurement === nothing
                u_local[observer_id] = zeros(6)
                uhat_local[observer_id] = zeros(6)
                U_local[observer_id] = zeros(6, 6)
                has_valid_local_measure[observer_id] = false
                continue
            end

            u_i, uhat_i, U_i, valid_i = _build_local_information_terms(
                model,
                x_prior,
                observer_pos,
                measurement.los_unit
            )
            u_local[observer_id] = u_i
            uhat_local[observer_id] = uhat_i
            U_local[observer_id] = U_i
            has_valid_local_measure[observer_id] = valid_i
        end

        consensus_participants = [
            observer_id for observer_id in observers_with_slot
            if model.filter_initialized[observer_id, slot_id] && get(has_valid_local_measure, observer_id, false)
        ]
        participant_set = Set(consensus_participants)

        has_consensus = !isempty(consensus_participants)
        z_cons = Dict{Int, Any}()
        zhat_cons = Dict{Int, Any}()
        S_cons = Dict{Int, Any}()
        if has_consensus
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
            model.last_no_measure_warning_t[slot_id] = NaN
        elseif !isfinite(model.last_no_measure_warning_t[slot_id]) || abs(t - model.last_no_measure_warning_t[slot_id]) > NAVIGATION_DT_TOL_SEC
            @warn "DEKF correction skipped: no valid LOS for slot $(slot_id) at t=$(round(t; digits=3)) s."
            model.last_no_measure_warning_t[slot_id] = t
        end

        # Apply correction/prediction for every initialized observer on this slot.
        for observer_id in observers_with_slot
            model.filter_initialized[observer_id, slot_id] || continue
            key = (observer_id, slot_id)
            P_prior = _is_finite_cov(P_prior_snap[observer_id, slot_id]) ? P_prior_snap[observer_id, slot_id] : model.covariance_pred[observer_id, slot_id]
            x_prior = _is_finite_state(x_prior_snap[observer_id, slot_id]) ? x_prior_snap[observer_id, slot_id] : model.state_pred[observer_id, slot_id]

            if has_consensus && (observer_id in participant_set)
                S_i = S_cons[observer_id]
                P_upd = _symmetrize_psd(inv(inv(P_prior) + S_i))
                x_upd = x_prior + P_upd * (z_cons[observer_id] - zhat_cons[observer_id])
                x_pred = _propagate_keplerian(SVector{6, Float64}(x_upd), model.μ, NAVIGATION_RATE_SEC)
                F = _compute_process_jacobian(SVector{6, Float64}(x_upd), model.μ, NAVIGATION_RATE_SEC)
                P_pred = _symmetrize_psd(F * P_upd * F' + num_agents * Matrix(Q))
                x_upd_buf[key] = x_upd
                P_upd_buf[key] = P_upd
                x_pred_buf[key] = x_pred
                P_pred_buf[key] = P_pred
            else
                x_upd = x_prior
                P_upd = P_prior
                x_pred = _propagate_keplerian(SVector{6, Float64}(x_upd), model.μ, NAVIGATION_RATE_SEC)
                F = _compute_process_jacobian(SVector{6, Float64}(x_upd), model.μ, NAVIGATION_RATE_SEC)
                P_pred = _symmetrize_psd(F * P_upd * F' + num_agents * Matrix(Q))
                x_upd_buf[key] = x_upd
                P_upd_buf[key] = P_upd
                x_pred_buf[key] = x_pred
                P_pred_buf[key] = P_pred
            end
            push!(touched, key)
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
            return planet.J2000_to_pci * r_moon_j2000
        catch
            if !_moon_fallback_warned[]
                _moon_fallback_warned[] = true
                @warn "MOON SPICE ephemeris is unavailable; using analytic lunar orbit fallback for third-body acceleration."
            end
            return planet.J2000_to_pci * moon_pos_analytic_j2000(et)
        end
    end

    r_body_j2000 = lock(RuntimeServices.SPICE_LOCK) do
        SVector{3, Float64}(spkpos(name, et, "J2000", "NONE", "EARTH")[1]) * 1e3
    end
    return planet.J2000_to_pci * r_body_j2000
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

observer_defs = (
    (a_m=observer_a_m, i_deg=observer_i_deg, M_deg=observer_mean_anomalies_deg[1]),
    (a_m=observer_a_m, i_deg=observer_i_deg, M_deg=observer_mean_anomalies_deg[2]),
    (a_m=observer_a_m, i_deg=observer_i_deg, M_deg=observer_mean_anomalies_deg[3])
)
target_defs = (
    (a_m=target_a_m, i_deg=target_primary_i_deg, M_deg=target_primary_mean_anomaly_deg),
    # symmetric inclination with respect to observer orbit
    (a_m=target_a_m, i_deg=target_secondary_i_deg, M_deg=target_secondary_mean_anomaly_deg)
    #(a_m=target_a_m, i_deg=85.5, M_deg=target_secondary_mean_anomaly_deg),
    #(a_m=target_a_m, i_deg=84.5, M_deg=target_secondary_mean_anomaly_deg)
)

# Create spacecraft models for all observers and targets.
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
    push!(
        spacecraft,
        make_translational_spacecraft(
            sat_id,
            tgt.a_m,
            e,
            tgt.i_deg,
            raan_deg,
            aop_deg,
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
    length(spacecraft)
)
fusion_model = DistributedFusionModel(
    observer_nav,
    observer_idxs,
    target_idxs,
    1,
    SIGMA_THETA_RAD,
    EKF_CONSENSUS_ITERS,
    EKF_PROCESS_Q_DIAG,
    planet.μ,
    length(spacecraft)
)

navigation_effectors = Any[optical_sensor]
push!(navigation_effectors, observer_nav)
push!(navigation_effectors, fusion_model)
navigation_rates = fill(NAVIGATION_RATE_SEC, length(navigation_effectors))

areas = [sc.root.ref_area for sc in spacecraft]
dynamic_effectors = (
    InverseSquaredGravityModel(),
    GravitationalHarmonicsModel(4, 0, EARTH_HARMONICS_FILE, planet), # up to J4 (zonal)
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
println("Mode: observer-centric noDA")

extra_save_fields = SaveField[]
obs_metrics = (
    ("obs_raw_now", () -> observer_nav.raw_count_now),
    ("obs_assigned_now", () -> observer_nav.assigned_count_now),
    ("obs_unassigned_now", () -> observer_nav.unassigned_count_now)
)
for (metric_name, metric_getter) in obs_metrics
    let name = metric_name, getter = metric_getter
        push!(
            extra_save_fields,
            SaveField(
                Symbol(name),
                (u, t, integrator) -> begin
                    values = getter()
                    return Dict("obs$(obs)" => values[obs] for obs in observer_idxs)
                end;
                per_satellite=false,
                column_prefix=name
            )
        )
    end
end
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
    end
end

save_fields = vcat(SimulationModel.default_save_fields(args), extra_save_fields)
t = @elapsed run_simulation(args; isolate_state=false, save_fields=save_fields)

println("COMPUTATIONAL TIME = $(t) s")
