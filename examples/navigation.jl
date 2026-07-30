include(joinpath(@__DIR__, "common.jl"))
if !isdefined(@__MODULE__, :NavigationPaths)
    include(joinpath(@__DIR__, "navigation", "paths.jl"))
end
using .SimulationModel
using Dates
using SPICE
using StaticArrays
using LinearAlgebra
using Random
using CSV
using DataFrames
using Statistics
const RuntimeServices = SpaceAGORA.RuntimeServices
const EARTH_HARMONICS_FILE = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
const COMPARISON_OUTPUT_ROOT = NavigationPaths.env_override_path(
    "SPACEAGORA_COMPARISON_OUTPUT",
    NavigationPaths.transient_navigation_path("runs", "single")
)
const DEFAULT_START_UTC = "2026-06-15T00:00:00"
const MISSION_TIME_SEC = parse(Float64, get(ENV, "SPACEAGORA_MISSION_TIME_SEC", "10000.0"))
const THRESHOLD_DISTANCE_KM = 500.0
const COMMUNICATION_RANGE_KM = 300.0
const SIGMA_THETA_RAD = parse(Float64, get(ENV, "SPACEAGORA_SIGMA_THETA_RAD", "0.0001"))
const NAVIGATION_RATE_SEC = 5.0
const NAVIGATION_DT_TOL_SEC = 0.25
const CONSENSUS_GROUPS_REFRESH_STEPS = 50
const EKF_CONSENSUS_ITERS = 10
const EKF_PROCESS_Q_DIAG = SVector{6, Float64}(5e-1, 5e-1, 5e-1, 5e-3, 5e-3, 5e-3)
const ORPHAN_ATTACH_MAX_DTHETA_RAD = 0.5
const BASELINE_ORPHAN_ATTACH_MAX_DTHETA_RAD = parse(Float64, get(ENV, "SPACEAGORA_BASELINE_ORPHAN_ATTACH_MAX_DTHETA_RAD", "0.01"))
const ORPHAN_ATTACH_MAX_LOS_RATE_DELTA_RADPS = 0.001
const LOCAL_INIT_MIN_MEASUREMENTS = 3
const IOD_NEIGHBOR_MISS_DISTANCE_MAX_M = 200.0
const IOD_MAX_POSITION_RMS_STD_M = parse(
    Float64,
    get(ENV, "SPACEAGORA_IOD_MAX_POSITION_RMS_STD_M", "10000.0")
)
const REQUEST_IOD_ONE_STEP_VALIDATION = lowercase(strip(get(ENV, "SPACEAGORA_ENABLE_IOD_ONE_STEP_VALIDATION", "true"))) in ("1", "true", "yes", "y")
const IOD_VALIDATION_MAHAL_MAX_D2 = parse(Float64, get(ENV, "SPACEAGORA_IOD_VALIDATION_MAHAL_MAX_D2", "13.82"))
const CONSENSUS_MATCH_MAHAL_IOD_MAX_D2 = 1.0
const CONSENSUS_MATCH_MAHAL_FILTER_MAX_D2 = 22.46
const CONSENSUS_MATCH_COV_REG_EPS = 1e-6
const MEAS_ASSOC_MAHAL_MAX_D2 = 13.82
const MEAS_ASSOC_DISAMBIG_RATIO_MIN = 1.2
const TRACK_ASSOC_DISAMBIG_RATIO_MIN = 1.2
const MEAS_ASSOC_A_INCLUDE_SELF_SCORE = true
const CONSENSUS_GROUP_FALLBACK_ENABLE = false
const DEKF_INIT_ALLOW_SINGLETON_IF_NO_GROUP = true
const DEKF_CONSENSUS_UPDATE_MODE = :icf # :information, :icf, or :kcf
const DEKF_KCF_EPSILON = 0.05
const TRACK_CLOSE_AFTER_MISSED_MEASUREMENTS = 20
const ENABLE_TRACK_CLOSE_AFTER_SOLO_TARGET_MEASURE_STREAK = true
const SOLO_TARGET_MEASURE_STREAK_CLOSE_STEPS = 20
const ENABLE_OBSERVER_OD_PERTURBATION = true
const OBSERVER_OD_POS_STD_M = 0.5
const OBSERVER_OD_VEL_STD_MPS = 0.005
const NAN_LOS = SVector{3, Float64}(NaN, NaN, NaN)
const NAN_STATE6 = SVector{6, Float64}(ntuple(_ -> NaN, 6))
@inline _nan_cov6() = fill(NaN, 6, 6)
const NAVIGATION_CASES = (
    :proposed,
    :no_da,
    :centralized_oracle,
    :independent_local_da,
    :distributed_oracle_da,
    :baseline_da
)
# Select one case per run with:
#   SPACEAGORA_NAV_CASE=proposed julia --project=. examples/navigation.jl
# Output is written under SPACEAGORA_COMPARISON_OUTPUT/<case>. Direct runs use
# the operating-system temporary directory unless an output override is set.
const NAVIGATION_CASE = Symbol(get(ENV, "SPACEAGORA_NAV_CASE", "proposed"))
NAVIGATION_CASE in NAVIGATION_CASES ||
    error("Invalid SPACEAGORA_NAV_CASE=$(NAVIGATION_CASE). Valid cases: $(NAVIGATION_CASES)")
const NAVIGATION_OUTPUT_LABEL = get(ENV, "SPACEAGORA_NAV_OUTPUT_LABEL", String(NAVIGATION_CASE))
const USE_NO_DA = NAVIGATION_CASE === :no_da
const USE_CENTRALIZED_ORACLE = NAVIGATION_CASE === :centralized_oracle
const USE_INDEPENDENT_LOCAL_DA = NAVIGATION_CASE === :independent_local_da
const USE_DISTRIBUTED_ORACLE_DA = NAVIGATION_CASE === :distributed_oracle_da
const USE_DIRECT_TARGET_IDS = USE_NO_DA || USE_DISTRIBUTED_ORACLE_DA
const USE_BASELINE_DA = NAVIGATION_CASE === :baseline_da
const USE_DISTRIBUTED_ESTIMATOR = !USE_CENTRALIZED_ORACLE
const ENABLE_IOD_ONE_STEP_VALIDATION =
    (
        NAVIGATION_CASE in
        (:proposed, :no_da, :independent_local_da, :distributed_oracle_da)
    ) &&
    REQUEST_IOD_ONE_STEP_VALIDATION
const SCENARIO_RNG_SEED = parse(Int, get(ENV, "SPACEAGORA_SCENARIO_SEED", "20260612"))
const SENSOR_RNG_SEED = parse(Int, get(ENV, "SPACEAGORA_SENSOR_SEED", string(SCENARIO_RNG_SEED + 1)))
const OBSERVER_OD_RNG_SEED = parse(Int, get(ENV, "SPACEAGORA_OBSERVER_OD_SEED", string(SCENARIO_RNG_SEED + 2)))
const MEASUREMENT_BIAS_RNG_SEED = parse(Int, get(ENV, "SPACEAGORA_MEASUREMENT_BIAS_SEED", string(SCENARIO_RNG_SEED + 3)))
const MISDETECTION_RNG_SEED = parse(Int, get(ENV, "SPACEAGORA_MISDETECTION_SEED", string(SCENARIO_RNG_SEED + 4)))
const FALSE_ALARM_COUNT_RNG_SEED = parse(Int, get(ENV, "SPACEAGORA_FALSE_ALARM_COUNT_SEED", string(SCENARIO_RNG_SEED + 5)))
const FALSE_ALARM_DIRECTION_RNG_SEED = parse(Int, get(ENV, "SPACEAGORA_FALSE_ALARM_DIRECTION_SEED", string(SCENARIO_RNG_SEED + 6)))
const TRACKING_SUCCESS_ERROR_MAX_M = parse(Float64, get(ENV, "SPACEAGORA_TRACKING_SUCCESS_ERROR_MAX_M", "1000.0"))
const TRACKING_POSSIBLE_MIN_JOINT_TICKS = parse(Int, get(ENV, "SPACEAGORA_TRACKING_MIN_JOINT_TICKS", "3"))
@inline _truthy(value::AbstractString)::Bool = lowercase(strip(value)) in ("1", "true", "yes", "y", "on")
const SAVE_SIMULATION_RESULTS = _truthy(get(ENV, "SPACEAGORA_SAVE_SIMULATION_RESULTS", "false"))
const SAVE_TARGET_ESTIMATE_FIELDS = _truthy(get(ENV, "SPACEAGORA_SAVE_TARGET_ESTIMATE_FIELDS", "false"))
const SAVE_COMPARISON_DETAILED_TABLES = _truthy(get(ENV, "SPACEAGORA_SAVE_COMPARISON_DETAILED_TABLES", "false"))
const SAVE_AUXILIARY_METRIC_TABLES = _truthy(get(ENV, "SPACEAGORA_SAVE_AUXILIARY_METRIC_TABLES", "true"))
const SAVE_IOD_EVENT_GEOMETRY = _truthy(get(ENV, "SPACEAGORA_SAVE_IOD_EVENT_GEOMETRY", "false"))
const SAVE_IOD_PAIRWISE_DIAGNOSTICS = _truthy(get(ENV, "SPACEAGORA_SAVE_IOD_PAIRWISE_DIAGNOSTICS", "false"))
const TRUTH_GEOMETRY_REPLAY_ONLY = _truthy(get(ENV, "SPACEAGORA_TRUTH_GEOMETRY_REPLAY_ONLY", "false"))
const TRUTH_GEOMETRY_REPLAY_INPUT = get(ENV, "SPACEAGORA_TRUTH_GEOMETRY_REPLAY_INPUT", "")
const TRUTH_GEOMETRY_REPLAY_OUTPUT = get(ENV, "SPACEAGORA_TRUTH_GEOMETRY_REPLAY_OUTPUT", "")
const LOG_NAV_EVENTS = _truthy(get(ENV, "SPACEAGORA_LOG_NAV_EVENTS", "false"))
const LOG_SIM_PROGRESS = _truthy(get(ENV, "SPACEAGORA_LOG_SIM_PROGRESS", "false"))
const SIM_PROGRESS_INTERVAL_SEC = parse(Float64, get(ENV, "SPACEAGORA_SIM_PROGRESS_INTERVAL_SEC", "1000.0"))
const FORCE_INCREMENTAL_GC = _truthy(get(ENV, "SPACEAGORA_FORCE_INCREMENTAL_GC", "false"))
const ENABLE_NAV_TIMING = (NAVIGATION_CASE === :proposed) && _truthy(get(ENV, "SPACEAGORA_ENABLE_NAV_TIMING", "true"))
const MISDETECTION_RATE = parse(Float64, get(ENV, "SPACEAGORA_MISDETECTION_RATE", "0.05"))
const FALSE_ALARM_RATE = parse(Float64, get(ENV, "SPACEAGORA_FALSE_ALARM_RATE", "0.02"))
const MEASUREMENT_BIAS_RAD = parse(Float64, get(ENV, "SPACEAGORA_MEASUREMENT_BIAS_RAD", "1e-5")) # 1-sigma per component.
const ENABLE_MISDETECTIONS = _truthy(get(ENV, "SPACEAGORA_ENABLE_MISDETECTIONS", "false"))
const ENABLE_FALSE_ALARMS = _truthy(get(ENV, "SPACEAGORA_ENABLE_FALSE_ALARMS", "false"))
const ENABLE_MEASUREMENT_BIAS = _truthy(get(ENV, "SPACEAGORA_ENABLE_MEASUREMENT_BIAS", "false"))
const USE_LEGACY_GLOBAL_SENSOR_RNG = _truthy(get(
    ENV,
    "SPACEAGORA_USE_LEGACY_GLOBAL_SENSOR_RNG",
    "false"
))
Random.seed!(SENSOR_RNG_SEED)

function _parse_start_utc(raw::AbstractString)::DateTime
    text = strip(raw)
    endswith(text, "Z") && (text = text[1:end - 1])
    return DateTime(text)
end

const START_UTC_STRING = get(ENV, "SPACEAGORA_START_UTC", DEFAULT_START_UTC)
const START_UTC = _parse_start_utc(START_UTC_STRING)

(DEKF_CONSENSUS_UPDATE_MODE == :information || DEKF_CONSENSUS_UPDATE_MODE == :icf || DEKF_CONSENSUS_UPDATE_MODE == :kcf) ||
    error("Invalid DEKF_CONSENSUS_UPDATE_MODE=$(DEKF_CONSENSUS_UPDATE_MODE). Use :information, :icf, or :kcf.")

Base.@kwdef mutable struct LOSMeasurement
    t::Float64
    observer::Int
    target::Int
    range_m::Float64
    los_unit::SVector{3, Float64}
    observer_pos::SVector{3, Float64} = NAN_LOS
end

Base.@kwdef mutable struct OpticalLOSSensorModel
    observer_idxs::Vector{Int}
    target_idxs::Vector{Int}
    detection_range_m::Float64
    sigma_theta_rad::Float64
    measurement_bias_enabled::Bool
    measurement_bias_std_rad::Float64
    observer_bias_rotation_vectors::Vector{SVector{3, Float64}}
    measurement_noise_rng::MersenneTwister
    misdetection_rng::MersenneTwister
    false_alarm_count_rng::MersenneTwister
    false_alarm_direction_rng::MersenneTwister
    counts::Matrix{Int}
    consecutive_counts::Matrix{Int}
    latest::Dict{Tuple{Int, Int}, LOSMeasurement}
    previous::Dict{Tuple{Int, Int}, LOSMeasurement}
    measurements_now::Vector{Vector{LOSMeasurement}}
    visible_now::BitMatrix
    visible_opportunity_total::Int
    true_detection_total::Int
    missed_detection_total::Int
    false_alarm_total::Int
    observer_epoch_total::Int
    false_alarm_nonzero_epoch_total::Int
    false_alarm_multiple_epoch_total::Int
    false_alarm_max_per_epoch::Int
    biased_measurement_total::Int
    next_false_alarm_id::Int
end

function OpticalLOSSensorModel(
    observer_idxs::Vector{Int},
    target_idxs::Vector{Int},
    detection_range_m::Float64,
    sigma_theta_rad::Float64,
    num_sats::Int;
    measurement_bias_enabled::Bool=false,
    measurement_bias_std_rad::Float64=0.0,
    bias_rng::AbstractRNG=Random.default_rng(),
    measurement_noise_rng::MersenneTwister=MersenneTwister(SENSOR_RNG_SEED),
    misdetection_rng::MersenneTwister=MersenneTwister(MISDETECTION_RNG_SEED),
    false_alarm_count_rng::MersenneTwister=MersenneTwister(FALSE_ALARM_COUNT_RNG_SEED),
    false_alarm_direction_rng::MersenneTwister=MersenneTwister(FALSE_ALARM_DIRECTION_RNG_SEED)
)
    measurement_bias_std_rad >= 0.0 ||
        throw(ArgumentError("measurement_bias_std_rad must be nonnegative"))
    observer_bias_rotation_vectors = [
        SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:num_sats
    ]
    if measurement_bias_enabled && measurement_bias_std_rad > 0.0
        for observer_idx in observer_idxs
            observer_bias_rotation_vectors[observer_idx] = _sample_rotation_vector(
                measurement_bias_std_rad,
                bias_rng
            )
        end
    end
    return OpticalLOSSensorModel(
        observer_idxs=observer_idxs,
        target_idxs=target_idxs,
        detection_range_m=detection_range_m,
        sigma_theta_rad=sigma_theta_rad,
        measurement_bias_enabled=measurement_bias_enabled,
        measurement_bias_std_rad=measurement_bias_std_rad,
        observer_bias_rotation_vectors=observer_bias_rotation_vectors,
        measurement_noise_rng=measurement_noise_rng,
        misdetection_rng=misdetection_rng,
        false_alarm_count_rng=false_alarm_count_rng,
        false_alarm_direction_rng=false_alarm_direction_rng,
        counts=zeros(Int, num_sats, num_sats),
        consecutive_counts=zeros(Int, num_sats, num_sats),
        latest=Dict{Tuple{Int, Int}, LOSMeasurement}(),
        previous=Dict{Tuple{Int, Int}, LOSMeasurement}(),
        measurements_now=[LOSMeasurement[] for _ in 1:num_sats],
        visible_now=falses(num_sats, num_sats),
        visible_opportunity_total=0,
        true_detection_total=0,
        missed_detection_total=0,
        false_alarm_total=0,
        observer_epoch_total=0,
        false_alarm_nonzero_epoch_total=0,
        false_alarm_multiple_epoch_total=0,
        false_alarm_max_per_epoch=0,
        biased_measurement_total=0,
        next_false_alarm_id=1
    )
end

@inline function _sample_rotation_vector(
    sigma_rad::Float64,
    rng::AbstractRNG=Random.default_rng()
)::SVector{3, Float64}
    sigma_rad <= 0.0 && return SVector{3, Float64}(0.0, 0.0, 0.0)
    return sigma_rad * SVector{3, Float64}(randn(rng), randn(rng), randn(rng))
end

@inline function _rotate_los(
    los_true::SVector{3, Float64},
    theta_vec::SVector{3, Float64}
)::SVector{3, Float64}
    theta_sq = dot(theta_vec, theta_vec)
    iszero(theta_sq) && return los_true

    theta = sqrt(theta_sq)
    theta_cross_los = cross(theta_vec, los_true)
    sinc_theta = sin(theta) / theta
    cosc_theta = 2.0 * (sin(0.5 * theta) / theta)^2
    return los_true +
           sinc_theta * theta_cross_los +
           cosc_theta * cross(theta_vec, theta_cross_los)
end

@inline function _perturb_los(
    los_true::SVector{3, Float64},
    sigma_theta_rad::Float64,
    rng::AbstractRNG=Random.default_rng()
)::SVector{3, Float64}
    return _rotate_los(los_true, _sample_rotation_vector(sigma_theta_rad, rng))
end

@inline function _apply_observer_los_bias(
    los_true::SVector{3, Float64},
    bias_rotation_vector::SVector{3, Float64}
)::SVector{3, Float64}
    return _rotate_los(los_true, bias_rotation_vector)
end

function _random_unit_los(rng::AbstractRNG=Random.default_rng())::SVector{3, Float64}
    v = SVector{3, Float64}(randn(rng), randn(rng), randn(rng))
    while norm(v) <= 1e-12
        v = SVector{3, Float64}(randn(rng), randn(rng), randn(rng))
    end
    return v / norm(v)
end

function _sample_poisson_inverse(lambda::Float64, rng::AbstractRNG=Random.default_rng())::Int
    lambda <= 0.0 && return 0
    u = rand(rng)
    k = 0
    probability = exp(-lambda)
    cumulative = probability
    while u > cumulative
        k += 1
        probability *= lambda / k
        cumulative += probability
    end
    return k
end

function SimulationModel.calcNavigationEffect!(
    model::OpticalLOSSensorModel,
    u,
    p,
    t::Float64,
    sat_idx::Int
)
    sat_idx in model.observer_idxs || return nothing # Only compute measurements for observer satellites.
    model.observer_epoch_total += 1
    empty!(model.measurements_now[sat_idx]) # Clear measurements from previous nav tick for this observer.
    model.visible_now[sat_idx, :] .= false
    observer_pos = SVector{3, Float64}(u.sc[sat_idx].pos)

    @inbounds for target_idx in model.target_idxs # Iterate over all target satellites, Skip if observer and target are the same satellite.
        target_idx == sat_idx && continue
        rel_pos = SVector{3, Float64}(u.sc[target_idx].pos) - observer_pos
        dist = norm(rel_pos)
        if dist <= model.detection_range_m
            model.visible_now[sat_idx, target_idx] = true
            model.visible_opportunity_total += 1
            # Draw the nominal LOS error for every visibility opportunity, also
            # when the measurement is later missed. This preserves the same
            # noise realization across paired misdetection-rate sweeps.
            measurement_noise_rng = USE_LEGACY_GLOBAL_SENSOR_RNG ?
                Random.default_rng() :
                model.measurement_noise_rng
            noise_rotation = _sample_rotation_vector(
                model.sigma_theta_rad,
                measurement_noise_rng
            )
            if ENABLE_MISDETECTIONS && rand(model.misdetection_rng) < MISDETECTION_RATE
                model.missed_detection_total += 1
                continue
            end
            los_true = rel_pos / dist
            if model.measurement_bias_enabled
                los_true = _apply_observer_los_bias(
                    los_true,
                    model.observer_bias_rotation_vectors[sat_idx]
                )
                model.biased_measurement_total += 1
            end
            los_meas = _rotate_los(los_true, noise_rotation)
            model.counts[sat_idx, target_idx] += 1
            model.true_detection_total += 1
            key = (sat_idx, target_idx)
            # Before updating the latest measurement, move the current latest to previous if it exists.
            prev_meas = get(model.latest, key, nothing)
            if prev_meas !== nothing
                model.previous[key] = prev_meas
                if _is_consecutive_measure_pair(t, prev_meas.t)
                    model.consecutive_counts[sat_idx, target_idx] += 1
                else
                    model.consecutive_counts[sat_idx, target_idx] = 1
                end
            else
                model.consecutive_counts[sat_idx, target_idx] = 1
            end
            meas = LOSMeasurement(
                t=t,
                observer=sat_idx,
                target=target_idx,
                range_m=dist,
                los_unit=los_meas,
                observer_pos=observer_pos
            )
            model.latest[key] = meas
            push!(model.measurements_now[sat_idx], meas)
        end
    end

    if ENABLE_FALSE_ALARMS
        # One inverse-CDF draw per observer epoch makes counts monotone under
        # common random numbers when lambda is increased across the sweep.
        n_false = _sample_poisson_inverse(FALSE_ALARM_RATE, model.false_alarm_count_rng)
        n_false > 0 && (model.false_alarm_nonzero_epoch_total += 1)
        n_false > 1 && (model.false_alarm_multiple_epoch_total += 1)
        model.false_alarm_max_per_epoch = max(model.false_alarm_max_per_epoch, n_false)
        for _ in 1:n_false
            false_target_id = -model.next_false_alarm_id
            model.next_false_alarm_id += 1
            false_meas = LOSMeasurement(
                t=t,
                observer=sat_idx,
                target=false_target_id,
                range_m=NaN,
                los_unit=_random_unit_los(model.false_alarm_direction_rng),
                observer_pos=observer_pos
            )
            model.false_alarm_total += 1
            push!(model.measurements_now[sat_idx], false_meas)
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
    return @SMatrix [
         0.0  -v[3]   v[2]
         v[3]  0.0   -v[1]
        -v[2]  v[1]   0.0
    ]
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
        row2_l = rho' * (Sl * SR + Sn1)
        row3_l = -(rho' * SRdot + eta' * SR)
        row3_ldot = -(rho' * SR)
        row4_l = rho' * (Sl * SRdot + Sn1dot + Sldot * SR) + eta' * (Sl * SR + Sn1)
        row4_ldot = rho' * (Sl * SR + Sn1)
        J = zeros(4, 6)
        J[1, 1:3] .= vec(row1_l)
        J[2, 1:3] .= vec(row2_l)
        J[3, 1:3] .= vec(row3_l)
        J[3, 4:6] .= vec(row3_ldot)
        J[4, 1:3] .= vec(row4_l)
        J[4, 4:6] .= vec(row4_ldot)

        P = Matrix(I, 3, 3) - l * l'
        l_prev = hasproperty(nd, :los_prev) ? nd.los_prev : l
        P_prev = Matrix(I, 3, 3) - l_prev * l_prev'
        cov_l = sigma_theta^2 * P
        # For l_dot=(l_k-l_{k-1})/dt with independent LOS errors.
        cov_l_dot = (sigma_theta^2 / dt^2) * (P + P_prev)
        cov_cross = (sigma_theta^2 / dt) * P
        Sigma_i = [cov_l cov_cross; cov_cross cov_l_dot]
        Cov_i = J * Sigma_i * J'
        Cov_e[row + 1:row + 4, row + 1:row + 4] .= Cov_i
    end
    singular_values = svdvals(H)
    isempty(singular_values) && return nothing
    rank_tol = max(size(H)...) * eps(Float64) * singular_values[1]
    singular_values[end] > rank_tol || return nothing

    K = pinv(H)
    covariance = K * Cov_e * K'
    return 0.5 * (covariance + covariance')
end

# Local track memory for each observer and slot, storing the latest two LOS measures and associated context for seeding and fusion.
Base.@kwdef mutable struct LocalTrack
    slot::Int
    last_meas::Union{Nothing, LOSMeasurement}
    prev_meas::Union{Nothing, LOSMeasurement}
    plot_target_id::Int
    status::Symbol
    has_measure_now::Bool
    consecutive_missed::Int
    last_update_t::Float64
    state_estimate_now::SVector{6, Float64}
    covariance_estimate_now::Matrix{Float64}
    observer_pos_now::SVector{3, Float64}
    observer_pos_prev::SVector{3, Float64}
    is_freshly_initialized::Bool
    created_t::Float64
    first_measurement_t::Float64
    filter_initialized_t::Float64
    initialization_position_error_m::Float64
    iod_group_same_target::Int
    iod_group_class::UInt8
    closed_t::Float64
    first_target_id::Int
    last_target_id::Int
    id_switch_count::Int
    measurement_update_count::Int
end

# Compact terminal record retained after a closed track is removed from the
# operational catalog.  These fields are sufficient to reconstruct all track
# lifecycle, initialization, identity-switch, duration, and fragmentation
# metrics without keeping measurement objects or covariance matrices alive.
Base.@kwdef struct ClosedTrackLifecycle
    slot::Int
    final_status::Symbol
    created_t::Float64
    first_measurement_t::Float64
    filter_initialized_t::Float64
    initialization_position_error_m::Float64
    iod_group_same_target::Int
    iod_group_class::UInt8
    closed_t::Float64
    first_target_id::Int
    last_target_id::Int
    id_switch_count::Int
    measurement_update_count::Int
end

@inline function ClosedTrackLifecycle(track::LocalTrack)::ClosedTrackLifecycle
    return ClosedTrackLifecycle(
        slot=track.slot,
        final_status=track.status,
        created_t=track.created_t,
        first_measurement_t=track.first_measurement_t,
        filter_initialized_t=track.filter_initialized_t,
        initialization_position_error_m=track.initialization_position_error_m,
        iod_group_same_target=track.iod_group_same_target,
        iod_group_class=track.iod_group_class,
        closed_t=track.closed_t,
        first_target_id=track.first_target_id,
        last_target_id=track.last_target_id,
        id_switch_count=track.id_switch_count,
        measurement_update_count=track.measurement_update_count
    )
end

Base.@kwdef mutable struct LocalMeasurementHypothesis
    id::Int
    measurements::Vector{LOSMeasurement}
    score::Float64 = 0.0
end

Base.@kwdef struct LocalPromotionCandidate
    measurements::Vector{LOSMeasurement}
    score::Float64
end

@inline function _track_oracle_id(track::LocalTrack)::Int
    track.plot_target_id > 0 && return track.plot_target_id
    track.last_meas !== nothing && return track.last_meas.target
    return 0
end

@inline function _same_target_measurements(measurements::Vector{LOSMeasurement})::Bool
    isempty(measurements) && return false
    gid = measurements[1].target
    return gid > 0 && all(measurement -> measurement.target == gid, measurements)
end

@inline function _same_target_iod_pair(
    seed_meas::LOSMeasurement,
    seed_prev::LOSMeasurement,
    neighbor_meas::LOSMeasurement,
    neighbor_prev::LOSMeasurement
)::Bool
    gid = seed_meas.target
    return gid > 0 &&
        seed_prev.target == gid &&
        neighbor_meas.target == gid &&
        neighbor_prev.target == gid
end

@inline function _pct(num::Float64, den::Float64)::Float64
    return (isfinite(num) && isfinite(den) && den > 0.0) ? 100.0 * num / den : NaN
end

Base.@kwdef mutable struct NavTimingStats
    local_epoch_t::Float64 = NaN
    local_epoch_ns::UInt64 = UInt64(0)
    local_epoch_calls::Int = 0
    local_da_epoch_ns::Vector{Float64} = Float64[]
    cross_da_epoch_ns::Vector{Float64} = Float64[]
    filter_epoch_ns::Vector{Float64} = Float64[]
    fusion_epoch_ns::Vector{Float64} = Float64[]
    total_epoch_ns::Vector{Float64} = Float64[]
end

@inline function _same_epoch(a::Float64, b::Float64)::Bool
    return isfinite(a) && abs(a - b) <= NAVIGATION_DT_TOL_SEC
end

function _flush_local_da_timing!(timing::NavTimingStats)::Float64
    timing.local_epoch_calls > 0 || return NaN
    elapsed_ns = Float64(timing.local_epoch_ns)
    push!(timing.local_da_epoch_ns, elapsed_ns)
    timing.local_epoch_t = NaN
    timing.local_epoch_ns = UInt64(0)
    timing.local_epoch_calls = 0
    return elapsed_ns
end

function _accumulate_local_da_timing!(timing::NavTimingStats, t::Float64, elapsed_ns::UInt64)::Nothing
    if timing.local_epoch_calls > 0 && !_same_epoch(timing.local_epoch_t, t)
        _flush_local_da_timing!(timing)
    end
    if timing.local_epoch_calls == 0
        timing.local_epoch_t = t
        timing.local_epoch_ns = UInt64(0)
    end
    timing.local_epoch_ns += elapsed_ns
    timing.local_epoch_calls += 1
    return nothing
end

function _timing_summary_ms(samples_ns::Vector{Float64})
    samples_ms = Float64[
        1e-6 * ns for ns in samples_ns if isfinite(ns) && ns >= 0.0
    ]
    if isempty(samples_ms)
        return (count=0, mean=NaN, p05=NaN, median=NaN, p95=NaN, max=NaN)
    end
    return (
        count=length(samples_ms),
        mean=sum(samples_ms) / length(samples_ms),
        p05=quantile(samples_ms, 0.05),
        median=median(samples_ms),
        p95=quantile(samples_ms, 0.95),
        max=maximum(samples_ms)
    )
end


Base.@kwdef mutable struct ObserverNavigationModel
    sensor::OpticalLOSSensorModel
    comms::InterAgentCommunicationModel
    observer_idxs::Vector{Int}
    local_tracks::Vector{Dict{Int, LocalTrack}}
    closed_track_lifecycle::Vector{Vector{ClosedTrackLifecycle}}
    local_measurement_hypotheses::Vector{Dict{Int, LocalMeasurementHypothesis}}
    maintained_t2t_neighbors::Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}
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
    meas_commit_false_alarm_total::Vector{Int}
    meas_true_opportunity_total::Vector{Int}
    b_conflict_events_total::Vector{Int}
    skipped_collision_B_total::Vector{Int}
    hyp_h1_created_total::Vector{Int}
    hyp_h1_to_h2_created_total::Vector{Int}
    hyp_h1_to_h2_same_target_total::Vector{Int}
    hyp_h2_to_h3_attempted_total::Vector{Int}
    hyp_h3_los_rate_pass_total::Vector{Int}
    hyp_h3_los_rate_pass_same_target_total::Vector{Int}
    hyp_h3_los_rate_fail_total::Vector{Int}
    hyp_h3_los_rate_fail_same_target_total::Vector{Int}
    hyp_promoted_total::Vector{Int}
    hyp_promoted_same_target_total::Vector{Int}
    tracks_created_with_nonreal_measurements_total::Vector{Int}
    xm2m_candidate_pair_total::Int
    xm2m_gate_pass_total::Int
    xm2m_gate_pass_same_target_total::Int
    xm2m_gate_fail_total::Int
    xm2m_gate_fail_same_target_total::Int
    xm2m_selected_pair_total::Int
    xm2m_selected_pair_same_target_total::Int
    timing::NavTimingStats
    next_local_slot::Vector{Int}       # counter for assigning new local slot IDs to orphan measurements.
    next_hypothesis_id::Vector{Int}
end

function ObserverNavigationModel(
    sensor::OpticalLOSSensorModel,
    comms::InterAgentCommunicationModel,
    observer_idxs::Vector{Int},
    initial_slot_ids::Vector{Int},
    num_sats::Int;
    od_perturb_enabled::Bool=ENABLE_OBSERVER_OD_PERTURBATION,
    od_pos_std_m::Float64=OBSERVER_OD_POS_STD_M,
    od_vel_std_mps::Float64=OBSERVER_OD_VEL_STD_MPS,
    od_rng_seed::Int=OBSERVER_OD_RNG_SEED
)
    local_slot_seed = max(maximum(initial_slot_ids), num_sats) + 1
    return ObserverNavigationModel(
        sensor=sensor,
        comms=comms,
        observer_idxs=observer_idxs,
        local_tracks=[Dict{Int, LocalTrack}() for _ in 1:num_sats],
        closed_track_lifecycle=[ClosedTrackLifecycle[] for _ in 1:num_sats],
        local_measurement_hypotheses=[Dict{Int, LocalMeasurementHypothesis}() for _ in 1:num_sats],
        maintained_t2t_neighbors=Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}(),
        od_perturb_enabled=od_perturb_enabled,
        od_pos_std_m=od_pos_std_m,
        od_vel_std_mps=od_vel_std_mps,
        od_rng=MersenneTwister(od_rng_seed),
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
        meas_commit_false_alarm_total=zeros(Int, num_sats),
        meas_true_opportunity_total=zeros(Int, num_sats),
        b_conflict_events_total=zeros(Int, num_sats),
        skipped_collision_B_total=zeros(Int, num_sats),
        hyp_h1_created_total=zeros(Int, num_sats),
        hyp_h1_to_h2_created_total=zeros(Int, num_sats),
        hyp_h1_to_h2_same_target_total=zeros(Int, num_sats),
        hyp_h2_to_h3_attempted_total=zeros(Int, num_sats),
        hyp_h3_los_rate_pass_total=zeros(Int, num_sats),
        hyp_h3_los_rate_pass_same_target_total=zeros(Int, num_sats),
        hyp_h3_los_rate_fail_total=zeros(Int, num_sats),
        hyp_h3_los_rate_fail_same_target_total=zeros(Int, num_sats),
        hyp_promoted_total=zeros(Int, num_sats),
        hyp_promoted_same_target_total=zeros(Int, num_sats),
        tracks_created_with_nonreal_measurements_total=zeros(Int, num_sats),
        xm2m_candidate_pair_total=0,
        xm2m_gate_pass_total=0,
        xm2m_gate_pass_same_target_total=0,
        xm2m_gate_fail_total=0,
        xm2m_gate_fail_same_target_total=0,
        xm2m_selected_pair_total=0,
        xm2m_selected_pair_same_target_total=0,
        timing=NavTimingStats(),
        next_local_slot=fill(local_slot_seed, num_sats),
        next_hypothesis_id=ones(Int, num_sats)
    )
end

@inline function _true_observer_state(u, observer_id::Int)::SVector{6, Float64}
    pos = SVector{3, Float64}(u.sc[observer_id].pos)
    vel = SVector{3, Float64}(u.sc[observer_id].vel)
    return SVector{6, Float64}(pos[1], pos[2], pos[3], vel[1], vel[2], vel[3])
end

@inline function _known_observer_state!(
    model,
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
    model,
    observer_id::Int,
    u,
    t::Float64
)::SVector{3, Float64}
    x = _known_observer_state!(model, observer_id, u, t)
    return SVector{3, Float64}(x[1], x[2], x[3])
end

@inline function _known_observer_vel!(
    model,
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
            # New local tracks stay inactive until promoted by measurement hypotheses.
            status=:inactive,
            has_measure_now=false,
            consecutive_missed=0,
            last_update_t=t,
            state_estimate_now=NAN_STATE6,
            covariance_estimate_now=_nan_cov6(),
            observer_pos_now=NAN_LOS,
            observer_pos_prev=NAN_LOS,
            is_freshly_initialized=false,
            created_t=t,
            first_measurement_t=NaN,
            filter_initialized_t=NaN,
            initialization_position_error_m=NaN,
            iod_group_same_target=-1,
            iod_group_class=IOD_CLASS_UNKNOWN,
            closed_t=NaN,
            first_target_id=0,
            last_target_id=0,
            id_switch_count=0,
            measurement_update_count=0,
        )
        LOG_NAV_EVENTS && println("Local track init | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s")
    end
    return tracks[slot_id]
end

@inline function _record_track_measurement_identity!(track::LocalTrack, measurement::LOSMeasurement)
    target_id = measurement.target
    target_id <= 0 && return nothing
    if track.first_target_id == 0
        track.first_target_id = target_id
        track.last_target_id = target_id
    elseif track.last_target_id != target_id
        track.id_switch_count += 1
        track.last_target_id = target_id
    end
    track.measurement_update_count += 1
    return nothing
end

# Assigns a new local slot ID for an orphan measurement for this observer, and increments the counter for the next assignment.
@inline function _next_local_slot!(model::ObserverNavigationModel, observer_id::Int)::Int
    slot_id = model.next_local_slot[observer_id]
    model.next_local_slot[observer_id] += 1
    return slot_id
end

@inline function _next_hypothesis_id!(model::ObserverNavigationModel, observer_id::Int)::Int
    hyp_id = model.next_hypothesis_id[observer_id]
    model.next_hypothesis_id[observer_id] += 1
    return hyp_id
end

@inline function _measurement_key(measurement::LOSMeasurement)::Tuple{Int, Int, Float64}
    return (measurement.observer, measurement.target, measurement.t)
end

@inline function _los_angle_rad(a::SVector{3, Float64}, b::SVector{3, Float64})::Float64
    au = _safe_unit(a)
    bu = _safe_unit(b)
    return acos(clamp(dot(au, bu), -1.0, 1.0))
end

@inline function _hypothesis_los_rate_error(
    m1::LOSMeasurement,
    m2::LOSMeasurement,
    m3::LOSMeasurement
)::Float64
    _is_consecutive_measure_pair(m2.t, m1.t) || return Inf
    _is_consecutive_measure_pair(m3.t, m2.t) || return Inf
    rate12 = (m2.los_unit - m1.los_unit) / (m2.t - m1.t)
    rate23 = (m3.los_unit - m2.los_unit) / (m3.t - m2.t)
    return norm(rate23 - rate12)
end

function _insert_hypothesis_dedup!(
    hypotheses::Dict{Int, LocalMeasurementHypothesis},
    hypothesis::LocalMeasurementHypothesis
)::Bool
    new_key = Tuple(_measurement_key(measurement) for measurement in hypothesis.measurements)
    duplicate_id = 0
    duplicate_score = Inf
    for (hyp_id, existing) in hypotheses
        Tuple(_measurement_key(measurement) for measurement in existing.measurements) == new_key || continue
        duplicate_id = hyp_id
        duplicate_score = existing.score
        break
    end

    if duplicate_id == 0
        hypotheses[hypothesis.id] = hypothesis
        return true
    elseif hypothesis.score < duplicate_score
        delete!(hypotheses, duplicate_id)
        hypotheses[hypothesis.id] = hypothesis
        return false
    end
    return false
end

function _compact_unlocked_recent_measurements(
    measurements::Vector{LOSMeasurement},
    locked_keys::Set{Tuple{Int, Int, Float64}},
    t::Float64
)::Vector{LOSMeasurement}
    # A hypothesis survives only if it keeps an unlocked suffix ending at the current nav tick.
    kept = [measurement for measurement in measurements if !(_measurement_key(measurement) in locked_keys)]
    isempty(kept) && return LOSMeasurement[]

    last = kept[end]
    abs(last.t - t) <= NAVIGATION_DT_TOL_SEC || return LOSMeasurement[]
    length(kept) == 1 && return LOSMeasurement[last]

    prev = kept[end - 1]
    if _is_consecutive_measure_pair(last.t, prev.t)
        return LOSMeasurement[prev, last]
    end
    return LOSMeasurement[last]
end

function _promote_measurement_hypothesis_to_seed_track!(
    model::ObserverNavigationModel,
    observer_id::Int,
    measurements::Vector{LOSMeasurement},
    t::Float64
)
    length(measurements) >= 2 || return nothing
    prev_measurement = measurements[end - 1]
    last_measurement = measurements[end]
    _is_consecutive_measure_pair(last_measurement.t, prev_measurement.t) || return nothing

    slot_id = _next_local_slot!(model, observer_id)
    track = _ensure_local_track!(model, observer_id, slot_id, t)
    track.prev_meas = prev_measurement
    track.last_meas = last_measurement
    track.has_measure_now = abs(last_measurement.t - t) <= NAVIGATION_DT_TOL_SEC
    track.consecutive_missed = 0
    track.last_update_t = t
    track.status = :seed_ready
    track.first_measurement_t = minimum(measurement.t for measurement in measurements)
    track.observer_pos_prev = prev_measurement.observer_pos
    track.observer_pos_now = last_measurement.observer_pos
    for measurement in measurements
        _record_track_measurement_identity!(track, measurement)
    end
    if any(measurement -> measurement.target <= 0, measurements)
        model.tracks_created_with_nonreal_measurements_total[observer_id] += 1
    end
    return nothing
end

@inline function _refresh_track_tick_context!(
    track::LocalTrack,
    observer_pos::SVector{3, Float64}
)
    track.has_measure_now = false
    track.observer_pos_prev = track.observer_pos_now
    track.observer_pos_now = observer_pos
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
    t::Float64
)::Tuple{Vector{Tuple{Int, Int}}, Int, Int}
    seed_meas = _latest_track_measure(nav_model, observer_id, slot_id, t)
    seed_meas === nothing && return (Tuple{Int, Int}[], 0, 0)
    seed_los = _safe_unit(seed_meas.los_unit)
    norm(seed_los) > 1e-12 || return (Tuple{Int, Int}[], 0, 0)
    seed_pos = seed_meas.observer_pos
    all(isfinite, seed_pos) || return (Tuple{Int, Int}[], 0, 0)
    seed_visibility = nav_model.sensor.detection_range_m

    seed_meas_prev = nav_model.local_tracks[observer_id][slot_id].prev_meas
    seed_meas_prev === nothing && return (Tuple{Int, Int}[], 0, 0)
    seed_los_prev = _safe_unit(seed_meas_prev.los_unit)
    norm(seed_los_prev) > 1e-12 || return (Tuple{Int, Int}[], 0, 0)
    seed_pos_prev = seed_meas_prev.observer_pos
    all(isfinite, seed_pos_prev) || return (Tuple{Int, Int}[], 0, 0)
    seed_gid = seed_meas.target

    valid_neighbors = Tuple{Int, Int}[]
    multi_pass_events = 0
    multi_pass_extra = 0
    for neighbor_id in neighbor_ids
        best_slot = 0
        best_miss = Inf
        best_same_target = false
        pass_count = 0
        for neighbor_slot in keys(nav_model.local_tracks[neighbor_id])
            neighbor_track = nav_model.local_tracks[neighbor_id][neighbor_slot]
            neighbor_track.last_meas === nothing && continue
            neighbor_track.prev_meas === nothing && continue
            neighbor_track.has_measure_now || continue
            neighbor_meas = neighbor_track.last_meas
            abs(t - neighbor_meas.t) <= NAVIGATION_DT_TOL_SEC || continue
            _is_consecutive_measure_pair(neighbor_meas.t, neighbor_track.prev_meas.t) || continue
            if USE_DIRECT_TARGET_IDS && neighbor_meas.target != seed_gid
                continue
            end
            neighbor_los = _safe_unit(neighbor_meas.los_unit)
            neighbor_los_prev = _safe_unit(neighbor_track.prev_meas.los_unit)
            norm(neighbor_los) > 1e-12 && norm(neighbor_los_prev) > 1e-12 || continue
            neighbor_pos = neighbor_meas.observer_pos
            neighbor_pos_prev = neighbor_track.prev_meas.observer_pos
            all(isfinite, neighbor_pos) && all(isfinite, neighbor_pos_prev) || continue
            miss_d = _los_ray_miss_distance(
                seed_pos,
                seed_los,
                seed_visibility,
                neighbor_pos,
                neighbor_los,
                nav_model.sensor.detection_range_m
            )
            miss_d_prev = _los_ray_miss_distance(
                seed_pos_prev,
                seed_los_prev,
                seed_visibility,
                neighbor_pos_prev,
                neighbor_los_prev,
                nav_model.sensor.detection_range_m
            )
            same_target_pair = _same_target_iod_pair(
                seed_meas,
                seed_meas_prev,
                neighbor_meas,
                neighbor_track.prev_meas
            )
            nav_model.xm2m_candidate_pair_total += 1
            passes_iod_gate = if USE_DIRECT_TARGET_IDS
                true
            elseif USE_BASELINE_DA
                miss_d <= IOD_NEIGHBOR_MISS_DISTANCE_MAX_M
            else
                miss_d <= IOD_NEIGHBOR_MISS_DISTANCE_MAX_M &&
                    miss_d_prev <= IOD_NEIGHBOR_MISS_DISTANCE_MAX_M
            end
            if passes_iod_gate
                nav_model.xm2m_gate_pass_total += 1
                same_target_pair && (nav_model.xm2m_gate_pass_same_target_total += 1)
                pass_count += 1
                combined_miss = USE_BASELINE_DA ? miss_d : max(miss_d, miss_d_prev)
                if combined_miss < best_miss
                    best_miss = combined_miss
                    best_slot = neighbor_slot
                    best_same_target = same_target_pair
                end
            else
                nav_model.xm2m_gate_fail_total += 1
                same_target_pair && (nav_model.xm2m_gate_fail_same_target_total += 1)
            end
        end
        if pass_count > 1
            multi_pass_events += 1
            multi_pass_extra += (pass_count - 1)
        end

        # The best slot for this neighbor is valid only if its miss distance passes the IOD gate.
        if best_slot != 0 && (
            USE_DIRECT_TARGET_IDS ||
            best_miss <= IOD_NEIGHBOR_MISS_DISTANCE_MAX_M
        )
            nav_model.xm2m_selected_pair_total += 1
            best_same_target && (nav_model.xm2m_selected_pair_same_target_total += 1)
            push!(valid_neighbors, (neighbor_id, best_slot))
        end
    end
    return valid_neighbors, multi_pass_events, multi_pass_extra
end

@inline function _same_target_iod_group(
    nav_model::ObserverNavigationModel,
    nodes_idx::Vector{Tuple{Int, Int}}
)::Bool
    return _iod_group_class(nav_model, nodes_idx) == IOD_CLASS_SAME_TARGET
end

const IOD_CLASS_UNKNOWN = UInt8(0)
const IOD_CLASS_SAME_TARGET = UInt8(1)
const IOD_CLASS_MIXED_REAL_TARGET = UInt8(2)
const IOD_CLASS_FALSE_CONTAMINATED = UInt8(3)
const IOD_CLASS_FALSE_ONLY = UInt8(4)

@inline function _iod_class_label(iod_class::UInt8)::String
    iod_class == IOD_CLASS_SAME_TARGET && return "same_target"
    iod_class == IOD_CLASS_MIXED_REAL_TARGET && return "mixed_real_target"
    iod_class == IOD_CLASS_FALSE_CONTAMINATED && return "false_contaminated"
    iod_class == IOD_CLASS_FALSE_ONLY && return "false_only"
    return "unknown"
end

function _iod_group_class(
    nav_model::ObserverNavigationModel,
    nodes_idx::Vector{Tuple{Int, Int}}
)::UInt8
    isempty(nodes_idx) && return IOD_CLASS_UNKNOWN
    target_ids = Int[]
    for (observer_id, slot_id) in nodes_idx
        track = get(nav_model.local_tracks[observer_id], slot_id, nothing)
        track === nothing && return IOD_CLASS_UNKNOWN
        track.last_meas === nothing && return IOD_CLASS_UNKNOWN
        track.prev_meas === nothing && return IOD_CLASS_UNKNOWN
        push!(target_ids, track.prev_meas.target, track.last_meas.target)
    end
    isempty(target_ids) && return IOD_CLASS_UNKNOWN
    has_false = any(<=(0), target_ids)
    has_real = any(>(0), target_ids)
    has_false && !has_real && return IOD_CLASS_FALSE_ONLY
    has_false && has_real && return IOD_CLASS_FALSE_CONTAMINATED
    return length(unique(target_ids)) == 1 ? IOD_CLASS_SAME_TARGET : IOD_CLASS_MIXED_REAL_TARGET
end

# Resolve which estimator mode is available for this local track.
@inline function _track_match_mode(
    model,
    observer_id::Int,
    slot_id::Int
)::Union{Nothing, Symbol}
    cap = size(model.filter_initialized, 2)
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
    cap = size(model.filter_initialized, 2)
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
# Gate on local Mahalanobis distance, resolve mutual degree-one pairs, then run Hungarian on the residual.
function _cross_track_match_for_consensus_pairwise(
    fusion_model,
    source_observer_id::Int,
    source_slot_ids::Vector{Int},
    candidate_observer_id::Int,
    candidate_slot_ids::Vector{Int}
)::Tuple{Dict{Int, Int}, Dict{Int, Float64}, Int, Int}
    matched_slot_by_source_slot = Dict{Int, Int}()
    match_d2_by_source_slot = Dict{Int, Float64}()
    n_rows = length(source_slot_ids)
    n_cols = length(candidate_slot_ids)
    n_rows == 0 && return (matched_slot_by_source_slot, match_d2_by_source_slot, 0, 0)

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
        source_gid = row_source_gid[r]
        if source_gid > 0 && any(candidate_slot_ids) do candidate_slot_id
            candidate_track = get(
                fusion_model.nav.local_tracks[candidate_observer_id],
                candidate_slot_id,
                nothing
            )
            candidate_track !== nothing &&
                candidate_track.last_meas !== nothing &&
                candidate_track.last_meas.target == source_gid &&
                _track_state_cov_for_matching(
                    fusion_model,
                    candidate_observer_id,
                    candidate_slot_id,
                    mode
                ) !== nothing
        end
            fusion_model.tt_true_opportunity_total += 1
        end
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
    if USE_BASELINE_DA && n_cols > 0
        row_best_col = zeros(Int, n_rows)
        for r in 1:n_rows
            row_attempted[r] || continue
            best_c = 0
            best_d2 = Inf
            for c in 1:n_cols
                if local_d2[r, c] < best_d2
                    best_d2 = local_d2[r, c]
                    best_c = c
                end
            end
            isfinite(best_d2) && (row_best_col[r] = best_c)
        end

        col_best_row = zeros(Int, n_cols)
        for c in 1:n_cols
            best_r = 0
            best_d2 = Inf
            for r in 1:n_rows
                row_attempted[r] || continue
                if local_d2[r, c] < best_d2
                    best_d2 = local_d2[r, c]
                    best_r = r
                end
            end
            isfinite(best_d2) && (col_best_row[c] = best_r)
        end

        for r in 1:n_rows
            c = row_best_col[r]
            c == 0 && continue
            col_best_row[c] == r || continue
            push!(assigned_pairs, (r, c))
        end
    elseif n_cols > 0
        active_rows = collect(row_attempted)
        active_cols = collect(trues(n_cols))
        dummy_orphans = collect(falses(n_cols))
        # First resolve mutual degree-one pairs iteratively, then run Hungarian on the residual.
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

    ratio_dropped_rows = Set{Int}()
    ratio_filtered_pairs = Tuple{Int, Int}[]
    for (r, c) in assigned_pairs
        if USE_BASELINE_DA
            push!(ratio_filtered_pairs, (r, c))
            continue
        end

        chosen_score = local_d2[r, c]
        alt_best = Inf

        # Same source track, alternative neighbor track (row alternatives).
        for cc in 1:n_cols
            cc == c && continue
            alt_best = min(alt_best, local_d2[r, cc])
        end

        # Same neighbor track, alternative source track (column alternatives).
        for rr in 1:n_rows
            rr == r && continue
            alt_best = min(alt_best, local_d2[rr, c])
        end

        if isfinite(chosen_score) && isfinite(alt_best)
            ratio = chosen_score <= 1e-12 ? (alt_best > 1e-12 ? Inf : 1.0) : (alt_best / chosen_score)
            if ratio < TRACK_ASSOC_DISAMBIG_RATIO_MIN
                fusion_model.tt_ratio_fail_total += 1
                push!(ratio_dropped_rows, r)
                continue
            end
        end

        push!(ratio_filtered_pairs, (r, c))
    end
    assigned_pairs = ratio_filtered_pairs

    matched_rows = Set{Int}()
    for (r, c) in assigned_pairs
        source_slot_id = source_slot_ids[r]
        candidate_slot_id = candidate_slot_ids[c]
        matched_slot_by_source_slot[source_slot_id] = candidate_slot_id
        match_d2_by_source_slot[source_slot_id] = local_d2[r, c]
        push!(matched_rows, r)
    end

    for r in 1:n_rows
        row_attempted[r] || continue
        r in matched_rows && continue

        fusion_model.tt_skipped_total += 1
        if r in ratio_dropped_rows
            # Counted as a ratio failure above; keep skip diagnostics distinct from gate failures.
        elseif pass_count_by_row[r] == 0
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

    return matched_slot_by_source_slot, match_d2_by_source_slot, multi_pass_events, multi_pass_extra
end

struct T2TMutualEdge
    first::Tuple{Int, Int}
    second::Tuple{Int, Int}
    weight::Float64
end

Base.@kwdef struct MatchGroup
    tracks::Vector{Tuple{Int, Int}} = Tuple{Int, Int}[]
    selected_tracks::Vector{Tuple{Int, Int}} = Tuple{Int, Int}[]
    retained_edges::Vector{T2TMutualEdge} = T2TMutualEdge[]
end

@inline function _canonical_t2t_edge(
    first::Tuple{Int, Int},
    second::Tuple{Int, Int},
    weight::Float64
)::T2TMutualEdge
    if isless(second, first)
        return T2TMutualEdge(second, first, weight)
    end
    return T2TMutualEdge(first, second, weight)
end

# Deterministic component-level grouping for a catalog of mutually confirmed
# T2T links. Every participant that holds the same catalog obtains the same
# result: links are processed by increasing weight, with canonical endpoint
# identifiers as a tie-break, and components are merged only when their
# observer sets are disjoint.
function _build_observer_unique_t2t_components(
    active_tracks::Vector{Tuple{Int, Int}},
    mutual_edges::Vector{T2TMutualEdge}
)::Tuple{Vector{Vector{Tuple{Int, Int}}}, Vector{T2TMutualEdge}, Int}
    active_unique = sort(unique(active_tracks))
    isempty(active_unique) && return (
        Vector{Vector{Tuple{Int, Int}}}(),
        T2TMutualEdge[],
        0
    )

    active_set = Set(active_unique)
    parent = Dict{Tuple{Int, Int}, Tuple{Int, Int}}(key => key for key in active_unique)
    observers_by_root = Dict{Tuple{Int, Int}, Set{Int}}(
        key => Set{Int}([key[1]]) for key in active_unique
    )

    function find_root(key::Tuple{Int, Int})::Tuple{Int, Int}
        root = key
        while parent[root] != root
            root = parent[root]
        end
        while parent[key] != key
            next_key = parent[key]
            parent[key] = root
            key = next_key
        end
        return root
    end

    ordered_edges = [
        edge for edge in mutual_edges
        if edge.first in active_set && edge.second in active_set &&
            edge.first != edge.second && isfinite(edge.weight)
    ]
    sort!(ordered_edges; by=edge -> (
        edge.weight,
        edge.first[1], edge.first[2],
        edge.second[1], edge.second[2]
    ))

    retained_edges = T2TMutualEdge[]
    conflict_rejected = 0
    for edge in ordered_edges
        first_root = find_root(edge.first)
        second_root = find_root(edge.second)
        first_root == second_root && continue

        first_observers = observers_by_root[first_root]
        second_observers = observers_by_root[second_root]
        has_observer_conflict = any(observer_id -> observer_id in second_observers, first_observers)
        if has_observer_conflict
            conflict_rejected += 1
            continue
        end

        # Keep the lexicographically smaller root so the internal union-find
        # representation is deterministic as well.
        keep_root, merge_root = isless(second_root, first_root) ?
            (second_root, first_root) : (first_root, second_root)
        parent[merge_root] = keep_root
        union!(observers_by_root[keep_root], observers_by_root[merge_root])
        delete!(observers_by_root, merge_root)
        push!(retained_edges, edge)
    end

    tracks_by_root = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}()
    for key in active_unique
        root = find_root(key)
        push!(get!(tracks_by_root, root, Tuple{Int, Int}[]), key)
    end

    components = collect(values(tracks_by_root))
    for component in components
        sort!(component)
        observer_ids = first.(component)
        @assert length(unique(observer_ids)) == length(observer_ids) "T2T consensus group contains multiple tracks from the same observer"
    end
    sort!(components; by=component -> first(component))
    return components, retained_edges, conflict_rejected
end

function _set_maintained_t2t_neighbors!(
    nav_model::ObserverNavigationModel,
    groups::Vector{MatchGroup}
)::Nothing
    empty!(nav_model.maintained_t2t_neighbors)
    for group in groups
        for edge in group.retained_edges
            push!(get!(nav_model.maintained_t2t_neighbors, edge.first, Tuple{Int, Int}[]), edge.second)
            push!(get!(nav_model.maintained_t2t_neighbors, edge.second, Tuple{Int, Int}[]), edge.first)
        end
    end
    for neighbors in values(nav_model.maintained_t2t_neighbors)
        sort!(neighbors)
        unique!(neighbors)
    end
    return nothing
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

    active_set = Set(active_unique)

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

    if USE_INDEPENDENT_LOCAL_DA
        for key in active_unique
            selected_tracks = Tuple{Int, Int}[key]
            _record_group_target_consistency!(selected_tracks)
            push!(groups, MatchGroup(tracks=selected_tracks, selected_tracks=selected_tracks))
        end
        return groups, fallback_count, mahal_multi_pass_events, mahal_multi_pass_extra
    end

    if USE_DIRECT_TARGET_IDS
        tracks_by_gid = Dict{Int, Vector{Tuple{Int, Int}}}()
        unknown_tracks = Tuple{Int, Int}[]
        for key in active_unique
            obs, slot = key
            track = get(model.nav.local_tracks[obs], slot, nothing)
            gid = track === nothing ? 0 : _track_oracle_id(track)
            if gid > 0
                push!(get!(tracks_by_gid, gid, Tuple{Int, Int}[]), key)
            else
                push!(unknown_tracks, key)
            end
        end

        for gid in sort(collect(keys(tracks_by_gid)))
            component_tracks = sort(tracks_by_gid[gid])
            selected_tracks = _first_track_per_observer(component_tracks)
            _record_group_target_consistency!(selected_tracks)
            push!(groups, MatchGroup(tracks=component_tracks, selected_tracks=selected_tracks))
        end
        for key in sort(unknown_tracks)
            selected_tracks = Tuple{Int, Int}[key]
            _record_group_target_consistency!(selected_tracks)
            push!(groups, MatchGroup(tracks=selected_tracks, selected_tracks=selected_tracks))
        end
        return groups, fallback_count, mahal_multi_pass_events, mahal_multi_pass_extra
    end

    active_slots_by_observer = Dict{Int, Vector{Int}}()
    for (observer_id, slot_id) in active_unique
        push!(get!(active_slots_by_observer, observer_id, Int[]), slot_id)
    end
    for slot_ids in values(active_slots_by_observer)
        sort!(slot_ids)
    end

    # For each directed observer-neighbor pair, solve TT association in one shot.
    # Retain the score with each proposal so mutually confirmed links can later
    # be ordered consistently at component level.
    directed_match_proposals = Dict{
        Tuple{Int, Int},
        Dict{Tuple{Int, Int}, Float64}
    }()
    for observer_id in sort(collect(keys(active_slots_by_observer)))
        source_slots = active_slots_by_observer[observer_id]
        neighbor_ids = sort(unique(get(neighbor_map, observer_id, Int[])))
        for candidate_observer in neighbor_ids
            candidate_observer == observer_id && continue
            haskey(active_slots_by_observer, candidate_observer) || continue
            candidate_slots = active_slots_by_observer[candidate_observer]

            matched_slot_by_source_slot, match_d2_by_source_slot, pass_events, pass_extra = _cross_track_match_for_consensus_pairwise(
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
                d2 = get(match_d2_by_source_slot, source_slot, Inf)
                isfinite(d2) || continue
                proposals = get!(directed_match_proposals, key, Dict{Tuple{Int, Int}, Float64}())
                proposals[matched_key] = d2
            end
        end
    end

    mutual_edge_weights = Dict{
        Tuple{Tuple{Int, Int}, Tuple{Int, Int}},
        Float64
    }()
    for (key, matched_scores) in directed_match_proposals
        for (matched_key, forward_d2) in matched_scores
            reverse_matches = get(
                directed_match_proposals,
                matched_key,
                Dict{Tuple{Int, Int}, Float64}()
            )
            reverse_d2 = get(reverse_matches, key, Inf)
            if !isfinite(reverse_d2)
                model.tt_skipped_total += 1
                model.tt_mutual_fail_total += 1
                continue
            end
            key in active_set || continue
            matched_key in active_set || continue

            observer_id, source_slot = key
            candidate_observer, matched_slot = matched_key
            source_track = get(model.nav.local_tracks[observer_id], source_slot, nothing)
            candidate_track = get(model.nav.local_tracks[candidate_observer], matched_slot, nothing)
            source_gid = (source_track === nothing || source_track.last_meas === nothing) ? 0 : source_track.last_meas.target
            candidate_gid = (candidate_track === nothing || candidate_track.last_meas === nothing) ? 0 : candidate_track.last_meas.target

            model.tt_commit_total += 1
            if source_gid > 0 && candidate_gid > 0
                if source_gid == candidate_gid
                    model.tt_commit_correct_total += 1
                else
                    model.tt_commit_wrong_total += 1
                end
            else
                model.tt_commit_unknown_total += 1
            end

            # A mutual link has one canonical record. The maximum of the two
            # directional scores is used as its conservative ordering weight.
            edge = _canonical_t2t_edge(key, matched_key, max(forward_d2, reverse_d2))
            edge_key = (edge.first, edge.second)
            mutual_edge_weights[edge_key] = min(
                get(mutual_edge_weights, edge_key, Inf),
                edge.weight
            )
        end
    end

    mutual_edges = T2TMutualEdge[
        T2TMutualEdge(edge_key[1], edge_key[2], weight)
        for (edge_key, weight) in mutual_edge_weights
    ]
    conflict_free_components, retained_edges, conflict_rejected = _build_observer_unique_t2t_components(
        active_unique,
        mutual_edges
    )
    model.tt_component_conflict_rejected_total += conflict_rejected

    # The centralized simulation directly evaluates the deterministic result
    # that all participants would obtain after disseminating the same mutual
    # edge catalog inside each raw T2T component.
    for component_tracks in conflict_free_components
        selected_tracks = component_tracks
        component_set = Set(component_tracks)
        component_edges = T2TMutualEdge[
            edge for edge in retained_edges
            if edge.first in component_set && edge.second in component_set
        ]

        # Optional global-id fallback: if enabled and this conflict-free component
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
                    gid_set = Set(gid_tracks)
                    gid_edges = T2TMutualEdge[
                        edge for edge in component_edges
                        if edge.first in gid_set && edge.second in gid_set
                    ]
                    _record_group_target_consistency!(gid_tracks)
                    push!(groups, MatchGroup(
                        tracks=gid_tracks,
                        selected_tracks=gid_tracks,
                        retained_edges=gid_edges
                    ))
                end
                unknown_tracks = sort([key for key in component_tracks if gids[key] == 0])
                if !isempty(unknown_tracks)
                    _record_group_target_consistency!(unknown_tracks)
                    push!(groups, MatchGroup(tracks=unknown_tracks, selected_tracks=unknown_tracks))
                end
                continue
            end
        end

        _record_group_target_consistency!(selected_tracks)
        push!(groups, MatchGroup(
            tracks=component_tracks,
            selected_tracks=selected_tracks,
            retained_edges=component_edges
        ))
    end
    return groups, fallback_count, mahal_multi_pass_events, mahal_multi_pass_extra
end

Base.@kwdef struct MeasAssocCandidate
    slot_id::Int
    d2::Float64
    id_match::Bool
end

@inline function _measurement_mahalanobis_context(
    model::ObserverNavigationModel,
    x::SVector{6, Float64},
    P::Matrix{Float64},
    observer_pos::SVector{3, Float64}
)
    (_is_finite_state(x) && _is_finite_cov(P)) || return nothing
    (size(P, 1) == 6 && size(P, 2) == 6) || return nothing

    all(isfinite, observer_pos) || return nothing
    h = _safe_unit(x[1:3] - observer_pos)
    norm(h) > 1e-12 || return nothing

    H = _measurement_jacobian(x, observer_pos)
    σ2 = model.sensor.sigma_theta_rad^2
    I3 = Matrix(I, 3, 3)
    h_vec = Vector(h)
    R = σ2 * (I3 - h_vec * h_vec')
    S = H * P * H' + R
    S_sym = 0.5 * (S + S')
    all(isfinite, S_sym) || return nothing

    S_pinv = try
        pinv(S_sym)
    catch
        return nothing
    end
    all(isfinite, S_pinv) || return nothing
    return (h=h, S_pinv=S_pinv)
end

@inline function _measurement_mahalanobis_d2_from_context(
    context,
    measurement::LOSMeasurement
)::Float64
    z = measurement.los_unit
    all(isfinite, z) || return Inf
    z_unit = _safe_unit(z)
    norm(z_unit) > 1e-12 || return Inf

    ν = Vector(z_unit - context.h)
    y = context.S_pinv * ν
    d2 = dot(ν, y)
    return isfinite(d2) ? d2 : Inf
end

@inline function _measurement_mahalanobis_d2_from_state_cov(
    model::ObserverNavigationModel,
    x::SVector{6, Float64},
    P::Matrix{Float64},
    observer_pos::SVector{3, Float64},
    measurement::LOSMeasurement
)::Float64
    context = _measurement_mahalanobis_context(model, x, P, observer_pos)
    context === nothing && return Inf
    return _measurement_mahalanobis_d2_from_context(context, measurement)
end

@inline function _collaborative_score_for_candidate(
    model::ObserverNavigationModel,
    observer_id::Int,
    observer_pos::SVector{3, Float64},
    measurement::LOSMeasurement,
    candidate::MeasAssocCandidate
)::Union{Nothing, Float64}
    if USE_INDEPENDENT_LOCAL_DA || USE_BASELINE_DA || USE_DIRECT_TARGET_IDS
        return candidate.d2
    end

    tracks = model.local_tracks[observer_id]
    haskey(tracks, candidate.slot_id) || return nothing
    track_key = (observer_id, candidate.slot_id)
    maintained_neighbors = get(
        model.maintained_t2t_neighbors,
        track_key,
        Tuple{Int, Int}[]
    )

    support_scores = Float64[]
    MEAS_ASSOC_A_INCLUDE_SELF_SCORE && push!(support_scores, candidate.d2)
    communication_neighbors = model.comms.neighbors[observer_id]
    for (neighbor_id, neighbor_slot) in maintained_neighbors
        neighbor_id in communication_neighbors || continue
        neighbor_track = get(model.local_tracks[neighbor_id], neighbor_slot, nothing)
        neighbor_track === nothing && continue
        neighbor_track.status == :filter_initialized || continue
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

    # First remove isolated active nodes. Removing them can expose new
    # mutual degree-one pairs on the next pass.
    for c in 1:n_cols
        active_cols[c] || continue
        deg = 0
        for r in 1:n_rows
            active_rows[r] || continue
            isfinite(local_d2[r, c]) || continue
            deg += 1
        end
        if deg == 0
            active_cols[c] = false
            orphan_cols[c] = true
            changed = true
        end
    end

    for r in 1:n_rows
        active_rows[r] || continue
        deg = 0
        for c in 1:n_cols
            active_cols[c] || continue
            isfinite(local_d2[r, c]) || continue
            deg += 1
        end
        if deg == 0
            active_rows[r] = false
            changed = true
        end
    end

    # Commit a pair only when it is the unique finite entry in both its
    # active row AND its active column.
    unique_col_by_row = zeros(Int, n_rows)
    for r in 1:n_rows
        active_rows[r] || continue
        candidate_col = 0
        for c in 1:n_cols
            active_cols[c] || continue
            isfinite(local_d2[r, c]) || continue
            if candidate_col != 0
                candidate_col = 0
                break
            end
            candidate_col = c
        end
        unique_col_by_row[r] = candidate_col
    end

    unique_row_by_col = zeros(Int, n_cols)
    for c in 1:n_cols
        active_cols[c] || continue
        candidate_row = 0
        for r in 1:n_rows
            active_rows[r] || continue
            isfinite(local_d2[r, c]) || continue
            if candidate_row != 0
                candidate_row = 0
                break
            end
            candidate_row = r
        end
        unique_row_by_col[c] = candidate_row
    end

    for r in 1:n_rows
        active_rows[r] || continue
        c = unique_col_by_row[r]
        c == 0 && continue
        unique_row_by_col[c] == r || continue
        push!(assigned_pairs, (r, c))
        active_rows[r] = false
        active_cols[c] = false
        changed = true
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

    # A true M2T opportunity is one real LOS for which an initialized local
    # track of the same physical target is available at this epoch. Count the
    # measurement once even if duplicate tracks of that target exist.
    for measurement in raw_measurements
        measurement.target > 0 || continue
        any(slot_id -> _track_oracle_id(tracks[slot_id]) == measurement.target, slot_ids) || continue
        model.meas_true_opportunity_total[observer_id] += 1
    end

    if USE_DIRECT_TARGET_IDS
        slots_by_target = Dict{Int, Vector{Int}}()
        for slot_id in slot_ids
            track = tracks[slot_id]
            gid = _track_oracle_id(track)
            gid > 0 || continue
            push!(get!(slots_by_target, gid, Int[]), slot_id)
        end
        for slots in values(slots_by_target)
            sort!(
                slots;
                by=slot_id -> (
                    isfinite(tracks[slot_id].last_update_t) ? -tracks[slot_id].last_update_t : Inf,
                    slot_id
                )
            )
        end

        used_slots = Set{Int}()
        for measurement in raw_measurements
            best_slot = 0
            for slot_id in get(slots_by_target, measurement.target, Int[])
                slot_id in used_slots && continue
                best_slot = slot_id
                break
            end

            if best_slot == 0
                push!(orphan_measurements, measurement)
                continue
            end

            assigned_by_slot[best_slot] = measurement
            push!(used_slots, best_slot)
            model.meas_commit_unique_total[observer_id] += 1
            model.meas_commit_unique_correct_total[observer_id] += 1
        end
        return assigned_by_slot, orphan_measurements, conflict_count, skipped_count
    end

    local_d2 = fill(Inf, n_rows, n_cols)
    local_id_match = falses(n_rows, n_cols)
    for r in 1:n_rows
        slot_id = slot_ids[r]
        track = tracks[slot_id]
        for c in 1:n_cols
            measurement = raw_measurements[c]
            d2 = _measurement_mahalanobis_d2_from_state_cov(
                model,
                track.state_estimate_now,
                track.covariance_estimate_now,
                track.observer_pos_now,
                measurement
            )
            if isfinite(d2) && d2 <= MEAS_ASSOC_MAHAL_MAX_D2
                local_d2[r, c] = d2
                local_id_match[r, c] = (track.last_meas !== nothing) && (track.last_meas.target == measurement.target)
            end
        end
    end

    orphan_cols = collect(falses(n_cols))
    assigned_pairs = Tuple{Int, Int}[]

    if USE_BASELINE_DA
        row_best_col = zeros(Int, n_rows)
        for r in 1:n_rows
            best_c = 0
            best_d2 = Inf
            for c in 1:n_cols
                if local_d2[r, c] < best_d2
                    best_d2 = local_d2[r, c]
                    best_c = c
                end
            end
            isfinite(best_d2) && (row_best_col[r] = best_c)
        end

        col_best_row = zeros(Int, n_cols)
        for c in 1:n_cols
            best_r = 0
            best_d2 = Inf
            for r in 1:n_rows
                if local_d2[r, c] < best_d2
                    best_d2 = local_d2[r, c]
                    best_r = r
                end
            end
            isfinite(best_d2) && (col_best_row[c] = best_r)
        end

        matched_cols = Set{Int}()
        for c in 1:n_cols
            r = col_best_row[c]
            if r != 0 && row_best_col[r] == c
                push!(assigned_pairs, (r, c))
                push!(matched_cols, c)
            else
                orphan_cols[c] = true
            end
        end
    else
        # First assign mutual degree-one pairs and compute the residual matrix.
        active_rows = collect(trues(n_rows))
        active_cols = collect(trues(n_cols))

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
    end

    # Check ratio gate for ambiguous assignments and update counters.
    for (r, c) in assigned_pairs
        measurement = raw_measurements[c]
        slot_id = slot_ids[r]

        candidate_rows = Int[]
        for rr in 1:n_rows
            isfinite(local_d2[rr, c]) && push!(candidate_rows, rr)
        end
        candidate_cols = Int[]
        for cc in 1:n_cols
            isfinite(local_d2[r, cc]) && push!(candidate_cols, cc)
        end

        chosen_candidate = MeasAssocCandidate(slot_id=slot_id, d2=local_d2[r, c], id_match=local_id_match[r, c])
        chosen = chosen_candidate
        ratio_failed = false
        ratio_used = false
        is_ambiguous = (length(candidate_rows) > 1) || (length(candidate_cols) > 1)

        if is_ambiguous && !USE_BASELINE_DA
            conflict_count += 1
            model.disambiguation_calls_total[observer_id] += 1

            chosen_score = _collaborative_score_for_candidate(
                model,
                observer_id,
                observer_pos,
                measurement,
                chosen_candidate
            )
            alt_best = Inf

            # Alternative tracks for the selected measurement (same column).
            for rr in candidate_rows
                rr == r && continue
                cand = MeasAssocCandidate(
                    slot_id=slot_ids[rr],
                    d2=local_d2[rr, c],
                    id_match=local_id_match[rr, c]
                )
                sc = _collaborative_score_for_candidate(model, observer_id, observer_pos, measurement, cand)
                sc === nothing && continue
                alt_best = min(alt_best, sc)
            end

            # Alternative measurements for the selected track (same row).
            for cc in candidate_cols
                cc == c && continue
                alternative_measurement = raw_measurements[cc]
                cand = MeasAssocCandidate(
                    slot_id=slot_id,
                    d2=local_d2[r, cc],
                    id_match=local_id_match[r, cc]
                )
                sc = _collaborative_score_for_candidate(
                    model,
                    observer_id,
                    observer_pos,
                    alternative_measurement,
                    cand
                )
                sc === nothing && continue
                alt_best = min(alt_best, sc)
            end

            if chosen_score !== nothing && isfinite(chosen_score) && isfinite(alt_best)
                ratio = chosen_score <= 1e-12 ? (alt_best > 1e-12 ? Inf : 1.0) : (alt_best / chosen_score)
                ratio_used = true
                if ratio < MEAS_ASSOC_DISAMBIG_RATIO_MIN
                    ratio_failed = true
                    chosen = nothing
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
        measurement.target <= 0 && (model.meas_commit_false_alarm_total[observer_id] += 1)
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
    t::Float64
)
    tracks = model.local_tracks[observer_id]
    for (slot_id, measurement) in selected_by_slot
        # By construction selected slots come from existing filter tracks.
        track = tracks[slot_id]
        prev_measurement = track.last_meas
        track.prev_meas = prev_measurement
        track.last_meas = measurement
        track.has_measure_now = true
        track.consecutive_missed = 0
        track.last_update_t = t
        track.status = :filter_initialized
        track.observer_pos_prev = prev_measurement === nothing ? NAN_LOS : prev_measurement.observer_pos
        track.observer_pos_now = measurement.observer_pos
        _record_track_measurement_identity!(track, measurement)
    end
end

function _process_orphan_measurements_theta_seed!(
    model::ObserverNavigationModel,
    observer_id::Int,
    orphan_measurements::Vector{LOSMeasurement},
    t::Float64;
    oracle_identity::Bool=false,
    theta_gate_rad::Float64=ORPHAN_ATTACH_MAX_DTHETA_RAD,
    promote_measure_count::Int=2
)
    remaining_orphans = LOSMeasurement[]
    assigned_seed_slots = Set{Int}()
    tracks = model.local_tracks[observer_id]

    for measurement in orphan_measurements
        best_slot = 0
        best_score = Inf
        for (slot_id, track) in tracks
            slot_id in assigned_seed_slots && continue
            (track.status in (:seed_ready, :iod_pending, :iod_initialized)) || continue
            track.last_meas === nothing && continue
            _is_consecutive_measure_pair(measurement.t, track.last_meas.t) || continue

            theta_error = _los_angle_rad(track.last_meas.los_unit, measurement.los_unit)
            passes_gate = oracle_identity ? (_track_oracle_id(track) == measurement.target) :
                          (theta_error <= theta_gate_rad)
            passes_gate || continue

            score = oracle_identity ? 0.0 : theta_error
            if score < best_score
                best_score = score
                best_slot = slot_id
            end
        end

        if best_slot == 0
            push!(remaining_orphans, measurement)
            continue
        end

        track = tracks[best_slot]
        prev_measurement = track.last_meas
        track.prev_meas = prev_measurement
        track.last_meas = measurement
        track.has_measure_now = true
        track.consecutive_missed = 0
        track.last_update_t = t
        track.observer_pos_prev = prev_measurement === nothing ? NAN_LOS : prev_measurement.observer_pos
        track.observer_pos_now = measurement.observer_pos
        _record_track_measurement_identity!(track, measurement)
        push!(assigned_seed_slots, best_slot)
    end

    current_hypotheses = model.local_measurement_hypotheses[observer_id]
    branched_hypotheses = Dict{Int, LocalMeasurementHypothesis}()
    promotion_candidates = LocalPromotionCandidate[]
    theta_scale = max(theta_gate_rad, eps(Float64))

    for hypothesis in values(current_hypotheses)
        isempty(hypothesis.measurements) && continue
        last_measurement = hypothesis.measurements[end]

        for measurement in remaining_orphans
            any(_measurement_key(m) == _measurement_key(measurement) for m in hypothesis.measurements) && continue
            _is_consecutive_measure_pair(measurement.t, last_measurement.t) || continue

            theta_error = _los_angle_rad(last_measurement.los_unit, measurement.los_unit)
            passes_gate = oracle_identity ? (last_measurement.target == measurement.target) :
                          (theta_error <= theta_gate_rad)
            passes_gate || continue

            score = oracle_identity ? 0.0 : hypothesis.score + theta_error / theta_scale
            if length(hypothesis.measurements) == 1
                candidate_h2 = LOSMeasurement[hypothesis.measurements[1], measurement]
                model.hyp_h1_to_h2_created_total[observer_id] += 1
                if _same_target_measurements(candidate_h2)
                    model.hyp_h1_to_h2_same_target_total[observer_id] += 1
                end
            end

            if length(hypothesis.measurements) + 1 >= promote_measure_count
                start_idx = max(1, length(hypothesis.measurements) - (promote_measure_count - 2))
                candidate_measurements = LOSMeasurement[]
                append!(candidate_measurements, hypothesis.measurements[start_idx:end])
                push!(candidate_measurements, measurement)
                if promote_measure_count >= 3
                    model.hyp_h2_to_h3_attempted_total[observer_id] += 1
                end
                push!(promotion_candidates, LocalPromotionCandidate(measurements=candidate_measurements, score=score))
            else
                next_measurements = copy(hypothesis.measurements)
                push!(next_measurements, measurement)
                new_hypothesis = LocalMeasurementHypothesis(
                    id=_next_hypothesis_id!(model, observer_id),
                    measurements=next_measurements,
                    score=score
                )
                _insert_hypothesis_dedup!(branched_hypotheses, new_hypothesis)
            end
        end
    end

    sort!(promotion_candidates; by=candidate -> candidate.score)
    locked_keys = Set{Tuple{Int, Int, Float64}}()
    rejected_candidates = LocalPromotionCandidate[]

    for candidate in promotion_candidates
        if any(_measurement_key(measurement) in locked_keys for measurement in candidate.measurements)
            push!(rejected_candidates, candidate)
            continue
        end

        for measurement in candidate.measurements
            push!(locked_keys, _measurement_key(measurement))
        end
        model.hyp_promoted_total[observer_id] += 1
        if _same_target_measurements(candidate.measurements)
            model.hyp_promoted_same_target_total[observer_id] += 1
        end
        _promote_measurement_hypothesis_to_seed_track!(model, observer_id, candidate.measurements, t)
    end

    next_hypotheses = Dict{Int, LocalMeasurementHypothesis}()
    for hypothesis in values(branched_hypotheses)
        measurements = _compact_unlocked_recent_measurements(hypothesis.measurements, locked_keys, t)
        isempty(measurements) && continue
        new_hypothesis = LocalMeasurementHypothesis(
            id=_next_hypothesis_id!(model, observer_id),
            measurements=measurements,
            score=hypothesis.score
        )
        if _insert_hypothesis_dedup!(next_hypotheses, new_hypothesis) && length(measurements) == 1
            model.hyp_h1_created_total[observer_id] += 1
        end
    end

    for candidate in rejected_candidates
        measurements = _compact_unlocked_recent_measurements(candidate.measurements, locked_keys, t)
        isempty(measurements) && continue
        new_hypothesis = LocalMeasurementHypothesis(
            id=_next_hypothesis_id!(model, observer_id),
            measurements=measurements,
            score=0.0
        )
        if _insert_hypothesis_dedup!(next_hypotheses, new_hypothesis) && length(measurements) == 1
            model.hyp_h1_created_total[observer_id] += 1
        end
    end

    covered_keys = Set{Tuple{Int, Int, Float64}}()
    for hypothesis in values(next_hypotheses)
        for measurement in hypothesis.measurements
            push!(covered_keys, _measurement_key(measurement))
        end
    end
    union!(covered_keys, locked_keys)

    for measurement in remaining_orphans
        key = _measurement_key(measurement)
        key in covered_keys && continue
        new_hypothesis = LocalMeasurementHypothesis(
            id=_next_hypothesis_id!(model, observer_id),
            measurements=LOSMeasurement[measurement],
            score=0.0
        )
        if _insert_hypothesis_dedup!(next_hypotheses, new_hypothesis)
            model.hyp_h1_created_total[observer_id] += 1
        end
        push!(covered_keys, key)
    end

    model.local_measurement_hypotheses[observer_id] = next_hypotheses
    return nothing
end

function _process_orphan_measurements_simple_greedy_seed!(
    model::ObserverNavigationModel,
    observer_id::Int,
    orphan_measurements::Vector{LOSMeasurement},
    t::Float64
)
    remaining_orphans = LOSMeasurement[]
    assigned_seed_slots = Set{Int}()
    tracks = model.local_tracks[observer_id]

    for measurement in orphan_measurements
        best_slot = 0
        best_theta = Inf
        for (slot_id, track) in tracks
            slot_id in assigned_seed_slots && continue
            (track.status in (:seed_ready, :iod_pending, :iod_initialized)) || continue
            track.last_meas === nothing && continue
            _is_consecutive_measure_pair(measurement.t, track.last_meas.t) || continue

            theta_error = _los_angle_rad(track.last_meas.los_unit, measurement.los_unit)
            theta_error <= BASELINE_ORPHAN_ATTACH_MAX_DTHETA_RAD || continue
            if theta_error < best_theta
                best_theta = theta_error
                best_slot = slot_id
            end
        end

        if best_slot == 0
            push!(remaining_orphans, measurement)
            continue
        end

        track = tracks[best_slot]
        prev_measurement = track.last_meas
        track.prev_meas = prev_measurement
        track.last_meas = measurement
        track.has_measure_now = true
        track.consecutive_missed = 0
        track.last_update_t = t
        track.observer_pos_prev = prev_measurement === nothing ? NAN_LOS : prev_measurement.observer_pos
        track.observer_pos_now = measurement.observer_pos
        _record_track_measurement_identity!(track, measurement)
        push!(assigned_seed_slots, best_slot)
    end

    current_hypotheses = model.local_measurement_hypotheses[observer_id]
    theta_scale = max(BASELINE_ORPHAN_ATTACH_MAX_DTHETA_RAD, eps(Float64))
    candidates = Tuple{Float64, Int, LOSMeasurement}[]

    for (hyp_id, hypothesis) in current_hypotheses
        isempty(hypothesis.measurements) && continue
        last_measurement = hypothesis.measurements[end]
        for measurement in remaining_orphans
            any(_measurement_key(m) == _measurement_key(measurement) for m in hypothesis.measurements) && continue
            _is_consecutive_measure_pair(measurement.t, last_measurement.t) || continue
            theta_error = _los_angle_rad(last_measurement.los_unit, measurement.los_unit)
            theta_error <= BASELINE_ORPHAN_ATTACH_MAX_DTHETA_RAD || continue
            push!(candidates, (theta_error, hyp_id, measurement))
        end
    end

    sort!(candidates; by=candidate -> candidate[1])
    used_hypotheses = Set{Int}()
    used_measurements = Set{Tuple{Int, Int, Float64}}()
    next_hypotheses = Dict{Int, LocalMeasurementHypothesis}()

    for (theta_error, hyp_id, measurement) in candidates
        hyp_id in used_hypotheses && continue
        key = _measurement_key(measurement)
        key in used_measurements && continue
        haskey(current_hypotheses, hyp_id) || continue

        hypothesis = current_hypotheses[hyp_id]
        next_measurements = copy(hypothesis.measurements)
        push!(next_measurements, measurement)
        score = hypothesis.score + theta_error / theta_scale

        push!(used_hypotheses, hyp_id)
        push!(used_measurements, key)

        if length(hypothesis.measurements) == 1
            model.hyp_h1_to_h2_created_total[observer_id] += 1
            if _same_target_measurements(next_measurements)
                model.hyp_h1_to_h2_same_target_total[observer_id] += 1
            end
        end

        if length(next_measurements) >= LOCAL_INIT_MIN_MEASUREMENTS
            model.hyp_h2_to_h3_attempted_total[observer_id] += 1
            model.hyp_promoted_total[observer_id] += 1
            if _same_target_measurements(next_measurements)
                model.hyp_promoted_same_target_total[observer_id] += 1
            end
            _promote_measurement_hypothesis_to_seed_track!(model, observer_id, next_measurements, t)
        else
            new_hyp_id = _next_hypothesis_id!(model, observer_id)
            next_hypotheses[new_hyp_id] = LocalMeasurementHypothesis(
                id=new_hyp_id,
                measurements=next_measurements,
                score=score
            )
        end
    end

    for measurement in remaining_orphans
        key = _measurement_key(measurement)
        key in used_measurements && continue
        hyp_id = _next_hypothesis_id!(model, observer_id)
        next_hypotheses[hyp_id] = LocalMeasurementHypothesis(
            id=hyp_id,
            measurements=LOSMeasurement[measurement],
            score=0.0
        )
        model.hyp_h1_created_total[observer_id] += 1
    end

    model.local_measurement_hypotheses[observer_id] = next_hypotheses
    return nothing
end

# For measurements not selected by filter-track association, update seed/IOD tracks or grow local measurement hypotheses.
function _process_orphan_measurements!(
    model::ObserverNavigationModel,
    observer_id::Int,
    orphan_measurements::Vector{LOSMeasurement},
    t::Float64
)
    if USE_BASELINE_DA
        _process_orphan_measurements_simple_greedy_seed!(
            model,
            observer_id,
            orphan_measurements,
            t
        )
        return nothing
    end

    if USE_DIRECT_TARGET_IDS
        _process_orphan_measurements_theta_seed!(
            model,
            observer_id,
            orphan_measurements,
            t;
            oracle_identity=true,
            theta_gate_rad=ORPHAN_ATTACH_MAX_DTHETA_RAD,
            promote_measure_count=LOCAL_INIT_MIN_MEASUREMENTS
        )
        return nothing
    end

    remaining_orphans = LOSMeasurement[]
    assigned_seed_slots = Set{Int}()
    tracks = model.local_tracks[observer_id]

    # Existing seed/IOD tracks still need fresh LOS samples while waiting for fusion/filter initialization.
    for measurement in orphan_measurements
        best_slot = 0
        best_score = Inf
        for (slot_id, track) in tracks
            slot_id in assigned_seed_slots && continue
            (track.status in (:seed_ready, :iod_pending, :iod_initialized)) || continue
            track.last_meas === nothing && continue
            track.prev_meas === nothing && continue
            _is_consecutive_measure_pair(measurement.t, track.last_meas.t) || continue
            theta_error = _los_angle_rad(track.last_meas.los_unit, measurement.los_unit)
            theta_error <= ORPHAN_ATTACH_MAX_DTHETA_RAD || continue
            los_rate_error = _hypothesis_los_rate_error(track.prev_meas, track.last_meas, measurement)
            isfinite(los_rate_error) || continue
            los_rate_error <= ORPHAN_ATTACH_MAX_LOS_RATE_DELTA_RADPS || continue
            score = theta_error
            if score < best_score
                best_score = score
                best_slot = slot_id
            end
        end

        if best_slot == 0
            push!(remaining_orphans, measurement)
            continue
        end

        track = tracks[best_slot]
        prev_measurement = track.last_meas
        track.prev_meas = prev_measurement
        track.last_meas = measurement
        track.has_measure_now = true
        track.consecutive_missed = 0
        track.last_update_t = t
        track.observer_pos_prev = prev_measurement === nothing ? NAN_LOS : prev_measurement.observer_pos
        track.observer_pos_now = measurement.observer_pos
        _record_track_measurement_identity!(track, measurement)
        push!(assigned_seed_slots, best_slot)
    end

    current_hypotheses = model.local_measurement_hypotheses[observer_id]
    branched_hypotheses = Dict{Int, LocalMeasurementHypothesis}()
    promotion_candidates = LocalPromotionCandidate[]

    rate_scale = max(ORPHAN_ATTACH_MAX_LOS_RATE_DELTA_RADPS, eps(Float64))

    # Grow with LOS angle only; LOS-rate is checked once a three-measurement candidate exists.
    for hypothesis in values(current_hypotheses)
        isempty(hypothesis.measurements) && continue
        last_measurement = hypothesis.measurements[end]

        for measurement in remaining_orphans
            any(_measurement_key(m) == _measurement_key(measurement) for m in hypothesis.measurements) && continue

            _is_consecutive_measure_pair(measurement.t, last_measurement.t) || continue
            theta_error = _los_angle_rad(last_measurement.los_unit, measurement.los_unit)
            theta_error <= ORPHAN_ATTACH_MAX_DTHETA_RAD || continue

            if length(hypothesis.measurements) == 1
                measurements = LOSMeasurement[hypothesis.measurements[1], measurement]
                new_hypothesis = LocalMeasurementHypothesis(
                    id=_next_hypothesis_id!(model, observer_id),
                    measurements=measurements,
                    score=0.0
                )
                if _insert_hypothesis_dedup!(branched_hypotheses, new_hypothesis)
                    model.hyp_h1_to_h2_created_total[observer_id] += 1
                    first_target = measurements[1].target
                    if all(m -> m.target == first_target, measurements)
                        model.hyp_h1_to_h2_same_target_total[observer_id] += 1
                    end
                end

            else
                m1 = hypothesis.measurements[end - 1]
                m2 = hypothesis.measurements[end]
                measurements = LOSMeasurement[m1, m2, measurement]
                los_rate_error = _hypothesis_los_rate_error(m1, m2, measurement)
                score = los_rate_error / rate_scale
                first_target = measurements[1].target
                same_target = all(m -> m.target == first_target, measurements)
                model.hyp_h2_to_h3_attempted_total[observer_id] += 1

                if isfinite(los_rate_error) && los_rate_error <= ORPHAN_ATTACH_MAX_LOS_RATE_DELTA_RADPS
                    model.hyp_h3_los_rate_pass_total[observer_id] += 1
                    if same_target
                        model.hyp_h3_los_rate_pass_same_target_total[observer_id] += 1
                    end
                    push!(
                        promotion_candidates,
                        LocalPromotionCandidate(measurements=measurements, score=score)
                    )
                else
                    model.hyp_h3_los_rate_fail_total[observer_id] += 1
                    if same_target
                        model.hyp_h3_los_rate_fail_same_target_total[observer_id] += 1
                    end
                    shifted = LOSMeasurement[m2, measurement]
                    shifted_hypothesis = LocalMeasurementHypothesis(
                        id=_next_hypothesis_id!(model, observer_id),
                        measurements=shifted,
                        score=0.0
                    )
                    _insert_hypothesis_dedup!(branched_hypotheses, shifted_hypothesis)
                end
            end
        end
    end

    sort!(promotion_candidates; by=candidate -> candidate.score)
    promoted_candidates = LocalPromotionCandidate[]
    rejected_candidates = LocalPromotionCandidate[]
    locked_keys = Set{Tuple{Int, Int, Float64}}()

    for candidate in promotion_candidates
        candidate_conflicts = false
        for measurement in candidate.measurements
            if _measurement_key(measurement) in locked_keys
                candidate_conflicts = true
                break
            end
        end

        if candidate_conflicts
            push!(rejected_candidates, candidate)
            continue
        end

        push!(promoted_candidates, candidate)
        for measurement in candidate.measurements
            push!(locked_keys, _measurement_key(measurement))
        end
    end

    for candidate in promoted_candidates
        model.hyp_promoted_total[observer_id] += 1
        first_target = candidate.measurements[1].target
        if all(measurement -> measurement.target == first_target, candidate.measurements)
            model.hyp_promoted_same_target_total[observer_id] += 1
        end
        _promote_measurement_hypothesis_to_seed_track!(
            model,
            observer_id,
            candidate.measurements,
            t
        )
    end

    next_hypotheses = Dict{Int, LocalMeasurementHypothesis}()
    for hypothesis in values(branched_hypotheses)
        measurements = _compact_unlocked_recent_measurements(hypothesis.measurements, locked_keys, t)
        isempty(measurements) && continue
        new_hypothesis = LocalMeasurementHypothesis(
            id=_next_hypothesis_id!(model, observer_id),
            measurements=measurements,
            score=0.0
        )
        if _insert_hypothesis_dedup!(next_hypotheses, new_hypothesis) && length(measurements) == 1
            model.hyp_h1_created_total[observer_id] += 1
        end
    end

    for candidate in rejected_candidates
        measurements = _compact_unlocked_recent_measurements(candidate.measurements, locked_keys, t)
        isempty(measurements) && continue
        new_hypothesis = LocalMeasurementHypothesis(
            id=_next_hypothesis_id!(model, observer_id),
            measurements=measurements,
            score=0.0
        )
        if _insert_hypothesis_dedup!(next_hypotheses, new_hypothesis) && length(measurements) == 1
            model.hyp_h1_created_total[observer_id] += 1
        end
    end

    covered_keys = Set{Tuple{Int, Int, Float64}}()
    for hypothesis in values(next_hypotheses)
        for measurement in hypothesis.measurements
            push!(covered_keys, _measurement_key(measurement))
        end
    end
    union!(covered_keys, locked_keys)

    for measurement in remaining_orphans
        key = _measurement_key(measurement)
        key in covered_keys && continue
        new_hypothesis = LocalMeasurementHypothesis(
            id=_next_hypothesis_id!(model, observer_id),
            measurements=LOSMeasurement[measurement],
            score=0.0
        )
        if _insert_hypothesis_dedup!(next_hypotheses, new_hypothesis)
            model.hyp_h1_created_total[observer_id] += 1
        end
        push!(covered_keys, key)
    end

    model.local_measurement_hypotheses[observer_id] = next_hypotheses
    return nothing
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
    timing_start_ns = ENABLE_NAV_TIMING ? time_ns() : UInt64(0)

    # Refresh per-tick local-track context for this observer.
    observer_tracks = model.local_tracks[observer_id]
    for track in values(observer_tracks)
        _refresh_track_tick_context!(track, observer_pos)
    end

    # Raw detections at this navigation tick.
    raw_measurements = model.sensor.measurements_now[observer_id]
    for measurement in raw_measurements
        measurement.observer_pos = observer_pos
    end

    # Local track update pipeline:
    # 1) Split raw measurements into assigned (max 1 per slot) and orphan.
    # 2) Commit assigned measurements to local tracks.
    # 3) Use orphan measurements to update seed/IOD tracks or grow local hypotheses.

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
        t
    )

    # Update seed/IOD tracks or grow local measurement hypotheses from orphan measurements.
    _process_orphan_measurements!(
        model,
        observer_id,
        orphan_measurements,
        t
    )

    if ENABLE_NAV_TIMING
        _accumulate_local_da_timing!(model.timing, t, time_ns() - timing_start_ns)
    end

    return nothing
end

Base.@kwdef mutable struct DistributedFusionModel
    nav::ObserverNavigationModel
    observer_idxs::Vector{Int}
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
    iod_pending::BitMatrix
    iod_pending_keys::Set{Tuple{Int, Int}}
    iod_pending_same_target::BitMatrix
    iod_pending_class::Matrix{UInt8}
    iod_initialized::BitMatrix
    state::Matrix{SVector{6, Float64}}
    covariance::Matrix{Matrix{Float64}}
    state_pred::Matrix{SVector{6, Float64}}
    covariance_pred::Matrix{Matrix{Float64}}
    filter_initialized::BitMatrix
    solo_target_measure_streak::Matrix{Int}
    last_update_t::Matrix{Float64}
    grouping_fallback_total::Int
    miss_multi_pass_total::Int
    miss_multi_pass_extra_total::Int
    mahal_multi_pass_total::Int
    mahal_multi_pass_extra_total::Int
    buffer_matching_groups::Vector{MatchGroup}
    tt_attempt_total::Int
    tt_true_opportunity_total::Int
    tt_commit_total::Int
    tt_commit_correct_total::Int
    tt_commit_wrong_total::Int
    tt_commit_unknown_total::Int
    tt_skipped_total::Int
    tt_no_candidate_total::Int
    tt_gate_fail_total::Int
    tt_ratio_fail_total::Int
    tt_mutual_fail_total::Int
    tt_component_conflict_rejected_total::Int
    tt_skipped_same_target_present_total::Int
    tt_skipped_no_same_target_total::Int
    tt_skipped_unknown_source_target_total::Int
    iod_group_total::Int
    iod_group_same_target_total::Int
    iod_group_mixed_target_total::Int
    iod_init_total::Int
    iod_init_same_target_total::Int
    iod_init_mixed_target_total::Int
    iod_position_cov_gate_evaluated_total::Int
    iod_position_cov_gate_rejected_total::Int
    iod_position_cov_gate_rejected_same_target_total::Int
    iod_position_cov_gate_rejected_mixed_target_total::Int
    iod_validation_attempted_total::Int
    iod_validation_confirmed_total::Int
    iod_validation_rejected_total::Int
    iod_validation_confirmed_same_target_total::Int
    iod_validation_confirmed_mixed_target_total::Int
    iod_validation_rejected_same_target_total::Int
    iod_validation_rejected_mixed_target_total::Int
    iod_validation_no_measure_total::Int
    iod_validation_no_measure_same_target_total::Int
    iod_validation_no_measure_mixed_target_total::Int
    iod_diagnostic_rows::Vector{Any}
    iod_pairwise_rows::Vector{Any}
    iod_geometry_cache_time_s::Float64
    iod_geometry_cache::Dict{Int, Any}
    consensus_group_total::Int
    consensus_group_same_target_total::Int
    consensus_group_mixed_target_total::Int
    consensus_group_unknown_total::Int
end

function DistributedFusionModel(
    nav::ObserverNavigationModel,
    observer_idxs::Vector{Int},
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
    n_observer_rows = maximum(observer_idxs)
    n_slot_cols = num_sats
    return DistributedFusionModel(
        nav=nav,
        observer_idxs=observer_idxs,
        min_neighbor_count=min_neighbor_count,
        sigma_theta_rad=sigma_theta_rad,
        num_consensus_iter=num_consensus_iter,
        process_q_diag=process_q_diag,
        μ=μ,
        J2=J2,
        R_ref_m=R_ref_m,
        iod_estimate_state=[NAN_STATE6 for _ in 1:n_observer_rows, _ in 1:n_slot_cols],
        iod_estimate_covariance=[nan_cov() for _ in 1:n_observer_rows, _ in 1:n_slot_cols],
        iod_estimate_time_s=fill(NaN, n_observer_rows, n_slot_cols),
        iod_used_neighbors=zeros(Int, n_observer_rows, n_slot_cols),
        iod_triangulation_ready=falses(n_observer_rows, n_slot_cols),
        iod_pending=falses(n_observer_rows, n_slot_cols),
        iod_pending_keys=Set{Tuple{Int, Int}}(),
        iod_pending_same_target=falses(n_observer_rows, n_slot_cols),
        iod_pending_class=fill(IOD_CLASS_UNKNOWN, n_observer_rows, n_slot_cols),
        iod_initialized=falses(n_observer_rows, n_slot_cols),
        state=[NAN_STATE6 for _ in 1:n_observer_rows, _ in 1:n_slot_cols],
        covariance=[nan_cov() for _ in 1:n_observer_rows, _ in 1:n_slot_cols],
        state_pred=[NAN_STATE6 for _ in 1:n_observer_rows, _ in 1:n_slot_cols],
        covariance_pred=[nan_cov() for _ in 1:n_observer_rows, _ in 1:n_slot_cols],
        filter_initialized=falses(n_observer_rows, n_slot_cols),
        solo_target_measure_streak=zeros(Int, n_observer_rows, n_slot_cols),
        last_update_t=fill(NaN, n_observer_rows, n_slot_cols),
        grouping_fallback_total=0,
        miss_multi_pass_total=0,
        miss_multi_pass_extra_total=0,
        mahal_multi_pass_total=0,
        mahal_multi_pass_extra_total=0,
        buffer_matching_groups=MatchGroup[],
        tt_attempt_total=0,
        tt_true_opportunity_total=0,
        tt_commit_total=0,
        tt_commit_correct_total=0,
        tt_commit_wrong_total=0,
        tt_commit_unknown_total=0,
        tt_skipped_total=0,
        tt_no_candidate_total=0,
        tt_gate_fail_total=0,
        tt_ratio_fail_total=0,
        tt_mutual_fail_total=0,
        tt_component_conflict_rejected_total=0,
        tt_skipped_same_target_present_total=0,
        tt_skipped_no_same_target_total=0,
        tt_skipped_unknown_source_target_total=0,
        iod_group_total=0,
        iod_group_same_target_total=0,
        iod_group_mixed_target_total=0,
        iod_init_total=0,
        iod_init_same_target_total=0,
        iod_init_mixed_target_total=0,
        iod_position_cov_gate_evaluated_total=0,
        iod_position_cov_gate_rejected_total=0,
        iod_position_cov_gate_rejected_same_target_total=0,
        iod_position_cov_gate_rejected_mixed_target_total=0,
        iod_validation_attempted_total=0,
        iod_validation_confirmed_total=0,
        iod_validation_rejected_total=0,
        iod_validation_confirmed_same_target_total=0,
        iod_validation_confirmed_mixed_target_total=0,
        iod_validation_rejected_same_target_total=0,
        iod_validation_rejected_mixed_target_total=0,
        iod_validation_no_measure_total=0,
        iod_validation_no_measure_same_target_total=0,
        iod_validation_no_measure_mixed_target_total=0,
        iod_diagnostic_rows=Any[],
        iod_pairwise_rows=Any[],
        iod_geometry_cache_time_s=NaN,
        iod_geometry_cache=Dict{Int, Any}(),
        consensus_group_total=0,
        consensus_group_same_target_total=0,
        consensus_group_mixed_target_total=0,
        consensus_group_unknown_total=0
    )
end

function _ensure_fusion_slot_capacity!(
    model::DistributedFusionModel,
    slot_id::Int
)
    slot_id <= size(model.filter_initialized, 2) && return nothing

    n_rows = size(model.filter_initialized, 1)
    old_cols = size(model.filter_initialized, 2)
    new_cols = max(slot_id, old_cols + 128)

    new_iod_state = fill(NAN_STATE6, n_rows, new_cols)
    new_iod_state[:, 1:old_cols] = model.iod_estimate_state
    model.iod_estimate_state = new_iod_state

    new_iod_cov = [_nan_cov6() for _ in 1:n_rows, _ in 1:new_cols]
    new_iod_cov[:, 1:old_cols] = model.iod_estimate_covariance
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

    new_iod_pending = falses(n_rows, new_cols)
    new_iod_pending[:, 1:old_cols] = model.iod_pending
    model.iod_pending = new_iod_pending

    new_iod_pending_same = falses(n_rows, new_cols)
    new_iod_pending_same[:, 1:old_cols] = model.iod_pending_same_target
    model.iod_pending_same_target = new_iod_pending_same

    new_iod_pending_class = fill(IOD_CLASS_UNKNOWN, n_rows, new_cols)
    new_iod_pending_class[:, 1:old_cols] = model.iod_pending_class
    model.iod_pending_class = new_iod_pending_class

    new_iod_init = falses(n_rows, new_cols)
    new_iod_init[:, 1:old_cols] = model.iod_initialized
    model.iod_initialized = new_iod_init

    new_state = fill(NAN_STATE6, n_rows, new_cols)
    new_state[:, 1:old_cols] = model.state
    model.state = new_state

    new_cov = [_nan_cov6() for _ in 1:n_rows, _ in 1:new_cols]
    new_cov[:, 1:old_cols] = model.covariance
    model.covariance = new_cov

    new_state_pred = fill(NAN_STATE6, n_rows, new_cols)
    new_state_pred[:, 1:old_cols] = model.state_pred
    model.state_pred = new_state_pred

    new_cov_pred = [_nan_cov6() for _ in 1:n_rows, _ in 1:new_cols]
    new_cov_pred[:, 1:old_cols] = model.covariance_pred
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

function _reset_pending_iod!(
    model::DistributedFusionModel,
    observer_id::Int,
    slot_id::Int,
    track::LocalTrack
)::Nothing
    model.iod_pending[observer_id, slot_id] = false
    delete!(model.iod_pending_keys, (observer_id, slot_id))
    model.iod_pending_same_target[observer_id, slot_id] = false
    model.iod_pending_class[observer_id, slot_id] = IOD_CLASS_UNKNOWN
    model.iod_initialized[observer_id, slot_id] = false
    model.iod_triangulation_ready[observer_id, slot_id] = false
    model.iod_used_neighbors[observer_id, slot_id] = 0
    model.iod_estimate_state[observer_id, slot_id] = NAN_STATE6
    model.iod_estimate_covariance[observer_id, slot_id] = _nan_cov6()
    model.iod_estimate_time_s[observer_id, slot_id] = NaN
    track.status = :seed_ready
    track.iod_group_same_target = -1
    track.iod_group_class = IOD_CLASS_UNKNOWN
    track.state_estimate_now = NAN_STATE6
    track.covariance_estimate_now = _nan_cov6()
    empty!(model.buffer_matching_groups)
    empty!(model.nav.maintained_t2t_neighbors)
    return nothing
end

@inline function _iod_position_error_m(u, target_id::Int, state)::Float64
    target_id > 0 || return NaN
    _is_finite_state(state) || return NaN
    true_position = SVector{3, Float64}(u.sc[target_id].pos)
    return norm(SVector{3, Float64}(state[1], state[2], state[3]) - true_position)
end

function _iod_group_member_labels(
    nav_model::ObserverNavigationModel,
    nodes_idx::Vector{Tuple{Int, Int}}
)
    SAVE_IOD_EVENT_GEOMETRY || return (group_observers="", group_target_ids="")
    observers = sort!(unique(first.(nodes_idx)))
    target_ids = Int[]
    for (observer_id, slot_id) in nodes_idx
        track = get(nav_model.local_tracks[observer_id], slot_id, nothing)
        track === nothing && continue
        track.prev_meas === nothing || push!(target_ids, track.prev_meas.target)
        track.last_meas === nothing || push!(target_ids, track.last_meas.target)
    end
    sort!(unique!(target_ids))
    return (
        group_observers=join(observers, ';'),
        group_target_ids=join(target_ids, ';')
    )
end

function _record_iod_pairwise_consistency!(
    model::DistributedFusionModel;
    time_s::Float64,
    reference_observer::Int,
    reference_slot::Int,
    nodes_idx::Vector{Tuple{Int, Int}},
    iod_class::UInt8,
    covariance_gate_passed::Bool
)::Nothing
    SAVE_IOD_PAIRWISE_DIAGNOSTICS || return nothing
    length(nodes_idx) >= 2 || return nothing

    pair_rows = Any[]
    visibility_radius_m = model.nav.sensor.detection_range_m
    gate_m = IOD_NEIGHBOR_MISS_DISTANCE_MAX_M
    reference_key = (reference_observer, reference_slot)

    for first_idx in 1:(length(nodes_idx) - 1)
        first_key = nodes_idx[first_idx]
        first_track = get(model.nav.local_tracks[first_key[1]], first_key[2], nothing)
        first_track === nothing && continue
        first_track.last_meas === nothing && continue
        first_track.prev_meas === nothing && continue

        for second_idx in (first_idx + 1):length(nodes_idx)
            second_key = nodes_idx[second_idx]
            second_track = get(model.nav.local_tracks[second_key[1]], second_key[2], nothing)
            second_track === nothing && continue
            second_track.last_meas === nothing && continue
            second_track.prev_meas === nothing && continue

            first_meas = first_track.last_meas
            first_prev = first_track.prev_meas
            second_meas = second_track.last_meas
            second_prev = second_track.prev_meas
            miss_current_m = _los_ray_miss_distance(
                first_meas.observer_pos,
                first_meas.los_unit,
                visibility_radius_m,
                second_meas.observer_pos,
                second_meas.los_unit,
                visibility_radius_m
            )
            miss_previous_m = _los_ray_miss_distance(
                first_prev.observer_pos,
                first_prev.los_unit,
                visibility_radius_m,
                second_prev.observer_pos,
                second_prev.los_unit,
                visibility_radius_m
            )
            pair_passes = isfinite(miss_current_m) &&
                isfinite(miss_previous_m) &&
                miss_current_m <= gate_m &&
                miss_previous_m <= gate_m

            push!(pair_rows, (
                first_observer=first_key[1],
                first_slot=first_key[2],
                first_previous_target=first_prev.target,
                first_current_target=first_meas.target,
                second_observer=second_key[1],
                second_slot=second_key[2],
                second_previous_target=second_prev.target,
                second_current_target=second_meas.target,
                is_reference_pair=(first_key == reference_key || second_key == reference_key),
                pair_same_target=_same_target_iod_pair(
                    first_meas,
                    first_prev,
                    second_meas,
                    second_prev
                ),
                miss_current_m=miss_current_m,
                miss_previous_m=miss_previous_m,
                maximum_miss_m=max(miss_current_m, miss_previous_m),
                pair_passes=pair_passes
            ))
        end
    end

    isempty(pair_rows) && return nothing
    all_pairs_pass = all(row -> row.pair_passes, pair_rows)
    for row in pair_rows
        push!(model.iod_pairwise_rows, merge((
            nav_case=String(NAVIGATION_CASE),
            time_s=time_s,
            reference_observer=reference_observer,
            reference_slot=reference_slot,
            composition=_iod_class_label(iod_class),
            group_size=length(nodes_idx),
            pair_count=length(pair_rows),
            covariance_gate_passed=covariance_gate_passed,
            all_pairs_pass=all_pairs_pass,
            gate_m=gate_m
        ), row))
    end
    return nothing
end

function _iod_event_geometry(
    model::DistributedFusionModel,
    u,
    time_s::Float64,
    observer::Int,
    reference_target::Int
)
    SAVE_IOD_EVENT_GEOMETRY || return (
        visible_target_count=-1,
        minimum_angular_separation_deg=NaN,
        closest_target_1=0,
        closest_target_2=0,
        reference_minimum_angular_separation_deg=NaN,
        reference_closest_target=0,
        visible_target_ids="",
        minimum_3nn_radius_deg=NaN,
        maximum_3nn_density_targets_per_sr=NaN,
        reference_3nn_radius_deg=NaN,
        reference_3nn_density_targets_per_sr=NaN
    )

    if !isfinite(model.iod_geometry_cache_time_s) ||
       abs(model.iod_geometry_cache_time_s - time_s) > NAVIGATION_DT_TOL_SEC
        empty!(model.iod_geometry_cache)
        model.iod_geometry_cache_time_s = time_s
    end

    base_geometry = get!(model.iod_geometry_cache, observer) do
        sensor = model.nav.sensor
        _truth_geometry_base(u, sensor.target_idxs, sensor.detection_range_m, observer)
    end
    return _truth_geometry_for_reference(base_geometry, reference_target)
end

function _record_iod_diagnostic!(
    model::DistributedFusionModel;
    time_s::Float64,
    observer::Int,
    slot::Int,
    stage::String,
    outcome::String,
    same_target::Bool,
    reference_target::Int,
    composition::String=(same_target ? "same_target" : "mixed_real_target"),
    group_observers::String="",
    group_target_ids::String="",
    position_error_m::Float64=NaN,
    position_rms_std_m::Float64=NaN,
    validation_d2::Float64=NaN,
    truth_state=nothing,
    visible_target_count::Int=-1,
    minimum_angular_separation_deg::Float64=NaN,
    closest_target_1::Int=0,
    closest_target_2::Int=0,
    reference_minimum_angular_separation_deg::Float64=NaN,
    reference_closest_target::Int=0,
    visible_target_ids::String="",
    minimum_3nn_radius_deg::Float64=NaN,
    maximum_3nn_density_targets_per_sr::Float64=NaN,
    reference_3nn_radius_deg::Float64=NaN,
    reference_3nn_density_targets_per_sr::Float64=NaN
)::Nothing
    if truth_state !== nothing
        geometry = _iod_event_geometry(model, truth_state, time_s, observer, reference_target)
        visible_target_count = geometry.visible_target_count
        minimum_angular_separation_deg = geometry.minimum_angular_separation_deg
        closest_target_1 = geometry.closest_target_1
        closest_target_2 = geometry.closest_target_2
        reference_minimum_angular_separation_deg = geometry.reference_minimum_angular_separation_deg
        reference_closest_target = geometry.reference_closest_target
        visible_target_ids = geometry.visible_target_ids
        minimum_3nn_radius_deg = geometry.minimum_3nn_radius_deg
        maximum_3nn_density_targets_per_sr = geometry.maximum_3nn_density_targets_per_sr
        reference_3nn_radius_deg = geometry.reference_3nn_radius_deg
        reference_3nn_density_targets_per_sr = geometry.reference_3nn_density_targets_per_sr
    end
    push!(model.iod_diagnostic_rows, (
        nav_case=String(NAVIGATION_CASE), time_s=time_s, observer=observer, slot=slot,
        stage=stage, outcome=outcome, same_target=same_target,
        composition=composition, reference_target=reference_target, position_error_m=position_error_m,
        position_rms_std_m=position_rms_std_m, validation_d2=validation_d2,
        group_observers=group_observers, group_target_ids=group_target_ids,
        visible_target_count=visible_target_count,
        minimum_angular_separation_deg=minimum_angular_separation_deg,
        closest_target_1=closest_target_1, closest_target_2=closest_target_2,
        reference_minimum_angular_separation_deg=reference_minimum_angular_separation_deg,
        reference_closest_target=reference_closest_target,
        visible_target_ids=visible_target_ids,
        minimum_3nn_radius_deg=minimum_3nn_radius_deg,
        maximum_3nn_density_targets_per_sr=maximum_3nn_density_targets_per_sr,
        reference_3nn_radius_deg=reference_3nn_radius_deg,
        reference_3nn_density_targets_per_sr=reference_3nn_density_targets_per_sr
    ))
    return nothing
end

function _truth_geometry_base(u, target_idxs::Vector{Int}, detection_range_m::Float64, observer::Int)
    observer_pos = SVector{3, Float64}(u.sc[observer].pos)
    visible_targets = Int[]
    true_los = SVector{3, Float64}[]
    for target_id in target_idxs
        relative_position = SVector{3, Float64}(u.sc[target_id].pos) - observer_pos
        range_m = norm(relative_position)
        (isfinite(range_m) && range_m > 0.0 && range_m <= detection_range_m) || continue
        push!(visible_targets, target_id)
        push!(true_los, relative_position / range_m)
    end

    minimum_separation = Inf
    closest_target_1 = 0
    closest_target_2 = 0
    neighbor_angles = [Float64[] for _ in visible_targets]
    for first_idx in 1:(length(true_los) - 1)
        for second_idx in (first_idx + 1):length(true_los)
            separation = acos(clamp(dot(true_los[first_idx], true_los[second_idx]), -1.0, 1.0))
            push!(neighbor_angles[first_idx], separation)
            push!(neighbor_angles[second_idx], separation)
            if separation < minimum_separation
                minimum_separation = separation
                closest_target_1 = visible_targets[first_idx]
                closest_target_2 = visible_targets[second_idx]
            end
        end
    end

    third_neighbor_radius = Dict{Int, Float64}()
    maximum_3nn_density = NaN
    minimum_3nn_radius = Inf
    for (target_id, angles) in zip(visible_targets, neighbor_angles)
        length(angles) >= 3 || continue
        partialsort!(angles, 1:3)
        radius = angles[3]
        third_neighbor_radius[target_id] = radius
        solid_angle = 2pi * (1.0 - cos(radius))
        density = solid_angle > 0.0 ? 3.0 / solid_angle : Inf
        maximum_3nn_density = isfinite(maximum_3nn_density) ? max(maximum_3nn_density, density) : density
        minimum_3nn_radius = min(minimum_3nn_radius, radius)
    end

    return (
        visible_target_count=length(visible_targets),
        visible_target_ids=join(visible_targets, ';'),
        minimum_separation=minimum_separation,
        closest_target_1=closest_target_1,
        closest_target_2=closest_target_2,
        true_los=true_los,
        visible_targets=visible_targets,
        third_neighbor_radius=third_neighbor_radius,
        maximum_3nn_density=maximum_3nn_density,
        minimum_3nn_radius=minimum_3nn_radius
    )
end

function _truth_geometry_for_reference(base_geometry, reference_target::Int)
    reference_index = findfirst(==(reference_target), base_geometry.visible_targets)
    reference_minimum_separation = Inf
    reference_closest_target = 0
    if reference_index !== nothing
        for other_index in eachindex(base_geometry.visible_targets)
            other_index == reference_index && continue
            separation = acos(clamp(
                dot(base_geometry.true_los[reference_index], base_geometry.true_los[other_index]),
                -1.0,
                1.0
            ))
            if separation < reference_minimum_separation
                reference_minimum_separation = separation
                reference_closest_target = base_geometry.visible_targets[other_index]
            end
        end
    end
    reference_3nn_radius = get(base_geometry.third_neighbor_radius, reference_target, NaN)
    reference_solid_angle = isfinite(reference_3nn_radius) ?
        2pi * (1.0 - cos(reference_3nn_radius)) : NaN
    reference_3nn_density = isfinite(reference_solid_angle) && reference_solid_angle > 0.0 ?
        3.0 / reference_solid_angle : NaN

    return (
        visible_target_count=base_geometry.visible_target_count,
        minimum_angular_separation_deg=isfinite(base_geometry.minimum_separation) ?
            rad2deg(base_geometry.minimum_separation) : NaN,
        closest_target_1=base_geometry.closest_target_1,
        closest_target_2=base_geometry.closest_target_2,
        reference_minimum_angular_separation_deg=isfinite(reference_minimum_separation) ?
            rad2deg(reference_minimum_separation) : NaN,
        reference_closest_target=reference_closest_target,
        visible_target_ids=base_geometry.visible_target_ids,
        minimum_3nn_radius_deg=isfinite(base_geometry.minimum_3nn_radius) ?
            rad2deg(base_geometry.minimum_3nn_radius) : NaN,
        maximum_3nn_density_targets_per_sr=base_geometry.maximum_3nn_density,
        reference_3nn_radius_deg=isfinite(reference_3nn_radius) ? rad2deg(reference_3nn_radius) : NaN,
        reference_3nn_density_targets_per_sr=reference_3nn_density
    )
end

Base.@kwdef mutable struct IODTruthGeometryReplayModel
    observer_idxs::Vector{Int}
    target_idxs::Vector{Int}
    detection_range_m::Float64
    events_by_time_ms::Dict{Int, Vector{NamedTuple}}
    rows::Vector{Any}
end

function IODTruthGeometryReplayModel(
    input_path::String,
    observer_idxs::Vector{Int},
    target_idxs::Vector{Int},
    detection_range_m::Float64
)
    isfile(input_path) || error("Truth-geometry replay input not found: $(input_path)")
    input = CSV.read(input_path, DataFrame)
    events_by_time_ms = Dict{Int, Vector{NamedTuple}}()
    for row in eachrow(input)
        event = (
            event_id=Int(row.event_id),
            time_s=Float64(row.time_s),
            observer=Int(row.observer),
            reference_target=Int(row.reference_target)
        )
        push!(
            get!(events_by_time_ms, round(Int, 1_000.0 * event.time_s), NamedTuple[]),
            event
        )
    end
    return IODTruthGeometryReplayModel(
        observer_idxs=observer_idxs,
        target_idxs=target_idxs,
        detection_range_m=detection_range_m,
        events_by_time_ms=events_by_time_ms,
        rows=Any[]
    )
end

function SimulationModel.calcNavigationEffect!(
    model::IODTruthGeometryReplayModel,
    u,
    p,
    t::Float64,
    sat_idx::Int
)
    sat_idx == model.observer_idxs[end] || return nothing
    events = get(model.events_by_time_ms, round(Int, 1_000.0 * t), nothing)
    events === nothing && return nothing
    geometry_by_observer = Dict{Int, Any}()
    for event in events
        base_geometry = get!(geometry_by_observer, event.observer) do
            _truth_geometry_base(u, model.target_idxs, model.detection_range_m, event.observer)
        end
        geometry = _truth_geometry_for_reference(base_geometry, event.reference_target)
        push!(model.rows, merge(event, geometry))
    end
    return nothing
end

function _truth_geometry_replay_dataframe(model::IODTruthGeometryReplayModel)::DataFrame
    isempty(model.rows) && return DataFrame(
        event_id=Int[], time_s=Float64[], observer=Int[], reference_target=Int[],
        visible_target_count=Int[], minimum_angular_separation_deg=Float64[],
        closest_target_1=Int[], closest_target_2=Int[],
        reference_minimum_angular_separation_deg=Float64[], reference_closest_target=Int[],
        visible_target_ids=String[], minimum_3nn_radius_deg=Float64[],
        maximum_3nn_density_targets_per_sr=Float64[], reference_3nn_radius_deg=Float64[],
        reference_3nn_density_targets_per_sr=Float64[]
    )
    result = DataFrame(model.rows)
    sort!(result, :event_id)
    return result
end

function _validate_pending_iods!(model::DistributedFusionModel, t::Float64, u)::Nothing
    ENABLE_IOD_ONE_STEP_VALIDATION || return nothing
    isempty(model.iod_pending_keys) && return nothing

    # Iterate only over tentative IODs instead of scanning every local track.
    # The copy permits confirmed/rejected keys to be deleted during iteration.
    for (observer_id, slot_id) in collect(model.iod_pending_keys)
        tracks = model.nav.local_tracks[observer_id]
        track = get(tracks, slot_id, nothing)
        if track === nothing || slot_id > size(model.iod_pending, 2) ||
           !model.iod_pending[observer_id, slot_id]
            delete!(model.iod_pending_keys, (observer_id, slot_id))
            continue
        end

            iod_t = model.iod_estimate_time_s[observer_id, slot_id]
            if !isfinite(iod_t)
                _reset_pending_iod!(model, observer_id, slot_id, track)
                continue
            end
            dt_validation = t - iod_t
            dt_validation > NAVIGATION_DT_TOL_SEC || continue
            pending_same_target = model.iod_pending_same_target[observer_id, slot_id]
            pending_class = model.iod_pending_class[observer_id, slot_id]
            pending_composition = _iod_class_label(pending_class)

            measurement = track.last_meas
            has_new_measurement = measurement !== nothing &&
                track.has_measure_now &&
                abs(t - measurement.t) <= NAVIGATION_DT_TOL_SEC &&
                measurement.t > iod_t + NAVIGATION_DT_TOL_SEC

            if !has_new_measurement
                if dt_validation >= NAVIGATION_RATE_SEC - NAVIGATION_DT_TOL_SEC
                    model.iod_validation_no_measure_total += 1
                    if pending_same_target
                        model.iod_validation_no_measure_same_target_total += 1
                    else
                        model.iod_validation_no_measure_mixed_target_total += 1
                    end
                    reference_target = measurement === nothing ? 0 : measurement.target
                    _record_iod_diagnostic!(
                        model; time_s=t, observer=observer_id, slot=slot_id,
                        stage="next_step_validation", outcome="no_measurement",
                        same_target=pending_same_target, reference_target=reference_target,
                        composition=pending_composition, truth_state=u
                    )
                    _reset_pending_iod!(model, observer_id, slot_id, track)
                end
                continue
            end

            x_iod = model.iod_estimate_state[observer_id, slot_id]
            P_iod = model.iod_estimate_covariance[observer_id, slot_id]
            model.iod_validation_attempted_total += 1

            if !(_is_finite_state(x_iod) && _is_finite_cov(P_iod))
                model.iod_validation_rejected_total += 1
                if pending_same_target
                    model.iod_validation_rejected_same_target_total += 1
                else
                    model.iod_validation_rejected_mixed_target_total += 1
                end
                _record_iod_diagnostic!(
                    model; time_s=t, observer=observer_id, slot=slot_id,
                    stage="next_step_validation", outcome="rejected_invalid_state",
                    same_target=pending_same_target, reference_target=measurement.target,
                    composition=pending_composition, truth_state=u
                )
                _reset_pending_iod!(model, observer_id, slot_id, track)
                continue
            end

            x_validation = _propagate_keplerian(
                x_iod,
                model.μ,
                dt_validation,
                model.J2,
                model.R_ref_m
            )
            F_validation = _compute_process_jacobian(
                x_iod,
                model.μ,
                dt_validation,
                model.J2,
                model.R_ref_m
            )
            q_scale = max(dt_validation / NAVIGATION_RATE_SEC, 1.0)
            Q_validation = Matrix(Diagonal(q_scale .* model.process_q_diag))
            P_validation_raw = F_validation * P_iod * F_validation' + Q_validation
            P_validation = 0.5 * (P_validation_raw + P_validation_raw')
            d2_validation = _measurement_mahalanobis_d2_from_state_cov(
                model.nav,
                x_validation,
                P_validation,
                measurement.observer_pos,
                measurement
            )
            validation_position_error_m = _iod_position_error_m(
                u,
                measurement.target,
                x_validation
            )

            if isfinite(d2_validation) && d2_validation <= IOD_VALIDATION_MAHAL_MAX_D2
                model.iod_validation_confirmed_total += 1
                if pending_same_target
                    model.iod_validation_confirmed_same_target_total += 1
                else
                    model.iod_validation_confirmed_mixed_target_total += 1
                end

                model.iod_pending[observer_id, slot_id] = false
                delete!(model.iod_pending_keys, (observer_id, slot_id))
                model.iod_pending_same_target[observer_id, slot_id] = false
                model.iod_pending_class[observer_id, slot_id] = IOD_CLASS_UNKNOWN
                model.iod_initialized[observer_id, slot_id] = true
                model.iod_triangulation_ready[observer_id, slot_id] = true
                model.iod_estimate_state[observer_id, slot_id] = x_validation
                model.iod_estimate_covariance[observer_id, slot_id] = copy(P_validation)
                model.iod_estimate_time_s[observer_id, slot_id] = t
                empty!(model.buffer_matching_groups)
                empty!(model.nav.maintained_t2t_neighbors)

                track.status = :iod_initialized
                track.iod_group_same_target = pending_same_target ? 1 : 0
                track.iod_group_class = pending_class
                if track.plot_target_id == 0
                    track.plot_target_id = measurement.target
                end
                track.state_estimate_now = x_validation
                track.covariance_estimate_now = copy(P_validation)

                model.iod_init_total += 1
                if pending_same_target
                    model.iod_init_same_target_total += 1
                else
                    model.iod_init_mixed_target_total += 1
                end
                _record_iod_diagnostic!(
                    model; time_s=t, observer=observer_id, slot=slot_id,
                    stage="next_step_validation", outcome="promoted",
                    same_target=pending_same_target, reference_target=measurement.target,
                    composition=pending_composition,
                    position_error_m=validation_position_error_m,
                    validation_d2=d2_validation, truth_state=u
                )
                LOG_NAV_EVENTS && println(
                    "IOD validation confirmed | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s | d2=$(round(d2_validation; digits=3))"
                )
            else
                model.iod_validation_rejected_total += 1
                if pending_same_target
                    model.iod_validation_rejected_same_target_total += 1
                else
                    model.iod_validation_rejected_mixed_target_total += 1
                end
                _record_iod_diagnostic!(
                    model; time_s=t, observer=observer_id, slot=slot_id,
                    stage="next_step_validation", outcome="rejected_score",
                    same_target=pending_same_target, reference_target=measurement.target,
                    composition=pending_composition,
                    position_error_m=validation_position_error_m,
                    validation_d2=d2_validation, truth_state=u
                )
                LOG_NAV_EVENTS && println(
                    "IOD validation rejected | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s | d2=$(round(d2_validation; digits=3))"
                )
                _reset_pending_iod!(model, observer_id, slot_id, track)
            end
    end
    return nothing
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
        if model.filter_initialized[observer_id, slot_id]
            tracks = model.nav.local_tracks[observer_id]
            track = tracks[slot_id]
            track.is_freshly_initialized = false
            continue
        end
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
        track.status = :filter_initialized
        track.is_freshly_initialized = true
        if !isfinite(track.filter_initialized_t)
            track.filter_initialized_t = t
            target_id = track.last_meas === nothing ? 0 : track.last_meas.target
            if target_id > 0
                true_pos = SVector{3, Float64}(u.sc[target_id].pos)
                track.initialization_position_error_m = norm(
                    SVector{3, Float64}(x_init[1], x_init[2], x_init[3]) - true_pos
                )
            end
        end
        if track.plot_target_id == 0 && (track.last_meas !== nothing)
            # Safety: set the plot label once if this slot had no prior IOD label.
            track.plot_target_id = track.last_meas.target
        end
        track.consecutive_missed = 0
        # Keep local track estimate aligned to x^- for next-tick local gating.
        track.state_estimate_now = model.state_pred[observer_id, slot_id]
        track.covariance_estimate_now = copy(model.covariance_pred[observer_id, slot_id])

        LOG_NAV_EVENTS && println("DEKF init (IOD+consensus) | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s")
    end
    return nothing
end

function _retire_stale_filter_tracks!(model::DistributedFusionModel, t::Float64)
    cap = size(model.filter_initialized, 2)
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
        closed_slots = Int[]
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
                            LOG_NAV_EVENTS && println(
                                "Solo-target-measure streak | target=$(target_id) | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s"
                            )
                        end
                        if streak >= SOLO_TARGET_MEASURE_STREAK_CLOSE_STEPS
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
                track.status = :closed
                track.closed_t = t
                push!(closed_slots, slot_id)
                if close_reason === :missed
                    LOG_NAV_EVENTS && println("Local track closed | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s | missed=$(track.consecutive_missed) | reason=missed_measurements_threshold")
                else
                    target_id = (track.last_meas === nothing) ? 0 : track.last_meas.target
                    LOG_NAV_EVENTS && println("Local track closed | target=$(target_id) | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s | solo_streak=$(close_streak) | reason=solo_target_measure_streak_threshold")
                end
            end
        end
        for slot_id in closed_slots
            track = tracks[slot_id]
            push!(model.nav.closed_track_lifecycle[observer_id], ClosedTrackLifecycle(track))
            delete!(tracks, slot_id)
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
    cap = size(model.filter_initialized, 2)
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
    cap = size(model.filter_initialized, 2)
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

@inline function check_trigger_consensus_groups_refresh(
    model::DistributedFusionModel,
    active_tracks::Vector{Tuple{Int, Int}}
)::Bool
    # Trigger rebuild when a local track lifecycle changed in this tick.
    for obs in model.observer_idxs
        tracks = model.nav.local_tracks[obs]
        for (_slot_id, track) in tracks
            if track.status == :filter_initialized && track.is_freshly_initialized
                # Consume one-shot init trigger.
                track.is_freshly_initialized = false
                return true
            end
        end
    end

    # Trigger rebuild when active-track membership changed vs buffered groups.
    active_set = Set(active_tracks)
    buffered_set = Set{Tuple{Int, Int}}()
    for group in model.buffer_matching_groups
        for key in group.tracks
            push!(buffered_set, key)
        end
    end
    active_set == buffered_set || return true

    # Trigger rebuild when a retained T2T link is no longer self-consistent.
    for group in model.buffer_matching_groups
        for edge in group.retained_edges
            (obs_i, slot_i) = edge.first
            mode_i = _track_match_mode(model, obs_i, slot_i)
            mode_i === nothing && return true

            (obs_j, slot_j) = edge.second
            mode_j = _track_match_mode(model, obs_j, slot_j)
            mode_j === nothing && return true

            # Keep the same semantics used in pairwise TT matching.
            mode = (mode_i === :iod || mode_j === :iod) ? :iod : :filter
            d2_gate = mode === :iod ? CONSENSUS_MATCH_MAHAL_IOD_MAX_D2 : CONSENSUS_MATCH_MAHAL_FILTER_MAX_D2

            state_cov_i = _track_state_cov_for_matching(model, obs_i, slot_i, mode)
            state_cov_j = _track_state_cov_for_matching(model, obs_j, slot_j, mode)
            (state_cov_i === nothing || state_cov_j === nothing) && return true

            x_i, P_i = state_cov_i
            x_j, P_j = state_cov_j
            d2 = _mahalanobis_distance_sq(x_i - x_j, P_i + P_j)
            if !isfinite(d2) || d2 > d2_gate
                return true
            end
        end
    end

    return false
end

@inline function _build_local_information_terms(
    model,
    x_prior::SVector{6, Float64},
    observer_pos::SVector{3, Float64},
    los_meas::SVector{3, Float64}
)::Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}}
    # Build local information-form measurement terms:
    #   u = H'R^-1 z
    #   uhat = H'R^-1 h(x)
    #   U = H'R^-1 H
    H = _measurement_jacobian(x_prior, observer_pos)
    zero_vec = zeros(6)
    zero_mat = zeros(6, 6)
    all(isfinite, los_meas) || return (zero_vec, zero_vec, zero_mat)
    los_unit = _safe_unit(los_meas)

    I3 = Matrix(I, 3, 3)
    R_i = (model.sigma_theta_rad^2) * (I3 - los_unit * los_unit')
    R_inv = pinv(R_i)
    all(isfinite, R_inv) || return (zero_vec, zero_vec, zero_mat)
    u_local = H' * (R_inv * los_meas)
    all(isfinite, u_local) || return (zero_vec, zero_vec, zero_mat)

    los_pred = _safe_unit(x_prior[1:3] - observer_pos)
    uhat_local = H' * (R_inv * los_pred)
    U_local = H' * (R_inv * H)
    if !(all(isfinite, uhat_local) && all(isfinite, U_local))
        return (zero_vec, zero_vec, zero_mat)
    end
    return (u_local, uhat_local, U_local)
end

Base.@kwdef mutable struct CentralizedOracleFusionModel
    sensor::OpticalLOSSensorModel
    observer_idxs::Vector{Int}
    target_idxs::Vector{Int}
    od_perturb_enabled::Bool
    od_pos_std_m::Float64
    od_vel_std_mps::Float64
    od_rng::MersenneTwister
    known_observer_state::Vector{SVector{6, Float64}}
    known_observer_state_t::Vector{Float64}
    sigma_theta_rad::Float64
    process_q_diag::SVector{6, Float64}
    μ::Float64
    J2::Float64
    R_ref_m::Float64
    state::Dict{Int, SVector{6, Float64}}
    covariance::Dict{Int, Matrix{Float64}}
    state_pred::Dict{Int, SVector{6, Float64}}
    covariance_pred::Dict{Int, Matrix{Float64}}
    initialized::Dict{Int, Bool}
    last_update_t::Dict{Int, Float64}
    first_measurement_t::Dict{Int, Float64}
    initialization_t::Dict{Int, Float64}
    initialization_position_error_m::Dict{Int, Float64}
end

function CentralizedOracleFusionModel(
    sensor::OpticalLOSSensorModel,
    observer_idxs::Vector{Int},
    target_idxs::Vector{Int},
    sigma_theta_rad::Float64,
    process_q_diag::SVector{6, Float64},
    μ::Float64,
    J2::Float64,
    R_ref_m::Float64;
    od_perturb_enabled::Bool=ENABLE_OBSERVER_OD_PERTURBATION,
    od_pos_std_m::Float64=OBSERVER_OD_POS_STD_M,
    od_vel_std_mps::Float64=OBSERVER_OD_VEL_STD_MPS,
    od_rng_seed::Int=OBSERVER_OD_RNG_SEED
)
    num_sats = maximum(vcat(observer_idxs, target_idxs))
    return CentralizedOracleFusionModel(
        sensor=sensor,
        observer_idxs=observer_idxs,
        target_idxs=target_idxs,
        od_perturb_enabled=od_perturb_enabled,
        od_pos_std_m=od_pos_std_m,
        od_vel_std_mps=od_vel_std_mps,
        od_rng=MersenneTwister(od_rng_seed),
        known_observer_state=fill(NAN_STATE6, num_sats),
        known_observer_state_t=fill(NaN, num_sats),
        sigma_theta_rad=sigma_theta_rad,
        process_q_diag=process_q_diag,
        μ=μ,
        J2=J2,
        R_ref_m=R_ref_m,
        state=Dict(target_id => NAN_STATE6 for target_id in target_idxs),
        covariance=Dict(target_id => _nan_cov6() for target_id in target_idxs),
        state_pred=Dict(target_id => NAN_STATE6 for target_id in target_idxs),
        covariance_pred=Dict(target_id => _nan_cov6() for target_id in target_idxs),
        initialized=Dict(target_id => false for target_id in target_idxs),
        last_update_t=Dict(target_id => NaN for target_id in target_idxs),
        first_measurement_t=Dict(target_id => NaN for target_id in target_idxs),
        initialization_t=Dict(target_id => NaN for target_id in target_idxs),
        initialization_position_error_m=Dict(target_id => NaN for target_id in target_idxs)
    )
end

function SimulationModel.calcNavigationEffect!(
    model::CentralizedOracleFusionModel,
    u,
    p,
    t::Float64,
    sat_idx::Int
)
    sat_idx == model.observer_idxs[end] || return nothing
    Q = Matrix(Diagonal(model.process_q_diag))

    for observer_id in model.observer_idxs
        _known_observer_state!(model, observer_id, u, t)
    end

    for target_id in model.target_idxs
        target_measurements = LOSMeasurement[]
        for observer_id in model.observer_idxs
            for measurement in model.sensor.measurements_now[observer_id]
                measurement.target == target_id || continue
                measurement.observer_pos = _known_observer_pos!(model, observer_id, u, t)
                push!(target_measurements, measurement)
            end
        end

        if !isempty(target_measurements) && !isfinite(model.first_measurement_t[target_id])
            model.first_measurement_t[target_id] = t
        end

        if !model.initialized[target_id]
            isempty(target_measurements) && continue

            nodes = NamedTuple[]
            A_rows = Matrix{Float64}(undef, 0, 3)
            b_rows = Float64[]
            dA_rows = Matrix{Float64}(undef, 0, 3)
            db_rows = Float64[]

            for measurement in target_measurements
                model.sensor.consecutive_counts[measurement.observer, target_id] >= LOCAL_INIT_MIN_MEASUREMENTS || continue
                prev_measurement = get(model.sensor.previous, (measurement.observer, target_id), nothing)
                prev_measurement === nothing && continue
                _is_consecutive_measure_pair(measurement.t, prev_measurement.t) || continue

                los = _safe_unit(measurement.los_unit)
                los_prev = _safe_unit(prev_measurement.los_unit)
                norm(los) > 1e-12 && norm(los_prev) > 1e-12 || continue
                los_rate = (los - los_prev) / (measurement.t - prev_measurement.t)
                observer_pos = measurement.observer_pos
                observer_vel = _known_observer_vel!(model, measurement.observer, u, t)
                all(isfinite, observer_pos) && all(isfinite, observer_vel) || continue

                eq = _build_agent_equations(observer_pos, observer_vel, los, los_rate)
                eq === nothing && continue
                A_rows = vcat(A_rows, eq.A)
                append!(b_rows, eq.b)
                dA_rows = vcat(dA_rows, eq.dA)
                append!(db_rows, eq.db)
                push!(nodes, (
                    r=observer_pos,
                    v=observer_vel,
                    los=los,
                    los_prev=los_prev,
                    los_rate=los_rate
                ))
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

            x0 = SVector{6, Float64}(z_est)
            model.state[target_id] = x0
            model.covariance[target_id] = copy(cov)
            model.state_pred[target_id] = x0
            model.covariance_pred[target_id] = copy(cov)
            model.initialized[target_id] = true
            model.initialization_t[target_id] = t
            true_position = SVector{3, Float64}(u.sc[target_id].pos)
            model.initialization_position_error_m[target_id] = norm(x_est - true_position)
        end

        x_prior = model.state_pred[target_id]
        P_prior = model.covariance_pred[target_id]
        if !(_is_finite_state(x_prior) && _is_finite_cov(P_prior))
            x_prior = model.state[target_id]
            P_prior = model.covariance[target_id]
        end
        (_is_finite_state(x_prior) && _is_finite_cov(P_prior)) || continue

        U_sum = zeros(6, 6)
        residual_info = zeros(6)
        for measurement in target_measurements
            observer_pos = measurement.observer_pos
            all(isfinite, observer_pos) || continue
            u_i, uhat_i, U_i = _build_local_information_terms(
                model,
                x_prior,
                observer_pos,
                measurement.los_unit
            )
            residual_info += (u_i - uhat_i)
            U_sum += U_i
        end

        if isempty(target_measurements)
            x_upd = x_prior
            P_upd = P_prior
        else
            Y_prior = pinv(P_prior)
            Y_post = Y_prior + U_sum
            P_upd = _symmetrize_psd(pinv(Y_post))
            dx = P_upd * residual_info
            x_upd = SVector{6, Float64}(Vector{Float64}(x_prior) + dx)
        end

        F = _compute_process_jacobian(x_upd, model.μ, NAVIGATION_RATE_SEC, model.J2, model.R_ref_m)
        x_pred = _propagate_keplerian(x_upd, model.μ, NAVIGATION_RATE_SEC, model.J2, model.R_ref_m)
        P_pred = _symmetrize_psd(F * P_upd * F' + Q)

        model.state[target_id] = x_upd
        model.covariance[target_id] = P_upd
        model.state_pred[target_id] = x_pred
        model.covariance_pred[target_id] = P_pred
        model.last_update_t[target_id] = t
    end
    return nothing
end

@inline function _fusion_slot_for_save(
    model::CentralizedOracleFusionModel,
    obs::Int,
    target_id::Int,
    t::Float64
)::Int
    get(model.initialized, target_id, false) || return 0
    t_last = get(model.last_update_t, target_id, NaN)
    (!isfinite(t_last) || abs(t - t_last) > NAVIGATION_DT_TOL_SEC) && return 0
    return target_id
end

@inline function _fusion_state_for_save(
    model::CentralizedOracleFusionModel,
    obs::Int,
    target_id::Int,
    t::Float64
)::SVector{6, Float64}
    _fusion_slot_for_save(model, obs, target_id, t) == 0 && return NAN_STATE6
    return model.state[target_id]
end

Base.@kwdef mutable struct TrackingWindowMetricsModel
    sensor::OpticalLOSSensorModel
    estimator
    observer_idxs::Vector{Int}
    target_idxs::Vector{Int}
    joint_visible_prev::BitMatrix
    visible_observer_count::Vector{Int}
    current_window_id::Matrix{Int}
    current_window_start_t::Matrix{Float64}
    current_window_ticks::Matrix{Int}
    current_window_tracked::BitMatrix
    current_error_sum_m::Matrix{Float64}
    current_error_sq_sum_m2::Matrix{Float64}
    current_error_count::Matrix{Int}
    current_error_max_m::Matrix{Float64}
    current_error_samples::Matrix{Vector{Tuple{Int, Float64}}}
    next_window_id::Int
    window_rows::Vector{Any}
    geometry_rows::Vector{Any}
end

function TrackingWindowMetricsModel(
    sensor::OpticalLOSSensorModel,
    estimator,
    observer_idxs::Vector{Int},
    target_idxs::Vector{Int},
    num_sats::Int
)
    return TrackingWindowMetricsModel(
        sensor=sensor,
        estimator=estimator,
        observer_idxs=observer_idxs,
        target_idxs=target_idxs,
        joint_visible_prev=falses(num_sats, num_sats),
        visible_observer_count=zeros(Int, num_sats),
        current_window_id=zeros(Int, num_sats, num_sats),
        current_window_start_t=fill(NaN, num_sats, num_sats),
        current_window_ticks=zeros(Int, num_sats, num_sats),
        current_window_tracked=falses(num_sats, num_sats),
        current_error_sum_m=zeros(Float64, num_sats, num_sats),
        current_error_sq_sum_m2=zeros(Float64, num_sats, num_sats),
        current_error_count=zeros(Int, num_sats, num_sats),
        current_error_max_m=zeros(Float64, num_sats, num_sats),
        current_error_samples=[Tuple{Int, Float64}[] for _ in 1:num_sats, _ in 1:num_sats],
        next_window_id=1,
        window_rows=Any[],
        geometry_rows=Any[]
    )
end

function _start_tracking_window!(model::TrackingWindowMetricsModel, obs::Int, target::Int, t::Float64)
    model.current_window_id[obs, target] = model.next_window_id
    model.next_window_id += 1
    model.current_window_start_t[obs, target] = t
    model.current_window_ticks[obs, target] = 0
    model.current_window_tracked[obs, target] = false
    model.current_error_sum_m[obs, target] = 0.0
    model.current_error_sq_sum_m2[obs, target] = 0.0
    model.current_error_count[obs, target] = 0
    model.current_error_max_m[obs, target] = 0.0
    empty!(model.current_error_samples[obs, target])
    return nothing
end

function _close_tracking_window!(model::TrackingWindowMetricsModel, obs::Int, target::Int)
    ticks = model.current_window_ticks[obs, target]
    if ticks >= TRACKING_POSSIBLE_MIN_JOINT_TICKS
        err_count = model.current_error_count[obs, target]
        mean_err = err_count > 0 ? model.current_error_sum_m[obs, target] / err_count : NaN
        rmse_err = err_count > 0 ? sqrt(model.current_error_sq_sum_m2[obs, target] / err_count) : NaN
        tracked = model.current_window_tracked[obs, target]
        samples = model.current_error_samples[obs, target]
        convergence_tick = max(1, ceil(Int, 0.2 * ticks))
        converged_errors = [error_m for (tick, error_m) in samples if tick >= convergence_tick]
        converged_count = length(converged_errors)
        converged_mean = converged_count > 0 ? mean(converged_errors) : NaN
        converged_rmse = converged_count > 0 ? sqrt(mean(abs2, converged_errors)) : NaN
        start_t = model.current_window_start_t[obs, target]
        end_t = isfinite(start_t) ? start_t + (ticks - 1) * NAVIGATION_RATE_SEC : NaN
        push!(
            model.window_rows,
            (
                nav_case=String(NAVIGATION_CASE),
                observer=obs,
                target=target,
                window_id=model.current_window_id[obs, target],
                window_start_t_s=start_t,
                window_end_t_s=end_t,
                joint_ticks=ticks,
                joint_duration_s=ticks * NAVIGATION_RATE_SEC,
                tracked=tracked,
                estimate_samples=err_count,
                estimate_duration_s=err_count * NAVIGATION_RATE_SEC,
                mean_error_m=mean_err,
                rmse_error_m=rmse_err,
                max_error_m=err_count > 0 ? model.current_error_max_m[obs, target] : NaN,
                first_estimate_error_m=isempty(samples) ? NaN : samples[1][2],
                converged_samples=converged_count,
                converged_mean_error_m=converged_mean,
                converged_rmse_error_m=converged_rmse,
                success_under_1km=tracked && err_count > 0 && mean_err <= TRACKING_SUCCESS_ERROR_MAX_M
            )
        )
    end

    model.current_window_id[obs, target] = 0
    model.current_window_start_t[obs, target] = NaN
    model.current_window_ticks[obs, target] = 0
    model.current_window_tracked[obs, target] = false
    model.current_error_sum_m[obs, target] = 0.0
    model.current_error_sq_sum_m2[obs, target] = 0.0
    model.current_error_count[obs, target] = 0
    model.current_error_max_m[obs, target] = 0.0
    empty!(model.current_error_samples[obs, target])
    return nothing
end

function SimulationModel.calcNavigationEffect!(
    model::TrackingWindowMetricsModel,
    u,
    p,
    t::Float64,
    sat_idx::Int
)
    sat_idx == model.observer_idxs[end] || return nothing

    # Reuse the pre-misdetection geometric visibility already evaluated by the
    # LOS sensor. A missed detection therefore cannot split a truth window, and
    # range calculations are not repeated solely for metric collection.
    visible = model.sensor.visible_now
    fill!(model.visible_observer_count, 0)
    for target_id in model.target_idxs
        count_visible = 0
        for observer_id in model.observer_idxs
            visible[observer_id, target_id] && (count_visible += 1)
        end
        model.visible_observer_count[target_id] = count_visible
    end

    # Detailed angular-crowding data are needed by the nominal analysis only.
    # Stress runs save compact aggregate metrics and skip this per-epoch geometry.
    if SAVE_AUXILIARY_METRIC_TABLES
        for observer_id in model.observer_idxs
            observer_pos = SVector{3, Float64}(u.sc[observer_id].pos)
            visible_targets = [
                target_id for target_id in model.target_idxs
                if visible[observer_id, target_id]
            ]
            true_los = SVector{3, Float64}[]
            sizehint!(true_los, length(visible_targets))
            for target_id in visible_targets
                relative_position = SVector{3, Float64}(u.sc[target_id].pos) - observer_pos
                range_m = norm(relative_position)
                push!(true_los, relative_position / range_m)
            end
            min_separation_rad = NaN
            closest_first = 0
            closest_second = 0
            if length(true_los) >= 2
                for first_idx in 1:(length(true_los) - 1)
                    for second_idx in (first_idx + 1):length(true_los)
                        separation = acos(clamp(dot(true_los[first_idx], true_los[second_idx]), -1.0, 1.0))
                        if !isfinite(min_separation_rad) || separation < min_separation_rad
                            min_separation_rad = separation
                            closest_first = visible_targets[first_idx]
                            closest_second = visible_targets[second_idx]
                        end
                    end
                end
            end
            push!(model.geometry_rows, (
                nav_case=String(NAVIGATION_CASE),
                time_s=t,
                observer=observer_id,
                simultaneously_visible_targets=length(visible_targets),
                minimum_angular_separation_deg=isfinite(min_separation_rad) ? rad2deg(min_separation_rad) : NaN,
                closest_target_1=closest_first,
                closest_target_2=closest_second
            ))
        end
    end

    for observer_id in model.observer_idxs
        for target_id in model.target_idxs
            jointly_visible = visible[observer_id, target_id] &&
                model.visible_observer_count[target_id] >= 2

            if jointly_visible
                if !model.joint_visible_prev[observer_id, target_id]
                    _start_tracking_window!(model, observer_id, target_id, t)
                end
                model.current_window_ticks[observer_id, target_id] += 1

                slot = _fusion_slot_for_save(model.estimator, observer_id, target_id, t)
                if slot != 0
                    x_est = _fusion_state_for_save(model.estimator, observer_id, target_id, t)
                    if _is_finite_state(x_est)
                        true_pos = SVector{3, Float64}(u.sc[target_id].pos)
                        err_m = norm(SVector{3, Float64}(x_est[1], x_est[2], x_est[3]) - true_pos)
                        if isfinite(err_m)
                            model.current_window_tracked[observer_id, target_id] = true
                            model.current_error_sum_m[observer_id, target_id] += err_m
                            model.current_error_sq_sum_m2[observer_id, target_id] += err_m^2
                            model.current_error_count[observer_id, target_id] += 1
                            push!(
                                model.current_error_samples[observer_id, target_id],
                                (model.current_window_ticks[observer_id, target_id], err_m)
                            )
                            model.current_error_max_m[observer_id, target_id] = max(
                                model.current_error_max_m[observer_id, target_id],
                                err_m
                            )
                        end
                    end
                end
            elseif model.joint_visible_prev[observer_id, target_id]
                _close_tracking_window!(model, observer_id, target_id)
            end

            model.joint_visible_prev[observer_id, target_id] = jointly_visible
        end
    end

    return nothing
end

function _finalize_tracking_windows!(model::TrackingWindowMetricsModel)
    for observer_id in model.observer_idxs
        for target_id in model.target_idxs
            model.joint_visible_prev[observer_id, target_id] || continue
            _close_tracking_window!(model, observer_id, target_id)
            model.joint_visible_prev[observer_id, target_id] = false
        end
    end
    return nothing
end

function _tracking_window_dataframe(model::TrackingWindowMetricsModel)::DataFrame
    if isempty(model.window_rows)
        return DataFrame(
            nav_case=String[],
            observer=Int[],
            target=Int[],
            window_id=Int[],
            window_start_t_s=Float64[],
            window_end_t_s=Float64[],
            joint_ticks=Int[],
            joint_duration_s=Float64[],
            tracked=Bool[],
            estimate_samples=Int[],
            estimate_duration_s=Float64[],
            mean_error_m=Float64[],
            rmse_error_m=Float64[],
            max_error_m=Float64[],
            first_estimate_error_m=Float64[],
            converged_samples=Int[],
            converged_mean_error_m=Float64[],
            converged_rmse_error_m=Float64[],
            success_under_1km=Bool[]
        )
    end
    return DataFrame(model.window_rows)
end

function _tracking_observer_summary(model::TrackingWindowMetricsModel)::DataFrame
    windows = _tracking_window_dataframe(model)
    rows = NamedTuple[]
    for observer_id in model.observer_idxs
        obs_windows = windows[windows.observer .== observer_id, :]
        possible = nrow(obs_windows)
        tracked = possible > 0 ? count(obs_windows.tracked) : 0
        successful = possible > 0 ? count(obs_windows.success_under_1km) : 0
        detected_unique_targets = count(
            target_id -> model.sensor.counts[observer_id, target_id] > 0,
            model.target_idxs
        )
        possible_unique_targets = length(unique([Int(row.target) for row in eachrow(obs_windows)]))
        tracked_unique_targets = length(unique([
            Int(row.target) for row in eachrow(obs_windows)
            if row.tracked
        ]))
        successful_unique_targets = length(unique([
            Int(row.target) for row in eachrow(obs_windows)
            if row.success_under_1km
        ]))
        error_values = [
            Float64(row.mean_error_m) for row in eachrow(obs_windows)
            if row.tracked && isfinite(Float64(row.mean_error_m))
        ]
        mean_error = isempty(error_values) ? NaN : sum(error_values) / length(error_values)
        sq_error_sum = 0.0
        sq_error_count = 0
        for row in eachrow(obs_windows)
            row.tracked || continue
            samples = Int(row.estimate_samples)
            rmse = Float64(row.rmse_error_m)
            if samples > 0 && isfinite(rmse)
                sq_error_sum += rmse^2 * samples
                sq_error_count += samples
            end
        end
        rmse_error = sq_error_count > 0 ? sqrt(sq_error_sum / sq_error_count) : NaN
        good_error_values = [
            Float64(row.mean_error_m) for row in eachrow(obs_windows)
            if row.success_under_1km && isfinite(Float64(row.mean_error_m))
        ]
        mean_good_error = isempty(good_error_values) ? NaN : sum(good_error_values) / length(good_error_values)
        push!(
            rows,
            (
                nav_case=String(NAVIGATION_CASE),
                observer=observer_id,
                possible_windows=Float64(possible),
                tracked_windows=Float64(tracked),
                tracking_coverage_pct=_pct(Float64(tracked), Float64(possible)),
                successful_windows_under_1km=Float64(successful),
                success_rate_tracked_pct=_pct(Float64(successful), Float64(tracked)),
                success_rate_possible_pct=_pct(Float64(successful), Float64(possible)),
                detected_unique_targets=Float64(detected_unique_targets),
                possible_unique_targets=Float64(possible_unique_targets),
                tracked_unique_targets=Float64(tracked_unique_targets),
                successful_unique_targets=Float64(successful_unique_targets),
                mean_error_tracked_windows_m=mean_error,
                rmse_error_tracked_windows_m=rmse_error,
                mean_error_successful_windows_m=mean_good_error
            )
        )
    end
    return DataFrame(rows)
end

@inline _track_final_status(track::LocalTrack)::Symbol = track.status
@inline _track_final_status(track::ClosedTrackLifecycle)::Symbol = track.final_status

function _track_lifecycle_row(
    observer_id::Int,
    track::Union{LocalTrack, ClosedTrackLifecycle},
    final_t::Float64
)
    end_t = isfinite(track.closed_t) ? track.closed_t : final_t
    duration_s = (isfinite(track.created_t) && isfinite(end_t)) ?
        max(0.0, end_t - track.created_t) : NaN
    filter_duration_s = (isfinite(track.filter_initialized_t) && isfinite(end_t)) ?
        max(0.0, end_t - track.filter_initialized_t) : NaN
    identity_class = if !isfinite(track.filter_initialized_t)
        "uninitialized"
    elseif track.iod_group_same_target == 0 || track.id_switch_count > 0
        "bad"
    elseif track.iod_group_same_target == 1
        "good"
    else
        "unknown"
    end
    return (
        nav_case=String(NAVIGATION_CASE),
        observer=observer_id,
        slot=track.slot,
        final_status=String(_track_final_status(track)),
        created_t_s=track.created_t,
        first_measurement_t_s=track.first_measurement_t,
        filter_initialized_t_s=track.filter_initialized_t,
        initialization_latency_s=(
            isfinite(track.first_measurement_t) && isfinite(track.filter_initialized_t)
        ) ? max(0.0, track.filter_initialized_t - track.first_measurement_t) : NaN,
        initialization_position_error_m=track.initialization_position_error_m,
        iod_group_same_target=track.iod_group_same_target,
        iod_group_composition=_iod_class_label(track.iod_group_class),
        identity_class=identity_class,
        closed_t_s=track.closed_t,
        lifecycle_end_t_s=end_t,
        duration_s=duration_s,
        filter_duration_s=filter_duration_s,
        first_target_id=track.first_target_id,
        last_target_id=track.last_target_id,
        id_switch_count=track.id_switch_count,
        measurement_update_count=track.measurement_update_count
    )
end

function _track_lifecycle_dataframe(model::ObserverNavigationModel, final_t::Float64)::DataFrame
    rows = NamedTuple[]
    for observer_id in model.observer_idxs
        for track in model.closed_track_lifecycle[observer_id]
            push!(rows, _track_lifecycle_row(observer_id, track, final_t))
        end
        tracks = model.local_tracks[observer_id]
        for slot_id in sort(collect(keys(tracks)))
            push!(rows, _track_lifecycle_row(observer_id, tracks[slot_id], final_t))
        end
    end

    if isempty(rows)
        return DataFrame(
            nav_case=String[],
            observer=Int[],
            slot=Int[],
            final_status=String[],
            created_t_s=Float64[],
            first_measurement_t_s=Float64[],
            filter_initialized_t_s=Float64[],
            initialization_latency_s=Float64[],
            initialization_position_error_m=Float64[],
            iod_group_same_target=Int[],
            iod_group_composition=String[],
            identity_class=String[],
            closed_t_s=Float64[],
            lifecycle_end_t_s=Float64[],
            duration_s=Float64[],
            filter_duration_s=Float64[],
            first_target_id=Int[],
            last_target_id=Int[],
            id_switch_count=Int[],
            measurement_update_count=Int[]
        )
    end
    return DataFrame(rows)
end

function _track_lifecycle_dataframe(
    model::CentralizedOracleFusionModel,
    final_t::Float64
)::DataFrame
    rows = NamedTuple[]
    for target_id in model.target_idxs
        first_t = model.first_measurement_t[target_id]
        initialized_t = model.initialization_t[target_id]
        is_initialized = model.initialized[target_id] && isfinite(initialized_t)
        (isfinite(first_t) || is_initialized) || continue
        lifecycle_start_t = isfinite(first_t) ? first_t : initialized_t
        push!(
            rows,
            (
                nav_case=String(NAVIGATION_CASE),
                observer=0,
                slot=target_id,
                final_status=is_initialized ? "filter_initialized" : "uninitialized",
                created_t_s=lifecycle_start_t,
                first_measurement_t_s=first_t,
                filter_initialized_t_s=initialized_t,
                initialization_latency_s=(
                    isfinite(first_t) && isfinite(initialized_t)
                ) ? max(0.0, initialized_t - first_t) : NaN,
                initialization_position_error_m=
                    model.initialization_position_error_m[target_id],
                iod_group_same_target=is_initialized ? 1 : -1,
                iod_group_composition=is_initialized ? "same_target" : "unknown",
                identity_class=is_initialized ? "good" : "uninitialized",
                closed_t_s=NaN,
                lifecycle_end_t_s=final_t,
                duration_s=isfinite(lifecycle_start_t) ?
                    max(0.0, final_t - lifecycle_start_t) : NaN,
                filter_duration_s=is_initialized ?
                    max(0.0, final_t - initialized_t) : NaN,
                first_target_id=target_id,
                last_target_id=target_id,
                id_switch_count=0,
                measurement_update_count=0
            )
        )
    end

    if isempty(rows)
        return DataFrame(
            nav_case=String[], observer=Int[], slot=Int[], final_status=String[],
            created_t_s=Float64[], first_measurement_t_s=Float64[],
            filter_initialized_t_s=Float64[], initialization_latency_s=Float64[],
            initialization_position_error_m=Float64[], iod_group_same_target=Int[],
            iod_group_composition=String[], identity_class=String[],
            closed_t_s=Float64[], lifecycle_end_t_s=Float64[], duration_s=Float64[],
            filter_duration_s=Float64[], first_target_id=Int[], last_target_id=Int[],
            id_switch_count=Int[], measurement_update_count=Int[]
        )
    end
    return DataFrame(rows)
end

function _geometry_difficulty_dataframe(model::TrackingWindowMetricsModel)::DataFrame
    if isempty(model.geometry_rows)
        return DataFrame(
            nav_case=String[], time_s=Float64[], observer=Int[],
            simultaneously_visible_targets=Int[], minimum_angular_separation_deg=Float64[],
            closest_target_1=Int[], closest_target_2=Int[]
        )
    end
    return DataFrame(model.geometry_rows)
end

function _iod_diagnostics_dataframe(model::DistributedFusionModel)::DataFrame
    if isempty(model.iod_diagnostic_rows)
        return DataFrame(
            nav_case=String[], time_s=Float64[], observer=Int[], slot=Int[],
            stage=String[], outcome=String[], same_target=Bool[], composition=String[], reference_target=Int[],
            position_error_m=Float64[], position_rms_std_m=Float64[], validation_d2=Float64[],
            group_observers=String[], group_target_ids=String[], visible_target_count=Int[],
            minimum_angular_separation_deg=Float64[], closest_target_1=Int[], closest_target_2=Int[],
            reference_minimum_angular_separation_deg=Float64[], reference_closest_target=Int[],
            visible_target_ids=String[], minimum_3nn_radius_deg=Float64[],
            maximum_3nn_density_targets_per_sr=Float64[], reference_3nn_radius_deg=Float64[],
            reference_3nn_density_targets_per_sr=Float64[]
        )
    end
    return DataFrame(model.iod_diagnostic_rows)
end

function _iod_pairwise_dataframe(model::DistributedFusionModel)::DataFrame
    if isempty(model.iod_pairwise_rows)
        return DataFrame(
            nav_case=String[], time_s=Float64[],
            reference_observer=Int[], reference_slot=Int[],
            composition=String[], group_size=Int[], pair_count=Int[],
            covariance_gate_passed=Bool[], all_pairs_pass=Bool[], gate_m=Float64[],
            first_observer=Int[], first_slot=Int[],
            first_previous_target=Int[], first_current_target=Int[],
            second_observer=Int[], second_slot=Int[],
            second_previous_target=Int[], second_current_target=Int[],
            is_reference_pair=Bool[], pair_same_target=Bool[],
            miss_current_m=Float64[], miss_previous_m=Float64[],
            maximum_miss_m=Float64[], pair_passes=Bool[]
        )
    end
    return DataFrame(model.iod_pairwise_rows)
end

function _tracking_fragmentation_metrics(
    windows::DataFrame,
    lifecycle::DataFrame
)
    window_count = nrow(windows)
    tracks_per_window = zeros(Int, window_count)
    fragment_excess_per_window = zeros(Int, window_count)

    # Fragmentation is defined observer-by-observer and target-by-target.  Index
    # the lifecycle table once instead of scanning every track for every
    # visibility window.  The previous implementation was mathematically
    # equivalent, but its O(N_window * N_track) DataFrame-row iteration caused
    # very large allocation/GC costs in stressed realizations.
    lifecycle_by_identity = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Float64, Float64}}}()
    for track in eachrow(lifecycle)
        observer_id = Int(track.observer)
        target_id = Int(track.first_target_id)
        target_id > 0 || continue
        initialization_t = Float64(track.filter_initialized_t_s)
        lifecycle_end_t = Float64(track.lifecycle_end_t_s)
        isfinite(initialization_t) && isfinite(lifecycle_end_t) || continue
        push!(
            get!(lifecycle_by_identity, (observer_id, target_id)) do
                Tuple{Int, Float64, Float64}[]
            end,
            (Int(track.slot), initialization_t, lifecycle_end_t)
        )
    end

    for (window_idx, window) in enumerate(eachrow(windows))
        observer_id = Int(window.observer)
        target_id = Int(window.target)
        start_t = Float64(window.window_start_t_s)
        end_t = Float64(window.window_end_t_s)
        slots = Set{Int}()
        for (slot, initialization_t, lifecycle_end_t) in
            get(lifecycle_by_identity, (observer_id, target_id), Tuple{Int, Float64, Float64}[])
            initialization_t <= end_t && lifecycle_end_t >= start_t || continue
            push!(slots, slot)
        end
        tracks_per_window[window_idx] = length(slots)
        fragment_excess_per_window[window_idx] = max(length(slots) - 1, 0)
    end
    return (
        tracks_per_window=tracks_per_window,
        fragment_excess_per_window=fragment_excess_per_window,
        fragmented_windows=count(>(0), fragment_excess_per_window),
        fragment_excess_total=sum(fragment_excess_per_window)
    )
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
    fusion_timing_start_ns = ENABLE_NAV_TIMING ? time_ns() : UInt64(0)

    _retire_stale_filter_tracks!(model, t)

    observer_set = Set(model.observer_idxs)
    neighbor_map = Dict{Int, Vector{Int}}()
    for sid in model.observer_idxs
        neighbor_map[sid] = [n for n in model.nav.comms.neighbors[sid] if n in observer_set]
    end

    Q = Diagonal(model.process_q_diag)

    # Validate tentative IOD states with the first LOS acquired after initialization.
    # Confirmed states are time-aligned to the current epoch before group bootstrap.
    _validate_pending_iods!(model, t, u)

    # Snapshot only initialized active tracks. A full deepcopy of covariance_pred
    # scales with all satellites/slots and dominates memory for large comparisons.
    x_prior_snap = Dict{Tuple{Int, Int}, SVector{6, Float64}}()
    P_prior_snap = Dict{Tuple{Int, Int}, Matrix{Float64}}()

    # Buffers for updated estimates and consensus; keyed by (observer_id, slot_id).
    x_upd_buf = Dict{Tuple{Int, Int}, SVector{6, Float64}}()
    P_upd_buf = Dict{Tuple{Int, Int}, Matrix{Float64}}()
    x_pred_buf = Dict{Tuple{Int, Int}, SVector{6, Float64}}()
    P_pred_buf = Dict{Tuple{Int, Int}, Matrix{Float64}}()
    touched = Set{Tuple{Int, Int}}()
    all_active_tracks = Tuple{Int, Int}[]

    for observer_id in model.observer_idxs
        for slot_id in keys(model.nav.local_tracks[observer_id])
            track = model.nav.local_tracks[observer_id][slot_id]
            track.status == :closed && continue
            can_seed_iod = (
                track.last_meas !== nothing &&
                track.prev_meas !== nothing &&
                track.has_measure_now &&
                abs(t - track.last_meas.t) <= NAVIGATION_DT_TOL_SEC &&
                _is_consecutive_measure_pair(track.last_meas.t, track.prev_meas.t)
            )
            is_initialized = (slot_id <= size(model.filter_initialized, 2)) &&
                             model.filter_initialized[observer_id, slot_id] &&
                             (track.status == :filter_initialized)
            (is_initialized || can_seed_iod) || continue
            _ensure_fusion_slot_capacity!(model, slot_id)
            push!(all_active_tracks, (observer_id, slot_id))
            if is_initialized
                key = (observer_id, slot_id)
                x_prior_snap[key] = model.state_pred[observer_id, slot_id]
                P_prior_snap[key] = model.covariance_pred[observer_id, slot_id]
            end
        end
    end

    # STEP 1/4: IOD attempt pass (observer-centric, all observers first).
    for (observer_id, slot_id) in all_active_tracks
        # Initialized tracks do not need another IOD attempt.
        (
            model.iod_pending[observer_id, slot_id] ||
            model.iod_initialized[observer_id, slot_id] ||
            model.filter_initialized[observer_id, slot_id]
        ) && continue

        tracks = model.nav.local_tracks[observer_id]
        track = tracks[slot_id]
        track.status = :seed_ready

        valid_neighbors, miss_pass_events, miss_pass_extra = _select_iod_neighbors(
            model.nav,
            observer_id,
            slot_id,
            get(neighbor_map, observer_id, Int[]),
            t
        )
        model.miss_multi_pass_total += miss_pass_events
        model.miss_multi_pass_extra_total += miss_pass_extra
        model.iod_used_neighbors[observer_id, slot_id] = length(valid_neighbors)
        if length(valid_neighbors) < model.min_neighbor_count
            model.iod_triangulation_ready[observer_id, slot_id] = false
            continue
        end

        nodes_idx = vcat([(observer_id, slot_id)], valid_neighbors)
        model.iod_group_total += 1
        iod_group_class = _iod_group_class(model.nav, nodes_idx)
        if iod_group_class == IOD_CLASS_SAME_TARGET
            model.iod_group_same_target_total += 1
        else
            model.iod_group_mixed_target_total += 1
        end
        group_reference_target = track.last_meas === nothing ? 0 : track.last_meas.target
        group_members = _iod_group_member_labels(model.nav, nodes_idx)
        _record_iod_diagnostic!(
            model; time_s=t, observer=observer_id, slot=slot_id,
            stage="grouping", outcome="formed",
            same_target=iod_group_class == IOD_CLASS_SAME_TARGET,
            reference_target=group_reference_target,
            composition=_iod_class_label(iod_group_class),
            group_observers=group_members.group_observers,
            group_target_ids=group_members.group_target_ids,
            truth_state=u
        )

        nodes = NamedTuple[]
        used_nodes_idx = Tuple{Int, Int}[]
        A_rows = Matrix{Float64}(undef, 0, 3)
        b_rows = Float64[]
        dA_rows = Matrix{Float64}(undef, 0, 3)
        db_rows = Float64[]

        for (idx, matched_slot) in nodes_idx
            meas = _latest_track_measure(model.nav, idx, matched_slot, t)
            meas === nothing && continue

            lrate = _track_los_rate(model.nav, idx, matched_slot)
            lrate === nothing && continue

            l_unit = _safe_unit(meas.los_unit)
            prev_meas = model.nav.local_tracks[idx][matched_slot].prev_meas
            prev_meas === nothing && continue
            l_prev_unit = _safe_unit(prev_meas.los_unit)
            norm(l_prev_unit) > 1e-12 || continue
            r = meas.observer_pos
            all(isfinite, r) || continue
            v = _known_observer_vel!(model.nav, idx, u, t)
            eq = _build_agent_equations(r, v, l_unit, lrate)
            eq === nothing && continue

            A_rows = vcat(A_rows, eq.A)
            append!(b_rows, eq.b)
            dA_rows = vcat(dA_rows, eq.dA)
            append!(db_rows, eq.db)
            push!(nodes, (r=r, v=v, los=l_unit, los_prev=l_prev_unit, los_rate=lrate))
            push!(used_nodes_idx, (idx, matched_slot))
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

        model.iod_position_cov_gate_evaluated_total += 1
        position_variance_trace = tr(@view cov[1:3, 1:3])
        position_rms_std = (
            isfinite(position_variance_trace) && position_variance_trace >= 0.0
        ) ? sqrt(position_variance_trace / 3.0) : Inf
        iod_class = _iod_group_class(model.nav, used_nodes_idx)
        iod_same_target = iod_class == IOD_CLASS_SAME_TARGET
        iod_composition = _iod_class_label(iod_class)
        reference_target = track.last_meas === nothing ? 0 : track.last_meas.target
        iod_group_members = _iod_group_member_labels(model.nav, used_nodes_idx)
        iod_position_error_m = _iod_position_error_m(
            u,
            reference_target,
            SVector{6, Float64}(z_est)
        )
        _record_iod_pairwise_consistency!(
            model;
            time_s=t,
            reference_observer=observer_id,
            reference_slot=slot_id,
            nodes_idx=used_nodes_idx,
            iod_class=iod_class,
            covariance_gate_passed=position_rms_std <= IOD_MAX_POSITION_RMS_STD_M
        )
        if position_rms_std > IOD_MAX_POSITION_RMS_STD_M
            model.iod_position_cov_gate_rejected_total += 1
            if iod_same_target
                model.iod_position_cov_gate_rejected_same_target_total += 1
            else
                model.iod_position_cov_gate_rejected_mixed_target_total += 1
            end
            _record_iod_diagnostic!(
                model; time_s=t, observer=observer_id, slot=slot_id,
                stage="covariance_gate", outcome="rejected",
                same_target=iod_same_target, reference_target=reference_target,
                composition=iod_composition,
                group_observers=iod_group_members.group_observers,
                group_target_ids=iod_group_members.group_target_ids,
                position_error_m=iod_position_error_m,
                position_rms_std_m=position_rms_std,
                truth_state=u
            )
            continue
        end

        _record_iod_diagnostic!(
            model; time_s=t, observer=observer_id, slot=slot_id,
            stage="covariance_gate", outcome="passed",
            same_target=iod_same_target, reference_target=reference_target,
            composition=iod_composition,
            group_observers=iod_group_members.group_observers,
            group_target_ids=iod_group_members.group_target_ids,
            position_error_m=iod_position_error_m,
            position_rms_std_m=position_rms_std,
            truth_state=u
        )

        model.iod_estimate_state[observer_id, slot_id] = SVector{6, Float64}(z_est)
        model.iod_estimate_covariance[observer_id, slot_id] = copy(cov)
        model.iod_estimate_time_s[observer_id, slot_id] = t
        model.iod_triangulation_ready[observer_id, slot_id] = true
        track.state_estimate_now = model.iod_estimate_state[observer_id, slot_id]
        track.covariance_estimate_now = copy(model.iod_estimate_covariance[observer_id, slot_id])

        if ENABLE_IOD_ONE_STEP_VALIDATION
            model.iod_pending[observer_id, slot_id] = true
            push!(model.iod_pending_keys, (observer_id, slot_id))
            model.iod_pending_same_target[observer_id, slot_id] = iod_same_target
            model.iod_pending_class[observer_id, slot_id] = iod_class
            model.iod_initialized[observer_id, slot_id] = false
            track.status = :iod_pending
            LOG_NAV_EVENTS && println("IOD pending validation | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s")
        else
            model.iod_pending[observer_id, slot_id] = false
            delete!(model.iod_pending_keys, (observer_id, slot_id))
            model.iod_pending_same_target[observer_id, slot_id] = false
            model.iod_pending_class[observer_id, slot_id] = IOD_CLASS_UNKNOWN
            model.iod_initialized[observer_id, slot_id] = true
            track.status = :iod_initialized
            track.iod_group_same_target = iod_same_target ? 1 : 0
            track.iod_group_class = iod_class
            if track.plot_target_id == 0 && (track.last_meas !== nothing)
                # Plot-only stable label: freeze at first successful IOD init.
                track.plot_target_id = track.last_meas.target
            end

            model.iod_init_total += 1
            if iod_same_target
                model.iod_init_same_target_total += 1
            else
                model.iod_init_mixed_target_total += 1
            end
            _record_iod_diagnostic!(
                model; time_s=t, observer=observer_id, slot=slot_id,
                stage="direct_promotion", outcome="promoted",
                same_target=iod_same_target, reference_target=reference_target,
                composition=iod_composition,
                group_observers=iod_group_members.group_observers,
                group_target_ids=iod_group_members.group_target_ids,
                position_error_m=iod_position_error_m,
                position_rms_std_m=position_rms_std,
                truth_state=u
            )
            LOG_NAV_EVENTS && println("IOD init | slot=$(slot_id) | obs=$(observer_id) | t=$(round(t; digits=3)) s")
        end
    end

    # Refresh groups on schedule, on trigger event, or when buffer is empty.
    refresh_by_time = t % (NAVIGATION_RATE_SEC * CONSENSUS_GROUPS_REFRESH_STEPS) < NAVIGATION_DT_TOL_SEC
    refresh_by_event = check_trigger_consensus_groups_refresh(model, all_active_tracks)
    must_refresh_groups = refresh_by_time || refresh_by_event || isempty(model.buffer_matching_groups)

    if must_refresh_groups
        match_groups, grouping_fallback_count, mahal_pass_events, mahal_pass_extra = _build_match_groups(model, all_active_tracks, neighbor_map)
        model.buffer_matching_groups = match_groups
        _set_maintained_t2t_neighbors!(model.nav, match_groups)
        model.grouping_fallback_total += grouping_fallback_count
        model.mahal_multi_pass_total += mahal_pass_events
        model.mahal_multi_pass_extra_total += mahal_pass_extra
    else
        match_groups = model.buffer_matching_groups
    end

    filter_timing_start_ns = ENABLE_NAV_TIMING ? time_ns() : UInt64(0)
    cross_elapsed_ns = ENABLE_NAV_TIMING ? (filter_timing_start_ns - fusion_timing_start_ns) : UInt64(0)

    for group in match_groups
        group_pairs = group.selected_tracks
        isempty(group_pairs) && continue

        group_slot_match = Dict{Int, Int}(observer_id => slot_id for (observer_id, slot_id) in group_pairs)
        _initialize_group_filter_from_iod_consensus!(model, group_slot_match, neighbor_map, t, u)

        u_local = Dict{Int, Any}()
        uhat_local = Dict{Int, Any}()
        U_local = Dict{Int, Any}()
        x_prior_by_observer = Dict{Int, SVector{6, Float64}}()
        P_prior_by_observer = Dict{Int, Matrix{Float64}}()

        for (observer_id, slot_id) in group_pairs
            u_local[observer_id] = zeros(6)
            uhat_local[observer_id] = zeros(6)
            U_local[observer_id] = zeros(6, 6)

            if !model.filter_initialized[observer_id, slot_id]
                continue
            end

            key = (observer_id, slot_id)
            x_prior = get(x_prior_snap, key, model.state_pred[observer_id, slot_id])
            P_prior = get(P_prior_snap, key, model.covariance_pred[observer_id, slot_id])
            if !(_is_finite_state(x_prior) && _is_finite_cov(P_prior))
                x_prior = model.state[observer_id, slot_id]
                P_prior = model.covariance[observer_id, slot_id]
            end
            if !(_is_finite_state(x_prior) && _is_finite_cov(P_prior))
                continue
            end
            x_prior_by_observer[observer_id] = x_prior
            P_prior_by_observer[observer_id] = P_prior

            measurement = _latest_track_measure(model.nav, observer_id, slot_id, t)
            measurement === nothing && continue
            observer_pos = measurement.observer_pos
            all(isfinite, observer_pos) || continue

            u_i, uhat_i, U_i = _build_local_information_terms(
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

        end

        # Apply correction/prediction for every initialized observer/local-track in this matched group.
        kcf_local_pairs = Tuple{Int, Int}[]
        kcf_x_upd_local = Dict{Tuple{Int, Int}, SVector{6, Float64}}()
        kcf_P_upd_local = Dict{Tuple{Int, Int}, Matrix{Float64}}()
        for (observer_id, slot_id) in group_pairs
            model.filter_initialized[observer_id, slot_id] || continue
            key = (observer_id, slot_id)
            x_prior = get(x_prior_snap, key, model.state_pred[observer_id, slot_id])
            P_prior = get(P_prior_snap, key, model.covariance_pred[observer_id, slot_id])
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
        if !isfinite(track.filter_initialized_t)
            track.filter_initialized_t = t
        end
    end

    if LOG_SIM_PROGRESS && SIM_PROGRESS_INTERVAL_SEC > 0.0 && t > 0.0 &&
       abs(t / SIM_PROGRESS_INTERVAL_SEC - round(t / SIM_PROGRESS_INTERVAL_SEC)) <= 1e-9
        gc_start_ns = time_ns()
        FORCE_INCREMENTAL_GC && GC.gc(false)
        gc_elapsed_s = (time_ns() - gc_start_ns) / 1e9
        active_filter_tracks = count(model.filter_initialized)
        local_track_count = sum(length(model.nav.local_tracks[obs]) for obs in model.observer_idxs)
        active_local_tracks = sum(
            count(track -> track.status != :closed, values(model.nav.local_tracks[obs]))
            for obs in model.observer_idxs
        )
        archived_track_count = sum(
            length(model.nav.closed_track_lifecycle[obs]) for obs in model.observer_idxs
        )
        local_hypothesis_count = sum(
            length(model.nav.local_measurement_hypotheses[obs]) for obs in model.observer_idxs
        )
        println(
            "SIM_PROGRESS | t_s=$(round(t; digits=1))" *
            " | active_filter_tracks=$(active_filter_tracks)" *
            " | active_local_tracks=$(active_local_tracks)" *
            " | retained_local_tracks=$(local_track_count)" *
            " | archived_tracks=$(archived_track_count)" *
            " | local_hypotheses=$(local_hypothesis_count)" *
            " | pending_iod=$(length(model.iod_pending_keys))" *
            " | max_rss_gib=$(round(Sys.maxrss() / 2.0^30; digits=3))" *
            " | incremental_gc_s=$(round(gc_elapsed_s; digits=3))"
        )
        flush(stdout)
        flush(stderr)
    end

    if ENABLE_NAV_TIMING
        filter_elapsed_ns = time_ns() - filter_timing_start_ns
        fusion_elapsed_ns = time_ns() - fusion_timing_start_ns
        local_elapsed_ns = _flush_local_da_timing!(model.nav.timing)
        local_elapsed_uint = isfinite(local_elapsed_ns) ? UInt64(round(Int, local_elapsed_ns)) : UInt64(0)
        push!(model.nav.timing.cross_da_epoch_ns, Float64(cross_elapsed_ns))
        push!(model.nav.timing.filter_epoch_ns, Float64(filter_elapsed_ns))
        push!(model.nav.timing.fusion_epoch_ns, Float64(fusion_elapsed_ns))
        push!(model.nav.timing.total_epoch_ns, Float64(local_elapsed_uint + fusion_elapsed_ns))
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


function _sample_truncated_gaussian(
    rng::AbstractRNG,
    mean_value::Float64,
    std_value::Float64,
    lower_bound::Float64,
    upper_bound::Float64;
    max_attempts::Int=100_000
)::Float64
    std_value > 0.0 || throw(ArgumentError("std_value must be positive"))
    lower_bound < upper_bound || throw(ArgumentError("lower_bound must be smaller than upper_bound"))
    for _ in 1:max_attempts
        sample = mean_value + std_value * randn(rng)
        lower_bound <= sample <= upper_bound && return sample
    end
    error("Unable to sample the requested truncated Gaussian after $(max_attempts) attempts")
end

@inline function _sample_wrapped_gaussian_deg(
    rng::AbstractRNG,
    mean_deg::Float64,
    std_deg::Float64
)::Float64
    return mod(mean_deg + std_deg * randn(rng), 360.0)
end

@inline function _sample_log_uniform(
    rng::AbstractRNG,
    lower_bound::Float64,
    upper_bound::Float64
)::Float64
    0.0 < lower_bound < upper_bound || throw(ArgumentError("log-uniform bounds must satisfy 0 < lower < upper"))
    log_lower = log10(lower_bound)
    return 10.0^(log_lower + (log10(upper_bound) - log_lower) * rand(rng))
end

function generate_observer_centered_debris_cluster(
    n_targets::Int,
    observer_a_m::Float64,
    observer_i_deg::Float64,
    observer_mean_anomalies_deg;
    observer_raan_deg::Float64=10.0,
    observer_aop_deg::Float64=14.0,
    earth_radius_m::Float64=6_378_136.3,
    minimum_periapsis_altitude_m::Float64=300_000.0,
    rng::AbstractRNG=MersenneTwister(SCENARIO_RNG_SEED)
)
    n_targets > 0 || throw(ArgumentError("n_targets must be positive"))
    earth_radius_m > 0.0 || throw(ArgumentError("earth_radius_m must be positive"))
    minimum_periapsis_altitude_m >= 0.0 || throw(ArgumentError("minimum_periapsis_altitude_m must be nonnegative"))

    # Table-II distributions.
    a_std_m = 200_000.0
    i_std_deg = 4.0
    raan_std_deg = 8.0
    aop_std_deg = 40.0
    mean_anomaly_std_deg = 1.5
    eccentricity_min = 1e-5
    eccentricity_max = 0.02

    m_center = sum(observer_mean_anomalies_deg) / length(observer_mean_anomalies_deg)
    target_defs = NamedTuple{(:a_m, :e, :i_deg, :raan_deg, :aop_deg, :M_deg), Tuple{Float64, Float64, Float64, Float64, Float64, Float64}}[]

    for _ in 1:n_targets
        # Draw eccentricity uniformly in log10(e). Conditional on that draw,
        # truncate the semi-major-axis Gaussian at the minimum value that
        # guarantees a periapsis altitude of at least 300 km.
        e_val = _sample_log_uniform(rng, eccentricity_min, eccentricity_max)
        a_min_m = (earth_radius_m + minimum_periapsis_altitude_m) / (1.0 - e_val)
        a_val_m = _sample_truncated_gaussian(
            rng,
            observer_a_m,
            a_std_m,
            a_min_m,
            Inf
        )

        # Inclination is a Gaussian truncated to its physical interval.
        i_val_deg = _sample_truncated_gaussian(
            rng,
            observer_i_deg,
            i_std_deg,
            0.0,
            180.0
        )

        push!(target_defs, (
            a_m = a_val_m,
            e = e_val,
            i_deg = i_val_deg,
            raan_deg = _sample_wrapped_gaussian_deg(rng, observer_raan_deg, raan_std_deg),
            aop_deg = _sample_wrapped_gaussian_deg(rng, observer_aop_deg, aop_std_deg),
            M_deg = _sample_wrapped_gaussian_deg(rng, m_center, mean_anomaly_std_deg)
        ))
    end

    return Tuple(target_defs)
end

planet = Earth("", SPICE_PATH)
isfile(EARTH_HARMONICS_FILE) || error("Earth harmonics file not found: $EARTH_HARMONICS_FILE")

e = 1e-4
raan_deg = 10.0
aop_deg = 14.0

observer_a_m = 6_963.0e3
target_a_m = 6_949.0e3

observer_i_deg = 70.0

observer_mean_anomaly_center_deg = 290.0
observer_raan_offsets_deg = (-1.6, -0.5, 0.5, 1.6)
observer_mean_anomaly_offsets_deg = (-4.0, -1.3, 1.3, 4.0)

# Previous compact 3-observer configuration.
# observer_i_deg = 85.0
# observer_mean_anomalies_deg = (288.0, 290.0, 292.0)
# observer_defs = (
#     (a_m=observer_a_m, i_deg=observer_i_deg, M_deg=observer_mean_anomalies_deg[1]),
#     (a_m=observer_a_m, i_deg=observer_i_deg, M_deg=observer_mean_anomalies_deg[2]),
#     (a_m=observer_a_m, i_deg=observer_i_deg, M_deg=observer_mean_anomalies_deg[3])
# )

# Sparse regional Walker-like observer arc: 4 nearby planes x 4 close-phased
# observers. The spacing keeps the graph local for R_com = 300 km while
# preserving repeated joint detections for R_det = 500 km.
observer_defs = Tuple([
    (
        a_m=observer_a_m,
        i_deg=observer_i_deg,
        raan_deg=mod(raan_deg + raan_offset_deg, 360.0),
        aop_deg=aop_deg,
        M_deg=observer_mean_anomaly_center_deg + m_offset_deg
    )
    for raan_offset_deg in observer_raan_offsets_deg
    for m_offset_deg in observer_mean_anomaly_offsets_deg
])
observer_mean_anomalies_deg = Tuple(obs.M_deg for obs in observer_defs)

const N_CLUSTER_TARGETS = parse(Int, get(ENV, "SPACEAGORA_N_CLUSTER_TARGETS", "300"))
target_defs = generate_observer_centered_debris_cluster(
    N_CLUSTER_TARGETS,
    observer_a_m,
    observer_i_deg,
    observer_mean_anomalies_deg;
    observer_raan_deg=raan_deg,
    observer_aop_deg=aop_deg,
    earth_radius_m=planet.Rp_e
)

target_population_df = DataFrame(
    target_number=collect(1:length(target_defs)),
    semi_major_axis_km=[target.a_m * 1e-3 for target in target_defs],
    eccentricity=[target.e for target in target_defs],
    log10_eccentricity=[log10(target.e) for target in target_defs],
    inclination_deg=[target.i_deg for target in target_defs],
    raan_deg=[target.raan_deg for target in target_defs],
    argument_of_perigee_deg=[target.aop_deg for target in target_defs],
    mean_anomaly_deg=[target.M_deg for target in target_defs],
    periapsis_altitude_km=[
        (target.a_m * (1.0 - target.e) - planet.Rp_e) * 1e-3
        for target in target_defs
    ]
)

# target_defs = (
#     (a_m=target_a_m , i_deg=observer_i_deg+3.0, M_deg=observer_mean_anomalies_deg[2]),
#     (a_m=target_a_m , i_deg=observer_i_deg+1.5, M_deg=observer_mean_anomalies_deg[2]),
#     (a_m=target_a_m , i_deg=observer_i_deg-1.5, M_deg=observer_mean_anomalies_deg[2]),
#     (a_m=target_a_m , i_deg=observer_i_deg-3.0, M_deg=observer_mean_anomalies_deg[2])
# )

spacecraft = SpacecraftModel[]
for (sat_id, obs) in enumerate(observer_defs)
    obs_raan_deg = hasproperty(obs, :raan_deg) ? obs.raan_deg : raan_deg
    obs_aop_deg = hasproperty(obs, :aop_deg) ? obs.aop_deg : aop_deg
    push!(
        spacecraft,
        make_translational_spacecraft(
            sat_id,
            obs.a_m,
            e,
            obs.i_deg,
            obs_raan_deg,
            obs_aop_deg,
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
    length(spacecraft);
    measurement_bias_enabled=ENABLE_MEASUREMENT_BIAS,
    measurement_bias_std_rad=MEASUREMENT_BIAS_RAD,
    bias_rng=MersenneTwister(MEASUREMENT_BIAS_RNG_SEED),
    measurement_noise_rng=MersenneTwister(SENSOR_RNG_SEED),
    misdetection_rng=MersenneTwister(MISDETECTION_RNG_SEED),
    false_alarm_count_rng=MersenneTwister(FALSE_ALARM_COUNT_RNG_SEED),
    false_alarm_direction_rng=MersenneTwister(FALSE_ALARM_DIRECTION_RNG_SEED)
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
    1,
    SIGMA_THETA_RAD,
    EKF_CONSENSUS_ITERS,
    EKF_PROCESS_Q_DIAG,
    planet.μ,
    fusion_j2,
    fusion_ref_radius_m,
    length(spacecraft)
)
central_oracle_model = CentralizedOracleFusionModel(
    optical_sensor,
    observer_idxs,
    target_idxs,
    SIGMA_THETA_RAD,
    EKF_PROCESS_Q_DIAG,
    planet.μ,
    fusion_j2,
    fusion_ref_radius_m
)
estimator_model = USE_CENTRALIZED_ORACLE ? central_oracle_model : fusion_model
tracking_metrics = TrackingWindowMetricsModel(
    optical_sensor,
    estimator_model,
    observer_idxs,
    target_idxs,
    length(spacecraft)
)

truth_geometry_replay_model = if TRUTH_GEOMETRY_REPLAY_ONLY
    isempty(TRUTH_GEOMETRY_REPLAY_INPUT) && error(
        "SPACEAGORA_TRUTH_GEOMETRY_REPLAY_INPUT is required in truth replay mode"
    )
    isempty(TRUTH_GEOMETRY_REPLAY_OUTPUT) && error(
        "SPACEAGORA_TRUTH_GEOMETRY_REPLAY_OUTPUT is required in truth replay mode"
    )
    IODTruthGeometryReplayModel(
        TRUTH_GEOMETRY_REPLAY_INPUT,
        observer_idxs,
        target_idxs,
        THRESHOLD_DISTANCE_KM * 1e3
    )
else
    nothing
end

navigation_effectors = if TRUTH_GEOMETRY_REPLAY_ONLY
    Any[truth_geometry_replay_model]
else
    effectors = Any[optical_sensor]
    if USE_CENTRALIZED_ORACLE
        push!(effectors, central_oracle_model)
    else
        push!(effectors, observer_nav)
        push!(effectors, fusion_model)
    end
    push!(effectors, tracking_metrics)
    effectors
end
navigation_rates = fill(NAVIGATION_RATE_SEC, length(navigation_effectors))

areas = [sc.root.ref_area for sc in spacecraft]
dynamic_effectors = (
    InverseSquaredGravityModel(),
    gh_model,
    SunMoonThirdBodyModel(planet),
    CannonballSRPModel(planet, areas, 1.3, 4.56e-6, 149_597_870_700.0)
)

now_utc = START_UTC
args = SimulationConfiguration(
    simulation_settings=SimulationSettings(
        results=SAVE_SIMULATION_RESULTS,
        verbose=false,
        generate_plots=false,
        results_directory=joinpath(COMPARISON_OUTPUT_ROOT, NAVIGATION_OUTPUT_LABEL),
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
println("Navigation comparison case: $(NAVIGATION_CASE)")
println(
    "Data association: ",
    USE_NO_DA ?
    "disabled; measurements are routed through their known target channel" :
    "enabled"
)
println("Scenario seed: $(SCENARIO_RNG_SEED), sensor seed: $(SENSOR_RNG_SEED), observer OD seed: $(OBSERVER_OD_RNG_SEED)")
println("Propagation time: $(MISSION_TIME_SEC / 3600.0) hours")
println("Observer idxs: $(observer_idxs)")
if length(target_idxs) <= 20
    println("Target idxs: $(target_idxs)")
else
    println("Target idxs: $(first(target_idxs)):$(last(target_idxs)) ($(length(target_idxs)) targets)")
end
flush(stdout)
flush(stderr)

extra_save_fields = SaveField[]
if SAVE_SIMULATION_RESULTS && SAVE_TARGET_ESTIMATE_FIELDS
    for target_idx in target_idxs
        let tgt = target_idx
            push!(
                extra_save_fields,
                SaveField(
                    Symbol("dekf_target$(tgt)_state"),
                    (u, t, integrator) -> begin
                        return Dict("obs$(obs)" => _fusion_state_for_save(estimator_model, obs, tgt, t) for obs in observer_idxs)
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
                        return Dict("obs$(obs)" => Float64(_fusion_slot_for_save(estimator_model, obs, tgt, t)) for obs in observer_idxs)
                    end;
                    per_satellite=false,
                    column_prefix="dekf_t$(tgt)_slot"
                )
            )
        end
    end
end

base_save_fields = SAVE_SIMULATION_RESULTS ? SimulationModel.default_save_fields(args) : SaveField[]
save_fields = vcat(base_save_fields, extra_save_fields)
println("RUN_SIMULATION_START")
flush(stdout)
flush(stderr)
simulation_timing = @timed run_simulation(args; isolate_state=false, save_fields=save_fields)
t = simulation_timing.time
simulation_gc_time_s = simulation_timing.gctime
simulation_allocated_gib = simulation_timing.bytes / 2.0^30
simulation_max_rss_gib = Sys.maxrss() / 2.0^30
println(
    "RUN_SIMULATION_DONE | elapsed_s=$(round(t; digits=3))" *
    " | gc_time_s=$(round(simulation_gc_time_s; digits=3))" *
    " | allocated_gib=$(round(simulation_allocated_gib; digits=3))" *
    " | max_rss_gib=$(round(simulation_max_rss_gib; digits=3))"
)
simulation_timing = nothing
flush(stdout)
flush(stderr)
if TRUTH_GEOMETRY_REPLAY_ONLY
    replay_df = _truth_geometry_replay_dataframe(truth_geometry_replay_model)
    mkpath(dirname(TRUTH_GEOMETRY_REPLAY_OUTPUT))
    CSV.write(TRUTH_GEOMETRY_REPLAY_OUTPUT, replay_df)
    println("TRUTH_GEOMETRY_REPLAY_DONE | events=$(nrow(replay_df)) | output=$(TRUTH_GEOMETRY_REPLAY_OUTPUT)")
    flush(stdout)
    exit(0)
end
println("POSTPROCESS_START")
flush(stdout)
flush(stderr)
postprocess_start_ns = time_ns()
function _postprocess_checkpoint(label::AbstractString)
    elapsed_s = (time_ns() - postprocess_start_ns) / 1e9
    println("POSTPROCESS_STAGE_DONE | stage=$(label) | elapsed_s=$(round(elapsed_s; digits=3))")
    flush(stdout)
    flush(stderr)
    return nothing
end
ENABLE_NAV_TIMING && _flush_local_da_timing!(observer_nav.timing)
_finalize_tracking_windows!(tracking_metrics)
_postprocess_checkpoint("finalize_windows")
tracking_window_df = _tracking_window_dataframe(tracking_metrics)
tracking_observer_summary_df = _tracking_observer_summary(tracking_metrics)
empty!(tracking_metrics.window_rows)
for samples in tracking_metrics.current_error_samples
    empty!(samples)
end
_postprocess_checkpoint("tracking_tables")
track_lifecycle_df = USE_CENTRALIZED_ORACLE ?
    _track_lifecycle_dataframe(central_oracle_model, MISSION_TIME_SEC) :
    _track_lifecycle_dataframe(observer_nav, MISSION_TIME_SEC)
for observer_id in observer_idxs
    empty!(observer_nav.local_tracks[observer_id])
    empty!(observer_nav.local_measurement_hypotheses[observer_id])
end
_postprocess_checkpoint("track_lifecycle")
geometry_difficulty_df = _geometry_difficulty_dataframe(tracking_metrics)
iod_diagnostics_df = _iod_diagnostics_dataframe(fusion_model)
iod_pairwise_df = _iod_pairwise_dataframe(fusion_model)
empty!(tracking_metrics.geometry_rows)
empty!(fusion_model.iod_diagnostic_rows)
empty!(fusion_model.iod_pairwise_rows)
GC.gc(true)
_postprocess_checkpoint("diagnostic_tables")
fragmentation = _tracking_fragmentation_metrics(tracking_window_df, track_lifecycle_df)
tracking_window_df.tracks_per_window = fragmentation.tracks_per_window
tracking_window_df.fragment_excess = fragmentation.fragment_excess_per_window
_postprocess_checkpoint("fragmentation")

# Simulation-side buffers are released as soon as their DataFrames exist.
# Collect only temporary allocations produced by fragmentation here.
GC.gc(false)
_postprocess_checkpoint("release_simulation_buffers")

runtime_local_da = _timing_summary_ms(observer_nav.timing.local_da_epoch_ns)
runtime_cross_da = _timing_summary_ms(observer_nav.timing.cross_da_epoch_ns)
runtime_association = _timing_summary_ms(Float64[
    observer_nav.timing.local_da_epoch_ns[index] +
    observer_nav.timing.cross_da_epoch_ns[index]
    for index in 1:min(
        length(observer_nav.timing.local_da_epoch_ns),
        length(observer_nav.timing.cross_da_epoch_ns)
    )
])
runtime_filter = _timing_summary_ms(observer_nav.timing.filter_epoch_ns)
runtime_fusion = _timing_summary_ms(observer_nav.timing.fusion_epoch_ns)
runtime_total = _timing_summary_ms(observer_nav.timing.total_epoch_ns)
sensor_visible_opportunities = optical_sensor.visible_opportunity_total
sensor_true_detections = optical_sensor.true_detection_total
sensor_missed_detections = optical_sensor.missed_detection_total
sensor_false_alarms = optical_sensor.false_alarm_total
sensor_observer_epochs = optical_sensor.observer_epoch_total
sensor_false_alarm_nonzero_epochs = optical_sensor.false_alarm_nonzero_epoch_total
sensor_false_alarm_multiple_epochs = optical_sensor.false_alarm_multiple_epoch_total
sensor_false_alarm_max_per_epoch = optical_sensor.false_alarm_max_per_epoch
sensor_biased_measurements = optical_sensor.biased_measurement_total
sensor_detection_rate_pct = sensor_visible_opportunities > 0 ?
    100.0 * sensor_true_detections / sensor_visible_opportunities :
    NaN
sensor_missed_rate_pct = sensor_visible_opportunities > 0 ?
    100.0 * sensor_missed_detections / sensor_visible_opportunities :
    NaN
sensor_false_alarm_mean_per_observer_epoch = sensor_observer_epochs > 0 ?
    sensor_false_alarms / sensor_observer_epochs : NaN
sensor_false_alarm_nonzero_epoch_pct = sensor_observer_epochs > 0 ?
    100.0 * sensor_false_alarm_nonzero_epochs / sensor_observer_epochs : NaN
sensor_false_alarm_multiple_epoch_pct = sensor_observer_epochs > 0 ?
    100.0 * sensor_false_alarm_multiple_epochs / sensor_observer_epochs : NaN
configured_false_alarm_nonzero_epoch_pct = ENABLE_FALSE_ALARMS ?
    100.0 * (1.0 - exp(-FALSE_ALARM_RATE)) : 0.0
configured_false_alarm_multiple_epoch_pct = ENABLE_FALSE_ALARMS ?
    100.0 * (1.0 - exp(-FALSE_ALARM_RATE) * (1.0 + FALSE_ALARM_RATE)) : 0.0
bias_components = Float64[]
for observer_id in observer_idxs
    append!(bias_components, optical_sensor.observer_bias_rotation_vectors[observer_id])
end
sensor_realized_bias_component_rms_rad = isempty(bias_components) ? NaN :
    sqrt(mean(abs2, bias_components))
sensor_realized_bias_norm_max_rad = isempty(observer_idxs) ? NaN : maximum(
    norm(optical_sensor.observer_bias_rotation_vectors[observer_id]) for observer_id in observer_idxs
)

println("LOS sensor error diagnostics:")
println("  enabled: missed_detections=$(ENABLE_MISDETECTIONS), false_alarms=$(ENABLE_FALSE_ALARMS), measurement_bias=$(ENABLE_MEASUREMENT_BIAS)")
println("  configured: missed_detection_rate=$(MISDETECTION_RATE), false_alarm_rate=$(FALSE_ALARM_RATE), measurement_bias_std_rad=$(MEASUREMENT_BIAS_RAD)")
println("  visible_opportunities=$(sensor_visible_opportunities), true_detections=$(sensor_true_detections), missed_detections=$(sensor_missed_detections)")
println("  realized_detection_rate_pct=$(round(sensor_detection_rate_pct; digits=3)), realized_missed_rate_pct=$(round(sensor_missed_rate_pct; digits=3))")
println("  false_alarms=$(sensor_false_alarms), mean_per_observer_epoch=$(sensor_false_alarm_mean_per_observer_epoch), nonzero_epochs=$(sensor_false_alarm_nonzero_epochs), multiple_epochs=$(sensor_false_alarm_multiple_epochs), max_per_epoch=$(sensor_false_alarm_max_per_epoch)")
println("  biased_measurements=$(sensor_biased_measurements), realized_bias_component_rms_rad=$(sensor_realized_bias_component_rms_rad), realized_bias_norm_max_rad=$(sensor_realized_bias_norm_max_rad)")
_postprocess_checkpoint("sensor_statistics")

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
println("  local track-to-track ambiguity (consensus Mahalanobis): events=$(fusion_model.mahal_multi_pass_total), extra_candidates=$(fusion_model.mahal_multi_pass_extra_total), ratio_fail=$(fusion_model.tt_ratio_fail_total), mutual_fail=$(fusion_model.tt_mutual_fail_total)")
println("  consensus-group fallback (global-id split): enabled=$(CONSENSUS_GROUP_FALLBACK_ENABLE), count=$(fusion_model.grouping_fallback_total)")
if ENABLE_NAV_TIMING
    println("Runtime diagnostics (proposed, navigation callback only):")
    println("  epochs=$(runtime_total.count)")
    println("  total_epoch_mean_ms=$(round(runtime_total.mean; digits=3)), p05_ms=$(round(runtime_total.p05; digits=3)), median_ms=$(round(runtime_total.median; digits=3)), p95_ms=$(round(runtime_total.p95; digits=3))")
    println("  association_mean_ms=$(round(runtime_association.mean; digits=3)), p05_ms=$(round(runtime_association.p05; digits=3)), median_ms=$(round(runtime_association.median; digits=3)), p95_ms=$(round(runtime_association.p95; digits=3))")
    println("  filter_mean_ms=$(round(runtime_filter.mean; digits=3)), p05_ms=$(round(runtime_filter.p05; digits=3)), median_ms=$(round(runtime_filter.median; digits=3)), p95_ms=$(round(runtime_filter.p95; digits=3))")
end

hyp_h1_created = sum(observer_nav.hyp_h1_created_total[obs] for obs in observer_idxs)
hyp_h1_to_h2_created = sum(observer_nav.hyp_h1_to_h2_created_total[obs] for obs in observer_idxs)
hyp_h1_to_h2_same_target = sum(observer_nav.hyp_h1_to_h2_same_target_total[obs] for obs in observer_idxs)
hyp_h1_to_h2_mixed_target = hyp_h1_to_h2_created - hyp_h1_to_h2_same_target
hyp_h2_to_h3_attempted = sum(observer_nav.hyp_h2_to_h3_attempted_total[obs] for obs in observer_idxs)
hyp_h3_los_rate_pass = sum(observer_nav.hyp_h3_los_rate_pass_total[obs] for obs in observer_idxs)
hyp_h3_los_rate_pass_same_target = sum(observer_nav.hyp_h3_los_rate_pass_same_target_total[obs] for obs in observer_idxs)
hyp_h3_los_rate_pass_mixed_target = hyp_h3_los_rate_pass - hyp_h3_los_rate_pass_same_target
hyp_h3_los_rate_fail = sum(observer_nav.hyp_h3_los_rate_fail_total[obs] for obs in observer_idxs)
hyp_h3_los_rate_fail_same_target = sum(observer_nav.hyp_h3_los_rate_fail_same_target_total[obs] for obs in observer_idxs)
hyp_h3_los_rate_fail_mixed_target = hyp_h3_los_rate_fail - hyp_h3_los_rate_fail_same_target
hyp_promoted = sum(observer_nav.hyp_promoted_total[obs] for obs in observer_idxs)
hyp_promoted_same_target = sum(observer_nav.hyp_promoted_same_target_total[obs] for obs in observer_idxs)
hyp_promoted_mixed_target = hyp_promoted - hyp_promoted_same_target
tracks_created_with_nonreal_measurements = sum(observer_nav.tracks_created_with_nonreal_measurements_total[obs] for obs in observer_idxs)
println("Local hypothesis diagnostics:")
println("  H1_created=$(hyp_h1_created)")
println("  H1_to_H2_created=$(hyp_h1_to_h2_created), same_target=$(hyp_h1_to_h2_same_target), mixed_target=$(hyp_h1_to_h2_mixed_target)")
println("  H2_to_H3_attempted=$(hyp_h2_to_h3_attempted)")
println("  H3_los_rate_pass=$(hyp_h3_los_rate_pass), same_target=$(hyp_h3_los_rate_pass_same_target), mixed_target=$(hyp_h3_los_rate_pass_mixed_target)")
println("  H3_los_rate_fail=$(hyp_h3_los_rate_fail), same_target=$(hyp_h3_los_rate_fail_same_target), mixed_target=$(hyp_h3_los_rate_fail_mixed_target)")
println("  promoted=$(hyp_promoted), same_target=$(hyp_promoted_same_target), mixed_target=$(hyp_promoted_mixed_target)")
println("  tracks_created_with_nonreal_measurements=$(tracks_created_with_nonreal_measurements)")
_postprocess_checkpoint("association_counters")

meas_unique_total = sum(observer_nav.meas_commit_unique_total[obs] for obs in observer_idxs)
meas_unique_correct = sum(observer_nav.meas_commit_unique_correct_total[obs] for obs in observer_idxs)
meas_unique_wrong = sum(observer_nav.meas_commit_unique_wrong_total[obs] for obs in observer_idxs)
meas_ambig_total = sum(observer_nav.meas_commit_ambig_total[obs] for obs in observer_idxs)
meas_ambig_correct = sum(observer_nav.meas_commit_ambig_correct_total[obs] for obs in observer_idxs)
meas_ambig_wrong = sum(observer_nav.meas_commit_ambig_wrong_total[obs] for obs in observer_idxs)
meas_ambig_dropped = sum(observer_nav.meas_commit_ambig_dropped_total[obs] for obs in observer_idxs)
meas_false_alarm_committed = sum(observer_nav.meas_commit_false_alarm_total[obs] for obs in observer_idxs)
meas_total_committed = meas_unique_total + meas_ambig_total
meas_total_correct = meas_unique_correct + meas_ambig_correct
meas_total_wrong = meas_unique_wrong + meas_ambig_wrong
meas_commit_acc_pct = meas_total_committed > 0 ? 100.0 * meas_total_correct / meas_total_committed : NaN
meas_true_opportunities = sum(observer_nav.meas_true_opportunity_total[obs] for obs in observer_idxs)
meas_true_not_committed = max(meas_true_opportunities - meas_total_correct, 0)
meas_recall_pct = _pct(Float64(meas_total_correct), Float64(meas_true_opportunities))
meas_unique_acc_pct = meas_unique_total > 0 ? 100.0 * meas_unique_correct / meas_unique_total : NaN
meas_ambig_acc_pct = meas_ambig_total > 0 ? 100.0 * meas_ambig_correct / meas_ambig_total : NaN

tt_commit_total = fusion_model.tt_commit_total
tt_commit_correct = fusion_model.tt_commit_correct_total
tt_commit_wrong = fusion_model.tt_commit_wrong_total
tt_commit_unknown = fusion_model.tt_commit_unknown_total
tt_skipped = fusion_model.tt_skipped_total
tt_no_candidate = fusion_model.tt_no_candidate_total
tt_gate_fail = fusion_model.tt_gate_fail_total
tt_ratio_fail = fusion_model.tt_ratio_fail_total
tt_mutual_fail = fusion_model.tt_mutual_fail_total
tt_component_conflict_rejected = fusion_model.tt_component_conflict_rejected_total
tt_skipped_same_target_present = fusion_model.tt_skipped_same_target_present_total
tt_skipped_no_same_target = fusion_model.tt_skipped_no_same_target_total
tt_skipped_unknown_source_target = fusion_model.tt_skipped_unknown_source_target_total
tt_attempt_total = fusion_model.tt_attempt_total
tt_true_opportunities = fusion_model.tt_true_opportunity_total
tt_true_not_committed = max(tt_true_opportunities - tt_commit_correct, 0)
tt_recall_pct = _pct(Float64(tt_commit_correct), Float64(tt_true_opportunities))
tt_known_total = tt_commit_correct + tt_commit_wrong
tt_acc_pct = tt_known_total > 0 ? 100.0 * tt_commit_correct / tt_known_total : NaN

xm2m_candidate_pairs = observer_nav.xm2m_candidate_pair_total
xm2m_gate_pass = observer_nav.xm2m_gate_pass_total
xm2m_gate_pass_same = observer_nav.xm2m_gate_pass_same_target_total
xm2m_gate_pass_mixed = xm2m_gate_pass - xm2m_gate_pass_same
xm2m_gate_fail = observer_nav.xm2m_gate_fail_total
xm2m_gate_fail_same = observer_nav.xm2m_gate_fail_same_target_total
xm2m_gate_fail_mixed = xm2m_gate_fail - xm2m_gate_fail_same
xm2m_selected = observer_nav.xm2m_selected_pair_total
xm2m_selected_same = observer_nav.xm2m_selected_pair_same_target_total
xm2m_selected_mixed = xm2m_selected - xm2m_selected_same
xm2m_gate_pass_same_pct = _pct(Float64(xm2m_gate_pass_same), Float64(xm2m_gate_pass))
xm2m_selected_same_pct = _pct(Float64(xm2m_selected_same), Float64(xm2m_selected))

iod_group_total = fusion_model.iod_group_total
iod_group_same = fusion_model.iod_group_same_target_total
iod_group_mixed = fusion_model.iod_group_mixed_target_total
iod_group_same_pct = _pct(Float64(iod_group_same), Float64(iod_group_total))
iod_init_total = fusion_model.iod_init_total
iod_init_same = fusion_model.iod_init_same_target_total
iod_init_mixed = fusion_model.iod_init_mixed_target_total
iod_init_same_pct = _pct(Float64(iod_init_same), Float64(iod_init_total))
iod_pos_cov_gate_evaluated = fusion_model.iod_position_cov_gate_evaluated_total
iod_pos_cov_gate_rejected = fusion_model.iod_position_cov_gate_rejected_total
iod_pos_cov_gate_rejected_same = fusion_model.iod_position_cov_gate_rejected_same_target_total
iod_pos_cov_gate_rejected_mixed = fusion_model.iod_position_cov_gate_rejected_mixed_target_total
iod_pos_cov_gate_rejected_pct = _pct(
    Float64(iod_pos_cov_gate_rejected),
    Float64(iod_pos_cov_gate_evaluated)
)
iod_validation_attempted = fusion_model.iod_validation_attempted_total
iod_validation_confirmed = fusion_model.iod_validation_confirmed_total
iod_validation_rejected = fusion_model.iod_validation_rejected_total
iod_validation_confirmed_same = fusion_model.iod_validation_confirmed_same_target_total
iod_validation_confirmed_mixed = fusion_model.iod_validation_confirmed_mixed_target_total
iod_validation_rejected_same = fusion_model.iod_validation_rejected_same_target_total
iod_validation_rejected_mixed = fusion_model.iod_validation_rejected_mixed_target_total
iod_validation_no_measure = fusion_model.iod_validation_no_measure_total
iod_validation_no_measure_same = fusion_model.iod_validation_no_measure_same_target_total
iod_validation_no_measure_mixed = fusion_model.iod_validation_no_measure_mixed_target_total
iod_validation_pending_end = count(fusion_model.iod_pending)
iod_validation_confirmed_pct = _pct(
    Float64(iod_validation_confirmed),
    Float64(iod_validation_attempted)
)

function _iod_error_values(stage::String, outcome::String, same_target::Bool)
    return Float64[
        Float64(row.position_error_m) for row in eachrow(iod_diagnostics_df)
        if String(row.stage) == stage && String(row.outcome) == outcome &&
            Bool(row.same_target) == same_target && isfinite(Float64(row.position_error_m))
    ]
end

function _error_stats(values::Vector{Float64})
    return isempty(values) ? (mean=NaN, median=NaN, p95=NaN, maximum=NaN) : (
        mean=mean(values), median=median(values), p95=quantile(values, 0.95), maximum=maximum(values)
    )
end

iod_cov_pass_same = nrow(iod_diagnostics_df) == 0 ? 0 : sum(
    (iod_diagnostics_df.stage .== "covariance_gate") .&
    (iod_diagnostics_df.outcome .== "passed") .& iod_diagnostics_df.same_target
)
iod_cov_pass_mixed = nrow(iod_diagnostics_df) == 0 ? 0 : sum(
    (iod_diagnostics_df.stage .== "covariance_gate") .&
    (iod_diagnostics_df.outcome .== "passed") .& .!iod_diagnostics_df.same_target
)
iod_cov_pass_same_error = _error_stats(_iod_error_values("covariance_gate", "passed", true))
iod_cov_pass_mixed_error = _error_stats(_iod_error_values("covariance_gate", "passed", false))
iod_promoted_same_error = _error_stats(vcat(
    _iod_error_values("next_step_validation", "promoted", true),
    _iod_error_values("direct_promotion", "promoted", true)
))
iod_promoted_mixed_error = _error_stats(vcat(
    _iod_error_values("next_step_validation", "promoted", false),
    _iod_error_values("direct_promotion", "promoted", false)
))
iod_rejected_same_error = _error_stats(_iod_error_values("next_step_validation", "rejected_score", true))
iod_rejected_mixed_error = _error_stats(_iod_error_values("next_step_validation", "rejected_score", false))
iod_cov_true_rejection_pct = _pct(
    Float64(iod_pos_cov_gate_rejected_same),
    Float64(iod_pos_cov_gate_rejected_same + iod_cov_pass_same)
)
iod_cov_wrong_rejection_pct = _pct(
    Float64(iod_pos_cov_gate_rejected_mixed),
    Float64(iod_pos_cov_gate_rejected_mixed + iod_cov_pass_mixed)
)
iod_validation_true_rejection_pct = _pct(
    Float64(iod_validation_rejected_same),
    Float64(iod_validation_rejected_same + iod_validation_confirmed_same)
)
iod_validation_wrong_rejection_pct = _pct(
    Float64(iod_validation_rejected_mixed),
    Float64(iod_validation_rejected_mixed + iod_validation_confirmed_mixed)
)

function _iod_composition_count(stage::String, outcome::String, composition::String)::Int
    nrow(iod_diagnostics_df) == 0 && return 0
    :composition in propertynames(iod_diagnostics_df) || return 0
    return sum(
        (iod_diagnostics_df.stage .== stage) .&
        (iod_diagnostics_df.outcome .== outcome) .&
        (iod_diagnostics_df.composition .== composition)
    )
end

iod_group_mixed_real = _iod_composition_count("grouping", "formed", "mixed_real_target")
iod_group_false_contaminated = _iod_composition_count("grouping", "formed", "false_contaminated")
iod_group_false_only = _iod_composition_count("grouping", "formed", "false_only")
iod_cov_rejected_false_contaminated = _iod_composition_count(
    "covariance_gate", "rejected", "false_contaminated"
)
iod_cov_rejected_false_only = _iod_composition_count(
    "covariance_gate", "rejected", "false_only"
)
iod_validation_rejected_false_contaminated = _iod_composition_count(
    "next_step_validation", "rejected_score", "false_contaminated"
) + _iod_composition_count(
    "next_step_validation", "rejected_invalid_state", "false_contaminated"
)
iod_validation_rejected_false_only = _iod_composition_count(
    "next_step_validation", "rejected_score", "false_only"
) + _iod_composition_count(
    "next_step_validation", "rejected_invalid_state", "false_only"
)
iod_initialized_false_contaminated = _iod_composition_count(
    "next_step_validation", "promoted", "false_contaminated"
) + _iod_composition_count(
    "direct_promotion", "promoted", "false_contaminated"
)
iod_initialized_false_only = _iod_composition_count(
    "next_step_validation", "promoted", "false_only"
) + _iod_composition_count(
    "direct_promotion", "promoted", "false_only"
)
fake_track_rows = nrow(track_lifecycle_df) == 0 || !(:iod_group_composition in propertynames(track_lifecycle_df)) ?
    falses(nrow(track_lifecycle_df)) :
    in.(track_lifecycle_df.iod_group_composition, Ref(("false_contaminated", "false_only"))) .&
    isfinite.(track_lifecycle_df.filter_initialized_t_s)
fake_track_count = count(fake_track_rows)
fake_track_duration_values = nrow(track_lifecycle_df) == 0 ? Float64[] : Float64[
    Float64(value) for value in track_lifecycle_df.filter_duration_s[fake_track_rows]
    if isfinite(Float64(value))
]
fake_track_duration_stats = _error_stats(fake_track_duration_values)

println("Cross-observer M2M / IOD diagnostics:")
println("  candidate_pairs=$(xm2m_candidate_pairs)")
println("  gate_pass=$(xm2m_gate_pass), same_target=$(xm2m_gate_pass_same), mixed_target=$(xm2m_gate_pass_mixed)")
println("  selected_pairs=$(xm2m_selected), same_target=$(xm2m_selected_same), mixed_target=$(xm2m_selected_mixed)")
println("  iod_groups=$(iod_group_total), same_target=$(iod_group_same), mixed_target=$(iod_group_mixed)")
println("  iod_position_cov_gate: threshold_rms_m=$(IOD_MAX_POSITION_RMS_STD_M), evaluated=$(iod_pos_cov_gate_evaluated), passed_same_target=$(iod_cov_pass_same), passed_mixed_target=$(iod_cov_pass_mixed), rejected_same_target=$(iod_pos_cov_gate_rejected_same), rejected_mixed_target=$(iod_pos_cov_gate_rejected_mixed)")
println("  iod_one_step_validation: enabled=$(ENABLE_IOD_ONE_STEP_VALIDATION), threshold_d2=$(IOD_VALIDATION_MAHAL_MAX_D2), attempted=$(iod_validation_attempted), confirmed=$(iod_validation_confirmed), rejected=$(iod_validation_rejected), no_measure=$(iod_validation_no_measure), pending_end=$(iod_validation_pending_end)")
println("    confirmed_same_target=$(iod_validation_confirmed_same), confirmed_mixed_target=$(iod_validation_confirmed_mixed), rejected_same_target=$(iod_validation_rejected_same), rejected_mixed_target=$(iod_validation_rejected_mixed)")
println("  iod_initialized=$(iod_init_total), same_target=$(iod_init_same), mixed_target=$(iod_init_mixed)")
println("  false-alarm IOD: groups_contaminated=$(iod_group_false_contaminated), groups_false_only=$(iod_group_false_only), initialized_contaminated=$(iod_initialized_false_contaminated), initialized_false_only=$(iod_initialized_false_only), fake_tracks=$(fake_track_count)")
println("  iod_error_after_cov_gate_mean_m: same_target=$(iod_cov_pass_same_error.mean), mixed_target=$(iod_cov_pass_mixed_error.mean)")
println("  iod_error_after_promotion_mean_m: same_target=$(iod_promoted_same_error.mean), mixed_target=$(iod_promoted_mixed_error.mean)")
_postprocess_checkpoint("iod_statistics")

group_total = fusion_model.consensus_group_total
group_same = fusion_model.consensus_group_same_target_total
group_mixed = fusion_model.consensus_group_mixed_target_total
group_unknown = fusion_model.consensus_group_unknown_total
group_known = group_same + group_mixed
group_same_pct_known = group_known > 0 ? 100.0 * group_same / group_known : NaN
group_same_pct_all = group_total > 0 ? 100.0 * group_same / group_total : NaN

tracking_possible_windows = nrow(tracking_window_df)
tracking_tracked_windows = tracking_possible_windows > 0 ? count(tracking_window_df.tracked) : 0
tracking_successful_windows = tracking_possible_windows > 0 ? count(tracking_window_df.success_under_1km) : 0
tracking_coverage_pct = _pct(Float64(tracking_tracked_windows), Float64(tracking_possible_windows))
tracking_success_rate_tracked_pct = _pct(Float64(tracking_successful_windows), Float64(tracking_tracked_windows))
tracking_success_rate_possible_pct = _pct(Float64(tracking_successful_windows), Float64(tracking_possible_windows))
tracking_error_values = [
    Float64(row.mean_error_m) for row in eachrow(tracking_window_df)
    if row.tracked && isfinite(Float64(row.mean_error_m))
]
tracking_mean_error_m = isempty(tracking_error_values) ? NaN : sum(tracking_error_values) / length(tracking_error_values)
tracking_good_error_values = [
    Float64(row.mean_error_m) for row in eachrow(tracking_window_df)
    if row.success_under_1km && isfinite(Float64(row.mean_error_m))
]
tracking_mean_good_error_m = isempty(tracking_good_error_values) ? NaN : sum(tracking_good_error_values) / length(tracking_good_error_values)
tracking_weighted_rows = [
    (
        samples=Int(row.estimate_samples),
        mean_error_m=Float64(row.mean_error_m),
        rmse_error_m=Float64(row.rmse_error_m)
    )
    for row in eachrow(tracking_window_df)
    if row.tracked &&
        Int(row.estimate_samples) > 0 &&
        isfinite(Float64(row.mean_error_m)) &&
        isfinite(Float64(row.rmse_error_m))
]
tracking_sample_count = isempty(tracking_weighted_rows) ? 0 : sum(row.samples for row in tracking_weighted_rows)
tracking_weighted_error_sum_m = isempty(tracking_weighted_rows) ? 0.0 : sum(row.mean_error_m * row.samples for row in tracking_weighted_rows)
tracking_weighted_sq_error_sum_m2 = isempty(tracking_weighted_rows) ? 0.0 : sum(row.rmse_error_m^2 * row.samples for row in tracking_weighted_rows)
tracking_sample_weighted_mean_error_m = tracking_sample_count > 0 ? tracking_weighted_error_sum_m / tracking_sample_count : NaN
tracking_sample_rmse_error_m = tracking_sample_count > 0 ? sqrt(tracking_weighted_sq_error_sum_m2 / tracking_sample_count) : NaN
estimate_duration_values = [
    Float64(row.estimate_duration_s) for row in eachrow(tracking_window_df)
    if row.tracked && isfinite(Float64(row.estimate_duration_s))
]
estimate_good_duration_values = [
    Float64(row.estimate_duration_s) for row in eachrow(tracking_window_df)
    if row.success_under_1km && isfinite(Float64(row.estimate_duration_s))
]
tracking_mean_estimate_duration_s = isempty(estimate_duration_values) ? NaN : sum(estimate_duration_values) / length(estimate_duration_values)
tracking_max_estimate_duration_s = isempty(estimate_duration_values) ? NaN : maximum(estimate_duration_values)
tracking_mean_good_estimate_duration_s = isempty(estimate_good_duration_values) ? NaN : sum(estimate_good_duration_values) / length(estimate_good_duration_values)
detected_unique_targets = count(
    target_id -> any(observer_id -> optical_sensor.counts[observer_id, target_id] > 0, observer_idxs),
    target_idxs
)
jointly_detected_unique_targets = length(unique([Int(row.target) for row in eachrow(tracking_window_df)]))
tracked_unique_targets = length(unique([
    Int(row.target) for row in eachrow(tracking_window_df)
    if row.tracked
]))
successful_tracked_unique_targets = length(unique([
    Int(row.target) for row in eachrow(tracking_window_df)
    if row.success_under_1km
]))
trackable_target_ids = Set(Int(row.target) for row in eachrow(tracking_window_df))
initialized_target_ids = Set(
    Int(row.first_target_id) for row in eachrow(track_lifecycle_df)
    if isfinite(Float64(row.filter_initialized_t_s)) &&
        Int(row.first_target_id) in target_idxs
)
correctly_initialized_target_ids = Set(
    Int(row.first_target_id) for row in eachrow(track_lifecycle_df)
    if String(row.identity_class) == "good" &&
        Int(row.first_target_id) in target_idxs
)
initialization_metrics_applicable = !USE_CENTRALIZED_ORACLE
initialized_trackable_target_ids = intersect(initialized_target_ids, trackable_target_ids)
correctly_initialized_trackable_target_ids = intersect(
    correctly_initialized_target_ids,
    trackable_target_ids
)
never_initialized_unique_targets = length(setdiff(trackable_target_ids, initialized_target_ids))
never_correctly_initialized_unique_targets = length(setdiff(
    trackable_target_ids,
    correctly_initialized_target_ids
))
initialization_coverage_pct = initialization_metrics_applicable ? _pct(
    Float64(length(initialized_trackable_target_ids)),
    Float64(length(trackable_target_ids))
) : NaN
initialization_success_pct = initialization_metrics_applicable ? _pct(
    Float64(length(correctly_initialized_trackable_target_ids)),
    Float64(length(trackable_target_ids))
) : NaN
track_count = nrow(track_lifecycle_df)
track_closed_count = track_count > 0 ? count(track_lifecycle_df.final_status .== "closed") : 0
track_id_switch_total = track_count > 0 ? sum(Int(row.id_switch_count) for row in eachrow(track_lifecycle_df)) : 0
tracks_with_id_switch = track_count > 0 ? count(track_lifecycle_df.id_switch_count .> 0) : 0
track_duration_values = [
    Float64(row.duration_s) for row in eachrow(track_lifecycle_df)
    if isfinite(Float64(row.duration_s))
]
track_filter_duration_values = [
    Float64(row.filter_duration_s) for row in eachrow(track_lifecycle_df)
    if isfinite(Float64(row.filter_duration_s))
]
track_mean_duration_s = isempty(track_duration_values) ? NaN : sum(track_duration_values) / length(track_duration_values)
track_max_duration_s = isempty(track_duration_values) ? NaN : maximum(track_duration_values)
track_mean_filter_duration_s = isempty(track_filter_duration_values) ? NaN : sum(track_filter_duration_values) / length(track_filter_duration_values)
track_max_filter_duration_s = isempty(track_filter_duration_values) ? NaN : maximum(track_filter_duration_values)
good_track_filter_duration_values = Float64[
    Float64(row.filter_duration_s) for row in eachrow(track_lifecycle_df)
    if String(row.identity_class) == "good" && isfinite(Float64(row.filter_duration_s))
]
bad_track_filter_duration_values = Float64[
    Float64(row.filter_duration_s) for row in eachrow(track_lifecycle_df)
    if String(row.identity_class) == "bad" && isfinite(Float64(row.filter_duration_s))
]
good_track_duration_stats = _error_stats(good_track_filter_duration_values)
bad_track_duration_stats = _error_stats(bad_track_filter_duration_values)
initialization_latency_values = [
    Float64(row.initialization_latency_s) for row in eachrow(track_lifecycle_df)
    if isfinite(Float64(row.initialization_latency_s))
]
initialization_error_values = [
    Float64(row.initialization_position_error_m) for row in eachrow(track_lifecycle_df)
    if isfinite(Float64(row.initialization_position_error_m))
]
converged_error_values = [
    Float64(row.converged_mean_error_m) for row in eachrow(tracking_window_df)
    if isfinite(Float64(row.converged_mean_error_m))
]
initialization_latency_mean_s = isempty(initialization_latency_values) ? NaN : mean(initialization_latency_values)
initialization_latency_median_s = isempty(initialization_latency_values) ? NaN : median(initialization_latency_values)
initialization_latency_p95_s = isempty(initialization_latency_values) ? NaN : quantile(initialization_latency_values, 0.95)
initialization_error_mean_m = isempty(initialization_error_values) ? NaN : mean(initialization_error_values)
initialization_error_median_m = isempty(initialization_error_values) ? NaN : median(initialization_error_values)
converged_error_mean_m = isempty(converged_error_values) ? NaN : mean(converged_error_values)
converged_error_median_m = isempty(converged_error_values) ? NaN : median(converged_error_values)

visible_count_values = nrow(geometry_difficulty_df) == 0 ? Float64[] :
    Float64.(geometry_difficulty_df.simultaneously_visible_targets)
angular_separation_values = nrow(geometry_difficulty_df) == 0 ? Float64[] : [
    Float64(value) for value in geometry_difficulty_df.minimum_angular_separation_deg
    if isfinite(Float64(value))
]
geometry_mean_visible = isempty(visible_count_values) ? NaN : mean(visible_count_values)
geometry_max_visible = isempty(visible_count_values) ? NaN : maximum(visible_count_values)
geometry_min_separation_deg = isempty(angular_separation_values) ? NaN : minimum(angular_separation_values)
geometry_p05_separation_deg = isempty(angular_separation_values) ? NaN : quantile(angular_separation_values, 0.05)
wrong_association_total = meas_total_wrong + tt_commit_wrong + hyp_promoted_mixed_target + iod_init_mixed + group_mixed
identity_anomaly_total = wrong_association_total + track_id_switch_total
println("POSTPROCESS_DONE")
flush(stdout)
flush(stderr)

summary_rows = DataFrame(
    section=String[],
    metric=String[],
    value=Float64[]
)
push!(summary_rows, ("run_config", "baseline_theta_gate_rad", USE_BASELINE_DA ? BASELINE_ORPHAN_ATTACH_MAX_DTHETA_RAD : NaN))
push!(summary_rows, ("run_config", "data_association_enabled", USE_NO_DA ? 0.0 : 1.0))
push!(summary_rows, ("run_config", "direct_target_ids", USE_DIRECT_TARGET_IDS ? 1.0 : 0.0))
push!(summary_rows, ("run_config", "save_iod_event_geometry", SAVE_IOD_EVENT_GEOMETRY ? 1.0 : 0.0))
push!(summary_rows, ("run_config", "save_iod_pairwise", SAVE_IOD_PAIRWISE_DIAGNOSTICS ? 1.0 : 0.0))
push!(summary_rows, ("run_config", "mc_metrics_schema_version", 5.0))
push!(summary_rows, ("run_config", "tracking_windows_truth_visibility", 1.0))
push!(summary_rows, ("run_config", "start_epoch_unix_s", datetime2unix(START_UTC)))
push!(summary_rows, ("run_config", "scenario_seed", Float64(SCENARIO_RNG_SEED)))
push!(summary_rows, ("run_config", "sensor_seed", Float64(SENSOR_RNG_SEED)))
push!(summary_rows, ("run_config", "observer_od_seed", Float64(OBSERVER_OD_RNG_SEED)))
push!(summary_rows, ("run_config", "measurement_bias_seed", Float64(MEASUREMENT_BIAS_RNG_SEED)))
push!(summary_rows, ("run_config", "misdetection_seed", Float64(MISDETECTION_RNG_SEED)))
push!(summary_rows, ("run_config", "false_alarm_count_seed", Float64(FALSE_ALARM_COUNT_RNG_SEED)))
push!(summary_rows, ("run_config", "false_alarm_direction_seed", Float64(FALSE_ALARM_DIRECTION_RNG_SEED)))
push!(summary_rows, ("run_config", "m2t_gate_d2", MEAS_ASSOC_MAHAL_MAX_D2))
push!(summary_rows, ("run_config", "m2t_ratio_min", MEAS_ASSOC_DISAMBIG_RATIO_MIN))
push!(summary_rows, ("run_config", "t2t_iod_gate_d2", CONSENSUS_MATCH_MAHAL_IOD_MAX_D2))
push!(summary_rows, ("run_config", "t2t_filter_gate_d2", CONSENSUS_MATCH_MAHAL_FILTER_MAX_D2))
push!(summary_rows, ("run_config", "t2t_ratio_min", TRACK_ASSOC_DISAMBIG_RATIO_MIN))
push!(summary_rows, ("run_config", "iod_position_rms_gate_m", IOD_MAX_POSITION_RMS_STD_M))
push!(summary_rows, ("run_config", "iod_validation_gate_d2", IOD_VALIDATION_MAHAL_MAX_D2))

push!(summary_rows, ("target_population", "target_count", Float64(length(target_defs))))
push!(summary_rows, ("target_population", "realized_a_mean_km", mean(target_population_df.semi_major_axis_km)))
push!(summary_rows, ("target_population", "realized_a_std_km", std(target_population_df.semi_major_axis_km)))
push!(summary_rows, ("target_population", "realized_log10_e_mean", mean(target_population_df.log10_eccentricity)))
push!(summary_rows, ("target_population", "realized_log10_e_std", std(target_population_df.log10_eccentricity)))
push!(summary_rows, ("target_population", "realized_i_mean_deg", mean(target_population_df.inclination_deg)))
push!(summary_rows, ("target_population", "realized_i_std_deg", std(target_population_df.inclination_deg)))
push!(summary_rows, ("target_population", "minimum_periapsis_altitude_km", minimum(target_population_df.periapsis_altitude_km)))

push!(summary_rows, ("sensor_errors", "enable_missed_detections", ENABLE_MISDETECTIONS ? 1.0 : 0.0))
push!(summary_rows, ("sensor_errors", "enable_false_alarms", ENABLE_FALSE_ALARMS ? 1.0 : 0.0))
push!(summary_rows, ("sensor_errors", "enable_measurement_bias", ENABLE_MEASUREMENT_BIAS ? 1.0 : 0.0))
push!(summary_rows, ("sensor_errors", "configured_sigma_theta_rad", SIGMA_THETA_RAD))
push!(summary_rows, ("sensor_errors", "configured_missed_detection_rate", MISDETECTION_RATE))
push!(summary_rows, ("sensor_errors", "configured_false_alarm_rate", FALSE_ALARM_RATE))
push!(summary_rows, ("sensor_errors", "configured_false_alarm_mean_per_observer_epoch", ENABLE_FALSE_ALARMS ? FALSE_ALARM_RATE : 0.0))
push!(summary_rows, ("sensor_errors", "configured_false_alarm_nonzero_epoch_probability_pct", configured_false_alarm_nonzero_epoch_pct))
push!(summary_rows, ("sensor_errors", "configured_false_alarm_multiple_epoch_probability_pct", configured_false_alarm_multiple_epoch_pct))
push!(summary_rows, ("sensor_errors", "configured_measurement_bias_rad", MEASUREMENT_BIAS_RAD))
push!(summary_rows, ("sensor_errors", "visible_opportunities", Float64(sensor_visible_opportunities)))
push!(summary_rows, ("sensor_errors", "true_detections", Float64(sensor_true_detections)))
push!(summary_rows, ("sensor_errors", "missed_detections", Float64(sensor_missed_detections)))
push!(summary_rows, ("sensor_errors", "false_alarms", Float64(sensor_false_alarms)))
push!(summary_rows, ("sensor_errors", "observer_epochs", Float64(sensor_observer_epochs)))
push!(summary_rows, ("sensor_errors", "false_alarm_nonzero_epochs", Float64(sensor_false_alarm_nonzero_epochs)))
push!(summary_rows, ("sensor_errors", "false_alarm_multiple_epochs", Float64(sensor_false_alarm_multiple_epochs)))
push!(summary_rows, ("sensor_errors", "false_alarm_max_per_epoch", Float64(sensor_false_alarm_max_per_epoch)))
push!(summary_rows, ("sensor_errors", "realized_false_alarm_mean_per_observer_epoch", sensor_false_alarm_mean_per_observer_epoch))
push!(summary_rows, ("sensor_errors", "realized_false_alarm_nonzero_epoch_pct", sensor_false_alarm_nonzero_epoch_pct))
push!(summary_rows, ("sensor_errors", "realized_false_alarm_multiple_epoch_pct", sensor_false_alarm_multiple_epoch_pct))
push!(summary_rows, ("sensor_errors", "biased_measurements", Float64(sensor_biased_measurements)))
push!(summary_rows, ("sensor_errors", "realized_bias_component_rms_rad", sensor_realized_bias_component_rms_rad))
push!(summary_rows, ("sensor_errors", "realized_bias_norm_max_rad", sensor_realized_bias_norm_max_rad))
push!(summary_rows, ("sensor_errors", "realized_detection_rate_pct", sensor_detection_rate_pct))
push!(summary_rows, ("sensor_errors", "realized_missed_rate_pct", sensor_missed_rate_pct))

push!(summary_rows, ("false_alarm", "generated", Float64(sensor_false_alarms)))
push!(summary_rows, ("false_alarm", "committed_to_m2t", Float64(meas_false_alarm_committed)))
push!(summary_rows, ("false_alarm", "seed_tracks_created", Float64(tracks_created_with_nonreal_measurements)))
push!(summary_rows, ("false_alarm", "iod_groups_mixed_real_target", Float64(iod_group_mixed_real)))
push!(summary_rows, ("false_alarm", "iod_groups_false_contaminated", Float64(iod_group_false_contaminated)))
push!(summary_rows, ("false_alarm", "iod_groups_false_only", Float64(iod_group_false_only)))
push!(summary_rows, ("false_alarm", "iod_cov_rejected_false_contaminated", Float64(iod_cov_rejected_false_contaminated)))
push!(summary_rows, ("false_alarm", "iod_cov_rejected_false_only", Float64(iod_cov_rejected_false_only)))
push!(summary_rows, ("false_alarm", "iod_validation_rejected_false_contaminated", Float64(iod_validation_rejected_false_contaminated)))
push!(summary_rows, ("false_alarm", "iod_validation_rejected_false_only", Float64(iod_validation_rejected_false_only)))
push!(summary_rows, ("false_alarm", "iod_initialized_false_contaminated", Float64(iod_initialized_false_contaminated)))
push!(summary_rows, ("false_alarm", "iod_initialized_false_only", Float64(iod_initialized_false_only)))
push!(summary_rows, ("false_alarm", "fake_tracks_initialized", Float64(fake_track_count)))
push!(summary_rows, ("false_alarm", "fake_track_duration_mean_s", fake_track_duration_stats.mean))
push!(summary_rows, ("false_alarm", "fake_track_duration_median_s", fake_track_duration_stats.median))
push!(summary_rows, ("false_alarm", "fake_track_duration_p95_s", fake_track_duration_stats.p95))

push!(summary_rows, ("meas_assoc", "committed_total", Float64(meas_total_committed)))
push!(summary_rows, ("meas_assoc", "committed_correct", Float64(meas_total_correct)))
push!(summary_rows, ("meas_assoc", "committed_wrong", Float64(meas_total_wrong)))
push!(summary_rows, ("meas_assoc", "false_alarm_committed", Float64(meas_false_alarm_committed)))
push!(summary_rows, ("meas_assoc", "commit_accuracy_pct", meas_commit_acc_pct))
push!(summary_rows, ("meas_assoc", "true_opportunities", Float64(meas_true_opportunities)))
push!(summary_rows, ("meas_assoc", "true_not_committed", Float64(meas_true_not_committed)))
push!(summary_rows, ("meas_assoc", "recall_pct", meas_recall_pct))
push!(summary_rows, ("meas_assoc", "unique_committed", Float64(meas_unique_total)))
push!(summary_rows, ("meas_assoc", "unique_correct", Float64(meas_unique_correct)))
push!(summary_rows, ("meas_assoc", "unique_wrong", Float64(meas_unique_wrong)))
push!(summary_rows, ("meas_assoc", "unique_accuracy_pct", meas_unique_acc_pct))
push!(summary_rows, ("meas_assoc", "ambig_committed", Float64(meas_ambig_total)))
push!(summary_rows, ("meas_assoc", "ambig_correct", Float64(meas_ambig_correct)))
push!(summary_rows, ("meas_assoc", "ambig_wrong", Float64(meas_ambig_wrong)))
push!(summary_rows, ("meas_assoc", "ambig_dropped", Float64(meas_ambig_dropped)))
push!(summary_rows, ("meas_assoc", "ambig_accuracy_pct", meas_ambig_acc_pct))

push!(summary_rows, ("local_hypothesis", "H1_created", Float64(hyp_h1_created)))
push!(summary_rows, ("local_hypothesis", "H1_to_H2_created", Float64(hyp_h1_to_h2_created)))
push!(summary_rows, ("local_hypothesis", "H1_to_H2_same_target", Float64(hyp_h1_to_h2_same_target)))
push!(summary_rows, ("local_hypothesis", "H1_to_H2_mixed_target", Float64(hyp_h1_to_h2_mixed_target)))
push!(summary_rows, ("local_hypothesis", "H2_to_H3_attempted", Float64(hyp_h2_to_h3_attempted)))
push!(summary_rows, ("local_hypothesis", "H3_los_rate_pass", Float64(hyp_h3_los_rate_pass)))
push!(summary_rows, ("local_hypothesis", "H3_los_rate_pass_same_target", Float64(hyp_h3_los_rate_pass_same_target)))
push!(summary_rows, ("local_hypothesis", "H3_los_rate_pass_mixed_target", Float64(hyp_h3_los_rate_pass_mixed_target)))
push!(summary_rows, ("local_hypothesis", "H3_los_rate_fail", Float64(hyp_h3_los_rate_fail)))
push!(summary_rows, ("local_hypothesis", "H3_los_rate_fail_same_target", Float64(hyp_h3_los_rate_fail_same_target)))
push!(summary_rows, ("local_hypothesis", "H3_los_rate_fail_mixed_target", Float64(hyp_h3_los_rate_fail_mixed_target)))
push!(summary_rows, ("local_hypothesis", "promoted", Float64(hyp_promoted)))
push!(summary_rows, ("local_hypothesis", "promoted_same_target", Float64(hyp_promoted_same_target)))
push!(summary_rows, ("local_hypothesis", "promoted_mixed_target", Float64(hyp_promoted_mixed_target)))
push!(summary_rows, ("local_hypothesis", "tracks_created_with_nonreal_measurements", Float64(tracks_created_with_nonreal_measurements)))

push!(summary_rows, ("cross_m2m", "candidate_pairs", Float64(xm2m_candidate_pairs)))
push!(summary_rows, ("cross_m2m", "gate_pass", Float64(xm2m_gate_pass)))
push!(summary_rows, ("cross_m2m", "gate_pass_same_target", Float64(xm2m_gate_pass_same)))
push!(summary_rows, ("cross_m2m", "gate_pass_mixed_target", Float64(xm2m_gate_pass_mixed)))
push!(summary_rows, ("cross_m2m", "gate_pass_same_target_pct", xm2m_gate_pass_same_pct))
push!(summary_rows, ("cross_m2m", "gate_fail", Float64(xm2m_gate_fail)))
push!(summary_rows, ("cross_m2m", "gate_fail_same_target", Float64(xm2m_gate_fail_same)))
push!(summary_rows, ("cross_m2m", "gate_fail_mixed_target", Float64(xm2m_gate_fail_mixed)))
push!(summary_rows, ("cross_m2m", "selected_pairs", Float64(xm2m_selected)))
push!(summary_rows, ("cross_m2m", "selected_pairs_same_target", Float64(xm2m_selected_same)))
push!(summary_rows, ("cross_m2m", "selected_pairs_mixed_target", Float64(xm2m_selected_mixed)))
push!(summary_rows, ("cross_m2m", "selected_pairs_same_target_pct", xm2m_selected_same_pct))
push!(summary_rows, ("cross_m2m", "iod_groups", Float64(iod_group_total)))
push!(summary_rows, ("cross_m2m", "iod_groups_same_target", Float64(iod_group_same)))
push!(summary_rows, ("cross_m2m", "iod_groups_mixed_target", Float64(iod_group_mixed)))
push!(summary_rows, ("cross_m2m", "iod_groups_same_target_pct", iod_group_same_pct))
push!(summary_rows, ("cross_m2m", "iod_position_cov_gate_threshold_rms_m", IOD_MAX_POSITION_RMS_STD_M))
push!(summary_rows, ("cross_m2m", "iod_position_cov_gate_evaluated", Float64(iod_pos_cov_gate_evaluated)))
push!(summary_rows, ("cross_m2m", "iod_position_cov_gate_rejected", Float64(iod_pos_cov_gate_rejected)))
push!(summary_rows, ("cross_m2m", "iod_position_cov_gate_rejected_same_target", Float64(iod_pos_cov_gate_rejected_same)))
push!(summary_rows, ("cross_m2m", "iod_position_cov_gate_rejected_mixed_target", Float64(iod_pos_cov_gate_rejected_mixed)))
push!(summary_rows, ("cross_m2m", "iod_position_cov_gate_rejected_pct", iod_pos_cov_gate_rejected_pct))
push!(summary_rows, ("cross_m2m", "iod_position_cov_gate_passed_same_target", Float64(iod_cov_pass_same)))
push!(summary_rows, ("cross_m2m", "iod_position_cov_gate_passed_mixed_target", Float64(iod_cov_pass_mixed)))
push!(summary_rows, ("cross_m2m", "iod_position_cov_gate_true_rejection_pct", iod_cov_true_rejection_pct))
push!(summary_rows, ("cross_m2m", "iod_position_cov_gate_wrong_rejection_pct", iod_cov_wrong_rejection_pct))
push!(summary_rows, ("cross_m2m", "iod_after_cov_gate_same_target_error_mean_m", iod_cov_pass_same_error.mean))
push!(summary_rows, ("cross_m2m", "iod_after_cov_gate_same_target_error_median_m", iod_cov_pass_same_error.median))
push!(summary_rows, ("cross_m2m", "iod_after_cov_gate_same_target_error_p95_m", iod_cov_pass_same_error.p95))
push!(summary_rows, ("cross_m2m", "iod_after_cov_gate_mixed_target_error_mean_m", iod_cov_pass_mixed_error.mean))
push!(summary_rows, ("cross_m2m", "iod_after_cov_gate_mixed_target_error_median_m", iod_cov_pass_mixed_error.median))
push!(summary_rows, ("cross_m2m", "iod_after_cov_gate_mixed_target_error_p95_m", iod_cov_pass_mixed_error.p95))
push!(summary_rows, ("cross_m2m", "iod_one_step_validation_enabled", ENABLE_IOD_ONE_STEP_VALIDATION ? 1.0 : 0.0))
push!(summary_rows, ("cross_m2m", "iod_validation_threshold_d2", IOD_VALIDATION_MAHAL_MAX_D2))
push!(summary_rows, ("cross_m2m", "iod_validation_attempted", Float64(iod_validation_attempted)))
push!(summary_rows, ("cross_m2m", "iod_validation_confirmed", Float64(iod_validation_confirmed)))
push!(summary_rows, ("cross_m2m", "iod_validation_rejected", Float64(iod_validation_rejected)))
push!(summary_rows, ("cross_m2m", "iod_validation_confirmed_pct", iod_validation_confirmed_pct))
push!(summary_rows, ("cross_m2m", "iod_validation_confirmed_same_target", Float64(iod_validation_confirmed_same)))
push!(summary_rows, ("cross_m2m", "iod_validation_confirmed_mixed_target", Float64(iod_validation_confirmed_mixed)))
push!(summary_rows, ("cross_m2m", "iod_validation_rejected_same_target", Float64(iod_validation_rejected_same)))
push!(summary_rows, ("cross_m2m", "iod_validation_rejected_mixed_target", Float64(iod_validation_rejected_mixed)))
push!(summary_rows, ("cross_m2m", "iod_validation_no_measure", Float64(iod_validation_no_measure)))
push!(summary_rows, ("cross_m2m", "iod_validation_no_measure_same_target", Float64(iod_validation_no_measure_same)))
push!(summary_rows, ("cross_m2m", "iod_validation_no_measure_mixed_target", Float64(iod_validation_no_measure_mixed)))
push!(summary_rows, ("cross_m2m", "iod_validation_true_rejection_pct", iod_validation_true_rejection_pct))
push!(summary_rows, ("cross_m2m", "iod_validation_wrong_rejection_pct", iod_validation_wrong_rejection_pct))
push!(summary_rows, ("cross_m2m", "iod_promoted_same_target_error_mean_m", iod_promoted_same_error.mean))
push!(summary_rows, ("cross_m2m", "iod_promoted_same_target_error_median_m", iod_promoted_same_error.median))
push!(summary_rows, ("cross_m2m", "iod_promoted_same_target_error_p95_m", iod_promoted_same_error.p95))
push!(summary_rows, ("cross_m2m", "iod_promoted_mixed_target_error_mean_m", iod_promoted_mixed_error.mean))
push!(summary_rows, ("cross_m2m", "iod_promoted_mixed_target_error_median_m", iod_promoted_mixed_error.median))
push!(summary_rows, ("cross_m2m", "iod_promoted_mixed_target_error_p95_m", iod_promoted_mixed_error.p95))
push!(summary_rows, ("cross_m2m", "iod_validation_rejected_same_target_error_mean_m", iod_rejected_same_error.mean))
push!(summary_rows, ("cross_m2m", "iod_validation_rejected_mixed_target_error_mean_m", iod_rejected_mixed_error.mean))
push!(summary_rows, ("cross_m2m", "iod_validation_pending_end", Float64(iod_validation_pending_end)))
push!(summary_rows, ("cross_m2m", "iod_initialized", Float64(iod_init_total)))
push!(summary_rows, ("cross_m2m", "iod_initialized_same_target", Float64(iod_init_same)))
push!(summary_rows, ("cross_m2m", "iod_initialized_mixed_target", Float64(iod_init_mixed)))
push!(summary_rows, ("cross_m2m", "iod_initialized_same_target_pct", iod_init_same_pct))

push!(summary_rows, ("track_assoc", "tt_attempt_total", Float64(tt_attempt_total)))
push!(summary_rows, ("track_assoc", "tt_committed_total", Float64(tt_commit_total)))
push!(summary_rows, ("track_assoc", "tt_committed_correct", Float64(tt_commit_correct)))
push!(summary_rows, ("track_assoc", "tt_committed_wrong", Float64(tt_commit_wrong)))
push!(summary_rows, ("track_assoc", "tt_committed_unknown", Float64(tt_commit_unknown)))
push!(summary_rows, ("track_assoc", "tt_skipped", Float64(tt_skipped)))
push!(summary_rows, ("track_assoc", "tt_skipped_no_candidate", Float64(tt_no_candidate)))
push!(summary_rows, ("track_assoc", "tt_skipped_gate_fail", Float64(tt_gate_fail)))
push!(summary_rows, ("track_assoc", "tt_skipped_ratio_fail", Float64(tt_ratio_fail)))
push!(summary_rows, ("track_assoc", "tt_skipped_mutual_fail", Float64(tt_mutual_fail)))
push!(summary_rows, ("track_assoc", "tt_component_conflict_rejected", Float64(tt_component_conflict_rejected)))
push!(summary_rows, ("track_assoc", "tt_skipped_same_target_present", Float64(tt_skipped_same_target_present)))
push!(summary_rows, ("track_assoc", "tt_skipped_no_same_target", Float64(tt_skipped_no_same_target)))
push!(summary_rows, ("track_assoc", "tt_skipped_unknown_source_target", Float64(tt_skipped_unknown_source_target)))
push!(summary_rows, ("track_assoc", "tt_accuracy_pct_known_only", tt_acc_pct))
push!(summary_rows, ("track_assoc", "tt_true_opportunities", Float64(tt_true_opportunities)))
push!(summary_rows, ("track_assoc", "tt_true_not_committed", Float64(tt_true_not_committed)))
push!(summary_rows, ("track_assoc", "tt_recall_pct", tt_recall_pct))
push!(summary_rows, ("track_assoc", "consensus_group_total", Float64(group_total)))
push!(summary_rows, ("track_assoc", "consensus_group_same_target", Float64(group_same)))
push!(summary_rows, ("track_assoc", "consensus_group_mixed_target", Float64(group_mixed)))
push!(summary_rows, ("track_assoc", "consensus_group_unknown_target", Float64(group_unknown)))
push!(summary_rows, ("track_assoc", "consensus_group_same_target_pct_known_only", group_same_pct_known))
push!(summary_rows, ("track_assoc", "consensus_group_same_target_pct_all", group_same_pct_all))

push!(summary_rows, ("tracking", "possible_windows", Float64(tracking_possible_windows)))
push!(summary_rows, ("tracking", "tracked_windows", Float64(tracking_tracked_windows)))
push!(summary_rows, ("tracking", "tracking_coverage_pct", tracking_coverage_pct))
push!(summary_rows, ("tracking", "successful_windows_under_1km", Float64(tracking_successful_windows)))
push!(summary_rows, ("tracking", "success_rate_tracked_pct", tracking_success_rate_tracked_pct))
push!(summary_rows, ("tracking", "success_rate_possible_pct", tracking_success_rate_possible_pct))
push!(summary_rows, ("tracking", "mean_error_tracked_windows_m", tracking_mean_error_m))
push!(summary_rows, ("tracking", "mean_error_successful_windows_m", tracking_mean_good_error_m))
push!(summary_rows, ("tracking", "sample_weighted_mean_error_m", tracking_sample_weighted_mean_error_m))
push!(summary_rows, ("tracking", "sample_rmse_error_m", tracking_sample_rmse_error_m))
push!(summary_rows, ("tracking", "success_error_threshold_m", TRACKING_SUCCESS_ERROR_MAX_M))
push!(summary_rows, ("tracking", "possible_min_joint_ticks", Float64(TRACKING_POSSIBLE_MIN_JOINT_TICKS)))
push!(summary_rows, ("tracking", "mean_estimate_duration_s", tracking_mean_estimate_duration_s))
push!(summary_rows, ("tracking", "max_estimate_duration_s", tracking_max_estimate_duration_s))
push!(summary_rows, ("tracking", "mean_successful_estimate_duration_s", tracking_mean_good_estimate_duration_s))
push!(summary_rows, ("tracking", "converged_window_fraction", 0.8))
push!(summary_rows, ("tracking", "converged_mean_error_m", converged_error_mean_m))
push!(summary_rows, ("tracking", "converged_median_error_m", converged_error_median_m))

push!(summary_rows, ("object_coverage", "detected_unique_targets", Float64(detected_unique_targets)))
push!(summary_rows, ("object_coverage", "jointly_detected_unique_targets", Float64(jointly_detected_unique_targets)))
push!(summary_rows, ("object_coverage", "tracked_unique_targets", Float64(tracked_unique_targets)))
push!(summary_rows, ("object_coverage", "successful_tracked_unique_targets", Float64(successful_tracked_unique_targets)))
push!(summary_rows, ("object_coverage", "tracked_unique_over_detected_pct", _pct(Float64(tracked_unique_targets), Float64(detected_unique_targets))))
push!(summary_rows, ("object_coverage", "successful_unique_over_detected_pct", _pct(Float64(successful_tracked_unique_targets), Float64(detected_unique_targets))))
push!(summary_rows, ("object_coverage", "tracked_unique_over_jointly_detected_pct", _pct(Float64(tracked_unique_targets), Float64(jointly_detected_unique_targets))))
push!(summary_rows, ("object_coverage", "successful_unique_over_jointly_detected_pct", _pct(Float64(successful_tracked_unique_targets), Float64(jointly_detected_unique_targets))))

if ENABLE_NAV_TIMING
    push!(summary_rows, ("runtime", "timed_epoch_count", Float64(runtime_total.count)))
    push!(summary_rows, ("runtime", "local_da_mean_ms", runtime_local_da.mean))
    push!(summary_rows, ("runtime", "local_da_p05_ms", runtime_local_da.p05))
    push!(summary_rows, ("runtime", "local_da_median_ms", runtime_local_da.median))
    push!(summary_rows, ("runtime", "local_da_p95_ms", runtime_local_da.p95))
    push!(summary_rows, ("runtime", "cross_da_mean_ms", runtime_cross_da.mean))
    push!(summary_rows, ("runtime", "cross_da_p05_ms", runtime_cross_da.p05))
    push!(summary_rows, ("runtime", "cross_da_median_ms", runtime_cross_da.median))
    push!(summary_rows, ("runtime", "cross_da_p95_ms", runtime_cross_da.p95))
    push!(summary_rows, ("runtime", "association_mean_ms", runtime_association.mean))
    push!(summary_rows, ("runtime", "association_p05_ms", runtime_association.p05))
    push!(summary_rows, ("runtime", "association_median_ms", runtime_association.median))
    push!(summary_rows, ("runtime", "association_p95_ms", runtime_association.p95))
    push!(summary_rows, ("runtime", "filter_mean_ms", runtime_filter.mean))
    push!(summary_rows, ("runtime", "filter_p05_ms", runtime_filter.p05))
    push!(summary_rows, ("runtime", "filter_median_ms", runtime_filter.median))
    push!(summary_rows, ("runtime", "filter_p95_ms", runtime_filter.p95))
    push!(summary_rows, ("runtime", "fusion_mean_ms", runtime_fusion.mean))
    push!(summary_rows, ("runtime", "fusion_p05_ms", runtime_fusion.p05))
    push!(summary_rows, ("runtime", "fusion_median_ms", runtime_fusion.median))
    push!(summary_rows, ("runtime", "fusion_p95_ms", runtime_fusion.p95))
    push!(summary_rows, ("runtime", "total_epoch_mean_ms", runtime_total.mean))
    push!(summary_rows, ("runtime", "total_epoch_p05_ms", runtime_total.p05))
    push!(summary_rows, ("runtime", "total_epoch_median_ms", runtime_total.median))
    push!(summary_rows, ("runtime", "total_epoch_p95_ms", runtime_total.p95))
    push!(summary_rows, ("runtime", "total_epoch_max_ms", runtime_total.max))
end

push!(summary_rows, ("track_lifecycle", "track_count", Float64(track_count)))
push!(summary_rows, ("track_lifecycle", "closed_count", Float64(track_closed_count)))
push!(summary_rows, ("track_lifecycle", "id_switch_total", Float64(track_id_switch_total)))
push!(summary_rows, ("track_lifecycle", "tracks_with_id_switch", Float64(tracks_with_id_switch)))
push!(summary_rows, ("track_lifecycle", "mean_duration_s", track_mean_duration_s))
push!(summary_rows, ("track_lifecycle", "max_duration_s", track_max_duration_s))
push!(summary_rows, ("track_lifecycle", "mean_filter_duration_s", track_mean_filter_duration_s))
push!(summary_rows, ("track_lifecycle", "max_filter_duration_s", track_max_filter_duration_s))
push!(summary_rows, ("track_lifecycle", "good_track_count", Float64(length(good_track_filter_duration_values))))
push!(summary_rows, ("track_lifecycle", "good_track_filter_duration_mean_s", good_track_duration_stats.mean))
push!(summary_rows, ("track_lifecycle", "good_track_filter_duration_median_s", good_track_duration_stats.median))
push!(summary_rows, ("track_lifecycle", "good_track_filter_duration_p95_s", good_track_duration_stats.p95))
push!(summary_rows, ("track_lifecycle", "bad_track_count", Float64(length(bad_track_filter_duration_values))))
push!(summary_rows, ("track_lifecycle", "bad_track_filter_duration_mean_s", bad_track_duration_stats.mean))
push!(summary_rows, ("track_lifecycle", "bad_track_filter_duration_median_s", bad_track_duration_stats.median))
push!(summary_rows, ("track_lifecycle", "bad_track_filter_duration_p95_s", bad_track_duration_stats.p95))
push!(summary_rows, ("track_lifecycle", "initialization_latency_count", Float64(length(initialization_latency_values))))
push!(summary_rows, (
    "track_lifecycle", "trackable_unique_targets",
    initialization_metrics_applicable ? Float64(length(trackable_target_ids)) : NaN
))
push!(summary_rows, (
    "track_lifecycle", "initialized_unique_targets",
    initialization_metrics_applicable ? Float64(length(initialized_trackable_target_ids)) : NaN
))
push!(summary_rows, (
    "track_lifecycle", "correctly_initialized_unique_targets",
    initialization_metrics_applicable ? Float64(length(correctly_initialized_trackable_target_ids)) : NaN
))
push!(summary_rows, (
    "track_lifecycle", "never_initialized_unique_targets",
    initialization_metrics_applicable ? Float64(never_initialized_unique_targets) : NaN
))
push!(summary_rows, (
    "track_lifecycle", "never_correctly_initialized_unique_targets",
    initialization_metrics_applicable ? Float64(never_correctly_initialized_unique_targets) : NaN
))
push!(summary_rows, ("track_lifecycle", "initialization_coverage_pct", initialization_coverage_pct))
push!(summary_rows, ("track_lifecycle", "initialization_success_pct", initialization_success_pct))
push!(summary_rows, ("track_lifecycle", "initialization_latency_mean_s", initialization_latency_mean_s))
push!(summary_rows, ("track_lifecycle", "initialization_latency_median_s", initialization_latency_median_s))
push!(summary_rows, ("track_lifecycle", "initialization_latency_p95_s", initialization_latency_p95_s))
push!(summary_rows, ("track_lifecycle", "initialization_position_error_mean_m", initialization_error_mean_m))
push!(summary_rows, ("track_lifecycle", "initialization_position_error_median_m", initialization_error_median_m))
push!(summary_rows, ("track_lifecycle", "fragmented_windows", Float64(fragmentation.fragmented_windows)))
push!(summary_rows, ("track_lifecycle", "fragment_excess_total", Float64(fragmentation.fragment_excess_total)))
push!(summary_rows, (
    "track_lifecycle",
    "fragmented_window_rate_pct",
    _pct(Float64(fragmentation.fragmented_windows), Float64(tracking_possible_windows))
))

push!(summary_rows, ("geometry", "observer_epoch_count", Float64(length(visible_count_values))))
push!(summary_rows, ("geometry", "mean_simultaneously_visible_targets", geometry_mean_visible))
push!(summary_rows, ("geometry", "max_simultaneously_visible_targets", geometry_max_visible))
push!(summary_rows, ("geometry", "minimum_angular_separation_deg", geometry_min_separation_deg))
push!(summary_rows, ("geometry", "p05_minimum_angular_separation_deg", geometry_p05_separation_deg))

push!(summary_rows, ("association_health", "wrong_association_total", Float64(wrong_association_total)))
push!(summary_rows, ("association_health", "id_switch_total", Float64(track_id_switch_total)))
push!(summary_rows, ("association_health", "identity_anomaly_total", Float64(identity_anomaly_total)))
push!(summary_rows, ("association_health", "m2t_wrong", Float64(meas_total_wrong)))
push!(summary_rows, ("association_health", "t2t_wrong", Float64(tt_commit_wrong)))
push!(summary_rows, ("association_health", "m2m_promoted_mixed_target", Float64(hyp_promoted_mixed_target)))
push!(summary_rows, ("association_health", "xm2m_iod_initialized_mixed_target", Float64(iod_init_mixed)))
push!(summary_rows, ("association_health", "consensus_mixed_target", Float64(group_mixed)))

println("Association quality table:")
# The complete table is written to CSV below. Avoid rendering all rows into
# every child log when stressed realizations are close to the heap hint.
println("  rows=$(nrow(summary_rows)), columns=$(ncol(summary_rows))")
_postprocess_checkpoint("summary_table")

results_dir = args.simulation_settings.results_directory
mkpath(results_dir)
println("CSV_WRITE_START")
flush(stdout)
flush(stderr)
# Serialize on the container-local temporary filesystem and only then copy the
# completed file into the (potentially bind-mounted) results directory.
function _atomic_csv_write(path::AbstractString, table)
    temporary_path = tempname()
    try
        CSV.write(temporary_path, table)
        cp(temporary_path, path; force=true)
    finally
        isfile(temporary_path) && rm(temporary_path; force=true)
    end
    return path
end

_atomic_csv_write(joinpath(results_dir, "association_quality_table.csv"), summary_rows)
_postprocess_checkpoint("csv_association_quality")
if SAVE_AUXILIARY_METRIC_TABLES
    _atomic_csv_write(joinpath(results_dir, "target_population_table.csv"), target_population_df)
    _postprocess_checkpoint("csv_target_population")
    _atomic_csv_write(joinpath(results_dir, "tracking_observer_summary.csv"), tracking_observer_summary_df)
    GC.gc(false)
    _postprocess_checkpoint("csv_observer_summary")

    # Keep the comparison output lean, while retaining the fields needed to explain
    # tracking duration and sample-weighted errors.
    tracking_window_plot_df = select(
        tracking_window_df,
        [
            :nav_case,
            :observer,
            :target,
            :window_id,
            :window_start_t_s,
            :window_end_t_s,
            :tracked,
            :estimate_samples,
            :estimate_duration_s,
            :success_under_1km,
            :mean_error_m,
            :rmse_error_m,
            :first_estimate_error_m,
            :converged_mean_error_m,
            :converged_rmse_error_m,
            :tracks_per_window,
            :fragment_excess
        ]; copycols=false
    )
    _atomic_csv_write(joinpath(results_dir, "tracking_window_table.csv"), tracking_window_plot_df)
    _postprocess_checkpoint("csv_tracking_windows")

    track_initialization_df = select(
        track_lifecycle_df,
        [
            :nav_case, :observer, :slot, :first_target_id, :last_target_id,
            :first_measurement_t_s, :filter_initialized_t_s, :initialization_latency_s,
            :initialization_position_error_m, :iod_group_same_target, :identity_class,
            :id_switch_count, :filter_duration_s
        ]; copycols=false
    )
    _atomic_csv_write(joinpath(results_dir, "track_initialization_table.csv"), track_initialization_df)
    _postprocess_checkpoint("csv_track_initializations")
    if SAVE_COMPARISON_DETAILED_TABLES
        _atomic_csv_write(joinpath(results_dir, "tracking_window_detail_table.csv"), tracking_window_df)
        _atomic_csv_write(joinpath(results_dir, "track_lifecycle_table.csv"), track_lifecycle_df)
        _atomic_csv_write(joinpath(results_dir, "geometry_difficulty_table.csv"), geometry_difficulty_df)
    end
else
    println("CSV_AUXILIARY_TABLES_SKIPPED")
end

_atomic_csv_write(joinpath(results_dir, "iod_diagnostics_table.csv"), iod_diagnostics_df)
_postprocess_checkpoint("csv_iod_diagnostics")
if SAVE_IOD_PAIRWISE_DIAGNOSTICS
    _atomic_csv_write(
        joinpath(results_dir, "iod_pairwise_consistency_table.csv"),
        iod_pairwise_df
    )
    _postprocess_checkpoint("csv_iod_pairwise_consistency")
end

println("CSV_WRITE_DONE")
flush(stdout)
flush(stderr)
println("COMPUTATIONAL TIME = $(t) s")
