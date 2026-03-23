if !isdefined(@__MODULE__, :SimulationModel)
    include("../simulation_model/SimulationModel.jl")
end
using .SimulationModel
using Dates
using SPICE
using StaticArrays
using LinearAlgebra
using Random
using CSV
using DataFrames

const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :run_simulation)
    include("../simulation/run_simulation.jl")
end
if !isdefined(@__MODULE__, :run_and_report)
    include("typed_example_utils.jl")
end

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
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

Base.@kwdef mutable struct DistributedIODCovModel
    sensor::OpticalLOSSensorModel
    comms::InterAgentCommunicationModel
    target_idx::Int
    observer_idxs::Vector{Int}
    min_neighbor_count::Int
    sigma_theta_rad::Float64
    estimate_state::Vector{SVector{6, Float64}}
    estimate_covariance::Vector{Matrix{Float64}}
    estimate_time_s::Vector{Float64}
    used_neighbors::Vector{Int}
    triangulation_ready::BitVector
    initialized::BitVector
end

function DistributedIODCovModel(
    sensor::OpticalLOSSensorModel,
    comms::InterAgentCommunicationModel,
    target_idx::Int,
    observer_idxs::Vector{Int},
    min_neighbor_count::Int,
    sigma_theta_rad::Float64,
    num_sats::Int
)
    return DistributedIODCovModel(
        sensor=sensor,
        comms=comms,
        target_idx=target_idx,
        observer_idxs=observer_idxs,
        min_neighbor_count=min_neighbor_count,
        sigma_theta_rad=sigma_theta_rad,
        estimate_state=[SVector{6, Float64}(ntuple(_ -> NaN, 6)) for _ in 1:num_sats],
        estimate_covariance=[fill(NaN, 6, 6) for _ in 1:num_sats],
        estimate_time_s=fill(NaN, num_sats),
        used_neighbors=fill(0, num_sats),
        triangulation_ready=falses(num_sats),
        initialized=falses(num_sats)
    )
end

function _print_iod_initialization(
    model::DistributedIODCovModel,
    sat_idx::Int,
    t::Float64,
    z_est::SVector{6, Float64},
    cov::Matrix{Float64},
    u
)
    target_true = SVector{6, Float64}(
        SVector{3, Float64}(u.sc[model.target_idx].pos)...,
        SVector{3, Float64}(u.sc[model.target_idx].vel)...
    )
    err = z_est - target_true
    println("IOD init | tgt=$(model.target_idx) | obs=$(sat_idx) | t=$(round(t; digits=3)) s")
    println("  true   = ", target_true)
    println("  est    = ", z_est)
    println("  err    = ", err)
    println("  |err_r|= $(round(norm(err[1:3]); digits=3)) m, |err_v|= $(round(norm(err[4:6]); digits=6)) m/s")
    println("  cov(6x6) = ")
    show(stdout, "text/plain", cov)
    println("\n")
    return nothing
end

function SimulationModel.calcNavigationEffect!(
    model::DistributedIODCovModel,
    u,
    p,
    t::Float64,
    sat_idx::Int
)
    sat_idx in model.observer_idxs || return nothing # Only compute for observer satellites.
    if sat_idx == model.observer_idxs[1]
        _update_neighbors!(model.comms, u, p)
    end
    model.initialized[sat_idx] && return nothing # Skip if this satellite has already initialized its state estimate.
    if !_has_two_measures(model.sensor, sat_idx, model.target_idx)
        model.triangulation_ready[sat_idx] = false
        model.used_neighbors[sat_idx] = 0
        return nothing
    end
    valid_neighbors = Int[]
    @inbounds for nbr in model.comms.neighbors[sat_idx]
        # Only consider neighbors that are also observers and have valid measurements for the target.
        if nbr in model.observer_idxs && _has_two_measures(model.sensor, nbr, model.target_idx)
            push!(valid_neighbors, nbr)
        end
    end
    model.used_neighbors[sat_idx] = length(valid_neighbors)
    # Require a minimum number of valid neighbors to ensure enough measurements for triangulation.
    if length(valid_neighbors) < model.min_neighbor_count
        model.triangulation_ready[sat_idx] = false
        return nothing
    end
    nodes_idx = vcat([sat_idx], valid_neighbors)
    nodes = NamedTuple[]
    A_rows = Matrix{Float64}(undef, 0, 3)
    b_rows = Float64[]
    dA_rows = Matrix{Float64}(undef, 0, 3)
    db_rows = Float64[]
    for idx in nodes_idx
        key = (idx, model.target_idx)
        haskey(model.sensor.latest, key) || continue
        haskey(model.sensor.previous, key) || continue
        l = model.sensor.latest[key].los_unit
        l_prev = model.sensor.previous[key].los_unit
        lrate = _los_rate(model.sensor, idx, model.target_idx)
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
    length(nodes) >= 2 || return nothing
    H = [A_rows zeros(size(A_rows, 1), 3); dA_rows A_rows]
    y = vcat(b_rows, db_rows)
    size(H, 1) >= 6 || return nothing
    z_est = pinv(H) * y
    x_est = SVector{3, Float64}(z_est[1], z_est[2], z_est[3])
    x_dot_est = SVector{3, Float64}(z_est[4], z_est[5], z_est[6])
    cov = _compute_state_covariance(nodes, x_est, x_dot_est, model.sigma_theta_rad, NAVIGATION_RATE_SEC)
    cov === nothing && return nothing
    model.estimate_state[sat_idx] = SVector{6, Float64}(z_est)
    model.estimate_covariance[sat_idx] = copy(cov)
    model.estimate_time_s[sat_idx] = t
    model.triangulation_ready[sat_idx] = true
    model.initialized[sat_idx] = true
    _print_iod_initialization(model, sat_idx, t, model.estimate_state[sat_idx], cov, u)
    return nothing
end

Base.@kwdef mutable struct DistributedEKFModel
    sensor::OpticalLOSSensorModel
    comms::InterAgentCommunicationModel
    iod::DistributedIODCovModel
    observer_idxs::Vector{Int}
    target_idx::Int
    sigma_theta_rad::Float64
    num_consensus_iter::Int
    process_q_diag::SVector{6, Float64}
    μ::Float64
    state::Vector{SVector{6, Float64}}
    covariance::Vector{Matrix{Float64}}
    state_pred::Vector{SVector{6, Float64}}
    covariance_pred::Vector{Matrix{Float64}}
    initialized::BitVector
    last_update_t::Vector{Float64}
    last_no_measure_warning_t::Float64
end

function DistributedEKFModel(
    sensor::OpticalLOSSensorModel,
    comms::InterAgentCommunicationModel,
    iod::DistributedIODCovModel,
    observer_idxs::Vector{Int},
    target_idx::Int,
    sigma_theta_rad::Float64,
    num_consensus_iter::Int,
    process_q_diag::SVector{6, Float64},
    μ::Float64,
    num_sats::Int
)
    return DistributedEKFModel(
        sensor=sensor,
        comms=comms,
        iod=iod,
        observer_idxs=observer_idxs,
        target_idx=target_idx,
        sigma_theta_rad=sigma_theta_rad,
        num_consensus_iter=num_consensus_iter,
        process_q_diag=process_q_diag,
        μ=μ,
        state=[SVector{6, Float64}(ntuple(_ -> NaN, 6)) for _ in 1:num_sats],
        covariance=[fill(NaN, 6, 6) for _ in 1:num_sats],
        state_pred=[SVector{6, Float64}(ntuple(_ -> NaN, 6)) for _ in 1:num_sats],
        covariance_pred=[fill(NaN, 6, 6) for _ in 1:num_sats],
        initialized=falses(num_sats),
        last_update_t=fill(NaN, num_sats),
        last_no_measure_warning_t=NaN
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

@inline function _max_degree_in_neighborhood(sid::Int, neighbor_map::Dict{Int, Vector{Int}})
    d_self = length(get(neighbor_map, sid, Int[]))
    d_max = d_self
    for n in get(neighbor_map, sid, Int[])
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
    consensus_values = Dict{Int, Any}(sid => deepcopy(local_values[sid]) for sid in observer_idxs)
    for _ in 1:num_iter
        new_values = Dict{Int, Any}()
        for sid in observer_idxs
            neighbors = get(neighbor_map, sid, Int[])
            degree_i = length(neighbors)
            max_deg = _max_degree_in_neighborhood(sid, neighbor_map)
            w_ii = 1.0 - degree_i / (max_deg + 1.0)
            sum_val = w_ii * consensus_values[sid]
            for nid in neighbors
                degree_j = length(get(neighbor_map, nid, Int[]))
                w_ij = 1.0 / (max(degree_i, degree_j) + 1.0)
                sum_val += w_ij * consensus_values[nid]
            end
            new_values[sid] = sum_val
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

function _initialize_dekf_from_iod_consensus!(
    model::DistributedEKFModel,
    neighbor_map::Dict{Int, Vector{Int}},
    t::Float64,
    u
)
    active = Int[]
    for sid in model.observer_idxs
        if model.iod.initialized[sid] &&
           _is_finite_state(model.iod.estimate_state[sid]) &&
           _is_finite_cov(model.iod.estimate_covariance[sid])
            push!(active, sid)
        end
    end
    isempty(active) && return nothing

    active_set = Set(active)
    active_neighbor_map = Dict{Int, Vector{Int}}(
        sid => [n for n in get(neighbor_map, sid, Int[]) if n in active_set] for sid in active
    )
    x_local = Dict{Int, Any}(sid => model.iod.estimate_state[sid] for sid in active)
    P_local = Dict{Int, Any}(sid => model.iod.estimate_covariance[sid] for sid in active)

    x_cons = _low_pass_consensus(x_local, active, active_neighbor_map, model.num_consensus_iter)
    P_cons = _low_pass_consensus(P_local, active, active_neighbor_map, model.num_consensus_iter)

    target_true = SVector{6, Float64}(
        SVector{3, Float64}(u.sc[model.target_idx].pos)...,
        SVector{3, Float64}(u.sc[model.target_idx].vel)...
    )
    for sid in active
        model.initialized[sid] && continue
        model.state[sid] = SVector{6, Float64}(x_cons[sid])
        model.covariance[sid] = _symmetrize_psd(Matrix(P_cons[sid]))
        model.state_pred[sid] = model.state[sid]
        model.covariance_pred[sid] = model.covariance[sid]
        model.initialized[sid] = true
        model.last_update_t[sid] = t

        err = model.state[sid] - target_true
        println("DEKF init (IOD+consensus) | tgt=$(model.target_idx) | obs=$(sid) | t=$(round(t; digits=3)) s")
        println("  state = ", model.state[sid])
        println("  cov(6x6) = ")
        show(stdout, "text/plain", model.covariance[sid])
        println("\n  |err_r|= $(round(norm(err[1:3]); digits=3)) m, |err_v|= $(round(norm(err[4:6]); digits=6)) m/s\n")
    end
    return nothing
end

# one different model for every different target satellite -> it'all the same but with dict to access the right target idx in the sensor measurements and in the IOD estimates. 
function SimulationModel.calcNavigationEffect!(
    model::DistributedEKFModel,
    u,
    p,
    t::Float64,
    sat_idx::Int
)
    # Run one distributed EKF pass per nav tick.
    sat_idx == model.observer_idxs[1] || return nothing
    _update_neighbors!(model.comms, u, p)

    observer_set = Set(model.observer_idxs) 
    neighbor_map = Dict{Int, Vector{Int}}() 
    for sid in model.observer_idxs
        neighbor_map[sid] = [n for n in model.comms.neighbors[sid] if n in observer_set]
    end

    _initialize_dekf_from_iod_consensus!(model, neighbor_map, t, u)

    I3 = Matrix(I, 3, 3)
    u_local = Dict{Int, Any}() # Maps from sat_idx to local measurement vector (3x1)
    uhat_local = Dict{Int, Any}() # Maps from sat_idx to local predicted measurement vector (3x1)
    U_local = Dict{Int, Any}() # Maps from sat_idx to local measurement information matrix (3x3)
    has_valid_local_measure = Dict{Int, Bool}()

    # Each observer forms its local measurement equations and information contributions.
    for sid in model.observer_idxs
        if !model.initialized[sid]
            u_local[sid] = zeros(6)
            uhat_local[sid] = zeros(6)
            U_local[sid] = zeros(6, 6)
            has_valid_local_measure[sid] = false
            continue
        end

        x_prior = model.state_pred[sid]
        r_agent = SVector{3, Float64}(u.sc[sid].pos)
        H = _measurement_jacobian(x_prior, r_agent)
        key = (sid, model.target_idx)
        meas = get(model.sensor.latest, key, nothing)
        if meas === nothing
            u_local[sid] = zeros(6)
            uhat_local[sid] = zeros(6)
            U_local[sid] = zeros(6, 6)
            has_valid_local_measure[sid] = false
            continue
        end
        if (t - meas.t) > (NAVIGATION_RATE_SEC + NAVIGATION_DT_TOL_SEC) || (meas.t - t) > NAVIGATION_DT_TOL_SEC
            u_local[sid] = zeros(6)
            uhat_local[sid] = zeros(6)
            U_local[sid] = zeros(6, 6)
            has_valid_local_measure[sid] = false
            continue
        end
        los_meas = meas.los_unit

        # Check for valid LOS measurement before proceeding with EKF update.
        if !all(isfinite, los_meas)
            u_local[sid] = zeros(6)
            uhat_local[sid] = zeros(6)
            U_local[sid] = zeros(6, 6)
            has_valid_local_measure[sid] = false
            continue
        end
        l_unit = _safe_unit(los_meas)
        if l_unit === nothing
            u_local[sid] = zeros(6)
            uhat_local[sid] = zeros(6)
            U_local[sid] = zeros(6, 6)
            has_valid_local_measure[sid] = false
            continue
        end

        # LOS covariance in Cartesian components:
        # R_i = σ_θ² (I - l lᵀ), rank-2 by construction.
        R_i = (model.sigma_theta_rad^2) * (I3 - l_unit * l_unit')
        R_inv = pinv(R_i)
        u_local[sid] = H' * (R_inv * los_meas)

        los_pred = _safe_unit(x_prior[1:3] - r_agent)
        if los_pred === nothing
            uhat_local[sid] = zeros(6)
            U_local[sid] = zeros(6, 6)
            continue
        end
        uhat_local[sid] = H' * (R_inv * los_pred)

        U_local[sid] = H' * (R_inv * H)
        has_valid_local_measure[sid] = true
    end

    # If no observer has a valid LOS for this target at this nav epoch,
    # skip correction and perform prediction-only propagation.
    num_valid_meas = sum(v ? 1 : 0 for v in values(has_valid_local_measure))
    num_agents = length(model.observer_idxs)
    Q = Diagonal(model.process_q_diag)
    if num_valid_meas == 0
        if !isfinite(model.last_no_measure_warning_t) || abs(t - model.last_no_measure_warning_t) > NAVIGATION_DT_TOL_SEC
            @warn "DEKF correction skipped: no valid LOS for target $(model.target_idx) at t=$(round(t; digits=3)) s."
            model.last_no_measure_warning_t = t
        end
        for sid in model.observer_idxs
            model.initialized[sid] || continue
            x_upd = model.state_pred[sid]
            P_upd = model.covariance_pred[sid]

            x_pred = _propagate_keplerian(
                SVector{6, Float64}(x_upd),
                model.μ,
                NAVIGATION_RATE_SEC
            )
            F = _compute_process_jacobian(
                SVector{6, Float64}(x_upd),
                model.μ,
                NAVIGATION_RATE_SEC
            )
            P_pred = F * P_upd * F' + num_agents * Matrix(Q)
            P_pred = _symmetrize_psd(P_pred)

            model.state[sid] = SVector{6, Float64}(x_upd)
            model.covariance[sid] = P_upd
            model.state_pred[sid] = x_pred
            model.covariance_pred[sid] = P_pred
            model.last_update_t[sid] = t
        end
        return nothing
    end
    model.last_no_measure_warning_t = NaN

    # Perform consensus to fuse information across observers.
    z_cons = _low_pass_consensus(u_local, model.observer_idxs, neighbor_map, model.num_consensus_iter)
    zhat_cons = _low_pass_consensus(uhat_local, model.observer_idxs, neighbor_map, model.num_consensus_iter)
    S_cons = _low_pass_consensus(U_local, model.observer_idxs, neighbor_map, model.num_consensus_iter)

    # Each observer performs an EKF update using the consensus information.
    for sid in model.observer_idxs
        model.initialized[sid] || continue
        P_prior = model.covariance_pred[sid]
        x_prior = model.state_pred[sid]
        S_i = S_cons[sid]
        P_upd = inv(inv(P_prior) + S_i)
        P_upd = _symmetrize_psd(P_upd)
        x_upd = x_prior + P_upd * (z_cons[sid] - zhat_cons[sid])

        x_pred = _propagate_keplerian(
            SVector{6, Float64}(x_upd),
            model.μ,
            NAVIGATION_RATE_SEC
        )
        F = _compute_process_jacobian(
            SVector{6, Float64}(x_upd),
            model.μ,
            NAVIGATION_RATE_SEC
        )
        P_pred = F * P_upd * F' + num_agents * Matrix(Q)
        P_pred = _symmetrize_psd(P_pred)

        model.state[sid] = SVector{6, Float64}(x_upd)
        model.covariance[sid] = P_upd
        model.state_pred[sid] = x_pred
        model.covariance_pred[sid] = P_pred
        model.last_update_t[sid] = t
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
            r_moon_j2000 = lock(SimulationModel.SPICE_LOCK) do
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

    r_body_j2000 = lock(SimulationModel.SPICE_LOCK) do
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
    p::ODEParams,
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
    p::ODEParams,
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
iod_models = [
    DistributedIODCovModel(
        optical_sensor,
        comms_model,
        target_idx,
        observer_idxs,
        1,
        SIGMA_THETA_RAD,
        length(spacecraft)
    ) for target_idx in target_idxs # one IOD model per target satellite
]
dekf_models = [
    DistributedEKFModel(
        optical_sensor,
        comms_model,
        iod_models[k], 
        observer_idxs,
        target_idxs[k],
        SIGMA_THETA_RAD,
        EKF_CONSENSUS_ITERS,
        EKF_PROCESS_Q_DIAG,
        planet.μ,
        length(spacecraft)
    ) for k in eachindex(target_idxs) # one DEKF model per target satellite
]
iod_by_target = Dict(target_idxs[k] => iod_models[k] for k in eachindex(target_idxs))   # to access the right IOD model for each target satellite
dekf_by_target = Dict(target_idxs[k] => dekf_models[k] for k in eachindex(target_idxs)) # to access the right DEKF model for each target satellite

navigation_effectors = Any[optical_sensor]
append!(navigation_effectors, iod_models)
append!(navigation_effectors, dekf_models)
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
        results_directory=joinpath(REPO_ROOT, "output"),
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
                    dekf = dekf_by_target[tgt]
                    return Dict("obs$(obs)" => dekf.state[obs] for obs in observer_idxs)
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
