using LoopVectorization
using AssociatedLegendrePolynomials
using LinearAlgebra
using SatelliteToolbox
using SatelliteToolboxGeomagneticField
using CSV
using DataFrames
using SpecialFunctions: loggamma
include(joinpath(@__DIR__, "..", "..", "core", "numerics", "quaternion_utils.jl"))
include(joinpath(@__DIR__, "..", "..", "environment", "ephemerides", "planet_data.jl"))
# import .config
# Define delta function
δ(x,y) = ==(x,y)
δ(x) = δ(x,0)

# Constants for the Tilted Dipole Model (Epoch 2020.0)
# Reference magnetic field strength at the equator on the Earth's surface.
const B0_2020 = 3.12e-5  # Tesla
# WGS84 Earth radius for the model.
const R_EARTH_MODEL = 6371.2e3 # meters
# North Magnetic Pole location (geocentric).
const POLE_LAT_2020 = deg2rad(80.7)  # Radians
const POLE_LON_2020 = deg2rad(-72.7) # Radians

# Pre-calculate the magnetic pole axis vector in ECEF for efficiency.
const M_HAT_ECEF = SVector{3, Float64}(
    cos(POLE_LAT_2020) * cos(POLE_LON_2020),
    cos(POLE_LAT_2020) * sin(POLE_LON_2020),
    sin(POLE_LAT_2020)
)
const sqrt_2 = sqrt(2.0)
const sqrt_3 = sqrt(3.0)
const _SPICE_BARYCENTER_BODIES = ("mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto")
const _THIRD_BODY_MU = Dict{String, Float64}(
    "sun" => 1.3271244002331e20,
    "mercury" => 2.2032e13,
    "venus" => 3.24858592e14,
    "earth" => 3.986004418e14,
    "moon" => 4.9028005821478e12,
    "mars" => 4.2828314e13,
    "jupiter" => 1.26686534e17,
    "saturn" => 3.7931187e16,
    "uranus" => 5.793939e15,
    "neptune" => 6.836529e15,
    "pluto" => 8.71e11,
    "titan" => 8.981e12
)

@inline _canonical_spice_name(name::String) = replace(lowercase(strip(name)), ' ' => '_')
@inline _mu_lookup_name(name::String) = replace(_canonical_spice_name(name), "_barycenter" => "")
@inline function _spice_query_name(name::String)
    key = _canonical_spice_name(name)
    return (endswith(key, "_barycenter") || !(key in _SPICE_BARYCENTER_BODIES)) ? key : key * "_barycenter"
end
@inline function _resolve_third_body_mu(name::String)::Float64
    key = _mu_lookup_name(name)
    μ = get(_THIRD_BODY_MU, key, NaN)
    isfinite(μ) || throw(ArgumentError("Unsupported third body '$name' for NBodyGravityModel. Add its GM to _THIRD_BODY_MU in Perturbations.jl."))
    return μ
end

# N-body gravity perturbation model
struct NBodyGravityModel{P <: AbstractPlanet, NP <: Tuple{Vararg{String}}, NM <: Tuple{Vararg{Float64}}} <: AbstractForceTorqueModel
    body_names::NP  # Names of celestial bodies to include
    body_mus::NM    # Gravitational parameters [m^3/s^2] for each third body
    primary_body_name::String # Name of the primary body
    planet::P # Planet data for primary body
end

@inline function _make_nbody_scratch_workspace(n_bodies::Int)::NBodyScratchWorkspace
    n_bodies >= 0 || throw(ArgumentError("NBody scratch workspace size must be >= 0, got $n_bodies"))
    n_workers = 1
    pos_primary_k_all = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n_bodies]
    thread_force = [MVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n_workers]
    return NBodyScratchWorkspace(pos_primary_k_all, thread_force)
end

@inline function _ensure_nbody_workspace_capacity!(
    workspace::NBodyScratchWorkspace,
    n_bodies::Int,
    n_workers::Int
)::NBodyScratchWorkspace
    if length(workspace.pos_primary_k_all) < n_bodies
        resize!(workspace.pos_primary_k_all, n_bodies)
    end
    if length(workspace.thread_force) < n_workers
        old_len = length(workspace.thread_force)
        resize!(workspace.thread_force, n_workers)
        @inbounds for worker_id in (old_len + 1):n_workers
            workspace.thread_force[worker_id] = MVector{3, Float64}(0.0, 0.0, 0.0)
        end
    end
    return workspace
end

@inline function _nbody_workspace_for_sat!(
    param::ODEParams,
    sat_idx::Int,
    n_bodies::Int,
    n_workers::Int
)::NBodyScratchWorkspace
    workspaces = param.shared_buffers.nbody_workspaces
    if sat_idx > length(workspaces)
        return _ensure_nbody_workspace_capacity!(_make_nbody_scratch_workspace(n_bodies), n_bodies, n_workers)
    end

    sat_entry = @inbounds workspaces[sat_idx]
    workspace = if sat_entry === nothing
        created = _make_nbody_scratch_workspace(n_bodies)
        @inbounds workspaces[sat_idx] = created
        created
    else
        sat_entry::NBodyScratchWorkspace
    end
    _ensure_nbody_workspace_capacity!(workspace, n_bodies, n_workers)
    return workspace
end

@inline function _canonical_harmonics_normalization(normalization)::Symbol
    norm = Symbol(normalization)
    if norm === :full || norm === :fully_normalized
        return :full
    elseif norm === :schmidt || norm === :schmidt_quasi_normalized || norm === :schmidt_quasinormalized
        return :schmidt
    elseif norm === :unnormalized || norm === :none
        return :unnormalized
    end
    throw(ArgumentError(
        "Unsupported gravitational harmonics normalization '$normalization'. " *
        "Supported values are :full (or :fully_normalized), :schmidt, and :unnormalized."
    ))
end

@inline function _fully_normalized_legendre_scale(l::Int, m::Int)::Float64
    l >= 0 || throw(ArgumentError("Degree must be >= 0, got $l."))
    0 <= m <= l || throw(ArgumentError("Order must satisfy 0 <= m <= l, got l=$l, m=$m."))
    δ0m = (m == 0) ? 1.0 : 0.0
    ratio = exp(0.5 * (loggamma(l - m + 1) - loggamma(l + m + 1)))
    return sqrt((2.0 - δ0m) * (2.0 * l + 1.0)) * ratio
end

@inline function _convert_harmonics_coefficients_to_full!(
    C::AbstractMatrix{Float64},
    S::AbstractMatrix{Float64},
    L::Int,
    M::Int,
    normalization::Symbol
)::Nothing
    norm = _canonical_harmonics_normalization(normalization)
    norm === :full && return nothing

    @inbounds for l in 0:L
        i = l + 1
        for m in 0:min(M, l)
            j = m + 1
            scale = if norm === :schmidt
                sqrt(2.0 * l + 1.0)
            else
                _fully_normalized_legendre_scale(l, m)
            end
            C[i, j] /= scale
            S[i, j] /= scale
        end
    end
    return nothing
end

struct GravitationalHarmonicsModel{P <: AbstractPlanet} <: AbstractForceTorqueModel
    L::Int64 # Degree
    M::Int64 # Order
    C::Matrix{Float64} # Cosine coefficients stored in fully normalized form
    S::Matrix{Float64} # Sine coefficients stored in fully normalized form
    A::Matrix{Float64} # Preallocated ALF array
    R::Vector{Float64} # Preallocated real terms powers
    I::Vector{Float64} # Preallocated imaginary terms vector
    VR01::Matrix{Float64} # Preallocated VR01 array
    VR11::Matrix{Float64} # Preallocated VR11 array
    N1::Matrix{Float64} # Preallocated N1 array
    N2::Matrix{Float64} # Preallocated N2 array
    sqrt_2n_plus_3::Vector{Float64} # Precalculated sqrt(2n+3) values
    coefficient_normalization::Symbol # Canonicalized source normalization keyword
    planet::P # Planet data for primary body
end

@inline function _make_harmonics_scratch_workspace(model::GravitationalHarmonicsModel)::HarmonicsScratchWorkspace
    L = model.L
    A = zeros(Float64, L + 4, L + 4)
    R = zeros(Float64, L + 4)
    I = zeros(Float64, L + 4)

    # Keep the same static diagonal initialization used by the legacy shared workspace.
    A[1, 1] = 1.0
    @inbounds for l = 1:(L + 2)
        i = l + 1
        A[i, i] = sqrt((2 * l + 1) / (2 * l)) * A[i - 1, i - 1]
    end
    return HarmonicsScratchWorkspace(A, R, I)
end

@inline function _harmonics_workspace_for_sat!(
    model::GravitationalHarmonicsModel,
    param::ODEParams,
    sat_idx::Int
)::HarmonicsScratchWorkspace
    workspaces = param.shared_buffers.harmonics_workspaces
    if sat_idx > length(workspaces)
        return _make_harmonics_scratch_workspace(model)
    end

    sat_entry = @inbounds workspaces[sat_idx]
    sat_map = if sat_entry === nothing
        created = Dict{UInt, HarmonicsScratchWorkspace}()
        @inbounds workspaces[sat_idx] = created
        created
    else
        sat_entry::Dict{UInt, HarmonicsScratchWorkspace}
    end

    key = objectid(model)
    workspace = get(sat_map, key, nothing)
    if workspace === nothing
        workspace = _make_harmonics_scratch_workspace(model)
        sat_map[key] = workspace
    end
    return workspace
end

struct SolarRadiationPressureModel <: AbstractForceTorqueModel
    Cr::Float64 # Reflectivity coefficient
    A::Float64  # Cross-sectional area in m^2
    AU_m::Float64 # Astronomical unit in meters
    direct::Bool
    albedo::Bool
    ir::Bool
    planet_albedo::Float64
    planet_ir_flux_w_m2::Float64
end

function SolarRadiationPressureModel(
    Cr::Real,
    A::Real,
    AU_m::Real=149_597_870_700.0;
    direct::Bool=true,
    albedo::Bool=false,
    ir::Bool=false,
    planet_albedo::Real=0.3,
    planet_ir_flux_w_m2::Real=237.0,
)
    Cr_f = Float64(Cr)
    A_f = Float64(A)
    AU_f = Float64(AU_m)
    planet_albedo_f = Float64(planet_albedo)
    planet_ir_flux_f = Float64(planet_ir_flux_w_m2)
    Cr_f >= 0.0 || throw(ArgumentError("SolarRadiationPressureModel.Cr must be >= 0, got $Cr_f."))
    A_f >= 0.0 || throw(ArgumentError("SolarRadiationPressureModel.A must be >= 0 m^2, got $A_f."))
    AU_f > 0.0 || throw(ArgumentError("SolarRadiationPressureModel.AU_m must be > 0 m, got $AU_f."))
    0.0 <= planet_albedo_f <= 1.0 || throw(ArgumentError("SolarRadiationPressureModel.planet_albedo must be in [0, 1], got $planet_albedo_f."))
    planet_ir_flux_f >= 0.0 || throw(ArgumentError("SolarRadiationPressureModel.planet_ir_flux_w_m2 must be >= 0, got $planet_ir_flux_f."))
    return SolarRadiationPressureModel(Cr_f, A_f, AU_f, direct, albedo, ir, planet_albedo_f, planet_ir_flux_f)
end

function NBodyGravityModel(;
    body_names::Tuple{Vararg{String}},
    primary_body_name::String="Earth",
    planet::P=Earth(""),
    body_mus::Union{Nothing, Tuple{Vararg{Float64}}}=nothing
) where {P <: AbstractPlanet}
    mus = body_mus === nothing ? Tuple(_resolve_third_body_mu(name) for name in body_names) : body_mus
    length(mus) == length(body_names) || throw(ArgumentError("NBodyGravityModel.body_mus length ($(length(mus))) must match body_names length ($(length(body_names)))."))
    return NBodyGravityModel{P, typeof(body_names), typeof(mus)}(body_names, mus, primary_body_name, planet)
end

# Constructor to get planet data
function NBodyGravityModel(body_names::Vector{String}, primary_body_name::String="Earth", spice_path::String="data/GRAMSuite.jl/GRAM Suite 2.0/SPICE")
    pname = lowercase(primary_body_name)
    planet = if pname == "earth"
        Earth("", spice_path)
    elseif pname == "mars"
        Mars("", spice_path)
    elseif pname == "venus"
        Venus("", spice_path)
    elseif pname == "titan"
        Titan("", spice_path)
    else
        throw(ArgumentError("Unsupported primary body '$primary_body_name' for NBodyGravityModel"))
    end
    return NBodyGravityModel(body_names=Tuple(body_names), primary_body_name=primary_body_name, planet=planet)
end

"""
    GravitationalHarmonicsModel(L, M, coefficients_file, planet; coefficient_normalization=:full)

Load degree/order-limited spherical-harmonics gravity coefficients from `coefficients_file`.
The CSV coefficients are converted at load time into the fully normalized convention expected
by the Pines recurrence used in `calcForceTorque`.

Accepted `coefficient_normalization` values are:
- `:full` or `:fully_normalized`
- `:schmidt`
- `:unnormalized`
"""
function GravitationalHarmonicsModel(
    L::Int64,
    M::Int64,
    coefficients_file::String,
    planet::P;
    coefficient_normalization::Symbol=:full
) where P <: AbstractPlanet
    if L < 0 || M < 0
        throw(ArgumentError("Gravitational harmonics degree/order must be non-negative, got L=$L, M=$M."))
    end
    if M > L
        throw(ArgumentError("Gravitational harmonics order must satisfy M <= L, got L=$L, M=$M."))
    end
    normalized_source = _canonical_harmonics_normalization(coefficient_normalization)

    harmonics_data = CSV.File(coefficients_file)
    total_data_size = size(harmonics_data, 1)
    degree = harmonics_data.degree[end] + 1
    if L + 1 > degree
        throw(ArgumentError(
            "Requested harmonics degree L=$L exceeds coefficients file support (max degree=$(degree - 1))."
        ))
    end
    if M + 1 > degree
        throw(ArgumentError(
            "Requested harmonics order M=$M exceeds coefficients file support (max order=$(degree - 1))."
        ))
    end

    C = zeros(Float64, degree, degree)
    S = zeros(Float64, degree, degree)
    for i=1:total_data_size
        l = harmonics_data.degree[i] + 1 # Get the degree, l, from the data and convert to an index (subtract 1 because the data starts at 2nd degree coefficient)
        m = harmonics_data.order[i] + 1 # Get the order, m, from the data and convert to an index (add 1 because the data starts at 0th order coefficient)
        C[l, m] = harmonics_data.C[i]
        S[l, m] = harmonics_data.S[i]
    end
    _convert_harmonics_coefficients_to_full!(C, S, L, M, normalized_source)

    N1 = zeros(Float64, L+4, L+4)
    N2 = zeros(Float64, L+4, L+4)
    VR01 = zeros(Float64, L+1, L+1)
    VR11 = zeros(Float64, L+1, L+1)
    @inbounds for m = 0:M+2
        j = m + 1
        @inbounds for l = m+2:L+2
            i = l + 1
            N1[i, j] = √((2*l+1)*(2*l-1)/(l+m)/(l-m))
            N2[i, j] = √((l+m-1)*(2*l+1)*(l-m-1)/(2*l-3)/(l+m)/(l-m))
        end
    end
    
    @inbounds for l = 0:L
        i = l + 1
        @inbounds for m = 0:min(M, l)
            j = m + 1
            divisor = (m == 0 ? sqrt_2 : 1)
            VR01[i, j] = sqrt((l-m)*(l+m+1)) / divisor
            VR11[i, j] = sqrt((2*l+1)*(l+m+2)*(l+m+1)/(2*l+3)) / divisor
        end
    end

    A = Matrix{Float64}(zeros(L+4, L+4))
    R = Vector{Float64}(zeros(L+4))
    I = Vector{Float64}(zeros(L+4))
    A[1, 1] = 1
    # Fill the diagonal elements of A
    @inbounds for l = 1:L+2
        i = l + 1
        A[i, i] = sqrt((2*l+1)/(2*l))*A[i-1, i-1]
    end

    # Precalculate the values of √(2n+3)
    sqrt_2n_plus_3 = Vector{Float64}(zeros(L+1))
    @inbounds for n = 1:L+1
        sqrt_2n_plus_3[n] = sqrt(2*n + 3)
    end

    return GravitationalHarmonicsModel(
        L,
        M,
        C[1:(L + 1), 1:(M + 1)],
        S[1:(L + 1), 1:(M + 1)],
        A,
        R,
        I,
        VR01,
        VR11,
        N1,
        N2,
        sqrt_2n_plus_3,
        normalized_source,
        planet
    )
end

"""
    calcForceTorque(model::NBodyGravityModel, x::AbstractVector{Float64}, ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
"""
function calcForceTorque(model::NBodyGravityModel, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_ii = hasproperty(x, :pos) ? SVector{3, Float64}(x.pos) : SVector{3, Float64}(x[1:3])
    mass = hasproperty(x, :mass) ? x.mass : x[7]
    primary_body_name = _spice_query_name(model.primary_body_name)
    et = param.shared_buffers.et_start[] + param.shared_buffers.current_time[]
    force_ii = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize force vector
    n_bodies = length(model.body_names)
    decision = _multibody_thread_decision(n_bodies; heavy_work=true)
    use_threads = decision.use_threads
    n_workers = use_threads ? ParallelPolicy.thread_worker_count(n_bodies, decision.allotment) : 1
    workspace = _nbody_workspace_for_sat!(param, i, n_bodies, n_workers)
    pos_primary_k_all = workspace.pos_primary_k_all
    nbody_cache_entry = param.shared_buffers.nbody_ephemeris_cache[]
    spice_rhs_memo_enabled = param.shared_buffers.spice_rhs_memo_enabled[]
    spice_rhs_memo = param.shared_buffers.spice_rhs_memo
    for (k, body_name) in pairs(model.body_names)
        body_name_spice = _spice_query_name(body_name)
        pos_primary_body_j2000_m = if nbody_cache_entry isa NBodyEphemerisCache
            cached = _nbody_body_position_from_cache_j2000_m(nbody_cache_entry, et, body_name_spice, primary_body_name)
            cached === nothing ? _nbody_body_position_from_spice_j2000_m(body_name_spice, et, primary_body_name, spice_rhs_memo_enabled, spice_rhs_memo, param.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls) : cached
        else
            _nbody_body_position_from_spice_j2000_m(body_name_spice, et, primary_body_name, spice_rhs_memo_enabled, spice_rhs_memo, param.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls)
        end
        pos_primary_k_all[k] = model.planet.J2000_to_pci * pos_primary_body_j2000_m
    end

    started_ns = time_ns()
    if use_threads
        thread_force = workspace.thread_force
        @inbounds for worker_id in 1:n_workers
            thread_force[worker_id] .= 0.0
        end
        ParallelPolicy.threaded_foreach_worker(n_bodies, decision.allotment) do worker_id, k
            pos_primary_k = pos_primary_k_all[k]
            pos_spacecraft_k = pos_primary_k - pos_ii
            pos_spacecraft_k_mag = norm(pos_spacecraft_k)
            thread_force[worker_id] .+= mass * model.body_mus[k] * (
                (pos_spacecraft_k / pos_spacecraft_k_mag^3) - (pos_primary_k / norm(pos_primary_k)^3)
            )
        end
        @inbounds for worker_id in 1:n_workers
            force_ii .+= thread_force[worker_id]
        end
    else
        @inbounds for k in eachindex(pos_primary_k_all)
            pos_primary_k = pos_primary_k_all[k]
            pos_spacecraft_k = pos_primary_k - pos_ii
            pos_spacecraft_k_mag = norm(pos_spacecraft_k)
            force_ii += mass * model.body_mus[k] * (
                (pos_spacecraft_k / pos_spacecraft_k_mag^3) - (pos_primary_k / norm(pos_primary_k)^3)
            )
        end
    end
    ParallelPolicy.record_policy_observation!(
        :multibody;
        mode=decision.mode,
        num_items=n_bodies,
        use_threads=use_threads,
        elapsed_ns=(time_ns() - started_ns)
    )

    return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline function _spice_rhs_memo_reset_if_stale!(
    memo::SpiceRhsMemo,
    et::Float64,
    primary_body_name::String
)
    if memo.et != et || memo.primary_body_name != primary_body_name
        memo.et = et
        memo.primary_body_name = primary_body_name
        empty!(memo.body_positions_j2000_m)
    end
    return nothing
end

@inline function _nbody_body_position_from_spice_j2000_m(
    body_name_spice::String,
    et::Float64,
    primary_body_name::String,
    memo_enabled::Bool,
    memo::SpiceRhsMemo,
    counter::Base.Threads.Atomic{Int64}
)::SVector{3, Float64}
    return memo_enabled ?
        _nbody_body_position_from_spice_memoized_j2000_m(body_name_spice, et, primary_body_name, memo, counter) :
        _nbody_body_position_from_spice_direct_j2000_m(body_name_spice, et, primary_body_name, counter)
end

@inline function _nbody_body_position_from_spice_direct_j2000_m(
    body_name_spice::String,
    et::Float64,
    primary_body_name::String,
    counter::Base.Threads.Atomic{Int64}
)::SVector{3, Float64}
    Base.Threads.atomic_add!(counter, 1)
    return spice_position_j2000_m(body_name_spice, et, primary_body_name)
end

@inline function _nbody_body_position_from_spice_memoized_j2000_m(
    body_name_spice::String,
    et::Float64,
    primary_body_name::String,
    memo::SpiceRhsMemo,
    counter::Base.Threads.Atomic{Int64}
)::SVector{3, Float64}
    return lock(memo.lock) do
        _spice_rhs_memo_reset_if_stale!(memo, et, primary_body_name)
        if haskey(memo.body_positions_j2000_m, body_name_spice)
            return memo.body_positions_j2000_m[body_name_spice]
        end
        Base.Threads.atomic_add!(counter, 1)
        position = spice_position_j2000_m(body_name_spice, et, primary_body_name)
        memo.body_positions_j2000_m[body_name_spice] = position
        return position
    end
end

@inline function _nbody_body_position_from_cache_j2000_m(
    cache::NBodyEphemerisCache,
    et::Float64,
    body_name_spice::String,
    primary_body_name::String
)::Union{Nothing, SVector{3, Float64}}
    cache.primary_body_name == primary_body_name || return nothing
    body_idx = get(cache.body_index_by_name, body_name_spice, 0)
    body_idx > 0 || return nothing

    ets = cache.ets
    n_samples = length(ets)
    n_samples >= 2 || return nothing
    if et < ets[1] || et > ets[n_samples]
        return nothing
    end

    idx = searchsortedlast(ets, et)
    if idx <= 0
        return nothing
    elseif idx >= n_samples
        return cache.positions_j2000_m[n_samples, body_idx]
    end

    et0 = ets[idx]
    et1 = ets[idx + 1]
    if et1 <= et0
        return cache.positions_j2000_m[idx, body_idx]
    end
    tau = (et - et0) / (et1 - et0)
    p1 = cache.positions_j2000_m[idx, body_idx]
    p2 = cache.positions_j2000_m[idx + 1, body_idx]
    if idx <= 1 || (idx + 2) > n_samples
        return _interp_vec3_linear(p1, p2, tau)
    end
    p0 = cache.positions_j2000_m[idx - 1, body_idx]
    p3 = cache.positions_j2000_m[idx + 2, body_idx]
    return _interp_vec3_catmull_rom(p0, p1, p2, p3, tau)
end

@inline function srp_cannonball_accel(
    pos_ii::SVector{3, Float64},
    pos_primary_sun::SVector{3, Float64},
    planet_radius_m::Float64,
    p_srp_unscaled::Float64,
    reflection_coefficient::Float64,
    reference_area_m2::Float64,
    mass_kg::Float64;
    AU_m::Float64=149_597_870_700.0
)::SVector{3, Float64}
    if !(isfinite(mass_kg) && mass_kg > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(reference_area_m2) && reference_area_m2 > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(reflection_coefficient) && reflection_coefficient >= 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    r_sun_to_sc = pos_ii - pos_primary_sun
    r_sun_to_sc_mag = norm(r_sun_to_sc)
    if !(isfinite(r_sun_to_sc_mag) && r_sun_to_sc_mag > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    eclipse_ratio = eclipse_area_calc(pos_ii, pos_primary_sun, planet_radius_m)
    P_srp = p_srp_unscaled * (AU_m / r_sun_to_sc_mag)^2
    return _cannonball_radiation_accel(
        r_sun_to_sc,
        P_srp * eclipse_ratio,
        reflection_coefficient,
        reference_area_m2,
        mass_kg,
    )
end

@inline function _cannonball_radiation_accel(
    direction_ii::SVector{3, Float64},
    pressure_n_m2::Float64,
    reflection_coefficient::Float64,
    reference_area_m2::Float64,
    mass_kg::Float64,
)::SVector{3, Float64}
    if !(isfinite(mass_kg) && mass_kg > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(reference_area_m2) && reference_area_m2 > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(reflection_coefficient) && reflection_coefficient >= 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(pressure_n_m2) && pressure_n_m2 > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    direction_mag = norm(direction_ii)
    if !(isfinite(direction_mag) && direction_mag > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    scale = reflection_coefficient * reference_area_m2 * pressure_n_m2 / mass_kg / direction_mag
    return scale * direction_ii
end

@inline function _lambert_phase_function(alpha_rad::Float64)::Float64
    if !isfinite(alpha_rad)
        return 0.0
    end
    α = clamp(alpha_rad, 0.0, π)
    return max(0.0, (sin(α) + (π - α) * cos(α)) / π)
end

@inline function planetary_albedo_accel(
    pos_ii::SVector{3, Float64},
    pos_primary_sun::SVector{3, Float64},
    planet_radius_m::Float64,
    p_srp_unscaled::Float64,
    reflection_coefficient::Float64,
    reference_area_m2::Float64,
    mass_kg::Float64,
    planet_albedo::Float64;
    AU_m::Float64=149_597_870_700.0
)::SVector{3, Float64}
    if !(isfinite(planet_radius_m) && planet_radius_m > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(planet_albedo) && planet_albedo > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    spacecraft_range_m = norm(pos_ii)
    sun_range_m = norm(pos_primary_sun)
    if !(isfinite(spacecraft_range_m) && spacecraft_range_m > planet_radius_m)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(sun_range_m) && sun_range_m > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    cos_alpha = clamp(dot(pos_primary_sun, pos_ii) / (sun_range_m * spacecraft_range_m), -1.0, 1.0)
    phase = _lambert_phase_function(acos(cos_alpha))
    if phase <= 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    P_solar_at_planet = p_srp_unscaled * (AU_m / sun_range_m)^2
    pressure = P_solar_at_planet * planet_albedo * phase * (planet_radius_m / spacecraft_range_m)^2
    return _cannonball_radiation_accel(
        pos_ii,
        pressure,
        reflection_coefficient,
        reference_area_m2,
        mass_kg,
    )
end

@inline function planetary_ir_accel(
    pos_ii::SVector{3, Float64},
    planet_radius_m::Float64,
    reflection_coefficient::Float64,
    reference_area_m2::Float64,
    mass_kg::Float64,
    planet_ir_flux_w_m2::Float64,
)::SVector{3, Float64}
    if !(isfinite(planet_radius_m) && planet_radius_m > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !(isfinite(planet_ir_flux_w_m2) && planet_ir_flux_w_m2 > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    spacecraft_range_m = norm(pos_ii)
    if !(isfinite(spacecraft_range_m) && spacecraft_range_m > planet_radius_m)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    pressure = (planet_ir_flux_w_m2 / 299_792_458.0) * (planet_radius_m / spacecraft_range_m)^2
    return _cannonball_radiation_accel(
        pos_ii,
        pressure,
        reflection_coefficient,
        reference_area_m2,
        mass_kg,
    )
end

function srp(
    planet,
    p_srp_unscaled::Float64,
    reflection_coefficient::Float64,
    reference_area_m2::Float64,
    mass_kg::Float64,
    pos_ii::AbstractVector,
    et::Float64;
    AU_m::Float64=149_597_870_700.0
)
    pos_ii_sv = SVector{3, Float64}(pos_ii)
    primary_body_name = _spice_query_name(planet.name)
    pos_primary_sun_j2000_m = spice_position_j2000_m("sun", et, primary_body_name)
    pos_primary_sun = SVector{3, Float64}(planet.J2000_to_pci * pos_primary_sun_j2000_m)
    return srp_cannonball_accel(
        pos_ii_sv,
        pos_primary_sun,
        planet.Rp_e,
        p_srp_unscaled,
        reflection_coefficient,
        reference_area_m2,
        mass_kg;
        AU_m=AU_m
    )
end

function calcForceTorque(model::SolarRadiationPressureModel, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    if model.A == 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    planet = param.args.environment_model.planet
    pos_ii = hasproperty(x, :pos) ? SVector{3, Float64}(x.pos) : SVector{3, Float64}(x[1:3])
    mass = hasproperty(x, :mass) ? x.mass : x[7]
    et = param.shared_buffers.et_start[] + param.shared_buffers.current_time[]
    primary_body_name = _spice_query_name(planet.name)
    spice_rhs_memo_enabled = param.shared_buffers.spice_rhs_memo_enabled[]
    spice_rhs_memo = param.shared_buffers.spice_rhs_memo

    pos_primary_sun = if model.direct || model.albedo
        cache_entry = param.shared_buffers.srp_sun_ephemeris_cache[]
        pos_primary_sun_j2000_m = if cache_entry isa SRPSunEphemerisCache
            cached = _srp_sun_position_from_cache_j2000_m(cache_entry, et)
            cached === nothing ? _srp_sun_position_from_spice_j2000_m(et, primary_body_name, spice_rhs_memo_enabled, spice_rhs_memo, param.shared_buffers.spice_runtime_counters.srp_spkpos_runtime_calls) : cached
        else
            _srp_sun_position_from_spice_j2000_m(et, primary_body_name, spice_rhs_memo_enabled, spice_rhs_memo, param.shared_buffers.spice_runtime_counters.srp_spkpos_runtime_calls)
        end
        SVector{3, Float64}(planet.J2000_to_pci * pos_primary_sun_j2000_m)
    else
        SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    accel_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    if model.direct
        accel_ii .+= srp_cannonball_accel(
            pos_ii,
            pos_primary_sun,
            planet.Rp_e,
            4.56e-6,
            model.Cr,
            model.A,
            mass;
            AU_m=model.AU_m
        )
    end
    if model.albedo
        accel_ii .+= planetary_albedo_accel(
            pos_ii,
            pos_primary_sun,
            planet.Rp_e,
            4.56e-6,
            model.Cr,
            model.A,
            mass,
            model.planet_albedo;
            AU_m=model.AU_m
        )
    end
    if model.ir
        accel_ii .+= planetary_ir_accel(
            pos_ii,
            planet.Rp_e,
            model.Cr,
            model.A,
            mass,
            model.planet_ir_flux_w_m2,
        )
    end

    force_ii = mass * SVector{3, Float64}(accel_ii)
    return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline function _srp_sun_position_from_spice_j2000_m(
    et::Float64,
    primary_body_name::String,
    memo_enabled::Bool,
    memo::SpiceRhsMemo,
    counter::Base.Threads.Atomic{Int64}
)::SVector{3, Float64}
    return memo_enabled ?
        _srp_sun_position_from_spice_memoized_j2000_m(et, primary_body_name, memo, counter) :
        _srp_sun_position_from_spice_direct_j2000_m(et, primary_body_name, counter)
end

@inline function _srp_sun_position_from_spice_direct_j2000_m(
    et::Float64,
    primary_body_name::String,
    counter::Base.Threads.Atomic{Int64}
)::SVector{3, Float64}
    Base.Threads.atomic_add!(counter, 1)
    return spice_position_j2000_m("sun", et, primary_body_name)
end

@inline function _srp_sun_position_from_spice_memoized_j2000_m(
    et::Float64,
    primary_body_name::String,
    memo::SpiceRhsMemo,
    counter::Base.Threads.Atomic{Int64}
)::SVector{3, Float64}
    sun_query_name = "sun"
    return lock(memo.lock) do
        _spice_rhs_memo_reset_if_stale!(memo, et, primary_body_name)
        if haskey(memo.body_positions_j2000_m, sun_query_name)
            return memo.body_positions_j2000_m[sun_query_name]
        end
        Base.Threads.atomic_add!(counter, 1)
        position = spice_position_j2000_m(sun_query_name, et, primary_body_name)
        memo.body_positions_j2000_m[sun_query_name] = position
        return position
    end
end

@inline function _srp_sun_position_from_cache_j2000_m(cache::SRPSunEphemerisCache, et::Float64)::Union{Nothing, SVector{3, Float64}}
    ets = cache.ets
    n_samples = length(ets)
    n_samples >= 2 || return nothing
    if et < ets[1] || et > ets[n_samples]
        return nothing
    end

    idx = searchsortedlast(ets, et)
    if idx <= 0
        return nothing
    elseif idx >= n_samples
        return cache.positions_j2000_m[n_samples]
    end

    et0 = ets[idx]
    et1 = ets[idx + 1]
    if et1 <= et0
        return cache.positions_j2000_m[idx]
    end
    tau = (et - et0) / (et1 - et0)
    p1 = cache.positions_j2000_m[idx]
    p2 = cache.positions_j2000_m[idx + 1]
    if idx <= 1 || (idx + 2) > n_samples
        return _interp_vec3_linear(p1, p2, tau)
    end
    p0 = cache.positions_j2000_m[idx - 1]
    p3 = cache.positions_j2000_m[idx + 2]
    return _interp_vec3_catmull_rom(p0, p1, p2, p3, tau)
end

@inline function _interp_vec3_linear(
    p0::SVector{3, Float64},
    p1::SVector{3, Float64},
    tau::Float64
)::SVector{3, Float64}
    return p0 + tau * (p1 - p0)
end

@inline function _interp_vec3_catmull_rom(
    p0::SVector{3, Float64},
    p1::SVector{3, Float64},
    p2::SVector{3, Float64},
    p3::SVector{3, Float64},
    tau::Float64
)::SVector{3, Float64}
    tau2 = tau * tau
    tau3 = tau2 * tau
    return 0.5 * (
        (2.0 * p1) +
        (-p0 + p2) * tau +
        (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * tau2 +
        (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * tau3
    )
end

"""
    calcForceTorque(model::GravitationalHarmonicsModel, x, param, i)

Compute the non-central spherical-harmonics gravity perturbation using coefficients stored in
the model's internal fully normalized convention. The inverse-square central term is handled
separately by the base gravity effector and is intentionally not included here.
"""
function calcForceTorque(model::GravitationalHarmonicsModel, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    # cnf = param.cnf
    
    rVec_cart = SVector{3, Float64}(x.pos)
    mass = x.mass               # Mass of the spacecraft, change to x.m if using StructArrays in Complete_passage

    workspace = _harmonics_workspace_for_sat!(model, param, i)
    A = workspace.A
    R = workspace.R
    I = workspace.I

    RE = param.args.environment_model.planet.Rp_e
    r = norm(rVec_cart)
    s,t,u=normalize(rVec_cart)
    L = model.L
    M = model.M
    @fastmath begin
        A[2, 1] = u * sqrt_3
        # Fill the off diagonal elements of A
        @inbounds @simd for n = 1:L+1
            i = n + 1
            A[i + 1, i] = u * model.sqrt_2n_plus_3[n] * A[i, i]
        end
        # Fill the rest of A
        @inbounds for m = 0:M+1
            j = m + 1
            @inbounds for l = m+2:L+1
                i = l + 1
                A[i, j] = u * model.N1[i, j] * A[i - 1, j] - model.N2[i, j] * A[i - 2, j]
            end
            if m == 0
                R[j] = 1.0
                I[j] = 0.0
            else
                R_term = R[j - 1]
                I_term = I[j - 1]
                R[j] = s * R_term - t * I_term
                I[j] = s * I_term + t * R_term
            end
        end

        ρ = RE/r
        ρ_np1 = -model.planet.μ/r * ρ
        a1 = a2 = a3 = a4 = 0.0
        @inbounds for l = 1:L
            i = l + 1
            ρ_np1 *= ρ
            sum1 = 0.0
            sum2 = 0.0
            sum3 = 0.0
            sum4 = 0.0
            @inbounds for m = 0:min(l, M)
                j = m + 1
                C = model.C[i, j]
                S = model.S[i, j]
                if m == 0
                    R_term = 0.0
                    I_term = 0.0
                else
                    R_term = R[j - 1]
                    I_term = I[j - 1]
                end
                D = (model.C[i, j] * R[j] + model.S[i, j] * I[j]) * sqrt_2
                E = ifelse(m == 0, 0.0, (C * R_term + S * I_term) * sqrt_2)
                F = ifelse(m == 0, 0.0, (S * R_term - C * I_term) * sqrt_2)

                # Avv00, Avv01, Avv11 = A[i, j], model.VR01[i, j] * A[i, j+1], model.VR11[i, j] * A[i+1, j+1]

                sum1 += m * A[i, j] * E
                sum2 += m * A[i, j] * F
                sum3 += model.VR01[i, j] * A[i, j + 1] * D
                sum4 += model.VR11[i, j] * A[i + 1, j + 1] * D
            end
            rr = ρ_np1/RE
            a1 += rr * sum1
            a2 += rr * sum2
            a3 += rr * sum3
            a4 -= rr * sum4
        end
    end

    force_ii = mass * SVector{3, Float64}(-a1 - s*a4, -a2 - t*a4, -a3 - u*a4) # Store gravity in config for other uses

    # cnf.gravity_harmonics_ii .= force_ii

    return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
end

"""
    get_magnetic_field_dipole(r_ecef::AbstractVector)

Calculates the Earth's magnetic field using a fast tilted dipole approximation.

This model is significantly faster than the full WMM and is suitable for use
inside performance-critical code like numerical integrators.

# Args

- `r_ecef`: The position vector of the spacecraft in an Earth-Centered,
            Earth-Fixed (ECEF) frame [meters].

# Returns

- A 3-element `SVector` representing the magnetic field `[Bx, By, Bz]` in the
  ECEF frame [Tesla].
"""
function get_magnetic_field_dipole(r_ecef::AbstractVector, L_PI::MMatrix{3, 3, Float64})
    r_norm = norm(r_ecef)
    r_hat = r_ecef / r_norm

    # Cosine of the magnetic colatitude
    cos_colat = dot(r_hat, M_HAT_ECEF)

    # Dipole field equation. This is a standard formulation.
    B_ecef = -B0_2020 * (R_EARTH_MODEL / r_norm)^3 * (M_HAT_ECEF - 3 * cos_colat * r_hat)

    return L_PI' * B_ecef
end

"""
    get_magnetic_field(date::DateTime, lat_deg::Number, lon_deg::Number, alt_m::Number)

Computes the Earth's magnetic field vector in the local North-East-Down (NED)
frame using the World Magnetic Model (WMM).

The function automatically uses the correct WMM version based on the input `date`.

# Args

- `date`: The `DateTime` of the measurement.
- `lat_deg`: The geodetic latitude of the observer [degrees].
- `lon_deg`: The longitude of the observer [degrees].
- `alt_m`: The altitude above the WGS84 ellipsoid [meters].

# Returns

- A 3-element `SVector` representing the magnetic field `[B_north, B_east, B_down]`
  in nanoTeslas [nT].
"""
function get_magnetic_field(date::DateTime, lat_rad::Number, lon_rad::Number, alt_m::Number, L_PI::MMatrix{3, 3, Float64})
    # println("Calculating magnetic field at lat: $lat_rad, lon: $lon_rad, alt: $alt_m, date: $date")
    # Calculate the magnetic field vector using the World Magnetic Model.
    # The result is in the NED frame and has units of nT.
    B_ned = igrf(yeardecimal(date), alt_m, lat_rad, lon_rad, Val(:geodetic))
    B_pp = ned_to_ecef(B_ned, lat_rad, lon_rad, alt_m)
    B_ii = L_PI' * B_pp
    return B_ii
end

"""
    calculate_magnetic_torque(m::AbstractVector, B::AbstractVector)

Calculates the magnetic torque exerted on a magnetic dipole by an external
magnetic field.

**τ** = **m** × **B**

The magnetic dipole moment `m` and the magnetic field `B` vectors **must** be
expressed in the same reference frame. The resulting torque vector `τ` will be
in that same frame.

# Args

- `m`: The magnetic dipole moment vector of the spacecraft [A·m²].
- `B`: The external magnetic field vector [Tesla].

# Returns

- A 3-element `SVector` representing the magnetic torque `[τ_x, τ_y, τ_z]` [N·m].
"""
function calculate_magnetic_torque(m::AbstractVector, B::AbstractVector)
    # Ensure inputs are StaticArrays for performance
    m_svector = SVector{3, Float64}(m)
    B_svector = SVector{3, Float64}(B)

    # Calculate the torque using the cross product
    # τ = m × B
    τ = cross(m_svector, B_svector)

    return τ
end

function eclipse_area_calc(r_sat::SVector{3, Float64}, r_sun::SVector{3, Float64}, rp::Float64)
    """
    Calculate the exposed area of the satellite. Translated from Python to Julia. 
    See equations and diagrams in Basilisk documentation: https://avslab.github.io/basilisk/Documentation/simulation/environment/eclipse/eclipse.html

    Parameters 
    ----------
    r_sat : Vector{Float64}
        Position vector of the satellite relative to the planet.
    r_sun : Vector{Float64}
        Position vector of the Sun relative to the planet.
    A : Float64
        Area of the satellite.
    rp : Float64
        Radius of the planet.


    Returns
    -------
    eclipse_ratio : Float64
        Eclipse ratio of the satellite.
    
    """
    rs = 695000e3 # Radius of the Sun in meters
    @inline _clamp_unit(x::Float64) = clamp(x, -1.0, 1.0)

    r_sun_norm = norm(r_sun)
    r_sat_norm = norm(r_sat)
    if !isfinite(r_sun_norm) || !isfinite(r_sat_norm) || r_sun_norm <= 0.0 || r_sat_norm <= 0.0
        return 1.0
    end

    if dot(r_sun, r_sat) >= 0 # check the cos of the angle between the satellite and the Sun. If positive (angle less than 90 degrees), the satellite is not in eclipse
        return 1.0 # If the satellite is not in eclipse, return 1.0
    end
    # Eclipse conditions
    f1 = asin(_clamp_unit((rs + rp) / r_sun_norm)) # Penumbra angle
    f2 = asin(_clamp_unit((rs - rp) / r_sun_norm)) # Umbra angle
    s0 = -dot(r_sat, r_sun) / r_sun_norm # Plane-axis intersection and planet center distance
    c1 = s0 + rp * sin(f1) # Distance from fundamental plane to cone vertex V1
    c2 = s0 - rp * sin(f2) # Distance from fundamental plane to cone vertex V2
    l1 = c1*tan(f1) # Radius of penumbra cone in fundamental plane
    l2 = c2*tan(f2) # Radius of umbra cone in fundamental plane
    l = √(max(r_sat_norm^2 - s0^2, 0.0)) # Distance from fundamental plane to satellite
    
    # Apparent radii of sun, planet, and apparent separation of sun and planet, respectively
    a = asin(_clamp_unit(rs / r_sun_norm)) # Apparent radius of the Sun
    b = asin(_clamp_unit(rp / r_sat_norm)) # Apparent radius of the planet
    c = acos(_clamp_unit(-dot(r_sun, r_sat) / (r_sun_norm * r_sat_norm))) # Apparent separation of the Sun and planet
    if c < b - a # Total eclipse condition
        return 0.0 # If the satellite is in total eclipse, return 0.0
    elseif c < a - b # Annular eclipse condition
        return 1.0 - b^2 / a^2 # Exposed sun fraction
    elseif c < a + b # Partial eclipse condition
        x = (c^2 + a^2 - b^2) / (2 * c)
        y = √(max(a^2 - x^2, 0.0))
        A = a^2 * acos(_clamp_unit(x / a)) + b^2 * acos(_clamp_unit((c - x) / b)) - c * y
        return 1 - A / (π * a^2) # Shadow fraction
    else # No eclipse condition
        return 1.0 # If the satellite is not in eclipse, return 1.0
    end
end

function srp!(model, root_index::Int64, sun_dir_ii::SVector{3, Float64}, body, P_srp::Float64, eclipse_ratio::Float64, orientation::Bool)
    """
    Calculate force on a body due to solar radiation pressure.

    Parameters
    ----------
    pos_ii : SVector{3, Float64}
        Position of the body in the inertial frame (J2000)
    sun_dir_ii : SVector{3, Float64}
        Unit vector in the direction of the Sun expressed in the inertial frame
    body : Body struct
        Struct containing physical information about the body
    r_sun_norm : Float64
        Magnitude of the spacecraft distance to the Sun
    P_srp : Float64
        Magnitude of the solar radiation pressure force at r_sun_norm meters from the Sun
    
    Returns
    -------
    F_srp : SVector{3, Float64}
        Force on the body in the inertial frame
    """
    rot_inertial = config.rotate_to_inertial(model, body, root_index)
    rot_body_to_inertial = rot(model.links[root_index].q)
    @inbounds for facet in body.SRP_facets
        rot_RF = rot_inertial * rot(facet.attitude)' # Rotation matrix from facet frame to inertial frame
        n = normalize(rot_RF * facet.normal_vector) # Normal vector of the facet in the inertial frame
        cos_α_srp = dot(n, sun_dir_ii) / norm(n) / norm(sun_dir_ii)

        if cos_α_srp > 0 && eclipse_ratio != 0.0 # If the facet is illuminated by the Sun
            F_SRP = -P_srp * facet.area * cos_α_srp * ((1 - facet.δ) * sun_dir_ii + 2 * (facet.ρ / 3 + facet.δ * cos_α_srp) * n) * eclipse_ratio
            body.net_force += F_SRP # Rotate F_SRP from body frame to inertial frame

            if orientation
                R_facet = rot_inertial*facet.cp + rot_body_to_inertial'*body.r # Vector from CoM of spacecraft to facet Cp in inertial frame
                # R_facet_body = config.rotate_to_body(body)*facet.cp + body.r
                body.net_torque += rot_body_to_inertial * cross(R_facet, F_SRP) # Calculate body frame net torque
            end
        end
    end
    # CSV.write("facet_forces.csv", df)
    # return F_SRP_tracker
    # println("Total F_SRP: $F_SRP_tracker")
end
