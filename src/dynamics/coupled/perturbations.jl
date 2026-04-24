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
const _J2_COMPARE_LOG_COUNT = Base.Threads.Atomic{Int}(0)
const _MARS_MU_M3S2 = 0.4282837285418775e5 * 1e9
const _SPICE_FORCE_BARYCENTER_BODIES = (
    "mars",
    "jupiter",
    "saturn",
    "uranus",
    "neptune",
    "pluto",
)
const _THIRD_BODY_MU = Dict{String, Float64}(
    "sun" => 1.32712440041279419e20,
    "mercury" => 2.2032e13,
    "venus" => 3.248585926e14,
    "earth" => 3.986004418e14,
    "moon" => 4.9028005821478e12,
    "mars" => _MARS_MU_M3S2,
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
    if endswith(key, "_barycenter")
        return key
    end
    return key in _SPICE_FORCE_BARYCENTER_BODIES ? key * "_barycenter" : key
end
@inline function _resolve_third_body_mu(name::String)::Float64
    key = _mu_lookup_name(name)
    μ = get(_THIRD_BODY_MU, key, NaN)
    isfinite(μ) || throw(ArgumentError("Unsupported third body '$name' for NBodyGravityModel. Add its GM to _THIRD_BODY_MU in Perturbations.jl."))
    return μ
end

@inline function _spice_lock()
    mod = @__MODULE__
    while true
        if isdefined(mod, :RuntimeServices)
            return getproperty(mod, :RuntimeServices).SPICE_LOCK
        end
        parent = parentmodule(mod)
        parent === mod && break
        mod = parent
    end
    error("RuntimeServices.SPICE_LOCK not found in module ancestry for perturbations.jl")
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
    active_orders_by_degree::Vector{Vector{Int}} # Nonzero tesseral/sectoral orders per degree (m >= 1)
    reference_radius_m::Float64 # Reference radius associated with the harmonics coefficients
    planet::P # Planet data for primary body
end

@inline function _infer_harmonics_reference_radius_m(coefficients_file::String, planet::AbstractPlanet)::Float64
    file_name = lowercase(basename(coefficients_file))
    if file_name == "lp165p.csv"
        # LP165P is conventionally distributed with a 1738.0 km lunar reference radius.
        return 1.7380e6
    end
    if planet.name == "Mars"
        # Mars gravity field files used by this project are referenced to 3396.0 km.
        return 3.3960e6
    end
    return planet.Rp_e
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
    elseif pname == "moon"
        Moon("", spice_path)
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
    coefficient_normalization::Symbol=:full,
    coefficients_normalized::Bool=true,
    j2_source::Symbol=:file_c20,
    reference_radius_m::Union{Nothing, Float64}=nothing
) where P <: AbstractPlanet
    if L < 0 || M < 0
        throw(ArgumentError("Gravitational harmonics degree/order must be non-negative, got L=$L, M=$M."))
    end
    if M > L
        throw(ArgumentError("Gravitational harmonics order must satisfy M <= L, got L=$L, M=$M."))
    end
    normalized_source = _canonical_harmonics_normalization(coefficient_normalization)

    # Use CSV.Rows for streaming/lazy reading instead of materializing entire file.
    # This is critical for large files like Venus MGNP180U.csv (16K+ rows).
    # We process rows on-the-fly and skip those beyond our requested degree/order.

    # First: Quick structure validation using CSV.File to detect headers (minimal overhead for small inspection)
    test_parse = CSV.File(coefficients_file; delim=',', stripwhitespace=true, ignorerepeated=true, limit=1)
    has_degree_header = hasproperty(test_parse, :degree)

    # Now do the expensive iteration with streaming to avoid materializing full columns
    C_dict = Dict{Tuple{Int,Int}, Float64}()  # Temporary storage - only entries we need
    S_dict = Dict{Tuple{Int,Int}, Float64}()
    max_degree_file = 0

    # Parse the entire file row-by-row using CSV.Rows (lazy evaluation)
    open(coefficients_file) do io
        rows = CSV.Rows(io; delim=',', stripwhitespace=true, ignorerepeated=true)
        for (row_num, row) in enumerate(rows)
            # Skip header row if it exists
            if row_num == 1 && has_degree_header
                continue
            end

            try
                # Extract values by column position (works for both header and headerless files)
                l = Int(floor(parse(Float64, string(row[1])))) + 1  # degree (0-based in file)
                m = Int(floor(parse(Float64, string(row[2])))) + 1  # order (0-based in file)
                c = parse(Float64, string(row[3]))
                s = parse(Float64, string(row[4]))

                # Track max degree found
                max_degree_file = max(max_degree_file, l)

                # Only store coefficients within requested range
                if l <= L + 1 && m <= M + 1
                    C_dict[(l, m)] = c
                    S_dict[(l, m)] = s
                end
            catch
                # Skip rows with parsing errors
                continue
            end
        end
    end

    if max_degree_file == 0
        throw(ArgumentError("No valid harmonic coefficients found in $coefficients_file"))
    end

    if L + 1 > max_degree_file
        throw(ArgumentError(
            "Requested harmonics degree L=$L exceeds coefficients file support (max degree=$(max_degree_file - 1))."
        ))
    end
    if M + 1 > max_degree_file
        throw(ArgumentError(
            "Requested harmonics order M=$M exceeds coefficients file support (max order=$(max_degree_file - 1))."
        ))
    end

    # Cap the degree to only what we need to reduce memory allocation.
    # This is critical for files like Venus MGNP180U.csv which has degree 180 but we often only use degree 50.
    degree = L + 1

    C = zeros(Float64, degree, degree)
    S = zeros(Float64, degree, degree)

    # Convert unnormalized coefficients to fully normalized form when requested.
    # This follows the GMAT HarmonicGravity V(n,m) recurrence used during file load.
    norm_factor = nothing
    if !coefficients_normalized
        norm_factor = zeros(Float64, L + 1, M + 1)
        @inbounds for n = 0:L
            i = n + 1
            norm_factor[i, 1] = sqrt(2n + 1)
            if M >= 1
                vnm = sqrt(2 * (2n + 1))
                @inbounds for m = 1:min(M, n)
                    vnm /= sqrt((n + m) * (n - m + 1))
                    norm_factor[i, m + 1] = vnm
                end
            end
        end
    end

    # Transfer coefficients from dict to matrices, applying normalization if needed
    for ((l, m), c) in pairs(C_dict)
        s = S_dict[(l, m)]
        if !coefficients_normalized
            vnm = norm_factor[l, m]
            if vnm != 0.0
                c /= vnm
                s /= vnm
            else
                c = 0.0
                s = 0.0
            end
        end
        C[l, m] = c
        S[l, m] = s
    end

    if !(j2_source in (:file_c20, :planet_j2))
        throw(ArgumentError("Unsupported j2_source=$(j2_source). Expected :file_c20 or :planet_j2."))
    end
    if j2_source == :planet_j2 && L >= 2
        C[3, 1] = -planet.J2 / sqrt(5.0)
        S[3, 1] = 0.0
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

    C_trunc = C[1:(L + 1), 1:(M + 1)]
    S_trunc = S[1:(L + 1), 1:(M + 1)]

    active_orders_by_degree = [Int[] for _ in 1:(L + 1)]
    @inbounds for l = 1:L
        i = l + 1
        active = active_orders_by_degree[i]
        @inbounds for m = 1:min(M, l)
            j = m + 1
            if C_trunc[i, j] != 0.0 || S_trunc[i, j] != 0.0
                push!(active, m)
            end
        end
    end

    return GravitationalHarmonicsModel(
        L,
        M,
        C_trunc,
        S_trunc,
        A,
        R,
        I,
        VR01,
        VR11,
        N1,
        N2,
        sqrt_2n_plus_3,
        normalized_source,
        active_orders_by_degree,
        reference_radius_m === nothing ? _infer_harmonics_reference_radius_m(coefficients_file, planet) : reference_radius_m,
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
        pos_primary_k_all[k] = pos_primary_body_j2000_m
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
            pos_primary_k_mag = norm(pos_primary_k)
            @fastmath thread_force[worker_id] .+= mass * model.body_mus[k] * (
                (pos_spacecraft_k / pos_spacecraft_k_mag^3) - (pos_primary_k / (pos_primary_k_mag * pos_primary_k_mag * pos_primary_k_mag))
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
            pos_primary_k_mag = norm(pos_primary_k)
            @fastmath force_ii += mass * model.body_mus[k] * (
                (pos_spacecraft_k / pos_spacecraft_k_mag^3) - (pos_primary_k / (pos_primary_k_mag * pos_primary_k_mag * pos_primary_k_mag))
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

@inline environment_requirements(model::NBodyGravityModel) = EffectorEnvironmentRequirements(third_body_names=model.body_names)

@inline gravity_backbone_kick_structure(::NBodyGravityModel) = :velocity_kick_explicit

@inline function _nbody_acceleration_ii(
    model::NBodyGravityModel,
    x::StateSample,
    third_bodies::ThirdBodyEphemerisSample,
)::SVector{3, Float64}
    accel_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for k in eachindex(third_bodies.positions_ii)
        pos_primary_k = third_bodies.positions_ii[k]
        pos_spacecraft_k = pos_primary_k - x.pos_ii
        pos_spacecraft_k_mag = norm(pos_spacecraft_k)
        accel_ii .+= model.body_mus[k] * (
            (pos_spacecraft_k / pos_spacecraft_k_mag^3) - (pos_primary_k / norm(pos_primary_k)^3)
        )
    end
    return SVector{3, Float64}(accel_ii)
end

@inline function wrench(
    model::NBodyGravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    third_bodies = env.third_bodies
    third_bodies === nothing && throw(ArgumentError("NBodyGravityModel wrench requires env.third_bodies."))
    accel_ii = _nbody_acceleration_ii(model, x, third_bodies)
    return x.mass_kg * accel_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline function gravity_backbone_kick_acceleration_ii(
    model::NBodyGravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::SVector{3, Float64}
    third_bodies = env.third_bodies
    third_bodies === nothing && throw(ArgumentError("NBodyGravityModel gravity_backbone kick requires env.third_bodies."))
    return _nbody_acceleration_ii(model, x, third_bodies)
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
    pos_primary_sun = SVector{3, Float64}(pos_primary_sun_j2000_m)
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
        SVector{3, Float64}(pos_primary_sun_j2000_m)
    else
        SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    accel_ii = _srp_total_acceleration_ii(model, planet, pos_ii, pos_primary_sun, mass)
    force_ii = mass * SVector{3, Float64}(accel_ii)
    return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline environment_requirements(model::SolarRadiationPressureModel) = EffectorEnvironmentRequirements(solar=(model.direct || model.albedo))
@inline gravity_backbone_kick_structure(::SolarRadiationPressureModel) = :velocity_kick_explicit

@inline function _srp_total_acceleration_ii(
    model::SolarRadiationPressureModel,
    planet,
    pos_ii::SVector{3, Float64},
    pos_primary_sun::SVector{3, Float64},
    mass_kg::Float64,
)::SVector{3, Float64}
    if model.A == 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0)
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
            mass_kg;
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
            mass_kg,
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
            mass_kg,
            model.planet_ir_flux_w_m2,
        )
    end
    return SVector{3, Float64}(accel_ii)
end

@inline function wrench(
    model::SolarRadiationPressureModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    if model.A == 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    pos_primary_sun = if model.direct || model.albedo
        solar = env.solar
        solar === nothing && throw(ArgumentError("SolarRadiationPressureModel wrench requires env.solar."))
        solar.sun_pos_ii
    else
        SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    accel_ii = _srp_total_acceleration_ii(model, env.planet, x.pos_ii, pos_primary_sun, x.mass_kg)
    return x.mass_kg * accel_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline function gravity_backbone_kick_acceleration_ii(
    model::SolarRadiationPressureModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::SVector{3, Float64}
    pos_primary_sun = if model.direct || model.albedo
        solar = env.solar
        solar === nothing && throw(ArgumentError("SolarRadiationPressureModel gravity_backbone kick requires env.solar."))
        solar.sun_pos_ii
    else
        SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    return _srp_total_acceleration_ii(model, env.planet, x.pos_ii, pos_primary_sun, x.mass_kg)
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
    et = param.shared_buffers.et_start[] + param.shared_buffers.current_time[]

    # Transform from J2000 inertial to planet-fixed using the configured ephemerides model.
    L_PI = planet_frame_lpi(model.planet, et, param.args.environment_model.ephemerides_model)
    rVec_cart = L_PI * SVector{3, Float64}(x.pos) # convert from inertial to planet-fixed frame for gravity calculation
    mass = x.mass

    workspace = _harmonics_workspace_for_sat!(model, param, i)
    A = workspace.A
    R = workspace.R
    I = workspace.I

    RE = model.reference_radius_m
    r = norm(rVec_cart)
    inv_r = 1.0 / r
    s = rVec_cart[1] * inv_r
    t = rVec_cart[2] * inv_r
    u = rVec_cart[3] * inv_r
    L = model.L
    M = model.M

    # NOTE: A[1,1] is already initialised to 1.0 in _make_harmonics_scratch_workspace.
    # We only need to refresh the sub-diagonal entries that depend on the current
    # position (u = z/r changes every call).
    A[2, 1] = u * sqrt_3
    @inbounds @fastmath for n = 1:L+1
        row = n + 1
        A[row + 1, row] = u * model.sqrt_2n_plus_3[n] * A[row, row]
    end

    # Longitude trig recurrence: R[j] = Re[(s+it)^(j-1)], I[j] = Im[(s+it)^(j-1)]
    R[1] = 1.0
    I[1] = 0.0
    Rn = 1.0
    In = 0.0
    @inbounds @fastmath for j = 2:(M + 2)
        Rn, In = s * Rn - t * In, s * In + t * Rn
        R[j] = Rn
        I[j] = In
    end

    ρ = RE * inv_r
    a1 = a2 = a3 = a4 = 0.0

    max_recur_row = 2
    ρ_np1 = -model.planet.μ * inv_r * ρ
    @inbounds @fastmath for l = 1:L
        row = l + 1

        if row > max_recur_row
            jmax = min(max(M, 1) + 1, l - 1)
            @inbounds for j = 1:jmax
                A[row, j] = u * model.N1[row, j] * A[row - 1, j] - model.N2[row, j] * A[row - 2, j]
            end
            max_recur_row = row
        end

        next_row = row + 1
        if next_row > max_recur_row
            jmax_next = min(max(M, 1) + 1, l)
            @inbounds for j = 1:jmax_next
                A[next_row, j] = u * model.N1[next_row, j] * A[next_row - 1, j] - model.N2[next_row, j] * A[next_row - 2, j]
            end
            max_recur_row = next_row
        end

        ρ_np1 *= ρ
        rr = ρ_np1 / RE
        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
        sum4 = 0.0

        C0 = model.C[row, 1]
        # I[1] == 0.0 always; R[1] == 1.0 always — exploit this to save two multiplies.
        D0 = C0 * sqrt_2
        sum3 += model.VR01[row, 1] * A[row, 2] * D0
        sum4 += model.VR11[row, 1] * A[row + 1, 2] * D0

        active_orders = model.active_orders_by_degree[row]
        @inbounds @fastmath for idx in eachindex(active_orders)
            m = active_orders[idx]
            j = m + 1
            C = model.C[row, j]
            S = model.S[row, j]
            R_term = R[j - 1]
            I_term = I[j - 1]
            Rj  = R[j]
            Ij  = I[j]
            D = (C * Rj  + S * Ij)  * sqrt_2
            E = (C * R_term + S * I_term) * sqrt_2
            F = (S * R_term - C * I_term) * sqrt_2

            mA = m * A[row, j]
            sum1 += mA * E
            sum2 += mA * F
            sum3 += model.VR01[row, j] * A[row, j + 1] * D
            sum4 += model.VR11[row, j] * A[row + 1, j + 1] * D
        end

        a1 += rr * sum1
        a2 += rr * sum2
        a3 += rr * sum3
        a4 -= rr * sum4
    end

    g_pp_generic = SVector{3, Float64}(-a1 - s*a4, -a2 - t*a4, -a3 - u*a4)
    force_ii = mass * L_PI' * g_pp_generic

    compare_j2 = lowercase(strip(get(ENV, "SPACEAGORA_DEBUG_COMPARE_J2", "0"))) in ("1", "true", "yes", "on")
    if compare_j2 && model.L == 2 && model.M == 0
        C20 = model.C[3, 1]
        if isfinite(C20) && C20 != 0.0
            x_pp, y_pp, z_pp = rVec_cart
            r2 = r * r
            r4 = r2 * r2
            z2 = z_pp * z_pp
            J2 = -sqrt(5.0) * C20

            common = 1.5 * model.planet.μ * J2 * RE^2 / r4
            g_pp_analytic = SVector{3, Float64}(
                common * (x_pp / r) * (5.0 * z2 / r2 - 1.0),
                common * (y_pp / r) * (5.0 * z2 / r2 - 1.0),
                common * (z_pp / r) * (5.0 * z2 / r2 - 3.0)
            )

            Δg_pp = g_pp_generic - g_pp_analytic
            rel = norm(Δg_pp) / max(norm(g_pp_analytic), 1e-30)
            prev = Base.Threads.atomic_add!(_J2_COMPARE_LOG_COUNT, 1)
            if prev < 20
                @info "J2 generic-vs-analytic comparison" sat=i et=et r_m=r c20=C20 j2=J2 g_generic_pp=norm(g_pp_generic) g_analytic_pp=norm(g_pp_analytic) delta_pp=norm(Δg_pp) rel_pp=rel
            end
        end
    end

    return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline environment_requirements(::GravitationalHarmonicsModel) = EffectorEnvironmentRequirements(planet_frame=true)

@inline function wrench(
    model::GravitationalHarmonicsModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    planet_frame = env.planet_frame
    planet_frame === nothing && throw(ArgumentError("GravitationalHarmonicsModel wrench requires env.planet_frame."))

    workspace = _make_harmonics_scratch_workspace(model)
    A = workspace.A
    R = workspace.R
    I = workspace.I

    rVec_cart = planet_frame.pos_pp
    RE = model.planet.Rp_e
    r = norm(rVec_cart)
    s, n, u = normalize(rVec_cart)
    L = model.L
    M = model.M
    @fastmath begin
        A[2, 1] = u * sqrt_3
        @inbounds @simd for degree = 1:L+1
            idx = degree + 1
            A[idx + 1, idx] = u * model.sqrt_2n_plus_3[degree] * A[idx, idx]
        end
        @inbounds for order = 0:M+1
            j = order + 1
            @inbounds for degree = order+2:L+1
                idx = degree + 1
                A[idx, j] = u * model.N1[idx, j] * A[idx - 1, j] - model.N2[idx, j] * A[idx - 2, j]
            end
            if order == 0
                R[j] = 1.0
                I[j] = 0.0
            else
                R_term = R[j - 1]
                I_term = I[j - 1]
                R[j] = s * R_term - n * I_term
                I[j] = s * I_term + n * R_term
            end
        end

        ρ = RE / r
        ρ_np1 = -model.planet.μ / r * ρ
        a1 = a2 = a3 = a4 = 0.0
        @inbounds for degree = 1:L
            idx = degree + 1
            ρ_np1 *= ρ
            sum1 = 0.0
            sum2 = 0.0
            sum3 = 0.0
            sum4 = 0.0
            @inbounds for order = 0:min(degree, M)
                j = order + 1
                C = model.C[idx, j]
                S = model.S[idx, j]
                if order == 0
                    R_term = 0.0
                    I_term = 0.0
                else
                    R_term = R[j - 1]
                    I_term = I[j - 1]
                end
                D = (C * R[j] + S * I[j]) * sqrt_2
                E = ifelse(order == 0, 0.0, (C * R_term + S * I_term) * sqrt_2)
                F = ifelse(order == 0, 0.0, (S * R_term - C * I_term) * sqrt_2)

                sum1 += order * A[idx, j] * E
                sum2 += order * A[idx, j] * F
                sum3 += model.VR01[idx, j] * A[idx, j + 1] * D
                sum4 += model.VR11[idx, j] * A[idx + 1, j + 1] * D
            end
            rr = ρ_np1 / RE
            a1 += rr * sum1
            a2 += rr * sum2
            a3 += rr * sum3
            a4 -= rr * sum4
        end
        force_pp = x.mass_kg * SVector{3, Float64}(-a1 - s * a4, -a2 - n * a4, -a3 - u * a4)
        force_ii = planet_frame.l_pi' * force_pp
        return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
    end
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

    # Precompute dot product and reciprocals used multiple times below.
    dot_rs_rsat = dot(r_sun, r_sat)
    if dot_rs_rsat >= 0 # check the cos of the angle between the satellite and the Sun. If positive (angle less than 90 degrees), the satellite is not in eclipse
        return 1.0 # If the satellite is not in eclipse, return 1.0
    end
    inv_r_sun_norm = 1.0 / r_sun_norm
    inv_r_sat_norm = 1.0 / r_sat_norm
    # Eclipse conditions
    f1 = asin(_clamp_unit((rs + rp) * inv_r_sun_norm)) # Penumbra angle
    f2 = asin(_clamp_unit((rs - rp) * inv_r_sun_norm)) # Umbra angle
    s0 = -dot_rs_rsat * inv_r_sun_norm # Plane-axis intersection and planet center distance
    c1 = s0 + rp * sin(f1) # Distance from fundamental plane to cone vertex V1
    c2 = s0 - rp * sin(f2) # Distance from fundamental plane to cone vertex V2
    l1 = c1*tan(f1) # Radius of penumbra cone in fundamental plane
    l2 = c2*tan(f2) # Radius of umbra cone in fundamental plane
    l = √(max(r_sat_norm^2 - s0^2, 0.0)) # Distance from fundamental plane to satellite

    # Apparent radii of sun, planet, and apparent separation of sun and planet, respectively
    a = asin(_clamp_unit(rs * inv_r_sun_norm)) # Apparent radius of the Sun
    b = asin(_clamp_unit(rp * inv_r_sat_norm)) # Apparent radius of the planet
    c = acos(_clamp_unit(-dot_rs_rsat * inv_r_sun_norm * inv_r_sat_norm)) # Apparent separation of the Sun and planet
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
