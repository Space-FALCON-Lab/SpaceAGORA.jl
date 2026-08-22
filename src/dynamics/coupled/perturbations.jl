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
    reference_radius_m::Float64 # Reference radius the harmonics coefficients were solved/normalized against
    gm_m3s2::Float64 # Gravitational parameter associated with this specific gravity solution
    include_central::Bool # Include the inverse-square primary-body term in this gravity model
    planet::P # Planet data for primary body
end

# Gravity coefficient files may carry a `#`-prefixed metadata header (reference
# radius, GM, source citation) documenting the specific gravity solution the
# coefficients were fit against -- this is *not* the same thing as the generic
# per-planet Rp_e/μ used elsewhere (atmosphere, orbit propagation, etc.), since
# different gravity solutions for the same body are commonly normalized against
# slightly different reference radii/GM (e.g. Mars-50c vs GMM-2B differ by
# ~0.08% in reference radius). `CSV.File`/`CSV.Rows` are called with
# `comment="#"` elsewhere in this file so these lines are transparently
# skipped by the coefficient parser; this function re-reads them separately to
# recover the values.
@inline function _read_harmonics_source_metadata(coefficients_file::String)::Dict{String, Float64}
    metadata = Dict{String, Float64}()
    isfile(coefficients_file) || return metadata
    open(coefficients_file) do io
        for line in eachline(io)
            startswith(line, "#") || break
            m = match(r"^#\s*(reference_radius_m|gm_m3s2)\s*[:=]\s*([-+0-9.eE]+)\s*$", line)
            m === nothing && continue
            metadata[m.captures[1]] = parse(Float64, m.captures[2])
        end
    end
    return metadata
end

@inline function _infer_harmonics_reference_radius_m(coefficients_file::String, planet::AbstractPlanet)::Float64
    metadata = _read_harmonics_source_metadata(coefficients_file)
    haskey(metadata, "reference_radius_m") && return metadata["reference_radius_m"]
    # No sourced metadata available for this file: fall back to the generic
    # per-planet radius. This is a poor substitute for the file's own
    # normalization radius (see the metadata header format above) and should
    # be treated as a last resort, not a validated value.
    return planet.Rp_e
end

@inline function _infer_harmonics_gm_m3s2(coefficients_file::String, planet::AbstractPlanet)::Float64
    metadata = _read_harmonics_source_metadata(coefficients_file)
    haskey(metadata, "gm_m3s2") && return metadata["gm_m3s2"]
    # No sourced metadata available for this file: fall back to the generic
    # per-planet GM. See _infer_harmonics_reference_radius_m.
    return planet.μ
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

# ---------------------------------------------------------------------------
# Batched harmonics workspace — satellite index is the first (fastest) dimension
# so that @turbo loops over the batch are unit-stride in memory.
# ---------------------------------------------------------------------------

struct HarmonicsBatchWorkspace
    A::Array{Float64, 3}    # [B, L+4, L+4]
    R::Matrix{Float64}       # [B, M+3]
    I::Matrix{Float64}       # [B, M+3]
    s_vec::Vector{Float64}   # x/r per satellite
    t_vec::Vector{Float64}   # y/r per satellite
    u_vec::Vector{Float64}   # z/r per satellite
    inv_r::Vector{Float64}
    mass::Vector{Float64}
    rVec::Matrix{Float64}    # [3, B] planet-frame positions
    ρ::Vector{Float64}
    ρ_np1::Vector{Float64}
    rr::Vector{Float64}
    a1::Vector{Float64}
    a2::Vector{Float64}
    a3::Vector{Float64}
    a4::Vector{Float64}
    sum1::Vector{Float64}
    sum2::Vector{Float64}
    sum3::Vector{Float64}
    sum4::Vector{Float64}
end

function _make_harmonics_batch_workspace(model::GravitationalHarmonicsModel, batch_size::Int)
    L = model.L
    M = model.M
    B = batch_size
    A = zeros(Float64, B, L + 4, L + 4)
    # Diagonal elements are position-independent; initialise once and never overwrite.
    A[:, 1, 1] .= 1.0
    @inbounds for l = 1:(L + 2)
        i = l + 1
        diag_val = sqrt((2 * l + 1) / (2 * l)) * A[1, i - 1, i - 1]
        A[:, i, i] .= diag_val
    end
    z = B -> zeros(Float64, B)
    return HarmonicsBatchWorkspace(
        A,
        zeros(Float64, B, M + 3),
        zeros(Float64, B, M + 3),
        z(B), z(B), z(B), z(B), z(B),
        zeros(Float64, 3, B),
        z(B), z(B), z(B), z(B), z(B), z(B), z(B), z(B), z(B), z(B), z(B),
    )
end

# Per-worker workspace pool keyed by objectid(model).
# Indexed by worker number (1..n_workers); resized on first use or constellation change.
const _HARMONICS_BATCH_POOL = Dict{UInt, Vector{HarmonicsBatchWorkspace}}()
const _HARMONICS_BATCH_POOL_LOCK = ReentrantLock()

function _get_harmonics_batch_pool(
    model::GravitationalHarmonicsModel,
    n_workers::Int,
    batch_size::Int,
)::Vector{HarmonicsBatchWorkspace}
    key = objectid(model)
    pool = get(_HARMONICS_BATCH_POOL, key, nothing)
    if pool !== nothing && length(pool) >= n_workers && size(pool[1].A, 1) >= batch_size
        return pool
    end
    lock(_HARMONICS_BATCH_POOL_LOCK) do
        pool = get(_HARMONICS_BATCH_POOL, key, nothing)
        if pool === nothing || length(pool) < n_workers || size(pool[1].A, 1) < batch_size
            _HARMONICS_BATCH_POOL[key] = [_make_harmonics_batch_workspace(model, batch_size) for _ in 1:n_workers]
        end
    end
    return _HARMONICS_BATCH_POOL[key]
end

# Fetch the batch pool from a per-solve cached ref to avoid the global Dict lookup
# on every RHS call. Falls through to _get_harmonics_batch_pool on first use or resize.
function _get_harmonics_batch_pool_cached!(
    pool_ref::Base.RefValue{Any},
    model::GravitationalHarmonicsModel,
    n_workers::Int,
    batch_size::Int,
)::Vector{HarmonicsBatchWorkspace}
    cached = pool_ref[]
    if cached isa Vector{HarmonicsBatchWorkspace} &&
       length(cached) >= n_workers &&
       !isempty(cached) && size(cached[1].A, 1) >= batch_size
        return cached
    end
    fresh = _get_harmonics_batch_pool(model, n_workers, batch_size)
    pool_ref[] = fresh
    return fresh
end

# Batched harmonics kernel: processes satellites item_start..item_end together.
# Coefficients (C, S, VR01, VR11) are loaded once per (degree, order) pair
# and broadcast across all satellites via @turbo SIMD over the batch dimension.
function _harmonics_flat_batch_kernel!(
    totals::Matrix{Float64},
    model::GravitationalHarmonicsModel,
    sc_state,
    work_items::Vector{Int},
    item_start::Int,
    item_end::Int,
    lpi::SMatrix{3, 3, Float64, 9},
    ws::HarmonicsBatchWorkspace,
)::Nothing
    B = item_end - item_start + 1
    B <= 0 && return nothing

    L = model.L
    M = model.M
    RE = model.reference_radius_m
    μ = model.gm_m3s2
    A = ws.A
    R = ws.R
    I = ws.I

    # Phase 1: gather per-satellite position components and initialise accumulators.
    @inbounds for b = 1:B
        sat_idx = work_items[item_start + b - 1]
        sc = sc_state[sat_idx]
        pos_ii = SVector{3, Float64}(sc[1], sc[2], sc[3])
        rVec = lpi * pos_ii
        r = norm(rVec)
        inv_r_b = 1.0 / r
        ws.rVec[1, b] = rVec[1]
        ws.rVec[2, b] = rVec[2]
        ws.rVec[3, b] = rVec[3]
        ws.s_vec[b]   = rVec[1] * inv_r_b
        ws.t_vec[b]   = rVec[2] * inv_r_b
        ws.u_vec[b]   = rVec[3] * inv_r_b
        ws.inv_r[b]   = inv_r_b
        ws.mass[b]    = sc[7]
        ws.ρ[b]       = RE * inv_r_b
        ws.ρ_np1[b]   = -μ * inv_r_b * ws.ρ[b]
        ws.a1[b] = 0.0
        ws.a2[b] = 0.0
        ws.a3[b] = 0.0
        ws.a4[b] = 0.0
    end

    # Phase 2: A sub-diagonal — A[b, row+1, row] = u[b] * sqrt_2n_plus_3[n] * A[b, row, row].
    # Diagonal A[b,i,i] was set once at workspace construction and is never overwritten.
    @turbo for b = 1:B
        A[b, 2, 1] = ws.u_vec[b] * sqrt_3
    end
    @inbounds for n = 1:(L + 1)
        row = n + 1
        s2n3 = model.sqrt_2n_plus_3[n]
        @turbo for b = 1:B
            A[b, row + 1, row] = ws.u_vec[b] * s2n3 * A[b, row, row]
        end
    end

    # Phase 3: longitude trig recurrence R[b,j] + i*I[b,j] = (s[b]+i*t[b])^(j-1).
    @turbo for b = 1:B
        R[b, 1] = 1.0
        I[b, 1] = 0.0
    end
    @inbounds for j = 2:(M + 2)
        @turbo for b = 1:B
            sv  = ws.s_vec[b]
            tv  = ws.t_vec[b]
            Rn  = R[b, j - 1]
            In  = I[b, j - 1]
            R[b, j] = sv * Rn - tv * In
            I[b, j] = sv * In + tv * Rn
        end
    end

    # Phase 4: main degree/order accumulation — coefficients loaded once per (l,m) pair,
    # broadcast to all B satellites via SIMD.
    max_recur_row = 2
    @inbounds for l = 1:L
        row = l + 1

        if row > max_recur_row
            jmax = min(max(M, 1) + 1, l - 1)
            for j = 1:jmax
                N1v = model.N1[row, j]
                N2v = model.N2[row, j]
                @turbo for b = 1:B
                    A[b, row, j] = ws.u_vec[b] * N1v * A[b, row - 1, j] - N2v * A[b, row - 2, j]
                end
            end
            max_recur_row = row
        end

        next_row = row + 1
        if next_row > max_recur_row
            jmax_next = min(max(M, 1) + 1, l)
            for j = 1:jmax_next
                N1v = model.N1[next_row, j]
                N2v = model.N2[next_row, j]
                @turbo for b = 1:B
                    A[b, next_row, j] = ws.u_vec[b] * N1v * A[b, next_row - 1, j] - N2v * A[b, next_row - 2, j]
                end
            end
            max_recur_row = next_row
        end

        @turbo for b = 1:B
            ws.ρ_np1[b] *= ws.ρ[b]
            ws.rr[b]   = ws.ρ_np1[b] / RE
            ws.sum1[b] = 0.0
            ws.sum2[b] = 0.0
            ws.sum3[b] = 0.0
            ws.sum4[b] = 0.0
        end

        # m=0 zonal term: I[b,1]==0 and R[b,1]==1 always, so D0 = C0*sqrt_2.
        C0      = model.C[row, 1]
        D0      = C0 * sqrt_2
        VR01_r1 = model.VR01[row, 1]
        VR11_r1 = model.VR11[row, 1]
        @turbo for b = 1:B
            ws.sum3[b] += VR01_r1 * A[b, row, 2] * D0
            ws.sum4[b] += VR11_r1 * A[b, row + 1, 2] * D0
        end

        active_orders = model.active_orders_by_degree[row]
        for idx in eachindex(active_orders)
            ord    = active_orders[idx]
            j      = ord + 1
            Cv     = model.C[row, j]
            Sv     = model.S[row, j]
            VR01v  = model.VR01[row, j]
            VR11v  = model.VR11[row, j]
            ordf   = Float64(ord)
            @turbo for b = 1:B
                R_prev = R[b, j - 1]
                I_prev = I[b, j - 1]
                Rj     = R[b, j]
                Ij     = I[b, j]
                D = (Cv * Rj    + Sv * Ij)    * sqrt_2
                E = (Cv * R_prev + Sv * I_prev) * sqrt_2
                F = (Sv * R_prev - Cv * I_prev) * sqrt_2
                mA = ordf * A[b, row, j]
                ws.sum1[b] += mA * E
                ws.sum2[b] += mA * F
                ws.sum3[b] += VR01v * A[b, row, j + 1] * D
                ws.sum4[b] += VR11v * A[b, row + 1, j + 1] * D
            end
        end

        @turbo for b = 1:B
            ws.a1[b] += ws.rr[b] * ws.sum1[b]
            ws.a2[b] += ws.rr[b] * ws.sum2[b]
            ws.a3[b] += ws.rr[b] * ws.sum3[b]
            ws.a4[b] -= ws.rr[b] * ws.sum4[b]
        end
    end

    # Phase 5: back-transform to inertial frame and scatter into totals.
    lpi_t = lpi'
    include_central = model.include_central
    @inbounds for b = 1:B
        sat_idx  = work_items[item_start + b - 1]
        sv       = ws.s_vec[b]
        tv       = ws.t_vec[b]
        uv       = ws.u_vec[b]
        inv_r_b  = ws.inv_r[b]
        mass_b   = ws.mass[b]
        rVec_b   = SVector{3, Float64}(ws.rVec[1, b], ws.rVec[2, b], ws.rVec[3, b])
        g_pp_generic = SVector{3, Float64}(
            -ws.a1[b] - sv * ws.a4[b],
            -ws.a2[b] - tv * ws.a4[b],
            -ws.a3[b] - uv * ws.a4[b],
        )
        g_pp = include_central ?
            g_pp_generic - μ * inv_r_b^3 * rVec_b :
            g_pp_generic
        force_ii = mass_b * (lpi_t * g_pp)
        totals[1, sat_idx] += force_ii[1]
        totals[2, sat_idx] += force_ii[2]
        totals[3, sat_idx] += force_ii[3]
    end
    return nothing
end

@inline function _harmonics_lpi_cache_key(model::GravitationalHarmonicsModel, param::ODEParams, et::Float64)
    return (
        model.planet.name,
        ephemerides_cache_key(param.args.environment_model.ephemerides_model),
        et,
    )
end

function _harmonics_lpi_at!(
    model::GravitationalHarmonicsModel,
    param::ODEParams,
    et::Float64
)::SMatrix{3, 3, Float64, 9}
    shared = param.shared_buffers
    key = _harmonics_lpi_cache_key(model, param, et)

    # The RHS evaluates all spacecraft at the same t. Avoid repeating the same
    # SPICE pxform inside every harmonics force evaluation.
    if shared.harmonics_lpi_key[] == key
        return shared.harmonics_lpi[]
    end

    return lock(shared.harmonics_lpi_lock) do
        if shared.harmonics_lpi_key[] != key
            ephemerides_model = param.args.environment_model.ephemerides_model
            if ephemerides_requires_spice(ephemerides_model)
                Base.Threads.atomic_add!(shared.spice_runtime_counters.planet_pxform_runtime_calls, 1)
            end
            shared.harmonics_lpi[] = planet_frame_lpi(model.planet, et, ephemerides_model)
            shared.harmonics_lpi_key[] = key
        end
        shared.harmonics_lpi[]
    end
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

const _HARMONICS_MODEL_CACHE_LOCK = ReentrantLock()
const _HARMONICS_MODEL_CACHE = Dict{Tuple{Int64, Int64, String, Symbol, Bool, Symbol, Union{Nothing, Float64}, Union{Nothing, Float64}, Bool, UInt}, Any}()

# `GravitationalHarmonicsModel(L, M, coefficients_file, planet)` streams and parses the
# whole coefficients file on every call (row-by-row CSV.Rows iteration, real I/O; some
# planet files here run 16K+ rows). Callers that build a fresh mission config per
# iteration (e.g. a Monte Carlo sample loop constructing thousands of trials in one
# process, same as the SPICE planet-construction case above) redo that parse from
# scratch every time even though the result only depends on the arguments. Cache it.
#
# Safe to share the whole returned model across concurrent callers: `A`, `R`, `I`,
# `N1`, `N2`, `VR01`, `VR11`, `sqrt_2n_plus_3`, `C`, `S` are precomputed, read-only
# recurrence/coefficient tables — nothing in `calcForceTorque`/`wrench` ever mutates
# them. The actual per-call mutable scratch space comes from a separate, properly
# pooled/per-satellite workspace (`_make_harmonics_scratch_workspace`,
# `_harmonics_workspace_for_sat!`), never from the model's own fields.
@inline function _harmonics_model_cache_key(
    L::Int64, M::Int64, coefficients_file::String, coefficient_normalization::Symbol,
    coefficients_normalized::Bool, j2_source::Symbol, reference_radius_m::Union{Nothing, Float64},
    gm_m3s2::Union{Nothing, Float64}, include_central::Bool, planet::AbstractPlanet,
)
    return (L, M, coefficients_file, coefficient_normalization, coefficients_normalized,
            j2_source, reference_radius_m, gm_m3s2, include_central, objectid(planet))
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

Set `include_central=false` to recover the legacy non-central perturbation-only
behavior.
"""
function GravitationalHarmonicsModel(
    L::Int64,
    M::Int64,
    coefficients_file::String,
    planet::P;
    coefficient_normalization::Symbol=:full,
    coefficients_normalized::Bool=true,
    j2_source::Symbol=:file_c20,
    reference_radius_m::Union{Nothing, Float64}=nothing,
    gm_m3s2::Union{Nothing, Float64}=nothing,
    include_central::Bool=true
) where P <: AbstractPlanet
    if L < 0 || M < 0
        throw(ArgumentError("Gravitational harmonics degree/order must be non-negative, got L=$L, M=$M."))
    end
    if M > L
        throw(ArgumentError("Gravitational harmonics order must satisfy M <= L, got L=$L, M=$M."))
    end
    cache_key = _harmonics_model_cache_key(
        L, M, coefficients_file, coefficient_normalization, coefficients_normalized,
        j2_source, reference_radius_m, gm_m3s2, include_central, planet,
    )
    cached = lock(_HARMONICS_MODEL_CACHE_LOCK) do
        get(_HARMONICS_MODEL_CACHE, cache_key, nothing)
    end
    cached !== nothing && return cached::GravitationalHarmonicsModel{P}

    normalized_source = _canonical_harmonics_normalization(coefficient_normalization)

    # Use CSV.Rows for streaming/lazy reading instead of materializing entire file.
    # This is critical for large files like Venus MGNP180U.csv (16K+ rows).
    # We process rows on-the-fly and skip those beyond our requested degree/order.
    # `comment="#"` skips any leading source-metadata header lines (reference
    # radius/GM/citation for this specific gravity solution -- see
    # _read_harmonics_source_metadata above) so they never reach the
    # degree/order/C/S column parser below.

    # First: Quick structure validation using CSV.File to detect headers (minimal overhead for small inspection)
    test_parse = CSV.File(coefficients_file; delim=',', stripwhitespace=true, ignorerepeated=true, comment="#", limit=1)
    has_degree_header = hasproperty(test_parse, :degree)

    # Now do the expensive iteration with streaming to avoid materializing full columns
    C_dict = Dict{Tuple{Int,Int}, Float64}()  # Temporary storage - only entries we need
    S_dict = Dict{Tuple{Int,Int}, Float64}()
    max_degree_file = 0

    # Parse the entire file row-by-row using CSV.Rows (lazy evaluation)
    open(coefficients_file) do io
        rows = CSV.Rows(io; delim=',', stripwhitespace=true, ignorerepeated=true, comment="#")
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

    model = GravitationalHarmonicsModel(
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
        gm_m3s2 === nothing ? _infer_harmonics_gm_m3s2(coefficients_file, planet) : gm_m3s2,
        include_central,
        planet
    )
    return lock(_HARMONICS_MODEL_CACHE_LOCK) do
        get!(_HARMONICS_MODEL_CACHE, cache_key, model)::GravitationalHarmonicsModel{P}
    end
end

"""
    calcForceTorque(model::NBodyGravityModel, x::AbstractVector{Float64}, ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
"""
function calcForceTorque(model::NBodyGravityModel, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_ii = SVector{3, Float64}(x[1], x[2], x[3])
    mass = Float64(x[7])
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
        ParallelPolicy.threaded_foreach_worker_persistent(:rhs_nbody, n_bodies, decision.allotment) do worker_id, k
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
    pos_ii = SVector{3, Float64}(x[1], x[2], x[3])
    mass = Float64(x[7])
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

Compute spherical-harmonics gravity using coefficients stored in the model's internal fully
normalized convention. By default this includes the inverse-square central term, so it should
not be paired with a separate `InverseSquaredGravityModel` for the same primary body.
"""
# Inner kernel: compute harmonics force/torque given a pre-computed inertial→planet-fixed
# rotation matrix. Called directly by the flat-batch parallel region (which already computed
# L_PI once serially) to avoid per-satellite cache-key allocation inside the parallel loop.
@inline function _harmonics_calcforcetorque_with_lpi(
    model::GravitationalHarmonicsModel,
    x::AbstractVector{Float64},
    param::ODEParams,
    i::Int64,
    L_PI::SMatrix{3, 3, Float64, 9},
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_ii = SVector{3, Float64}(x[1], x[2], x[3])
    rVec_cart = L_PI * pos_ii # convert from inertial to planet-fixed frame for gravity calculation
    mass = Float64(x[7])

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
    ρ_np1 = -model.gm_m3s2 * inv_r * ρ
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
    g_pp = if model.include_central
        g_pp_generic - model.gm_m3s2 * inv_r^3 * rVec_cart
    else
        g_pp_generic
    end
    force_ii = mass * L_PI' * g_pp

    if _DEBUG_COMPARE_J2[] && model.L == 2 && model.M == 0
        C20 = model.C[3, 1]
        if isfinite(C20) && C20 != 0.0
            x_pp, y_pp, z_pp = rVec_cart
            r2 = r * r
            r4 = r2 * r2
            z2 = z_pp * z_pp
            J2 = -sqrt(5.0) * C20

            common = 1.5 * model.gm_m3s2 * J2 * RE^2 / r4
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

function calcForceTorque(model::GravitationalHarmonicsModel, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    et = param.shared_buffers.et_start[] + param.shared_buffers.current_time[]
    L_PI = _harmonics_lpi_at!(model, param, et)
    return _harmonics_calcforcetorque_with_lpi(model, x, param, i, L_PI)
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
    RE = model.reference_radius_m
    r = norm(rVec_cart)
    inv_r = 1.0 / r
    s = rVec_cart[1] * inv_r
    n = rVec_cart[2] * inv_r
    u = rVec_cart[3] * inv_r
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
        ρ_np1 = -model.gm_m3s2 / r * ρ
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
        g_pp_generic = SVector{3, Float64}(-a1 - s * a4, -a2 - n * a4, -a3 - u * a4)
        g_pp = if model.include_central
            g_pp_generic - model.gm_m3s2 * inv_r^3 * rVec_cart
        else
            g_pp_generic
        end
        force_pp = x.mass_kg * g_pp
        force_ii = planet_frame.l_pi' * force_pp
        return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
    end
end

@inline gravity_backbone_structure(::GravitationalHarmonicsModel) = :position_only_static_gravity

@inline function gravity_backbone_acceleration_ii(
    model::GravitationalHarmonicsModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::SVector{3, Float64}
    force_ii, _ = wrench(model, x, env, t)
    return force_ii / x.mass_kg
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

    # Dipole field equation with Earth's dipole MOMENT pointing toward the
    # geomagnetic SOUTH pole (m_hat = -M_HAT_ECEF):
    #   B = B0 (R/r)^3 [3 (m_hat.r_hat) r_hat - m_hat]
    #     = B0 (R/r)^3 [M_HAT_ECEF - 3 cos_colat r_hat]
    # Checks: equator (cos_colat = 0) gives B = B0 * M_HAT_ECEF (north, B0);
    # north pole gives -2 B0 M_HAT_ECEF (down). The pre-2026-07 sign returned
    # the ANTIPARALLEL field (170 deg from IGRF at LEO).
    B_ecef = B0_2020 * (R_EARTH_MODEL / r_norm)^3 * (M_HAT_ECEF - 3 * cos_colat * r_hat)

    return L_PI' * B_ecef
end

"""
    get_magnetic_field(date::DateTime, lat_rad::Number, lon_rad::Number, alt_m::Number, L_PI::MMatrix{3, 3, Float64})

Computes the Earth's magnetic field vector in the inertial frame using the
International Geomagnetic Reference Field (IGRF).

# Args

- `date`: The `DateTime` of the measurement (sets the IGRF epoch).
- `lat_rad`: The geodetic latitude of the observer [radians].
- `lon_rad`: The longitude of the observer [radians].
- `alt_m`: The altitude above the WGS84 ellipsoid [meters].
- `L_PI`: The inertial-to-planet-fixed rotation matrix.

# Returns

- A 3-element vector representing the magnetic field in the inertial frame in
  nanoTeslas [nT]. Note the unit: [`get_magnetic_field_dipole`](@ref) returns
  Tesla, so the two are NOT drop-in interchangeable.
"""
function get_magnetic_field(date::DateTime, lat_rad::Number, lon_rad::Number, alt_m::Number, L_PI::MMatrix{3, 3, Float64})
    # The IGRF evaluation returns NED components in nT.
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

"""
    MagneticTorqueRodModel <: AbstractForceTorqueModel

Dynamic effector for spacecraft magnetic torque rods / magnetorquers.

Sums the body-frame dipole moment `m` over every [`Magnet`](@ref) attached to
any link of the spacecraft, samples the Earth magnetic field at the
spacecraft's current position, and returns the resulting body-frame torque
`τ = m × B` ([`calculate_magnetic_torque`](@ref)). Produces zero force and
zero torque when the spacecraft carries no magnets or when attitude state is
unavailable (`orientation_sim=false`).

Field source (`field_model`):

- `:dipole` (default): fast tilted-dipole approximation
  ([`get_magnetic_field_dipole`](@ref)), epoch-2020 constants, typically good
  to 10-20% at LEO. Bit-identical to the pre-option behavior.
- `:igrf`: the International Geomagnetic Reference Field evaluated at the
  spacecraft's geodetic position. Requires `igrf_year` (decimal year, e.g.
  `2025.4`) — the secular field drift over a single simulation is negligible,
  so one epoch per model is deliberate and keeps the integrator hot loop free
  of calendar conversions.

Both sources are Earth models; neither is meaningful at another central body.
"""
struct MagneticTorqueRodModel <: AbstractForceTorqueModel
    field_model::Symbol
    igrf_year::Float64

    function MagneticTorqueRodModel(; field_model::Symbol=:dipole, igrf_year::Real=NaN)
        field_model in (:dipole, :igrf) ||
            throw(ArgumentError("field_model must be :dipole or :igrf, got $(repr(field_model))"))
        # SatelliteToolboxGeomagneticField's IGRF hard-rejects epochs outside
        # [1900, 2035) (and warns about reduced accuracy past 2030), so reject
        # unsupported epochs here at configuration time instead of at the
        # first wrench evaluation.
        if field_model === :igrf && !(isfinite(igrf_year) && 1900.0 <= igrf_year < 2035.0)
            throw(ArgumentError(
                "field_model=:igrf requires igrf_year in [1900, 2035) (decimal year, e.g. 2025.4)"))
        end
        return new(field_model, Float64(igrf_year))
    end
end

"""
    _magnetic_field_inertial(model, l_pi, pos_pp, lat_rad, lon_rad, alt_m)

Inertial-frame magnetic field [Tesla] for the model's configured field source.
The IGRF branch converts from the library's nT to Tesla here so both branches
share one unit contract.
"""
@inline function _magnetic_field_inertial(
    model::MagneticTorqueRodModel,
    l_pi::SMatrix{3, 3, Float64, 9},
    pos_pp::SVector{3, Float64},
    lat_rad::Float64,
    lon_rad::Float64,
    alt_m::Float64,
)::SVector{3, Float64}
    if model.field_model === :igrf
        B_ned_nT = igrf(model.igrf_year, alt_m, lat_rad, lon_rad, Val(:geodetic))
        B_pp_nT = ned_to_ecef(B_ned_nT, lat_rad, lon_rad, alt_m)
        return SVector{3, Float64}(l_pi' * B_pp_nT) .* 1e-9
    end
    return get_magnetic_field_dipole(pos_pp, MMatrix{3, 3, Float64}(l_pi))
end

@inline function _total_dipole_moment_body(spacecraft)::SVector{3, Float64}
    m_total = SVector{3, Float64}(0.0, 0.0, 0.0)
    for link in spacecraft.links
        for magnet in link.magnets
            m_total = m_total + SVector{3, Float64}(magnet.m)
        end
    end
    return m_total
end

function calcForceTorque(model::MagneticTorqueRodModel, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    zero_wrench = (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
    if !param.args.mission_configuration.orientation_sim
        return zero_wrench
    end
    if i < 1 || i > length(param.args.dynamics_model.spacecraft)
        return zero_wrench
    end
    if !hasproperty(x, :q)
        return zero_wrench
    end

    spacecraft = param.args.dynamics_model.spacecraft[i]
    m_body = _total_dipole_moment_body(spacecraft)
    iszero(m_body) && return zero_wrench

    pos_ii = SVector{3, Float64}(x[1], x[2], x[3])
    et = param.shared_buffers.et_start[] + param.shared_buffers.current_time[]
    planet = param.args.environment_model.planet
    l_pi = planet_frame_lpi(planet, et, param.args.environment_model.ephemerides_model)
    pos_pp = l_pi * pos_ii
    if model.field_model === :igrf
        alt_m, lat_rad, lon_rad = rtolatlong(pos_pp, planet)
        B_ii = _magnetic_field_inertial(model, SMatrix{3, 3, Float64, 9}(l_pi), pos_pp, lat_rad, lon_rad, alt_m)
    else
        B_ii = get_magnetic_field_dipole(pos_pp, MMatrix{3, 3, Float64}(l_pi))
    end

    q_ib = getproperty(x, :q)
    q_body = SVector{4, Float64}(Float64(q_ib[1]), Float64(q_ib[2]), Float64(q_ib[3]), Float64(q_ib[4]))
    B_body = rot(q_body) * B_ii
    torque_body = calculate_magnetic_torque(m_body, B_body)
    return SVector{3, Float64}(0.0, 0.0, 0.0), torque_body
end

@inline environment_requirements(::MagneticTorqueRodModel) = EffectorEnvironmentRequirements(planet_frame=true)

@inline function wrench(
    model::MagneticTorqueRodModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    zero_wrench = (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
    if x.q_ib === nothing || x.spacecraft === nothing
        return zero_wrench
    end

    m_body = _total_dipole_moment_body(x.spacecraft)
    iszero(m_body) && return zero_wrench

    planet_frame = env.planet_frame
    planet_frame === nothing && throw(ArgumentError("MagneticTorqueRodModel wrench requires env.planet_frame."))

    B_ii = _magnetic_field_inertial(
        model,
        planet_frame.l_pi,
        planet_frame.pos_pp,
        planet_frame.lat_rad,
        planet_frame.lon_rad,
        planet_frame.alt_m,
    )
    B_body = rot(x.q_ib) * B_ii
    torque_body = calculate_magnetic_torque(m_body, B_body)
    return SVector{3, Float64}(0.0, 0.0, 0.0), torque_body
end

"""
    LVLHCascadeAttitudeControlModel <: AbstractForceTorqueModel

Closed-loop LVLH-pointing attitude controller as a dynamic effector: an outer
proportional attitude loop commands a rate-limited body rate toward a fixed
LVLH-relative attitude setpoint, and an inner proportional rate loop converts
the rate error into a body torque with a per-axis saturation. This cascade
(rate-command) architecture is the common small-satellite ADCS structure
(e.g. Blue Canyon XACT-class flight software).

    omega_cmd = clamp.(-k_out .* theta_err, -w_max, w_max)
    tau       = clamp.(k_rate .* (omega_cmd - omega_rel), -tau_max, tau_max)

`theta_err` is the small-angle rotation vector (body coordinates) from the
commanded attitude to the current attitude, extracted from
`R_err = R_bl * R_cmd'`; `omega_rel` is the body rate relative to the rotating
LVLH frame. The LVLH triad is +Z nadir, +Y = nadir x velocity (negative orbit
normal), +X completing (near along-track), and `q_cmd_lb` is the commanded
LVLH-to-body quaternion (scalar-last, identity = LVLH-aligned).

The torque is applied directly to the body (zero force): the actuator stack
(wheel allocation, momentum management) is abstracted into `tau_max`. The
control law is evaluated continuously at the integrator rate — a
flight-software discrete update (a few Hz) is well inside the closed-loop
bandwidth this law can express, so the continuous approximation is benign for
tracking studies. Signs follow the engine's own conventions
([`quaternion_derivative`](@ref) gives `theta_dot = +omega_body` near
identity), which the probe suite pins with an end-to-end convergence test.
Returns zero wrench when attitude state is unavailable.
"""
struct LVLHCascadeAttitudeControlModel <: AbstractForceTorqueModel
    q_cmd_lb::SVector{4, Float64}
    k_out::SVector{3, Float64}
    w_max::Float64
    k_rate::SVector{3, Float64}
    tau_max::Float64
    tau_ff::SVector{3, Float64}

    function LVLHCascadeAttitudeControlModel(;
        q_cmd_lb::AbstractVector{<:Real}=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
        k_out::AbstractVector{<:Real},
        w_max::Real,
        k_rate::AbstractVector{<:Real},
        tau_max::Real,
        tau_ff::AbstractVector{<:Real}=SVector{3, Float64}(0.0, 0.0, 0.0),
    )
        q = SVector{4, Float64}(q_cmd_lb...)
        qn = norm(q)
        (isfinite(qn) && qn > 0) ||
            throw(ArgumentError("q_cmd_lb must be a finite, nonzero quaternion"))
        ko = SVector{3, Float64}(k_out...)
        kr = SVector{3, Float64}(k_rate...)
        all(isfinite, ko) && all(>=(0.0), ko) ||
            throw(ArgumentError("k_out must be finite and nonnegative [1/s]"))
        all(isfinite, kr) && all(>=(0.0), kr) ||
            throw(ArgumentError("k_rate must be finite and nonnegative [N m s]"))
        (isfinite(w_max) && w_max > 0) ||
            throw(ArgumentError("w_max must be finite and positive [rad/s]"))
        (isfinite(tau_max) && tau_max > 0) ||
            throw(ArgumentError("tau_max must be finite and positive [N m]"))
        tf = SVector{3, Float64}(tau_ff...)
        all(isfinite, tf) ||
            throw(ArgumentError("tau_ff must be finite [N m]"))
        return new(q / qn, ko, Float64(w_max), kr, Float64(tau_max), tf)
    end
end

"""
    _lvlh_cascade_torque(model, pos_ii, vel_ii, q_ib, w_body)

Body-frame control torque of [`LVLHCascadeAttitudeControlModel`](@ref) for the
given inertial position/velocity and attitude state. Returns zero for
degenerate orbits (near-zero radius or angular momentum).
"""
@inline function _lvlh_cascade_torque(
    model::LVLHCascadeAttitudeControlModel,
    pos_ii::SVector{3, Float64},
    vel_ii::SVector{3, Float64},
    q_ib::SVector{4, Float64},
    w_body::SVector{3, Float64},
)::SVector{3, Float64}
    r2 = dot(pos_ii, pos_ii)
    h_ii = cross(pos_ii, vel_ii)
    if !(r2 > 0.0) || !(dot(h_ii, h_ii) > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    z_l = -pos_ii / sqrt(r2)
    y_l = cross(z_l, vel_ii)
    y_l = y_l / norm(y_l)
    x_l = cross(y_l, z_l)
    R_li = SMatrix{3, 3, Float64, 9}(
        x_l[1], y_l[1], z_l[1],
        x_l[2], y_l[2], z_l[2],
        x_l[3], y_l[3], z_l[3],
    )                                    # inertial -> LVLH (rows = triad)
    R_bi = rot(q_ib)                     # inertial -> body
    R_bl = R_bi * R_li'                  # LVLH -> body
    R_e = R_bl * rot(model.q_cmd_lb)'    # commanded -> current, body coords
    # Rotation-log error extraction. The plain vee-map (0.5*(R_e[2,3]-R_e[3,2]),
    # ...) returns sin(theta)*axis, which vanishes at a half-turn — an
    # undesired zero-command equilibrium at 180 deg and weak authority near it
    # (Codex review). The quaternion log 2*atan2(|qv|, qs)*qv/|qv| equals the
    # vee-map to second order at small angles (2*qs*qv == sin(theta)*axis)
    # but grows monotonically to pi*axis at the antipode. Index order/signs
    # follow the engine conventions pinned empirically in the probe suite:
    # rotating the body by +theta yields rot(q) ≈ I - S(theta); the opposite
    # sign turns the attitude loop into positive feedback.
    tr_e = R_e[1, 1] + R_e[2, 2] + R_e[3, 3]
    q_e = if tr_e > 0.0
        s = sqrt(tr_e + 1.0) * 2
        SVector{4, Float64}((R_e[2, 3] - R_e[3, 2]) / s, (R_e[3, 1] - R_e[1, 3]) / s,
                            (R_e[1, 2] - R_e[2, 1]) / s, s / 4)
    elseif R_e[1, 1] > R_e[2, 2] && R_e[1, 1] > R_e[3, 3]
        s = sqrt(1.0 + R_e[1, 1] - R_e[2, 2] - R_e[3, 3]) * 2
        SVector{4, Float64}(s / 4, (R_e[1, 2] + R_e[2, 1]) / s,
                            (R_e[1, 3] + R_e[3, 1]) / s, (R_e[2, 3] - R_e[3, 2]) / s)
    elseif R_e[2, 2] > R_e[3, 3]
        s = sqrt(1.0 + R_e[2, 2] - R_e[1, 1] - R_e[3, 3]) * 2
        SVector{4, Float64}((R_e[1, 2] + R_e[2, 1]) / s, s / 4,
                            (R_e[2, 3] + R_e[3, 2]) / s, (R_e[3, 1] - R_e[1, 3]) / s)
    else
        s = sqrt(1.0 + R_e[3, 3] - R_e[1, 1] - R_e[2, 2]) * 2
        SVector{4, Float64}((R_e[1, 3] + R_e[3, 1]) / s, (R_e[2, 3] + R_e[3, 2]) / s,
                            s / 4, (R_e[1, 2] - R_e[2, 1]) / s)
    end
    q_e = q_e[4] < 0.0 ? -q_e : q_e      # shortest rotation (qs >= 0)
    qv = SVector{3, Float64}(q_e[1], q_e[2], q_e[3])
    sv = norm(qv)
    theta_err = sv > 1e-12 ? (2.0 * atan(sv, q_e[4]) / sv) * qv : 2.0 * qv
    # Body-frame LVLH feed-forward: with the body-rate quaternion kinematics
    # (see quaternion_derivative), the state that tracks the LVLH frame is
    # its rotation rate expressed in BODY coordinates — pinned by the
    # torque-free tracking and closed-loop convergence probes.
    w_lvlh_body = R_bi * (h_ii / r2)
    w_rel = w_body - w_lvlh_body
    w_cmd = clamp.(-model.k_out .* theta_err, -model.w_max, model.w_max)
    # tau_ff: constant body-frame feedforward torque, the stateless equivalent
    # of converged integral action against a quasi-steady disturbance (a pure
    # PD cascade holds a steady offset error of |tau_dist| / (k_rate*k_out)
    # under a constant torque; the feedforward nulls it). True discrete-rate
    # integral action belongs in a stateful control effector - a continuous
    # wrench cannot hold accumulator state correctly under adaptive stepping.
    return clamp.(model.k_rate .* (w_cmd - w_rel) + model.tau_ff, -model.tau_max, model.tau_max)
end

function calcForceTorque(model::LVLHCascadeAttitudeControlModel, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    zero_wrench = (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
    if !param.args.mission_configuration.orientation_sim
        return zero_wrench
    end
    if !hasproperty(x, :q) || !hasproperty(x, :ω)
        return zero_wrench
    end
    pos_ii = SVector{3, Float64}(x[1], x[2], x[3])
    vel_ii = SVector{3, Float64}(x[4], x[5], x[6])
    q_ib = SVector{4, Float64}(Float64.(getproperty(x, :q))...)
    w_body = SVector{3, Float64}(Float64.(getproperty(x, :ω))...)
    torque = _lvlh_cascade_torque(model, pos_ii, vel_ii, q_ib, w_body)
    return SVector{3, Float64}(0.0, 0.0, 0.0), torque
end

@inline environment_requirements(::LVLHCascadeAttitudeControlModel) = EffectorEnvironmentRequirements()

@inline function wrench(
    model::LVLHCascadeAttitudeControlModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    zero_wrench = (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
    if x.q_ib === nothing || x.ω_body === nothing
        return zero_wrench
    end
    torque = _lvlh_cascade_torque(model, x.pos_ii, x.vel_ii, x.q_ib, x.ω_body)
    return SVector{3, Float64}(0.0, 0.0, 0.0), torque
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
