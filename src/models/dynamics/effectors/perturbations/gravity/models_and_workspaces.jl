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

@inline function _spice_query_name_uncached(name::String)::String
    key = _canonical_spice_name(name)
    if endswith(key, "_barycenter")
        return key
    end
    return key in _SPICE_FORCE_BARYCENTER_BODIES ? key * "_barycenter" : key
end

# Canonical SPICE body names are resolved on the RHS hot path: the per-satellite
# environment sample asks for one per third body per spacecraft per derivative
# evaluation (simulation/engine/effector_sampling.jl), and `strip`, `lowercase`
# and `replace` each allocate a fresh short string every time. At 256
# spacecraft, degree-50 harmonics and a Sun/Moon third body, that name
# resolution was 741 of the 7606 profile samples taken inside the constellation
# RHS -- ~10% of the derivative evaluation spent lowercasing "Sun" and "Moon"
# once per spacecraft -- and 706 of the 1457 sampled allocation events.
#
# The names come from a handful of model and planet fields and never change
# during a run, so the resolved form is interned and the hot path becomes one
# dictionary read with no allocation.
#
# Publication is copy-on-write behind a lock: a reader takes the table by an
# atomic load and that table is never mutated afterwards, so a reader can never
# observe a rehash in progress even though the RHS resolves names from many
# threads. The interned value is exactly the String the uncached expression
# produces, so every downstream comparison, dictionary key and SPICE query is
# unchanged.
mutable struct _SpiceQueryNameIntern
    @atomic table::Dict{String, String}
end

const _SPICE_QUERY_NAME_INTERN = _SpiceQueryNameIntern(Dict{String, String}())
const _SPICE_QUERY_NAME_INTERN_LOCK = ReentrantLock()

@noinline function _intern_spice_query_name(name::String)::String
    return lock(_SPICE_QUERY_NAME_INTERN_LOCK) do
        current = @atomic :acquire _SPICE_QUERY_NAME_INTERN.table
        existing = get(current, name, nothing)
        existing === nothing || return existing
        resolved = _spice_query_name_uncached(name)
        # Copy-on-write: publish a table that no reader is walking.
        updated = copy(current)
        updated[name] = resolved
        @atomic :release _SPICE_QUERY_NAME_INTERN.table = updated
        return resolved
    end
end

@inline function _spice_query_name(name::String)::String
    table = @atomic :acquire _SPICE_QUERY_NAME_INTERN.table
    cached = get(table, name, nothing)
    cached === nothing || return cached
    return _intern_spice_query_name(name)
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
    pos_primary_k_all = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n_bodies]
    body_force_ii = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:n_bodies]
    return NBodyScratchWorkspace(pos_primary_k_all, body_force_ii)
end

@inline function _ensure_nbody_workspace_capacity!(
    workspace::NBodyScratchWorkspace,
    n_bodies::Int,
    n_workers::Int
)::NBodyScratchWorkspace
    if length(workspace.pos_primary_k_all) < n_bodies
        resize!(workspace.pos_primary_k_all, n_bodies)
    end
    if length(workspace.body_force_ii) < n_bodies
        resize!(workspace.body_force_ii, n_bodies)
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
# so that @inbounds loops over the batch are unit-stride in memory.
# ---------------------------------------------------------------------------

# Row r of the ALF triangle lives in slot ((r-1) % 3) + 1 of the rolling window.
# Degree l touches rows l, l+1 and l+2, which are three consecutive rows and so
# always three distinct slots.
@inline _harmonics_row_slot(row::Int)::Int = ((row - 1) % 3) + 1

struct HarmonicsBatchWorkspace
    # Rolling 3-row window over the ALF triangle, not the whole [B, L+4, L+4]
    # triangle. Degree l reads rows l and l+1 and writes row l+2, so exactly
    # three rows are live at any point; row r lives in slot ((r-1) % 3) + 1.
    #
    # The full triangle was 91 MB at B=4096, L=50 -- past this box's 64 MB L3,
    # which is what made the kernel bandwidth-bound at large B (see the
    # register-blocking result: strip-mining regressed 0.87x median because it
    # turned contiguous column walks into strided gathers). The window is
    # 3/(L+4) of that: 5.1 MB at the same point, which stays resident.
    A::Array{Float64, 3}    # [B, 3, L+4]
    # Position-independent diagonal A[r, r], previously stored in the triangle
    # itself and never overwritten. A rolling slot cannot hold it across the
    # rows that share the slot, so it is precomputed once per model instead.
    diag::Vector{Float64}    # [L+4]
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
    A = zeros(Float64, B, 3, L + 4)
    # Diagonal elements are position-independent. Same recurrence, same order,
    # same values as the in-triangle version this replaces.
    diag = zeros(Float64, L + 4)
    diag[1] = 1.0
    @inbounds for l = 1:(L + 2)
        i = l + 1
        diag[i] = sqrt((2 * l + 1) / (2 * l)) * diag[i - 1]
    end
    z = B -> zeros(Float64, B)
    return HarmonicsBatchWorkspace(
        A,
        diag,
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
            # Grow, never shrink either dimension: the per-region width search
            # alternates widths, and rebuilding on every switch would allocate
            # the whole pool each time. ponytail: keeps the widest x tallest
            # pool for the model's lifetime; a per-width cache if memory matters.
            if pool !== nothing
                n_workers = max(n_workers, length(pool))
                batch_size = max(batch_size, size(pool[1].A, 1))
            end
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
# and reused across the batch, whose loops LLVM is free to vectorise.
#
# BIT-IDENTICAL to `_harmonics_scalar_force_ii` by construction. Two properties
# make that true and both are load-bearing:
#
#   1. No `@turbo` or `@fastmath`, and `@simd ivdep` only on the three
#      per-degree batch loops (`b = 1:B`), never across `l` or `j`. `@turbo`
#      and `@fastmath` license reassociation and FMA contraction, which is
#      what made the previous batched kernel round differently from the
#      scalar one; `@simd ivdep` licenses neither on its own, it only tells
#      LLVM the iterations don't alias. That is true here regardless of batch
#      size: each `b` writes only its own workspace slot and reads no other
#      iteration's, so there is no loop-carried dependency to reorder and no
#      reduction for the annotation to reassociate — the per-satellite
#      floating-point operations stay exactly where the scalar kernel puts
#      them. Plain `@inbounds` on the other loops leaves LLVM to prove
#      vectorisability on its own, which it does above a batch-size threshold
#      only; `@simd ivdep` makes the same vectorisation unconditional (see
#      `git log -1 cd212833aa`).
#   2. The nesting is degree, then order, then batch. For any one satellite the
#      sum1..sum4 accumulations therefore run in exactly the scalar kernel's
#      sequence. Hoisting the batch loop outwards, or reducing across the batch,
#      would break it.
function _harmonics_flat_batch_kernel!(
    slots::Array{Float64, 3},
    eff_idx::Int,
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
        rVec, inv_r_b, s_b, t_b, u_b = _harmonics_frame_terms(lpi, pos_ii)
        ws.rVec[1, b] = rVec[1]
        ws.rVec[2, b] = rVec[2]
        ws.rVec[3, b] = rVec[3]
        ws.s_vec[b]   = s_b
        ws.t_vec[b]   = t_b
        ws.u_vec[b]   = u_b
        ws.inv_r[b]   = inv_r_b
        ws.mass[b]    = sc[7]
        ws.ρ[b]       = RE * inv_r_b
        ws.ρ_np1[b]   = -μ * inv_r_b * ws.ρ[b]
        ws.a1[b] = 0.0
        ws.a2[b] = 0.0
        ws.a3[b] = 0.0
        ws.a4[b] = 0.0
    end

    # Phase 2: seed the two base rows of the recursion. Every later row is built
    # inside Phase 4, into its rolling slot, immediately before the degree that
    # reads it — the sub-diagonal pre-pass over all L+1 rows is gone with the
    # full triangle.
    #
    # `zcol` here and in Phase 4 is the one column per row that the kernel reads
    # but never writes, and which therefore has to hold zero: column r+1 while
    # r-1 <= M (the triangle's zero above the diagonal), column M+2 past that
    # (the order truncation — mathematically A[l, M+1] is nonzero for M+1 <= l,
    # but the recursion stops at M+1 columns, so the full-triangle version read
    # its allocation zero there and this must too). Every other column read for
    # row r is written by the recursion, the sub-diagonal or the diagonal.
    # Zeroing it is what makes a slot safe to reuse — three rows later within
    # this call, and on the next call, which maps different rows onto the same
    # three slots and would otherwise leave a late row's columns visible to an
    # early one.
    slot_1 = _harmonics_row_slot(1)
    slot_2 = _harmonics_row_slot(2)
    diag_2 = ws.diag[2]
    @inbounds for b = 1:B
        A[b, slot_1, 1] = 1.0
        A[b, slot_1, 2] = 0.0
        A[b, slot_2, 1] = ws.u_vec[b] * sqrt_3
        A[b, slot_2, 2] = diag_2
        A[b, slot_2, 3] = 0.0
    end

    # Phase 3: longitude trig recurrence R[b,j] + i*I[b,j] = (s[b]+i*t[b])^(j-1).
    @inbounds for b = 1:B
        R[b, 1] = 1.0
        I[b, 1] = 0.0
    end
    @inbounds for j = 2:(M + 2)
        @inbounds for b = 1:B
            sv  = ws.s_vec[b]
            tv  = ws.t_vec[b]
            Rn  = R[b, j - 1]
            In  = I[b, j - 1]
            R[b, j] = sv * Rn - tv * In
            I[b, j] = sv * In + tv * Rn
        end
    end

    # Phase 4: main degree/order accumulation — coefficients loaded once per (l,m)
    # pair and reused across the batch. The nesting (degree, then order, then
    # batch) is what makes this bit-identical to the scalar kernel: for a fixed
    # satellite the sums are accumulated in exactly the scalar sequence.
    @inbounds for l = 1:L
        row = l + 1
        next_row = row + 1
        slot_prev = _harmonics_row_slot(row - 1)
        slot_row = _harmonics_row_slot(row)
        slot_next = _harmonics_row_slot(next_row)

        # Build row l+2, the only row this degree creates. The full-triangle
        # version wrapped this in `next_row > max_recur_row` beside a companion
        # `row > max_recur_row` branch; the companion was dead and this guard
        # always true, because max_recur_row is l+2 after iteration l while row
        # is l+1. One row per degree is exactly what the three live slots hold.
        jmax_next = min(max(M, 1) + 1, l)
        for j = 1:jmax_next
            N1v = model.N1[next_row, j]
            N2v = model.N2[next_row, j]
            @inbounds for b = 1:B
                A[b, slot_next, j] = ws.u_vec[b] * N1v * A[b, slot_row, j] - N2v * A[b, slot_prev, j]
            end
        end
        # m_cap, not M: the recursion's column bound is min(max(M, 1) + 1, l), so
        # at M = 0 it still writes column 2 and zeroing M + 2 = 2 here would
        # clobber it. Every M >= 1 is unaffected, which is why only the
        # zonal-only configuration saw it.
        m_cap = max(M, 1)
        zcol_next = (next_row - 1) <= m_cap ? (next_row + 1) : (m_cap + 2)
        if zcol_next <= L + 4
            @inbounds for b = 1:B
                A[b, slot_next, zcol_next] = 0.0
            end
        end
        # Sub-diagonal and diagonal, written after the zeroing so that a zcol
        # landing on either is restored rather than lost.
        sub_j = next_row - 1
        s2n3 = model.sqrt_2n_plus_3[next_row - 2]
        sub_diag_src = ws.diag[sub_j]
        diag_next = ws.diag[next_row]
        @inbounds for b = 1:B
            A[b, slot_next, sub_j] = ws.u_vec[b] * s2n3 * sub_diag_src
            A[b, slot_next, next_row] = diag_next
        end

        @inbounds for b = 1:B
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
        @inbounds @simd ivdep for b = 1:B
            ws.sum3[b] += VR01_r1 * A[b, slot_row, 2] * D0
            ws.sum4[b] += VR11_r1 * A[b, slot_next, 2] * D0
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
            @inbounds @simd ivdep for b = 1:B
                R_prev = R[b, j - 1]
                I_prev = I[b, j - 1]
                Rj     = R[b, j]
                Ij     = I[b, j]
                D = (Cv * Rj    + Sv * Ij)    * sqrt_2
                E = (Cv * R_prev + Sv * I_prev) * sqrt_2
                F = (Sv * R_prev - Cv * I_prev) * sqrt_2
                mA = ordf * A[b, slot_row, j]
                ws.sum1[b] += mA * E
                ws.sum2[b] += mA * F
                ws.sum3[b] += VR01v * A[b, slot_row, j + 1] * D
                ws.sum4[b] += VR11v * A[b, slot_next, j + 1] * D
            end
        end

        @inbounds @simd ivdep for b = 1:B
            ws.a1[b] += ws.rr[b] * ws.sum1[b]
            ws.a2[b] += ws.rr[b] * ws.sum2[b]
            ws.a3[b] += ws.rr[b] * ws.sum3[b]
            ws.a4[b] -= ws.rr[b] * ws.sum4[b]
        end
    end

    # Phase 5: back-transform to inertial frame and scatter into totals.
    include_central = model.include_central
    @inbounds for b = 1:B
        sat_idx  = work_items[item_start + b - 1]
        sv       = ws.s_vec[b]
        tv       = ws.t_vec[b]
        uv       = ws.u_vec[b]
        inv_r_b  = ws.inv_r[b]
        mass_b   = ws.mass[b]
        rVec_b   = SVector{3, Float64}(ws.rVec[1, b], ws.rVec[2, b], ws.rVec[3, b])
        force_ii, _ = _harmonics_back_transform(
            ws.a1[b], ws.a2[b], ws.a3[b], ws.a4[b], sv, tv, uv, inv_r_b,
            mass_b, rVec_b, lpi, μ, include_central,
        )
        slots[1, eff_idx, sat_idx] = force_ii[1]
        slots[2, eff_idx, sat_idx] = force_ii[2]
        slots[3, eff_idx, sat_idx] = force_ii[3]
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
            shared.harmonics_lpi[] = ephemerides_frame_lpi(ephemerides_model, model.planet, et)
            shared.harmonics_lpi_key[] = key
        end
        shared.harmonics_lpi[]
    end
end
