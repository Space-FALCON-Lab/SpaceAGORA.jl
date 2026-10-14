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
    # Collect-then-sum: each body's contribution goes into its own slot (in
    # parallel or not), then the slots are summed in body order on this thread.
    # Serial and threaded evaluation are bit-identical by construction.
    body_force_ii = workspace.body_force_ii
    if use_threads
        ParallelPolicy.threaded_collect_persistent!(:rhs_nbody, body_force_ii, n_bodies, decision.allotment) do k
            @inbounds _nbody_body_force_ii(pos_primary_k_all[k], pos_ii, mass, model.body_mus[k])
        end
    else
        @inbounds for k in 1:n_bodies
            body_force_ii[k] = _nbody_body_force_ii(pos_primary_k_all[k], pos_ii, mass, model.body_mus[k])
        end
    end
    @inbounds for k in 1:n_bodies
        force_ii .+= body_force_ii[k]
    end
    # Per satellite, and under satellite_batch on a Polyester worker task: pass
    # the run's context rather than relying on task-local lookup.
    ParallelPolicy.record_policy_observation!(
        :multibody;
        mode=decision.mode,
        num_items=n_bodies,
        use_threads=use_threads,
        elapsed_ns=(time_ns() - started_ns),
        ctx=ParallelPolicy.policy_context_hint(param)
    )

    return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline environment_requirements(model::NBodyGravityModel) = EffectorEnvironmentRequirements(third_body_names=model.body_names)

@inline gravity_backbone_kick_structure(::NBodyGravityModel) = :velocity_kick_explicit

# The one n-body acceleration kernel. The wrench path and the flat batch
# pre-pass both call it, so the two routes perform the same arithmetic in the
# same order and agree bit for bit.
@inline function _nbody_acceleration_from_positions(
    model::NBodyGravityModel,
    pos_ii::SVector{3, Float64},
    positions_ii,
)::SVector{3, Float64}
    accel_ii = MVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for k in eachindex(positions_ii)
        pos_primary_k = positions_ii[k]
        pos_spacecraft_k = pos_primary_k - pos_ii
        pos_spacecraft_k_mag = norm(pos_spacecraft_k)
        accel_ii .+= model.body_mus[k] * (
            (pos_spacecraft_k / pos_spacecraft_k_mag^3) - (pos_primary_k / norm(pos_primary_k)^3)
        )
    end
    return SVector{3, Float64}(accel_ii)
end

@noinline function _nbody_force_ii(
    model::NBodyGravityModel,
    pos_ii::SVector{3, Float64},
    mass::Float64,
    positions_ii,
)::SVector{3, Float64}
    return mass * _nbody_acceleration_from_positions(model, pos_ii, positions_ii)
end

@inline function _nbody_acceleration_ii(
    model::NBodyGravityModel,
    x::StateSample,
    third_bodies::ThirdBodyEphemerisSample,
)::SVector{3, Float64}
    return _nbody_acceleration_from_positions(model, x.pos_ii, third_bodies.positions_ii)
end

@inline function wrench(
    model::NBodyGravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    third_bodies = env.third_bodies
    third_bodies === nothing && throw(ArgumentError("NBodyGravityModel wrench requires env.third_bodies."))
    return _nbody_force_ii(model, x.pos_ii, x.mass_kg, third_bodies.positions_ii), SVector{3, Float64}(0.0, 0.0, 0.0)
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

@noinline function _nbody_body_position_from_cache_j2000_m(
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


"""Return configured third-body acceleration at an arbitrary ephemeris epoch."""
function nbody_acceleration_ii_at_epoch(
    model::NBodyGravityModel,
    pos_ii::SVector{3, Float64},
    param::ODEParams,
    et::Float64,
)::SVector{3, Float64}
    primary_body_name = _spice_query_name(model.primary_body_name)
    nbody_cache_entry = param.shared_buffers.nbody_ephemeris_cache[]
    spice_rhs_memo_enabled = param.shared_buffers.spice_rhs_memo_enabled[]
    spice_rhs_memo = param.shared_buffers.spice_rhs_memo
    counter = param.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls
    positions_ii = SVector{3, Float64}[]

    for (k, body_name) in pairs(model.body_names)
        body_name_spice = _spice_query_name(body_name)
        pos_primary_body_j2000_m = if nbody_cache_entry isa NBodyEphemerisCache
            cached = _nbody_body_position_from_cache_j2000_m(
                nbody_cache_entry,
                et,
                body_name_spice,
                primary_body_name,
            )
            cached === nothing ? _nbody_body_position_from_spice_j2000_m(
                body_name_spice,
                et,
                primary_body_name,
                spice_rhs_memo_enabled,
                spice_rhs_memo,
                counter,
            ) : cached
        else
            _nbody_body_position_from_spice_j2000_m(
                body_name_spice,
                et,
                primary_body_name,
                spice_rhs_memo_enabled,
                spice_rhs_memo,
                counter,
            )
        end
        push!(positions_ii, pos_primary_body_j2000_m)
    end
    return _nbody_acceleration_from_positions(model, pos_ii, positions_ii)
end
