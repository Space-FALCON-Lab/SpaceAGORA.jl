# Coordinator side of the distributed density service (see
# src/parallel/process/density_service.jl for the pool and the worker side).
#
# This file owns two things: installing the Ref hooks ParallelProcess needs to
# build and evaluate a GRAM model without naming EnvironmentModels (it is loaded
# before SimulationModel), and the batched entry point that mirrors
# `_gram_isolated_pool_batch_eval!` so both can sit behind one call site.

# SimulationCallbacks is a nested module with no `using` of its own, and nested
# modules do not inherit their parent's bindings, so ParallelProcess has to be
# named explicitly. Same parentmodule walk environment/ephemerides/planets.jl
# uses for SPICE_LOCK. Safe at include time because src/SpaceAGORA.jl includes
# the parallel routing/process layer before SimulationModel.
# Resolved tolerantly rather than asserted at load time.
#
# The obvious `parentmodule(parentmodule(@__MODULE__)).ParallelProcess` breaks
# every test that includes src/core/simulation_model.jl straight into Main
# (test/probes/*, test/gates/ci_architecture_contract_gate.jl): there the parent
# chain ends at Main, which has no ParallelProcess, and the file fails to load at
# all. The density service is optional, so an unavailable module must disable the
# service, not break the module that merely mentions it.
const ParallelProcess = let root = parentmodule(parentmodule(@__MODULE__))
    isdefined(root, :ParallelProcess) ? getfield(root, :ParallelProcess) : nothing
end

"""
    _install_density_service_hooks!()

Wire the model constructor and evaluator into `ParallelProcess`.

`ParallelProcess` is included ahead of `SimulationModel` in `src/SpaceAGORA.jl`,
so it cannot reference `GRAMAtmosphereModel` or `_gram_core_density_state`
directly. Same Ref-injection pattern `density_models.jl` uses to reach the
GRAMSuite extension.
"""
function _install_density_service_hooks!()
    ParallelProcess === nothing && return nothing
    ParallelProcess._DENSITY_SERVICE_BUILD_MODEL_FN[] =
        (planet) -> EnvironmentModels.GRAMAtmosphereModel(planet_name=String(planet))
    ParallelProcess._DENSITY_SERVICE_EVAL_FN[] =
        (model, h, lat, lon, el, wind, vacT) -> EnvironmentModels._gram_core_density_state(
            model.core, Float64(h), Float64(lat), Float64(lon), Float64(el),
            wind, model.instance_lock, Float64(vacT),
        )
    return nothing
end

"""
    _gram_process_pool_mode() -> Symbol

`off` (default), `on`, or `auto`. `auto` engages only above the batch-size
threshold, since below it the round-trip cost is not amortised.
"""
@inline function _gram_process_pool_mode()::Symbol
    return ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_GRAM_PROCESS_POOL"; default="off")
end

@inline function _gram_process_pool_workers()::Int
    raw = strip(get(ENV, "SPACEAGORA_GRAM_PROCESS_POOL_WORKERS", ""))
    isempty(raw) && return max(1, Sys.CPU_THREADS ÷ 2)
    parsed = tryparse(Int, raw)
    return parsed === nothing ? max(1, Sys.CPU_THREADS ÷ 2) : max(1, parsed)
end

# Below this many satellites in a batch, the per-round-trip cost is a larger
# share of the batch than the GRAM work it displaces. 64 is a placeholder pending
# the measurement described in the file header of density_service.jl; it is
# deliberately conservative, because engaging too eagerly is a regression while
# engaging too late merely leaves the old behaviour in place.
@inline function _gram_process_pool_threshold()::Int
    raw = strip(get(ENV, "SPACEAGORA_GRAM_PROCESS_POOL_THRESHOLD", ""))
    isempty(raw) && return 64
    parsed = tryparse(Int, raw)
    return parsed === nothing ? 64 : max(1, parsed)
end

@inline function _gram_process_pool_enabled(num_items::Int)::Bool
    # No ParallelProcess in this module graph (standalone include): the service
    # cannot run, so it declines rather than erroring.
    ParallelProcess === nothing && return false
    mode = _gram_process_pool_mode()
    mode === :off && return false
    # Never nest inside an outer process split: those workers already own the
    # cores this pool would spawn onto, and the two pools would oversubscribe
    # each other. Same arbitration SPACEAGORA_OUTER_PARALLEL_ACTIVE performs for
    # the inner thread budget.
    ParallelPolicy.outer_parallel_active() && return false
    mode === :on && return true
    return num_items >= _gram_process_pool_threshold()
end

"""
    _gram_process_pool_batch_eval!(rhos, Ts, winds, density_model,
                                   hs, lats, lons, el_time, wind, p) -> Bool

Answer a whole batch of density queries from the distributed density service.

Returns `true` if the batch was served and the output arrays are filled, `false`
if the service declined (disabled, below threshold, nested inside an outer
process split, no live workers, or a non-native-GRAM model), in which case the
caller must fall through to its existing path. Never throws on a service-level
failure: a density service that cannot answer must degrade to the in-process
path, not fail the solve.

Signature deliberately mirrors `_gram_isolated_pool_batch_eval!` so the two can
sit behind a single call site.
"""
function _gram_process_pool_batch_eval!(
    rhos::AbstractVector{Float64},
    Ts::AbstractVector{Float64},
    winds::AbstractVector{SVector{3, Float64}},
    density_model,
    hs::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    lons::AbstractVector{<:Real},
    el_time::Union{Float64, AbstractVector{<:Real}},
    wind::Bool,
    p,
)::Bool
    density_model isa EnvironmentModels.GRAMAtmosphereModel || return false
    n = length(hs)
    _gram_process_pool_enabled(n) || return false

    planet = _density_service_planet_name(density_model)
    isempty(planet) && return false

    workers = try
        ParallelProcess.ensure_density_workers!(
            _gram_process_pool_workers(); planet_name=planet
        )
    catch err
        @warn "Density service unavailable; falling back to the in-process path." exception=(err, catch_backtrace())
        return false
    end
    isempty(workers) && return false

    hs_f = _density_service_f64(hs, n)
    lats_f = _density_service_f64(lats, n)
    lons_f = _density_service_f64(lons, n)
    els_f = el_time isa Float64 ? fill(el_time, n) : _density_service_f64(el_time, n)

    ranges = ParallelProcess.density_service_partition(n, length(workers))
    isempty(ranges) && return false

    results = ParallelProcess.density_service_dispatch(
        workers, ranges, hs_f, lats_f, lons_f, els_f, wind,
        p.args.environment_model.planet.T_ref,
    )
    results === nothing && return false

    @inbounds for (k, rng) in enumerate(ranges)
        (rho_k, T_k, wx, wy, wz) = results[k]
        for (j, i) in enumerate(rng)
            rhos[i] = rho_k[j]
            Ts[i] = T_k[j]
            winds[i] = SVector{3, Float64}(wx[j], wy[j], wz[j])
        end
    end
    return true
end

@inline function _density_service_f64(v::AbstractVector{<:Real}, n::Int)::Vector{Float64}
    out = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        out[i] = Float64(v[i])
    end
    return out
end

# The worker needs a planet name, not the model: a GRAMAtmosphereModel wraps a
# live native handle that cannot be serialised across a process boundary. Each
# worker builds its own from the name instead.
@inline function _density_service_planet_name(model)::String
    for field in (:planet_name, :planet)
        if hasproperty(model, field)
            value = getproperty(model, field)
            value isa AbstractString && return lowercase(strip(String(value)))
        end
    end
    return ""
end

"""
    _rhs_density_service_candidate(p, num_sats) -> Bool

Whether the RHS atmosphere prefill may be served by the distributed density
service for this solve.

Checked once per prefill pass, before the planet-frame loop, so the decision is
uniform across satellites within a derivative evaluation.

The track-cache guard matters: `_density_state_from_kinematics!` keeps
per-satellite interpolation state when the GRAM track cache is on, and a batched
query has no way to advance it. Track cache defaults to `off`, so the common
native-GRAM configuration is eligible.
"""
@inline function _rhs_density_service_candidate(p, num_sats::Int)::Bool
    num_sats > 0 || return false
    _gram_process_pool_enabled(num_sats) || return false
    model = p.args.environment_model.density_model
    model isa EnvironmentModels.GRAMAtmosphereModel || return false
    cb_env = _callback_env_config(p)
    _gram_track_cache_enabled(cb_env.gram_track_cache, model) && return false
    # The vacuum-predicted cache keeps a per-satellite spline that a batched
    # query cannot advance, same reason the track cache disqualifies a run.
    cb_env.vacuum_gram_cache_enabled && return false
    # Per-satellite density model instances mean per-satellite native handles;
    # the service builds one instance per worker from a planet name and cannot
    # reproduce a heterogeneous set.
    isempty(p.shared_buffers.density_models) || return false
    return true
end

"""
    _rhs_density_service_fill!(p, t, num_sats, alts, lats, lons) -> Bool

Fill the shared density/temperature/wind buffers for every active satellite from
the distributed density service, using planet-frame values the caller has
already computed. Returns `false` if the service declined or failed, leaving the
buffers untouched so the caller can fall back.

**This is exact, not an approximation.** The queries carry the same
`(alt, lat, lon, t)` the per-satellite path would have passed for this same
derivative evaluation, so the density each satellite sees is identical; only
*where* the native call runs changes. That is the whole reason this hook sits in
the RHS prefill rather than behind `density_freeze_per_step`, which reuses a
once-per-accepted-step sample and does change the answer.

Inactive satellites are excluded from the dispatch rather than sent with stale
buffer contents: their planet frames were never computed this pass, so their
`alt/lat/lon` slots hold whatever the previous step left, and feeding those to
native GRAM invites a domain error for a result nothing will read.
"""
function _rhs_density_service_fill!(
    p, t::Float64, num_sats::Int,
    alts::Vector{Float64}, lats::Vector{Float64}, lons::Vector{Float64},
)::Bool
    # Replicates _gram_isolated_pool_density_state's per-satellite gate exactly.
    # That gate decides, before any native call, whether a satellite gets a
    # native GRAM query at all: above 2000 km the model returns vacuum, and above
    # the entry interface on a non-keplerian run it returns the analytic
    # density_polyfit instead. Dispatching those to the service would not merely
    # waste a round trip -- it would return *different numbers* than the
    # in-process path, which is the difference between a speedup and a wrong
    # answer. Only the satellites that would really have called GRAM are sent.
    planet = p.args.environment_model.planet
    T_ref = planet.T_ref
    EI = p.args.environment_model.EI * 1e3
    keplerian = p.args.mission_configuration.keplerian
    zero_wind = SVector{3, Float64}(0.0, 0.0, 0.0)

    gram_idx = Int[]
    sizehint!(gram_idx, num_sats)
    @inbounds for i in 1:num_sats
        p.is_active[i] || continue
        h = alts[i]
        isfinite(h) || return false
        if h > 2000.0e3
            _write_density_buffers!(p, i, 0.0, T_ref, zero_wind, t)
        elseif (h - EI > 0.0) && !keplerian
            rho, T, w = EnvironmentModels.density_polyfit(h, p)
            _write_density_buffers!(p, i, rho, T, w, t)
        else
            push!(gram_idx, i)
        end
    end
    isempty(gram_idx) && return true

    n = length(gram_idx)
    q_alt = Vector{Float64}(undef, n)
    q_lat = Vector{Float64}(undef, n)
    q_lon = Vector{Float64}(undef, n)
    @inbounds for (j, i) in enumerate(gram_idx)
        q_alt[j] = alts[i]
        q_lat[j] = lats[i]
        q_lon[j] = lons[i]
    end

    rhos = Vector{Float64}(undef, n)
    Ts = Vector{Float64}(undef, n)
    winds = Vector{SVector{3, Float64}}(undef, n)
    served = _gram_process_pool_batch_eval!(
        rhos, Ts, winds,
        p.args.environment_model.density_model,
        q_alt, q_lat, q_lon, t, true, p,
    )
    # Partial credit is not available: satellites handled locally above already
    # have their buffers written, but those writes are identical to what the
    # fallback would produce, so re-running the per-satellite pass over all of
    # them is correct, just redundant for a few.
    served || return false

    @inbounds for (j, i) in enumerate(gram_idx)
        _write_density_buffers!(p, i, rhos[j], Ts[j], winds[j], t)
    end
    return true
end
