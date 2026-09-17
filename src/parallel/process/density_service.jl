# Distributed density service: a pool of worker processes that each own one
# native GRAM instance and answer batched density queries for a single
# propagation loop.
#
# Why processes rather than more threads. Native GRAM serializes every call on a
# process-wide lock that GRAMSuite refuses to release while `Threads.nthreads() > 1`
# (`gram_use_global_lock`, which returns `true` regardless once threaded). So N
# threads sharing one instance convoy, and N *instances* on N threads still
# convoy, because the lock is below the instance. A worker process started with
# `--threads=1` takes that lock uncontended, and separate processes have separate
# copies of it -- along with separate CSPICE global state, which is the other
# thing that makes concurrent GRAM unsafe in one address space.
#
# Why this only works batched. A native GRAM point call costs ~245 us (measured:
# `atmo256_gram_live_10min`, 37.67 s serial over ~154k calls). A Distributed
# round trip is the same order, so dispatching one call per satellite per RHS
# evaluation would trade a lock for a socket and gain nothing. The service is
# therefore built around one round trip per worker per *batch*: at N=256 split
# over 12 workers that is 12 round trips carrying ~21 queries each, and the
# ~8 KB of query data amortizes against ~5 ms of GRAM work per worker.
#
# Consequently this is only reachable from a call site that has every satellite's
# state in hand at once -- the density callback with
# `SPACEAGORA_DENSITY_FREEZE_PER_STEP=1`, not the per-satellite RHS path (see
# `SimulationCallbacks._gram_process_pool_batch_eval!`).

const _DENSITY_PROCESS_POOL = ProcessPool(Base.active_project())

"""
    density_process_pool() -> ProcessPool

The process-global pool backing the distributed density service. Kept separate
from [`campaign_process_pool`](@ref) on purpose: that pool's workers each run a
whole simulation, while these answer density queries for one propagation loop,
and a run can legitimately want both. Sharing one pool would let outer campaign
work and inner density work evict each other.
"""
density_process_pool()::ProcessPool = _DENSITY_PROCESS_POOL

# Filled in by SimulationCallbacks once SimulationModel is loaded. ParallelProcess
# is included *before* SimulationModel in src/SpaceAGORA.jl, so it cannot name
# EnvironmentModels at load time; this is the same Ref-injection pattern
# density_models.jl uses for the GRAMSuite extension.
const _DENSITY_SERVICE_BUILD_MODEL_FN = Ref{Function}(
    (_planet) -> error("Density service model constructor not installed.")
)
const _DENSITY_SERVICE_EVAL_FN = Ref{Function}(
    (_model, _h, _lat, _lon, _el, _wind, _vacT) -> error("Density service evaluator not installed.")
)

# Worker-side state. One instance per worker process, created on first use and
# reused for the life of the worker: construction is the expensive part (it loads
# MERRA2 data and takes GRAM's one-time init branch through CSPICE), and it is
# also the only part that is unsafe to do concurrently.
const _WORKER_DENSITY_MODEL = Ref{Any}(nothing)
const _WORKER_DENSITY_PLANET = Ref{String}("")

"""
    install_worker_density_model!(planet_name) -> Bool

Create this worker's GRAM instance for `planet_name` if it does not already have
one, and drive a single throwaway evaluation through it.

The warm-up is not an optimisation. A freshly constructed GRAM atmosphere takes a
one-time initialisation branch on its first `update()` that reaches CSPICE, and
CSPICE shares a global call-trace stack; taking that branch for the first time
inside a batch would put it in the same window as other work. Here it runs once,
alone, before the instance is ever published for queries.
"""
function install_worker_density_model!(planet_name::AbstractString)::Bool
    planet = String(planet_name)
    if _WORKER_DENSITY_MODEL[] !== nothing && _WORKER_DENSITY_PLANET[] == planet
        return true
    end
    model = _DENSITY_SERVICE_BUILD_MODEL_FN[](planet)
    try
        _DENSITY_SERVICE_EVAL_FN[](model, 1.0e5, 0.0, 0.0, 0.0, true, 200.0)
    catch
        # Non-fatal: a real failure here recurs on the first real query, where it
        # is reported with the caller's context instead of as a bare warm-up error.
    end
    _WORKER_DENSITY_MODEL[] = model
    _WORKER_DENSITY_PLANET[] = planet
    return true
end

"""
    density_batch_remote(hs, lats, lons, els, wind) -> (rhos, Ts, winds)

Evaluate a batch of density queries on this worker's own GRAM instance.

Runs *on the worker*. Takes and returns plain arrays only -- never the solver
parameter object, which is large, full of live runtime state, and in general not
serialisable. The whole payload for a 256-satellite batch is ~8 KB out and
~10 KB back.
"""
function density_batch_remote(
    hs::Vector{Float64},
    lats::Vector{Float64},
    lons::Vector{Float64},
    els::Vector{Float64},
    wind::Bool,
    vacuum_temperature::Float64,
)
    model = _WORKER_DENSITY_MODEL[]
    model === nothing && error(
        "Density service worker has no GRAM instance; install_worker_density_model! " *
        "was not called on this worker before dispatch."
    )
    n = length(hs)
    rhos = Vector{Float64}(undef, n)
    Ts = Vector{Float64}(undef, n)
    wind_x = Vector{Float64}(undef, n)
    wind_y = Vector{Float64}(undef, n)
    wind_z = Vector{Float64}(undef, n)
    eval_fn = _DENSITY_SERVICE_EVAL_FN[]
    @inbounds for i in 1:n
        rho, T, w = eval_fn(model, hs[i], lats[i], lons[i], els[i], wind, vacuum_temperature)
        rhos[i] = rho
        Ts[i] = T
        wind_x[i] = w[1]
        wind_y[i] = w[2]
        wind_z[i] = w[3]
    end
    # Wind is returned as three flat vectors rather than a Vector{SVector}: the
    # coordinator may hold a different StaticArrays version resolution than the
    # worker in a mixed environment, and flat Float64 vectors serialize without
    # depending on any of that.
    return (rhos, Ts, wind_x, wind_y, wind_z)
end

"""
    ensure_density_workers!(n; planet_name) -> Vector{Int}

Grow the density pool to `n` workers, each carrying its own GRAM instance for
`planet_name`, and return their ids.

Worker startup is expensive (SpaceAGORA + GRAMSuite + SPICE kernels, then the
GRAM instance itself: tens of seconds each), so this is worth calling once at
solve setup and never inside a timed region. Workers that fail to build an
instance are dropped from the returned list rather than left to fail on their
first query -- a short pool is degraded, but a pool with a dead member is wrong.
"""
function ensure_density_workers!(n::Int; planet_name::AbstractString)::Vector{Int}
    n <= 0 && return Int[]
    _ensure_density_atexit!()
    workers_all = ensure_process_workers!(_DENSITY_PROCESS_POOL, n)
    planet = String(planet_name)
    ready = Int[]
    for w in workers_all
        ok = try
            remotecall_fetch(install_worker_density_model!, w, planet)
        catch err
            @warn "Density service worker could not build its GRAM instance; excluding it from the pool." worker=w exception=(err, catch_backtrace())
            false
        end
        ok && push!(ready, w)
    end
    return ready
end

"""
    shutdown_density_workers!() -> Nothing

Remove every worker in the density pool.

Registered with `atexit` the first time the pool is grown, because without it a
coordinator that exits -- normally, or killed part-way through a benchmark --
leaves its workers running and reparented to init. Each one holds a full GRAM
instance and its atmosphere tables, measured at roughly 1.5 GB resident, so a
handful of runs is enough to exhaust a workstation. This was observed: twelve
orphans holding 18.3 GB survived a session of interrupted runs.

Safe to call more than once; a pool that is already empty is a no-op.
"""
function shutdown_density_workers!()::Nothing
    try
        shutdown_process_pool!(_DENSITY_PROCESS_POOL)
    catch
        # Teardown runs at process exit, where Distributed may already be tearing
        # down its own transport; a failure here must not mask the real exit.
    end
    return nothing
end

const _DENSITY_ATEXIT_REGISTERED = Ref(false)

@inline function _ensure_density_atexit!()
    _DENSITY_ATEXIT_REGISTERED[] && return nothing
    _DENSITY_ATEXIT_REGISTERED[] = true
    atexit(shutdown_density_workers!)
    return nothing
end

"""
    density_service_partition(n_items, n_workers) -> Vector{UnitRange{Int}}

Split `n_items` queries into at most `n_workers` contiguous ranges.

Contiguous and index-stable, not round-robin, and that is the point: each GRAM
instance carries per-track state, and the coordinator's own per-satellite caches
are keyed by satellite index. A satellite must therefore land on the same worker
on every step of a run, which a deterministic contiguous split gives for free as
long as the batch size is constant.
"""
function density_service_partition(n_items::Int, n_workers::Int)::Vector{UnitRange{Int}}
    (n_items <= 0 || n_workers <= 0) && return UnitRange{Int}[]
    w = min(n_workers, n_items)
    base, extra = divrem(n_items, w)
    ranges = Vector{UnitRange{Int}}(undef, w)
    start = 1
    @inbounds for i in 1:w
        len = base + (i <= extra ? 1 : 0)
        ranges[i] = start:(start + len - 1)
        start += len
    end
    return ranges
end

"""
    density_service_dispatch(workers, ranges, hs, lats, lons, els, wind)
        -> Union{Nothing, Vector{Tuple}}

Scatter one batch across `workers` and gather the per-range results, or `nothing`
if any worker failed.

Lives here rather than in the callback layer so `Distributed` stays confined to
`src/parallel/process/`, which owns process concerns; the caller works in plain
arrays and never imports Distributed.

All-or-nothing on purpose. A partial result would leave some satellites with
stale density and produce a silently wrong trajectory, which is worse than
falling back to the in-process path and being slow.
"""
function density_service_dispatch(
    workers::Vector{Int},
    ranges::Vector{UnitRange{Int}},
    hs::Vector{Float64},
    lats::Vector{Float64},
    lons::Vector{Float64},
    els::Vector{Float64},
    wind::Bool,
    vacuum_temperature::Float64,
)
    n_r = length(ranges)
    (n_r == 0 || length(workers) < n_r) && return nothing
    results = Vector{Union{Nothing, Tuple}}(nothing, n_r)
    failed = fill(false, n_r)
    @sync for k in 1:n_r
        rng = ranges[k]
        w = workers[k]
        Threads.@spawn begin
            try
                results[k] = remotecall_fetch(
                    density_batch_remote, w,
                    hs[rng], lats[rng], lons[rng], els[rng], wind, vacuum_temperature,
                )
            catch err
                @warn "Density service worker failed on a batch; falling back to the in-process path." worker=w exception=(err, catch_backtrace())
                failed[k] = true
            end
        end
    end
    any(failed) && return nothing
    any(r -> r === nothing, results) && return nothing
    return Tuple[r for r in results]
end
