"""
    ProcessPool

A cached set of `Distributed` worker processes bootstrapped with `SpaceAGORA`
(and, best-effort, `GRAMSuite`, so the GRAM package extension is active on
every worker for campaigns that need it), the default SPICE kernel set
(leapseconds + DE44x planetary ephemerides), and a warmed-up GRAM native
wrapper module (avoids a world-age `MethodError` on the first GRAM-density
sample dispatched to a fresh worker), so campaigns using the standard
vendored kernels and a GRAM density model don't need any extra per-worker
setup of their own. A campaign on a non-default kernel directory still needs
to furnish those kernels on the pool's workers itself.

Each worker is started with `--threads=1`: process workers do not share the
coordinator's Julia thread pool, so inner thread-based parallelism inside a
worker's own `run_simulation` call is unaffected by (and does not contend
with) how many process workers are active.
"""
mutable struct ProcessPool
    workers::Vector{Int}
    project_path::String
    lock::ReentrantLock
end

ProcessPool(project_path::AbstractString) = ProcessPool(Int[], String(project_path), ReentrantLock())

const _CAMPAIGN_PROCESS_POOL = ProcessPool(Base.active_project())

"""
    campaign_process_pool() -> ProcessPool

Return the process-global [`ProcessPool`](@ref) that campaign runners consult
and grow when routed to `:process`. Shared across repeated campaigns in the
same session so worker processes (and their one-time SpaceAGORA/GRAMSuite
precompilation cost) are reused rather than respawned per campaign.
"""
campaign_process_pool()::ProcessPool = _CAMPAIGN_PROCESS_POOL

@inline function _process_worker_exeflags(project_path::AbstractString)::Cmd
    return Cmd(["--threads=1", "--startup-file=no", "--project=$(project_path)"])
end

# `Distributed.remotecall_eval` (the plain function `@everywhere` itself
# expands to -- not the macro, which expands to a `:toplevel` Expr and is
# therefore only legal as a direct top-level statement, never inside a
# function body like this one) instead of a `remotecall_wait` closure: any
# closure defined inside the SpaceAGORA package itself has
# SpaceAGORA.ParallelProcess as its *home module*, and Distributed needs that
# module to already exist on the worker to deserialize the closure at all --
# a chicken-and-egg problem on a bare worker that hasn't loaded SpaceAGORA
# yet. `remotecall_eval` sidesteps this by shipping unevaluated source (an
# `Expr`) and `eval`ing it directly in the worker's Main, with no
# closure/module resolution step. Once SpaceAGORA itself is loaded this way,
# subsequent `remotecall_wait` closures defined inside the package work
# normally, since SpaceAGORA.ParallelProcess now exists on the worker to
# resolve against.
function _bootstrap_process_worker!(worker::Int, project_path::String)::Nothing
    Distributed.remotecall_eval(Main, [worker], :(using SpaceAGORA))
    # GRAMSuite is a weak/optional dependency (only in [weakdeps]/[extras]),
    # so plain `using GRAMSuite` doesn't resolve even on the coordinator
    # without first pushing its vendored path onto LOAD_PATH (see
    # examples/common.jl's `ensure_gramsuite_loaded!`, which this mirrors for
    # the worker side). `import GRAMSuite`, not `using GRAMSuite`: SpaceAGORA
    # and GRAMSuite both export a same-named `InitialTime` (different types),
    # so `using` both in the same Main makes any unqualified reference to it
    # ambiguous -- `import` only binds the module name itself, matching
    # `ensure_gramsuite_loaded!`'s own deliberate choice for the same reason.
    # Best-effort: a campaign that never touches a GRAM density model doesn't
    # need this to succeed.
    try
        repo_root = dirname(project_path)
        vendored_gramsuite = joinpath(repo_root, "data", "GRAMSuite.jl")
        remotecall_wait(worker, vendored_gramsuite) do path
            if Base.find_package("GRAMSuite") === nothing && isdir(path)
                pushfirst!(LOAD_PATH, path)
            end
            nothing
        end
        Distributed.remotecall_eval(Main, [worker], :(import GRAMSuite))
    catch err
        @warn "Process worker $(worker) could not load GRAMSuite; campaigns using a GRAM density model will fail on this worker." exception=(err, catch_backtrace())
    end
    _furnish_default_spice_kernels!(worker)
    return nothing
end

# SPICE kernel state (leapseconds, planetary ephemerides, orientation frames)
# is per-OS-process, not part of any value that crosses a Distributed
# boundary: a `SimulationConfiguration`'s `AbstractPlanet` deserializes fine
# as plain data (radius, mass, ...), but the *side effect* of furnishing
# kernels that happened on the coordinator when it was originally
# constructed never replays on the worker, since deserialization only
# rebuilds field values, not the original constructor call. Reconstructing a
# default Earth() on each worker furnishes the same shared kernel set
# (leapseconds + DE44x planetary ephemerides) that the vast majority of
# SpaceAGORA campaigns already rely on -- covering the common case
# automatically. A campaign built on a non-default SPICE kernel directory,
# or a non-Earth-primary mission whose required kernels aren't in that
# shared set, must still furnish its own kernels on the pool's workers
# (e.g. via `remotecall_wait` on `ensure_process_workers!`'s return value)
# before dispatching -- best-effort here, not a full solution.
function _furnish_default_spice_kernels!(worker::Int)::Nothing
    try
        # remotecall_eval, not a remotecall_wait closure referencing
        # `SpaceAGORA` by bare name: this function is itself defined inside
        # SpaceAGORA.ParallelProcess, a submodule that does not automatically
        # see its own parent's name (ordinary Julia scoping, unrelated to
        # Distributed) -- `..SpaceAGORA` isn't resolvable at ParallelProcess's
        # own include time either, since ParallelProcess is included before
        # SpaceAGORA finishes defining its other submodules. Shipping the
        # quoted expression and letting Main resolve `SpaceAGORA` on the
        # worker (already bound there by the `using SpaceAGORA` bootstrap
        # step above) sidesteps both problems at once.
        Distributed.remotecall_eval(Main, [worker], :(SpaceAGORA.SimulationModel.Earth("")))
    catch err
        @warn "Process worker $(worker) could not furnish default SPICE kernels; campaigns needing SPICE ephemerides (native GRAM, N-body, SRP, orientation) may fail on this worker unless the caller furnishes kernels itself." exception=(err, catch_backtrace())
    end
    _warm_gram_wrapper!(worker)
    return nothing
end

# GRAMSuite dynamically `Base.include`s its native wrapper module (GRAM.jl) on
# first-ever GRAMAtmosphereModel construction in a process (see
# _load_gram_wrapper! in data/GRAMSuite.jl), which defines new types/methods
# (e.g. GRAMSuite.GRAM.EphemerisStateC) at a *later* Julia world age than
# whatever code triggered the construction. If that first construction happens
# inside a single already-compiled closure `f` dispatched to this worker (the
# ordinary case: a campaign's per-sample closure builds its own density model),
# `f`'s own later code -- e.g. run_simulation's GRAM ephemeris-bypass density
# callback, which calls EphemerisStateC(...) without its own invokelatest --
# still executes under the pre-dynamic-load world age and throws
# `MethodError(..., world_age=...)`, even though GRAMAtmosphereModel's own
# constructor safely self-recurses via invokelatest for its own construction.
# Forcing this dynamic load here, in its own separate remote call (so it gets
# its own fresh world age before any real campaign closure ever runs), avoids
# the whole class of issue -- the same "serial warmup first" workaround this
# session's own GRAM-concurrency testing already relies on. Best-effort: a
# campaign that never uses a GRAM density model doesn't need this to succeed.
#
# Must exercise density SAMPLING too, not just construction: a bare
# GRAMAtmosphereModel() warmup alone was not sufficient in testing -- a later
# real campaign sample still hit the same world-age MethodError, this time
# from inside the GRAM ephemeris-bypass hook (`_gram_apply_user_ephemeris_state!`
# in ext/SpaceAGORAGRAMSuiteExt.jl, which calls `GRAMSuite.GRAM.EphemerisStateC`
# without its own invokelatest). That code path is only reached once a density
# *query* actually happens, not at model construction. Running a real
# `point_density_state` call here, in this same warmup remote call, forces
# every world-age-sensitive code path a later campaign sample will hit to
# settle before any real work is dispatched to this worker.
function _warm_gram_wrapper!(worker::Int)::Nothing
    try
        Distributed.remotecall_eval(Main, [worker], quote
            _spaceagora_parallel_process_gram_warmup_model = SpaceAGORA.SimulationModel.GRAMAtmosphereModel(planet_name="earth")
            GRAMSuite.point_density_state(
                _spaceagora_parallel_process_gram_warmup_model.core, 150.0e3, 0.0, 0.0, 0.0, false
            )
            nothing
        end)
    catch err
        @warn "Process worker $(worker) could not warm up the GRAM native wrapper; the first GRAM-density campaign sample dispatched here may fail with a world-age MethodError." exception=(err, catch_backtrace())
    end
    return nothing
end

"""
    ensure_process_workers!(pool::ProcessPool, n::Int; warmup_fn=nothing) -> Vector{Int}

Grow `pool` to at least `n` bootstrapped workers (spawning new `addprocs`
workers only for the shortfall) and return its current worker ids.

New workers are started against `pool.project_path` and bootstrapped with
`SpaceAGORA`/`GRAMSuite`/default SPICE kernels before being added to the pool,
so any worker id this function returns is immediately ready to accept
campaign work using the standard vendored kernel set.

`warmup_fn`, if given, is a zero-argument closure `remotecall_fetch`'d
(best-effort) on each *newly added* worker only, after the generic bootstrap
above and before the worker is added to the pool. This exists to close a gap
the generic bootstrap can't: `_warm_gram_wrapper!` above only exercises GRAM's
native wrapper module, not the much larger JIT/specialization surface of a
full `run_simulation` call, which is specialized per exact
`SimulationConfiguration` type combination -- something this generic,
caller-agnostic bootstrap has no way to predict. Measured directly: a worker's
first-ever full campaign-sample call can take ~70s cold vs. a fraction of a
second once warm, on an identical workload. A caller that knows its own
representative call (e.g. `() -> f(first(seeds))`, reusing the exact closure
it's about to dispatch for real) can pay that cost here, once per new worker,
instead of it silently landing inside the caller's own first real (measured)
dispatch. Everything `warmup_fn` references must already be resolvable on the
worker (ordinary Distributed closure-shipping rule -- see the campaign
dispatch closures this same pool already ships for real work).
"""
function ensure_process_workers!(pool::ProcessPool, n::Int; warmup_fn=nothing)::Vector{Int}
    desired = max(1, n)
    lock(pool.lock) do
        shortfall = desired - length(pool.workers)
        if shortfall > 0
            new_workers = addprocs(shortfall; exeflags=_process_worker_exeflags(pool.project_path))
            for w in new_workers
                _bootstrap_process_worker!(w, pool.project_path)
            end
            if warmup_fn !== nothing
                for w in new_workers
                    try
                        remotecall_fetch(warmup_fn, w)
                    catch err
                        @warn "Process worker $(w) warmup_fn failed; its first real dispatch may be slow." exception=(err, catch_backtrace())
                    end
                end
            end
            append!(pool.workers, new_workers)
        end
        return copy(pool.workers)
    end
end

"""
    shutdown_process_pool!(pool::ProcessPool) -> Nothing

Remove every worker currently in `pool` via `rmprocs` and clear it. Mainly
useful for tests; campaign code leaves the process-global pool warm across
calls by design.
"""
function shutdown_process_pool!(pool::ProcessPool)::Nothing
    lock(pool.lock) do
        isempty(pool.workers) && return nothing
        rmprocs(pool.workers)
        empty!(pool.workers)
        return nothing
    end
end
