"""
    MonteCarloSpec(; seeds, threads=1, fail_fast=false)

Configuration for [`run_monte_carlo`](@ref).

`seeds` is any finite iterable of per-sample identifiers. `threads` caps the
number of Monte Carlo worker tasks used by the runner; Julia must already have
been started with at least that many threads, for example
`julia --threads=8 --project=. script.jl`. `fail_fast=true` rethrows the first
sample failure after already-running threaded work reaches a scheduling point.
"""
struct MonteCarloSpec{S}
    seeds::S
    threads::Int
    fail_fast::Bool
end

function MonteCarloSpec(; seeds, threads::Integer=1, fail_fast::Bool=false)
    threads > 0 || throw(ArgumentError("MonteCarloSpec.threads must be > 0; got $(threads)."))
    return MonteCarloSpec(seeds, Int(threads), fail_fast)
end

"""
    MonteCarloSampleResult

Result for one Monte Carlo sample.

Fields include `index`, `seed`, `success`, `elapsed_s`, `value`, `error`, and
`backtrace`. Successful samples store the user function's return value in
`value`; failed samples store the thrown exception and captured stack trace.
"""
Base.@kwdef struct MonteCarloSampleResult
    index::Int
    seed::Any
    success::Bool
    elapsed_s::Float64
    value::Any = nothing
    error::Any = nothing
    backtrace::Any = nothing
end

"""
    MonteCarloResult

Aggregate result returned by [`run_monte_carlo`](@ref).

`samples` preserves seed order for all samples that ran. `successful` and
`failed` are convenience subsets of `samples`. `elapsed_s` is the total campaign
wall time, and `threads` is the worker-task count used by the run.
"""
struct MonteCarloResult
    samples::Vector{MonteCarloSampleResult}
    successful::Vector{MonteCarloSampleResult}
    failed::Vector{MonteCarloSampleResult}
    elapsed_s::Float64
    threads::Int
end

function MonteCarloResult(samples::Vector{MonteCarloSampleResult}, elapsed_s::Real, threads::Integer)
    successful = MonteCarloSampleResult[s for s in samples if s.success]
    failed = MonteCarloSampleResult[s for s in samples if !s.success]
    return MonteCarloResult(samples, successful, failed, Float64(elapsed_s), Int(threads))
end

function _validate_monte_carlo_threads(threads::Int)
    available = Base.Threads.nthreads()
    if threads > available
        throw(ArgumentError(
            "Monte Carlo requested threads=$(threads), but Julia was started with only $(available) thread(s). " *
            "Restart Julia with `--threads=$(threads)` or request fewer Monte Carlo threads."
        ))
    end
    return threads
end

function _run_monte_carlo_sample(f, index::Int, seed)
    start_ns = time_ns()
    try
        value = f(seed)
        elapsed_s = (time_ns() - start_ns) / 1.0e9
        return MonteCarloSampleResult(
            index=index,
            seed=seed,
            success=true,
            elapsed_s=elapsed_s,
            value=value
        )
    catch err
        bt = stacktrace(catch_backtrace())
        elapsed_s = (time_ns() - start_ns) / 1.0e9
        return MonteCarloSampleResult(
            index=index,
            seed=seed,
            success=false,
            elapsed_s=elapsed_s,
            error=err,
            backtrace=bt
        )
    end
end

function _throw_first_monte_carlo_failure(samples::Vector{MonteCarloSampleResult})
    idx = findfirst(s -> !s.success, samples)
    idx === nothing && return nothing
    sample = samples[idx]
    throw(ErrorException("Monte Carlo sample $(sample.index) with seed $(sample.seed) failed: $(sample.error)"))
end

function _run_monte_carlo_serial(f, seeds::Vector, spec::MonteCarloSpec)
    samples = MonteCarloSampleResult[]
    for (index, seed) in enumerate(seeds)
        sample = _run_monte_carlo_sample(f, index, seed)
        push!(samples, sample)
        if spec.fail_fast && !sample.success
            _throw_first_monte_carlo_failure(samples)
        end
    end
    return samples
end

function _run_monte_carlo_threaded(f, seeds::Vector, spec::MonteCarloSpec, worker_count::Int)
    jobs = Channel{Tuple{Int, Any}}(length(seeds))
    for (index, seed) in enumerate(seeds)
        put!(jobs, (index, seed))
    end
    close(jobs)

    samples = Vector{Union{Nothing, MonteCarloSampleResult}}(nothing, length(seeds))
    stop_requested = Base.Threads.Atomic{Bool}(false)

    Base.@sync begin
        for _ in 1:worker_count
            Base.Threads.@spawn begin
                for (index, seed) in jobs
                    spec.fail_fast && stop_requested[] && break
                    sample = _run_monte_carlo_sample(f, index, seed)
                    samples[index] = sample
                    if spec.fail_fast && !sample.success
                        Base.Threads.atomic_xchg!(stop_requested, true)
                        break
                    end
                end
            end
        end
    end

    # samples[index] writes plus the in-order filter above already leave
    # `completed` sorted by index; no re-sort needed.
    completed = MonteCarloSampleResult[s for s in samples if s !== nothing]
    if spec.fail_fast
        _throw_first_monte_carlo_failure(completed)
    end
    return completed
end

# Process-backend counterpart of `_run_monte_carlo_threaded`: identical
# job-queue structure, but each sample runs via `remotecall_fetch` on one of
# `worker_ids` instead of a local `Threads.@spawn` task. The remote callable is
# a closure over `f` dispatched through a `CachingPool`, so the (potentially
# large, configuration-capturing) closure is serialized to each worker once and
# referenced by identity afterwards, instead of being re-serialized with every
# sample. The dispatch loop itself uses `@async`/`@sync` (not `Threads.@spawn`):
# each task just blocks on IPC waiting for a worker's reply, so it should not
# occupy an OS thread the way genuinely CPU-bound work would.
function _run_monte_carlo_process(f, seeds::Vector, spec::MonteCarloSpec, worker_ids::Vector{Int})
    jobs = Channel{Tuple{Int, Any}}(length(seeds))
    for (index, seed) in enumerate(seeds)
        put!(jobs, (index, seed))
    end
    close(jobs)

    samples = Vector{Union{Nothing, MonteCarloSampleResult}}(nothing, length(seeds))
    stop_requested = Base.Threads.Atomic{Bool}(false)

    pool = CachingPool(worker_ids)
    run_sample = (index, seed) -> _run_monte_carlo_sample(f, index, seed)
    try
        Base.@sync begin
            for _ in worker_ids
                Base.@async begin
                    for (index, seed) in jobs
                        spec.fail_fast && stop_requested[] && break
                        sample = remotecall_fetch(run_sample, pool, index, seed)
                        samples[index] = sample
                        if spec.fail_fast && !sample.success
                            Base.Threads.atomic_xchg!(stop_requested, true)
                            break
                        end
                    end
                end
            end
        end
    finally
        # Drop the cached closure on the workers; it can capture large
        # configuration state that should not outlive the campaign.
        Distributed.clear!(pool)
    end

    # samples[index] writes plus the in-order filter above already leave
    # `completed` sorted by index; no re-sort needed.
    completed = MonteCarloSampleResult[s for s in samples if s !== nothing]
    if spec.fail_fast
        _throw_first_monte_carlo_failure(completed)
    end
    return completed
end

"""
    run_monte_carlo(f, seeds; threads=1, fail_fast=false,
                    route_features=nothing, route_state=nothing,
                    route_tuning=nothing) -> MonteCarloResult
    run_monte_carlo(f, spec::MonteCarloSpec) -> MonteCarloResult

Run `f(seed)` for each Monte Carlo seed and return ordered sample results.

`threads` controls the number of outer Monte Carlo worker tasks, not the number
of Julia threads created at runtime. Start Julia with enough threads before
calling this function:

```julia
result = run_monte_carlo(1:100; threads=8) do seed
    args = make_config_for_seed(seed)
    run_simulation(args; return_solution=true)
end
```

The runner records per-sample exceptions instead of throwing by default, so a
long campaign can finish and report all failures. Set `fail_fast=true` to throw
on the first failed sample.

# Adaptive mode

Pass `threads=:auto` to let the outer-route bandit pick serial or threaded
execution from empirical runtime history instead of a fixed worker count. The
runner builds [`OuterRouteFeatures`](@ref) from the campaign shape (pass
`route_features` from [`campaign_route_features`](@ref) to describe per-sample
satellite count, density-model family, and mission length; the sample count is
always filled in from `seeds`), consults [`select_outer_route!`](@ref), and
after the campaign records per-sample success and amortized wall-clock feedback
via [`record_outer_route_feedback!`](@ref), so repeated campaigns with the same
shape converge to the fastest allocation. History lives in
[`campaign_outer_route_state`](@ref) unless an isolated
[`OuterRouteState`](@ref) is passed as `route_state`; `route_tuning` overrides
the [`OuterRouteTuning`](@ref).

While adaptive threaded workers are active the runner sets
`SPACEAGORA_OUTER_PARALLEL_ACTIVE=1` and, unless one is already set, an
`SPACEAGORA_INNER_THREAD_BUDGET` of `nthreads() ÷ workers` so inner and outer
parallelism split the thread pool instead of oversubscribing it. Conversely,
when `SPACEAGORA_OUTER_PARALLEL_ACTIVE` is already set — a nested adaptive
campaign inside another campaign's worker — the runner yields to the enclosing
split: it executes serially and records no feedback, so contended timings never
poison the shared route statistics.
"""
function run_monte_carlo(f, spec::MonteCarloSpec)
    seeds = collect(spec.seeds)
    isempty(seeds) && return MonteCarloResult(MonteCarloSampleResult[], 0.0, 0)

    requested_threads = _validate_monte_carlo_threads(spec.threads)
    worker_count = min(requested_threads, length(seeds))

    start_ns = time_ns()
    samples = if worker_count == 1
        _run_monte_carlo_serial(f, seeds, spec)
    else
        _run_monte_carlo_threaded(f, seeds, spec, worker_count)
    end
    elapsed_s = (time_ns() - start_ns) / 1.0e9
    return MonteCarloResult(samples, elapsed_s, worker_count)
end

function run_monte_carlo(
    f,
    seeds;
    threads::Union{Integer, Symbol}=1,
    fail_fast::Bool=false,
    route_features::Union{Nothing, OuterRouteFeatures}=nothing,
    route_state::Union{Nothing, OuterRouteState}=nothing,
    route_tuning::Union{Nothing, OuterRouteTuning}=nothing
)
    if threads isa Symbol
        threads === :auto || throw(ArgumentError(
            "run_monte_carlo threads must be a positive integer or :auto; got :$(threads)."
        ))
        return _run_monte_carlo_adaptive(
            f,
            seeds;
            fail_fast=fail_fast,
            route_features=route_features,
            route_state=route_state,
            route_tuning=route_tuning
        )
    end
    if !(route_features === nothing && route_state === nothing && route_tuning === nothing)
        throw(ArgumentError(
            "run_monte_carlo route_features/route_state/route_tuning are only consulted with threads=:auto."
        ))
    end
    return run_monte_carlo(f, MonteCarloSpec(seeds=seeds, threads=threads, fail_fast=fail_fast))
end
