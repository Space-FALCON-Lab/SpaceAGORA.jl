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

    completed = MonteCarloSampleResult[s for s in samples if s !== nothing]
    sort!(completed; by=s -> s.index)
    if spec.fail_fast
        _throw_first_monte_carlo_failure(completed)
    end
    return completed
end

"""
    run_monte_carlo(f, seeds; threads=1, fail_fast=false) -> MonteCarloResult
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

function run_monte_carlo(f, seeds; threads::Integer=1, fail_fast::Bool=false)
    return run_monte_carlo(f, MonteCarloSpec(seeds=seeds, threads=threads, fail_fast=fail_fast))
end
