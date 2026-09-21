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
    # Coordinator clock (`time_ns`) at which this sample's result was in hand,
    # stamped by the dispatcher that collected it (`_stamp_finished`). NaN when
    # the sample was constructed elsewhere and never collected -- on a process
    # worker, before the round trip home. The route bandit reads the spread of
    # these to credit a route with its steady per-sample cost rather than the
    # campaign's mean (see `steady_per_sample_s`).
    finished_ns::Float64 = NaN
end

@inline function _stamp_finished(s::MonteCarloSampleResult)::MonteCarloSampleResult
    return MonteCarloSampleResult(index=s.index, seed=s.seed, success=s.success,
                                  elapsed_s=s.elapsed_s, value=s.value,
                                  error=s.error, backtrace=s.backtrace,
                                  finished_ns=Float64(time_ns()))
end

"""
    MonteCarloResult

Aggregate result returned by [`run_monte_carlo`](@ref).

`samples` preserves seed order for all samples that ran. `successful` and
`failed` are convenience subsets of `samples`. `elapsed_s` is the total campaign
wall time, and `threads` is the number of samples the run kept in flight at
once. `route` is the outer route that ran them (`:none`, `:threads` or
`:process`), and `local_slots` how many of those in-flight samples the
coordinator ran on its own threads beside the process pool (mixed dispatch;
zero for every other route).
"""
struct MonteCarloResult
    samples::Vector{MonteCarloSampleResult}
    successful::Vector{MonteCarloSampleResult}
    failed::Vector{MonteCarloSampleResult}
    elapsed_s::Float64
    threads::Int
    route::Symbol
    local_slots::Int
end

function MonteCarloResult(samples::Vector{MonteCarloSampleResult}, elapsed_s::Real, threads::Integer;
                          route::Symbol = (threads > 1 ? :threads : :none), local_slots::Integer = 0)
    successful = MonteCarloSampleResult[s for s in samples if s.success]
    failed = MonteCarloSampleResult[s for s in samples if !s.success]
    return MonteCarloResult(samples, successful, failed, Float64(elapsed_s), Int(threads), route, Int(local_slots))
end

"""
    steady_per_sample_s(result::MonteCarloResult) -> Float64

The campaign's per-sample cost once it was running steadily: the wall between
the median completion and the last one, divided by the samples that completed
in that window. Falls back to `elapsed_s / n` when fewer than four samples
carry a completion stamp.

`elapsed_s / n` is what a campaign cost; this is what the next one will. The
difference is everything the first campaign of a route pays once -- the
process pool spinning up, a worker JIT-compiling the sample closure, the
coordinator compiling the dispatcher -- and it is not small: measured on the
TRX50's L12 (independent_1sat_1hr, 64 samples, 24 workers) the process route's
first campaign took 3.16 s against 0.20 s at steady state. Credited with the
mean, the route bandit rated the pool at 49 ms per sample against the threads
route's 16 ms, chose threads, and had no reason ever to re-try the pool it had
mis-measured; the static process route ran the same shape 3x faster.
"""
function steady_per_sample_s(result::MonteCarloResult)::Float64
    n = length(result.samples)
    n <= 0 && return 0.0
    mean_s = result.elapsed_s / n
    stamps = Float64[s.finished_ns for s in result.samples if isfinite(s.finished_ns)]
    length(stamps) >= 4 || return mean_s
    sort!(stamps)
    half = length(stamps) ÷ 2
    span_s = (stamps[end] - stamps[half]) / 1.0e9
    tail = length(stamps) - half
    (isfinite(span_s) && span_s > 0.0 && tail > 0) || return mean_s
    return span_s / tail
end

"""
    _warm_campaign_dispatchers()

Run every Monte Carlo dispatcher once on a trivial sample: the serial loop, the
mixed dispatcher on local slots, and -- given a second thread -- the threaded
dispatcher and the mixed dispatcher at width two, then the steady-cost
estimator over the result. This is the body of the package's precompile
workload (see `@compile_workload` in `SpaceAGORA.jl`), kept as a function so
it is one piece of code compiled into the pkgimage at precompile time and
exercised by the test suite at run time; a `@compile_workload` block's own
lines never execute in a test process.
"""
function _warm_campaign_dispatchers()::Nothing
    sample = seed -> seed * 2
    seeds = collect(1:4)
    spec1 = MonteCarloSpec(seeds = seeds, threads = 1)
    serial = _run_monte_carlo_serial(sample, seeds, spec1)
    _run_monte_carlo_mixed(sample, seeds, spec1, Int[], 1)
    if Base.Threads.nthreads() > 1
        spec2 = MonteCarloSpec(seeds = seeds, threads = 2)
        _run_monte_carlo_threaded(sample, seeds, spec2, 2)
        _run_monte_carlo_mixed(sample, seeds, spec2, Int[], 2)
    end
    steady_per_sample_s(MonteCarloResult(serial, 0.01, 1))
    return nothing
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
        sample = _stamp_finished(_run_monte_carlo_sample(f, index, seed))
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
                    sample = _stamp_finished(_run_monte_carlo_sample(f, index, seed))
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
"""
    DispatchAdmission(workers, locals)

Per-consumer admission for [`_run_monte_carlo_mixed`](@ref): one flag per pool
worker and one per coordinator local slot, each of which that consumer checks
BEFORE taking its next job.

This is how a campaign changes its own shape without a barrier. The alternative
-- stopping the dispatch, deciding, and dispatching the remainder -- costs a
whole extra dispatch, and a dispatch's fixed cost is paid per dispatch rather
than per sample: measured on this repo's workstation, a guard round of eleven
samples cost 0.24-0.36 s against 0.24-0.26 s for the remaining fifty-three.
Closing a consumer costs one atomic store.

Closing is one-way and never interrupts work in flight: a closed consumer
finishes the sample it holds, declines to take another, and drops out of the
dispatch. Closing every local slot is "reduce L to zero"; closing the pool
class is "the rest of this campaign runs on the coordinator's threads".
Reducing L to `k > 0` closes the `L - k` highest-numbered local slots, so which
consumers stop is deterministic.

The flags can only be closed, never reopened, which is what keeps a guard from
becoming an optimizer that widens a running campaign.
"""
struct DispatchAdmission
    worker_open::Vector{Base.Threads.Atomic{Bool}}
    local_open::Vector{Base.Threads.Atomic{Bool}}
end

function DispatchAdmission(workers::Integer, locals::Integer)
    workers >= 0 || throw(ArgumentError("DispatchAdmission workers must be >= 0; got $(workers)."))
    locals >= 0 || throw(ArgumentError("DispatchAdmission locals must be >= 0; got $(locals)."))
    return DispatchAdmission(
        [Base.Threads.Atomic{Bool}(true) for _ in 1:Int(workers)],
        [Base.Threads.Atomic{Bool}(true) for _ in 1:Int(locals)],
    )
end

@inline function _dispatch_admission_flags(a::DispatchAdmission, class::Symbol)
    class === :worker && return a.worker_open
    class === :local && return a.local_open
    throw(ArgumentError("DispatchAdmission class must be :worker or :local; got :$(class)."))
end

"""
    dispatch_admits(admission, class, ordinal) -> Bool

Whether that consumer may take another job. An out-of-range ordinal admits, so
a dispatcher running more consumers than the admission was built for is never
silently starved.
"""
@inline function dispatch_admits(a::DispatchAdmission, class::Symbol, ordinal::Integer)::Bool
    flags = _dispatch_admission_flags(a, class)
    (1 <= ordinal <= length(flags)) || return true
    return flags[ordinal][]
end

"""
    close_dispatch_consumer!(admission, class, ordinal) -> Bool

Stop one consumer taking further work. Returns whether this call was the one
that closed it. Refuses to close the last open consumer of the whole dispatch:
something has to drain the queue.
"""
function close_dispatch_consumer!(a::DispatchAdmission, class::Symbol, ordinal::Integer)::Bool
    flags = _dispatch_admission_flags(a, class)
    (1 <= ordinal <= length(flags)) || return false
    dispatch_open_count(a) > 1 || return false
    return Base.Threads.atomic_xchg!(flags[ordinal], false)
end

"""
    close_dispatch_class!(admission, class) -> Int

Stop every consumer of one class, and return how many were closed. Refuses if
that class is all that is left open.
"""
function close_dispatch_class!(a::DispatchAdmission, class::Symbol)::Int
    flags = _dispatch_admission_flags(a, class)
    dispatch_open_count(a) > count(f -> f[], flags) || return 0
    closed = 0
    for i in eachindex(flags)
        Base.Threads.atomic_xchg!(flags[i], false) && (closed += 1)
    end
    return closed
end

"""
    dispatch_open_count(admission[, class]) -> Int

Consumers still taking work, in one class or in total.
"""
dispatch_open_count(a::DispatchAdmission, class::Symbol)::Int =
    count(f -> f[], _dispatch_admission_flags(a, class))
dispatch_open_count(a::DispatchAdmission)::Int =
    count(f -> f[], a.worker_open) + count(f -> f[], a.local_open)

"""
    _run_monte_carlo_mixed(f, seeds, spec, worker_ids, local_slots;
                           class_sink=nothing, take_sink=nothing,
                           admission=nothing, on_complete=nothing,
                           worker_runner=nothing) -> Vector{MonteCarloSampleResult}

One job queue, two kinds of consumer: one `@async` feeder per pool worker,
each blocking on a `remotecall_fetch` of a single sample, and `local_slots`
`Threads.@spawn` tasks running samples in this process. Whichever finishes a
sample first takes the next, so a slow slot (a 1-thread worker on a heavy
sample) is balanced against a fast one without any static partition.

`class_sink`, when given, is a `Vector{Symbol}` as long as `seeds` into which
each sample's consumer class (`:worker` or `:local`) is written at its own
index, the same way `samples[index]` is. Nothing in the result itself records
which kind of consumer ran a sample -- a worker's `elapsed_s` is measured on
the worker and comes home indistinguishable from a local slot's -- and the R7
planner's guard needs exactly that split to compare the two classes against
what it predicted.

`take_sink` is the companion: the coordinator's `time_ns()` at the moment that
consumer TOOK the job, written at the same index. With `finished_ns` (stamped
when its result was in hand) it brackets one consumer's occupancy for one
sample -- the work plus that sample's own share of the round trip and
serialization, and nothing else. Measured from a single dispatch-wide start
instead, the number carries the pool's shared setup and the queueing of
whoever was served earlier, which on an 8-worker round reads 965 ms for a
38 ms sample and is not a per-sample cost at all.

`on_complete(sample, class, ordinal)` is called by the consumer that finished
a sample, on that consumer's task, immediately after the sample is recorded.
With `admission` it is how a campaign re-plans itself while it runs: the hook
watches the classes and closes consumers (see [`DispatchAdmission`](@ref)),
and the consumers it closes drop out at their next take. The hook runs on the
dispatch's hot path, so it must be cheap and it must not throw.

`worker_runner` replaces `remotecall_fetch` for the worker class. It exists so
the consumer bookkeeping above -- admission, class tags, take stamps, the
completion hook -- can be tested against both classes in one process, without
a `Distributed` pool. It changes nothing when `nothing`.

Pool workers are `--threads=1`, so with W workers on a T-thread coordinator a
process-only campaign uses W cores and leaves T idle; the local slots are how
the process route fills them (see `ParallelProfiles.mixed_local_slots` for
how many). The caller sets the local samples' inner budget with
`outer_split_env_pairs(local_slots)` around this call; the workers' own pool
is their share. `local_slots = 0` is the process-only dispatch.
"""
function _run_monte_carlo_mixed(
    f, seeds::Vector, spec::MonteCarloSpec, worker_ids::Vector{Int}, local_slots::Int;
    class_sink::Union{Nothing, Vector{Symbol}} = nothing,
    take_sink::Union{Nothing, Vector{Float64}} = nothing,
    admission::Union{Nothing, DispatchAdmission} = nothing,
    on_complete = nothing,
    worker_runner = nothing
)
    (isempty(worker_ids) && local_slots < 1) && throw(ArgumentError(
        "_run_monte_carlo_mixed needs at least one pool worker or one local slot."))
    class_sink === nothing || length(class_sink) == length(seeds) || throw(ArgumentError(
        "_run_monte_carlo_mixed class_sink must have one entry per seed; got " *
        "$(length(class_sink)) for $(length(seeds)) seeds."))
    take_sink === nothing || length(take_sink) == length(seeds) || throw(ArgumentError(
        "_run_monte_carlo_mixed take_sink must have one entry per seed; got " *
        "$(length(take_sink)) for $(length(seeds)) seeds."))
    jobs = Channel{Tuple{Int, Any}}(length(seeds))
    for (index, seed) in enumerate(seeds)
        put!(jobs, (index, seed))
    end
    close(jobs)

    samples = Vector{Union{Nothing, MonteCarloSampleResult}}(nothing, length(seeds))
    stop_requested = Base.Threads.Atomic{Bool}(false)
    run_sample = (index, seed) -> _run_monte_carlo_sample(f, index, seed)
    consume = (run, class, ordinal) -> begin
        while true
            spec.fail_fast && stop_requested[] && break
            # Admission is checked BEFORE the take, not after: a consumer that
            # took a job and then withdrew would drop that sample on the floor.
            (admission === nothing || dispatch_admits(admission, class, ordinal)) || break
            # The channel is filled and closed before any consumer starts, so
            # `take!` either returns a job or throws because the queue is
            # drained. That is what ends the loop; there is no other producer.
            job = try
                take!(jobs)
            catch
                break
            end
            index, seed = job
            take_sink === nothing || (take_sink[index] = Float64(time_ns()))
            sample = _stamp_finished(run(index, seed))
            samples[index] = sample
            class_sink === nothing || (class_sink[index] = class)
            on_complete === nothing || on_complete(sample, class, ordinal)
            if spec.fail_fast && !sample.success
                Base.Threads.atomic_xchg!(stop_requested, true)
                break
            end
        end
    end

    pool = worker_runner === nothing && !isempty(worker_ids) ? CachingPool(worker_ids) : nothing
    try
        Base.@sync begin
            for (ordinal, _) in enumerate(worker_ids)
                run = worker_runner === nothing ?
                    ((index, seed) -> remotecall_fetch(run_sample, pool, index, seed)) :
                    ((index, seed) -> worker_runner(f, index, seed))
                Base.@async consume(run, :worker, ordinal)
            end
            for ordinal in 1:local_slots
                Base.Threads.@spawn consume(run_sample, :local, ordinal)
            end
        end
    finally
        # Drop the cached closure on the workers; it can capture large
        # configuration state that should not outlive the campaign.
        pool === nothing || Distributed.clear!(pool)
    end

    # samples[index] writes plus the in-order filter above already leave
    # `completed` sorted by index; no re-sort needed.
    completed = MonteCarloSampleResult[s for s in samples if s !== nothing]
    if spec.fail_fast
        _throw_first_monte_carlo_failure(completed)
    end
    return completed
end

# Process-backend dispatch: the mixed dispatcher with no local slots.
function _run_monte_carlo_process(f, seeds::Vector, spec::MonteCarloSpec, worker_ids::Vector{Int};
                                  class_sink::Union{Nothing, Vector{Symbol}} = nothing,
                                  take_sink::Union{Nothing, Vector{Float64}} = nothing,
                                  admission::Union{Nothing, DispatchAdmission} = nothing,
                                  on_complete = nothing)
    return _run_monte_carlo_mixed(f, seeds, spec, worker_ids, 0; class_sink = class_sink,
                                  take_sink = take_sink, admission = admission,
                                  on_complete = on_complete)
end

"""
    outer_split_env_pairs(worker_count) -> Vector{Pair{String,String}}

Environment a threaded outer split must run its samples under: the split is
declared active, and each sample's inner thread budget is at most its share of
the pool, `fld(Threads.nthreads(), worker_count)`.

Without the budget a sample resolves `effective_inner_thread_budget()` to the
whole pool. The shipped inner policy (R4/R5) survives that because its AIMD
controller sees the contention and backs off; the V2 static width rule takes
`min(items, budget)` and holds it, so every concurrent sample threads its
callbacks at full pool width. Measured on B15 mcgrid_16sat_8mc at 12 threads,
8 concurrent samples: 10.33 s per sample with the pool advertised, 3.07 s with
the share advertised, 4.56 s for R5 under the same overstatement.

The adaptive campaign runner already did this for the `threads=:auto` path
(adaptive_routing.jl); `run_monte_carlo(f, seeds; threads=N)` did not, so the
integer-threads API paid the full cost.

An inherited budget is a CEILING the split may lower, not one it must honor.

It used to outrank the share outright -- "an explicit user budget always
wins" -- and that is the wrong way round whenever the inherited value is the
wider of the two. Whoever set it did not know how many samples would run
beside each other; this function is the only place that does. Honoring a
whole-pool budget under a W-wide split hands every one of W concurrent
samples the whole pool, which is the overstatement the measurement above
prices, and nothing downstream can detect it: the width and the route are
both still correct, so the campaign looks healthy from the outside.

Measured on this repo's 24-logical-core workstation, mcgrid_8sat_16mc at
32 threads, 16 samples over a warm store, R6, idle box: 1.57 s per campaign
and 43.9 GB summed sample allocation with the share declared, against 2.50 s
and 69.9 GB with an inherited `SPACEAGORA_INNER_THREAD_BUDGET=32` honored --
1.59x on both. `policy_threads_enabled_total` separates the two regimes
cleanly (0 against 874), which is what makes the overstatement identifiable
in a row after the fact.

A narrower explicit budget still wins, since it can only reduce concurrency.
"""
function outer_split_env_pairs(worker_count::Int)::Vector{Pair{String, String}}
    pairs = Pair{String, String}["SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1"]
    share = max(1, fld(Base.Threads.nthreads(), max(1, worker_count)))
    budget = capped_inner_thread_budget(share)
    budget === nothing || push!(pairs, "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget))
    return pairs
end

"""
    capped_inner_thread_budget(share) -> Union{Nothing, Int}

The per-sample inner thread budget an outer split of this `share` must
declare, or `nothing` when the environment already says something the split
is happy with.

`nothing` means "leave `SPACEAGORA_INNER_THREAD_BUDGET` alone": either the
inherited value is already at or below `share`, or it is unparseable, in
which case the split declines to rewrite a setting it cannot read rather than
silently discarding it. Anything else is the value to declare -- `share` when
nothing was inherited, and `share` again when the inherited budget is wider
than the split can afford. See `outer_split_env_pairs` for why a wider
inherited budget does not win.
"""
function capped_inner_thread_budget(share::Int)::Union{Nothing, Int}
    floor_share = max(1, share)
    raw = strip(get(ENV, "SPACEAGORA_INNER_THREAD_BUDGET", ""))
    isempty(raw) && return floor_share
    inherited = tryparse(Int, raw)
    inherited === nothing && return nothing
    # A non-positive inherited budget is how the policy layer spells "the whole
    # pool" (effective_inner_thread_budget), so it is the widest value there is
    # and the split must cap it.
    inherited <= 0 && return floor_share
    return inherited > floor_share ? floor_share : nothing
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
        # Every sample that runs beside others must be told its share of the
        # pool. See outer_split_env_pairs: without this each of the concurrent
        # solves believes it owns every thread, and under the V2 static width
        # rule that belief is acted on for the whole solve.
        withenv(outer_split_env_pairs(worker_count)...) do
            _run_monte_carlo_threaded(f, seeds, spec, worker_count)
        end
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
