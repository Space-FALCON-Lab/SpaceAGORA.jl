# The predictive campaign planner (routing profile R7).
#
# R6 learns the campaign plan: a bandit over :none | :threads | :process, a
# split race over widths inside the first campaign of an unseen shape, and
# persisted state so the next session exploits what this one paid for. That
# works when a shape is run many times and costs a measurable amount when it is
# not -- an exploratory campaign on the wrong arm is a whole campaign, and on
# the TRX50's P1-P5 grid three points ended up more than 10% behind the best
# pinned static route because of what the learner spent or mis-measured.
#
# R7 computes the plan instead. The inputs are the campaign shape (sample
# count, how many pool workers are affordable, how many threads the
# coordinator has) and the machine's calibrated scalability constants; the
# output is one plan, chosen once, before any sample runs. Three rules keep it
# honest:
#
#   1. The default is the plan a PINNED STATIC ROUTE would run. Deviating from
#      it requires a predicted gain larger than a margin, so a model error
#      inside the margin cannot cost anything the static route would not also
#      have cost.
#   2. Nothing is fitted in the campaign. There is no exploration, no racing
#      and no persisted learning state on this path.
#   3. A guard after the first round of samples re-plans the remainder when the
#      observation contradicts the prediction, and it only ever moves TOWARD
#      the static-equivalent plan.
#
# The switch is `SPACEAGORA_CAMPAIGN_PLANNER`; `bandit` (the default) is R6
# unchanged, bit for bit.

using ..SimulationModel: ParallelCost

"""
    campaign_planner_mode() -> Symbol

Which campaign planner `run_monte_carlo(threads=:auto)` uses: `:bandit` (the
R6 learner, the default) or `:predictive` (the R7 planner in this file).

Read from `SPACEAGORA_CAMPAIGN_PLANNER` on every call rather than snapshotted,
so a `withenv` block around one campaign selects the planner for that campaign
only -- which is how the unit tests and the benchmark harness drive it. Values
are case-insensitive; an unrecognized value throws rather than falling back to
a default, because silently running the other planner is exactly the failure a
profile-level switch must not have.
"""
function campaign_planner_mode()::Symbol
    raw = lowercase(strip(get(ENV, "SPACEAGORA_CAMPAIGN_PLANNER", "")))
    isempty(raw) && return :bandit
    raw == "bandit" && return :bandit
    raw == "predictive" && return :predictive
    throw(ArgumentError(
        "SPACEAGORA_CAMPAIGN_PLANNER must be \"bandit\" or \"predictive\"; got \"$(raw)\"."
    ))
end

@inline function _predictive_env_float(name::AbstractString, fallback::Float64)::Float64
    raw = strip(get(ENV, String(name), ""))
    isempty(raw) && return fallback
    parsed = tryparse(Float64, raw)
    parsed === nothing && throw(ArgumentError("$(name) must be a real number; got \"$(raw)\"."))
    return parsed
end

@inline function _predictive_env_int(name::AbstractString, fallback::Int)::Int
    raw = strip(get(ENV, String(name), ""))
    isempty(raw) && return fallback
    parsed = tryparse(Int, raw)
    parsed === nothing && throw(ArgumentError("$(name) must be an integer; got \"$(raw)\"."))
    return parsed
end

# Round-trip cost of handing one sample to a pool worker and taking the result
# back, relative to one uncontended sample time. The default of zero is an
# ASSUMPTION with a stated domain -- a `remotecall_fetch` round trip is
# milliseconds against samples of 0.03-1 s -- and the reduced-scale run on this
# repo's workstation refutes it at the bottom of that domain: on
# independent_1sat_1hr, 64 samples of ~38 ms over 8 pool workers against 8
# coordinator threads, the pool campaign ran 0.722 s against the threads
# route's 0.427 s with the same per-sample work, which is a worker class
# roughly 1.6x a thread rather than 1.0x. The default stays zero because one
# point on one machine is not a value to hard-code; it is now a named field, so
# a machine that has measured its own can declare it.
const PREDICTIVE_REMOTE_OVERHEAD = 0.0

"""
    PredictivePlannerConfig(; margin, guard_factor, local_slots_max)

The R7 planner's three knobs. The no-argument constructor reads each one from
the environment; an explicit keyword overrides the environment (which is how
the unit tests drive a plan space the host machine does not have).

- `margin` (`SPACEAGORA_PREDICTIVE_MARGIN`, default `0.15`): the relative
  predicted gain a non-static plan must show over the best static-equivalent
  plan before the planner will leave the static plan. ASSUMED -- see
  `docs/architecture/predictive_routing_r7.md`.
- `guard_factor` (`SPACEAGORA_PREDICTIVE_GUARD_FACTOR`, default `1.5`): the
  observed/predicted per-class sample-time ratio beyond which the guard
  re-plans the remainder of the campaign. ASSUMED.
- `local_slots_max` (`SPACEAGORA_PREDICTIVE_LOCAL_SLOTS_MAX`, default
  `Threads.nthreads() - 1`): the cap on coordinator local slots under mixed
  dispatch. DERIVED from R6's measured practice: `mixed_local_slots` keeps
  thread 1 free for the `@async` feeders that keep the pool supplied, so the
  most slots the coordinator can ever offer is `T - 1`.
- `remote_overhead` (`SPACEAGORA_PREDICTIVE_REMOTE_OVERHEAD`, default `0.0`):
  the `remotecall_fetch` round trip relative to one uncontended sample time,
  i.e. how much more a pool worker costs per sample than a coordinator thread.
  ASSUMED at zero, and measurably wrong on short samples -- see the constants
  table in `docs/architecture/predictive_routing_r7.md` for the measurement and
  why the default is still zero.
"""
struct PredictivePlannerConfig
    margin::Float64
    guard_factor::Float64
    local_slots_max::Int
    remote_overhead::Float64
end

function PredictivePlannerConfig(;
    margin::Real = _predictive_env_float("SPACEAGORA_PREDICTIVE_MARGIN", 0.15),
    guard_factor::Real = _predictive_env_float("SPACEAGORA_PREDICTIVE_GUARD_FACTOR", 1.5),
    local_slots_max::Integer = _predictive_env_int(
        "SPACEAGORA_PREDICTIVE_LOCAL_SLOTS_MAX", max(0, Base.Threads.nthreads() - 1)),
    remote_overhead::Real = _predictive_env_float(
        "SPACEAGORA_PREDICTIVE_REMOTE_OVERHEAD", PREDICTIVE_REMOTE_OVERHEAD),
)
    margin >= 0.0 || throw(ArgumentError("PredictivePlannerConfig margin must be >= 0; got $(margin)."))
    guard_factor >= 1.0 ||
        throw(ArgumentError("PredictivePlannerConfig guard_factor must be >= 1; got $(guard_factor)."))
    local_slots_max >= 0 ||
        throw(ArgumentError("PredictivePlannerConfig local_slots_max must be >= 0; got $(local_slots_max)."))
    remote_overhead >= 0.0 ||
        throw(ArgumentError("PredictivePlannerConfig remote_overhead must be >= 0; got $(remote_overhead)."))
    return PredictivePlannerConfig(Float64(margin), Float64(guard_factor), Int(local_slots_max),
                                   Float64(remote_overhead))
end


"""
    PredictivePlan

One candidate campaign plan and its predicted cost.

- `route` -- `:none`, `:threads` or `:process`, as `MonteCarloResult.route`.
- `workers` -- concurrent thread tasks for `:threads`, pool worker processes
  for `:process`, `1` for `:none`.
- `local_slots` -- coordinator samples running beside the pool (mixed
  dispatch); zero for every route but `:process`.
- `consumers` -- `workers + local_slots`, i.e. how many samples are in flight.
- `inner_thread_budget` -- threads one sample may use. Always `1` when
  `consumers > 1` (see the module docs for why v1 never splits inward).
- `static_equivalent` -- true when this plan is what a pinned static route
  would have run: `:threads@min(n,T)`, `:process@L=0`, or `:none`.
- `makespan` -- predicted campaign time in units of one uncontended sample.
- `worker_slowdown` / `heap_slowdown` -- the per-class sample-time multipliers
  the makespan was built from. The guard compares these against what the first
  round actually did.
"""
struct PredictivePlan
    route::Symbol
    workers::Int
    local_slots::Int
    consumers::Int
    inner_thread_budget::Int
    static_equivalent::Bool
    makespan::Float64
    worker_slowdown::Float64
    heap_slowdown::Float64
end

"""
    PredictivePlanning

The planner's answer: every candidate it ranked, the one it chose, and why.
`gain` is the chosen plan's relative predicted improvement over the best
static-equivalent plan (`0.0` when the chosen plan is itself that plan).
"""
struct PredictivePlanning
    plans::Vector{PredictivePlan}
    chosen::PredictivePlan
    best_static::Union{Nothing, PredictivePlan}
    gain::Float64
    reason::Symbol
    constants_loaded::Bool
end

"""
    predictive_heap_slowdown(constants, k) -> Float64

How much longer one sample takes when `k` of them share this process's heap,
relative to running alone: `k / S(k)` with `S` the Universal Scalability Law at
the machine's calibrated `(usl_alpha_base, usl_beta_alloc)`.

The parameters are SOURCED -- `scripts/calibrate_machine.jl` fits them on this
machine -- but their application here is ASSUMED: they were fitted on the
allocation kernel, and what is claimed is only that a Monte Carlo sample's
contention with its neighbors on one heap has the same SHAPE as that kernel's,
not the same magnitude. That is why the margin rule sits on top of the model
rather than the model deciding alone.

`nothing` constants give `1.0` at every width: an uncalibrated machine has no
contention model, so the planner predicts no contention, every plan's makespan
is its round count, and the margin is what carries the safety.
"""
function predictive_heap_slowdown(constants::Union{Nothing, ParallelCost.MachineConstants}, k::Integer)::Float64
    k <= 1 && return 1.0
    constants === nothing && return 1.0
    alpha = max(0.0, constants.usl_alpha_base)
    beta = max(0.0, constants.usl_beta_alloc)
    (alpha > 0.0 || beta > 0.0) || return 1.0
    return Float64(k) / ParallelCost.usl_speedup(alpha, beta, Int(k))
end

"""
    predictive_makespan(n_samples, slowdowns) -> Float64

Greedy list-schedule makespan for `n_samples` identical unit samples over
consumers whose per-sample times are `slowdowns`, in units of one uncontended
sample.

Deterministic by construction: every consumer takes the next queued sample the
moment it is free, ties going to the lower index, which is exactly what the
one-queue dispatchers in `monte_carlo.jl` do. The samples are modeled as
identical because the planner has no per-sample information -- a Monte Carlo
campaign's seeds are draws from one distribution -- and treating them as equal
is the only assumption available that does not invent a spread.
"""
function predictive_makespan(n_samples::Integer, slowdowns::Vector{Float64})::Float64
    n = Int(n_samples)
    n <= 0 && return 0.0
    isempty(slowdowns) && return Inf
    # All consumers equal (the uncalibrated case, and every pure-route plan) is
    # a closed form; going through the loop would give the same answer more
    # slowly.
    first_s = slowdowns[1]
    if all(s -> s == first_s, slowdowns)
        return cld(n, length(slowdowns)) * first_s
    end
    free = zeros(Float64, length(slowdowns))
    for _ in 1:n
        j = argmin(free)
        free[j] += slowdowns[j]
    end
    return maximum(free)
end

function _predictive_slowdowns(route::Symbol, workers::Int, local_slots::Int,
                               constants::Union{Nothing, ParallelCost.MachineConstants},
                               remote_overhead::Float64)
    worker_s = 1.0 + remote_overhead
    if route === :process
        heap_s = predictive_heap_slowdown(constants, local_slots)
        return (vcat(fill(worker_s, workers), fill(heap_s, local_slots)), worker_s, heap_s)
    elseif route === :threads
        # Threads-route tasks and mixed-dispatch local slots are the same
        # thing: samples sharing one process's heap and allocator. That is why
        # both are priced with `predictive_heap_slowdown`.
        heap_s = predictive_heap_slowdown(constants, workers)
        return (fill(heap_s, workers), 1.0, heap_s)
    end
    return ([1.0], 1.0, 1.0)
end

function _predictive_plan(route::Symbol, workers::Int, local_slots::Int, n_samples::Int,
                          static_equivalent::Bool,
                          constants::Union{Nothing, ParallelCost.MachineConstants},
                          remote_overhead::Float64 = PREDICTIVE_REMOTE_OVERHEAD)::PredictivePlan
    slowdowns, worker_s, heap_s = _predictive_slowdowns(route, workers, local_slots, constants,
                                                        remote_overhead)
    consumers = route === :process ? workers + local_slots : workers
    return PredictivePlan(
        route, workers, local_slots, consumers,
        consumers > 1 ? 1 : 0,
        static_equivalent,
        predictive_makespan(n_samples, slowdowns),
        worker_s, heap_s,
    )
end

# `inner_thread_budget = 0` above is the planner's spelling of "this plan does
# not declare a budget" for the serial plan, which must leave the sample the
# whole pool. Everything concurrent declares 1.
@inline function _predictive_declared_budget(plan::PredictivePlan)::Int
    return plan.consumers > 1 ? 1 : max(1, plan.inner_thread_budget)
end

"""
    predictive_plan_candidates(; n_samples, threads, process_workers,
                               threads_candidate, local_slots_cap, constants,
                               config) -> Vector{PredictivePlan}

The v1 plan space, deliberately small:

- `:none` -- a single sample, or a single-threaded coordinator with no pool.
  One consumer with the whole thread budget.
- `:threads` at `W = min(n, T)`, budget 1. This is the static `outer_threads`
  route.
- `:process` with `W_p = min(process_workers, n)` pool workers and `L` local
  slots, `L` from `0` up to `min(local_slots_max, n - W_p, local_slots_cap)`.
  `L = 0` is the static `outer_process` route.

Narrower widths are not enumerated. A narrower split runs the same samples on
less of the machine: for the process route the workers are one thread each, so
a narrower pool only idles cores, and for the threads route v1 hands every
concurrent sample a budget of 1 regardless of width, so a narrower split buys
nothing to trade against the lost concurrency. Widths become a real decision
again only when the inner speedup is measurable a priori, which is the v2 item
named in the docs.
"""
function predictive_plan_candidates(;
    n_samples::Integer,
    threads::Integer,
    process_workers::Integer,
    threads_candidate::Bool,
    local_slots_cap::Integer,
    constants::Union{Nothing, ParallelCost.MachineConstants},
    config::PredictivePlannerConfig,
)::Vector{PredictivePlan}
    n = max(0, Int(n_samples))
    T = max(1, Int(threads))
    pool = Int(process_workers) >= 2 ? min(Int(process_workers), n) : 0
    plans = PredictivePlan[]
    if n <= 1 || (T <= 1 && pool == 0)
        push!(plans, _predictive_plan(:none, 1, 0, n, true, constants, config.remote_overhead))
    end
    if n > 1 && T > 1 && threads_candidate
        W = min(n, T)
        W > 1 && push!(plans, _predictive_plan(:threads, W, 0, n, true, constants, config.remote_overhead))
    end
    if n > 1 && pool >= 2
        lmax = max(0, min(config.local_slots_max, n - pool, Int(local_slots_cap)))
        for L in 0:lmax
            push!(plans, _predictive_plan(:process, pool, L, n, L == 0, constants, config.remote_overhead))
        end
    end
    # A shape with no parallel route at all still has to run.
    isempty(plans) && push!(plans, _predictive_plan(:none, 1, 0, n, true, constants, config.remote_overhead))
    return plans
end

# Ranking order. Makespan first; then static-equivalent plans ahead of
# deviations; then, among equals, the process route ahead of the threads route.
#
# That last preference is SOURCED, and it is a TIE-BREAK, not an override: on
# the TRX50 cold-11 measurements the one-heap threads route never beat the pool
# by more than 11% and lost to it by 5-60% at W >= 16 on P3/P4, so when the two
# are predicted to take the same number of rounds the pool is the safer of the
# two. When the threads route is predicted FASTER (P5 mcgrid_8sat_16mc at
# W_p=4, T=8: two rounds against four) it wins on the makespan and this rule
# never fires.
@inline function _predictive_route_rank(route::Symbol)::Int
    route === :process && return 0
    route === :threads && return 1
    return 2
end

@inline function _predictive_sort_key(p::PredictivePlan)
    return (round(p.makespan; sigdigits = 12), p.static_equivalent ? 0 : 1,
            _predictive_route_rank(p.route), p.local_slots)
end

"""
    predictive_plan(; n_samples, threads, process_workers, threads_candidate,
                    local_slots_cap, constants, config) -> PredictivePlanning

Rank the plan space and choose one.

The chosen plan is the best static-equivalent plan unless some other plan is
predicted to beat it by at least `config.margin` relative. The asymmetry is the
point: the static routes are what the paper's baselines measure and what a user
who pins a route gets, so a planner that deviates on a small predicted edge can
only lose ground to them when the model is wrong, while one that deviates on a
large predicted edge keeps R6's measured 25-40% mixed-dispatch wins at the
P3/P4 mid budgets, where the edge is a whole dispatch round.
"""
function predictive_plan(;
    n_samples::Integer,
    threads::Integer = Base.Threads.nthreads(),
    process_workers::Integer = 0,
    threads_candidate::Bool = Base.Threads.nthreads() > 1,
    local_slots_cap::Integer = 0,
    constants::Union{Nothing, ParallelCost.MachineConstants} = nothing,
    config::PredictivePlannerConfig = PredictivePlannerConfig(),
)::PredictivePlanning
    plans = predictive_plan_candidates(
        n_samples = n_samples, threads = threads, process_workers = process_workers,
        threads_candidate = threads_candidate, local_slots_cap = local_slots_cap,
        constants = constants, config = config)
    sort!(plans; by = _predictive_sort_key)
    best = first(plans)
    static_idx = findfirst(p -> p.static_equivalent, plans)
    best_static = static_idx === nothing ? nothing : plans[static_idx]
    if best.static_equivalent
        return PredictivePlanning(plans, best, best_static, 0.0, :static_equivalent_best,
                                  constants !== nothing)
    end
    if best_static === nothing
        return PredictivePlanning(plans, best, nothing, 0.0, :no_static_candidate,
                                  constants !== nothing)
    end
    gain = best_static.makespan > 0.0 ?
        (best_static.makespan - best.makespan) / best_static.makespan : 0.0
    if gain >= config.margin
        return PredictivePlanning(plans, best, best_static, gain, :predicted_gain,
                                  constants !== nothing)
    end
    return PredictivePlanning(plans, best_static, best_static, gain, :margin_not_met,
                              constants !== nothing)
end

"""
    predictive_guard_verdict(plan, config; worker_mean_s, local_mean_s,
                             failures) -> NamedTuple

What to do with the rest of the campaign after its first round.

The prediction under test is the ratio between the two consumer classes:
`heap_slowdown / worker_slowdown`, i.e. how much slower a coordinator local
slot should be than a pool worker. The first round measures the same ratio
directly, since the pool workers are the uncontended reference (`s_w = 1`).
`ratio` below is observed over predicted; `1.0` is the model being right.

The verdict only ever moves toward a static-equivalent plan -- fewer local
slots, never more, never a different route. A guard that could widen would be
an optimizer running inside the campaign, which is the thing R7 exists not to
be. A failed sample is treated as the strongest signal available and drops the
local slots outright: the plan that remains is `:process@L=0`, the static
route, which is the conservative answer when something about the deviation is
not understood.
"""
function predictive_guard_verdict(
    plan::PredictivePlan,
    config::PredictivePlannerConfig;
    worker_mean_s::Float64,
    local_mean_s::Float64,
    failures::Integer = 0,
)
    if plan.route !== :process || plan.local_slots <= 0
        return (replan = false, local_slots = plan.local_slots, ratio = NaN,
                reason = :nothing_to_reduce)
    end
    if failures > 0
        return (replan = true, local_slots = 0, ratio = NaN, reason = :sample_failure)
    end
    ok = isfinite(worker_mean_s) && worker_mean_s > 0.0 &&
         isfinite(local_mean_s) && local_mean_s > 0.0
    ok || return (replan = false, local_slots = plan.local_slots, ratio = NaN,
                  reason = :not_observed)
    predicted = plan.heap_slowdown / plan.worker_slowdown
    predicted > 0.0 || return (replan = false, local_slots = plan.local_slots, ratio = NaN,
                               reason = :not_observed)
    ratio = (local_mean_s / worker_mean_s) / predicted
    if ratio > 2.0 * config.guard_factor
        return (replan = true, local_slots = 0, ratio = ratio, reason = :local_slots_far_slower)
    elseif ratio > config.guard_factor
        return (replan = true, local_slots = plan.local_slots ÷ 2, ratio = ratio,
                reason = :local_slots_slower)
    end
    return (replan = false, local_slots = plan.local_slots, ratio = ratio,
            reason = :observation_matches)
end

"""
    predictive_replan(plan, local_slots, n_samples, constants) -> PredictivePlan

The same route and pool with a reduced local-slot count, re-priced so the
remainder's makespan and the guard's own record describe what is about to run.
"""
function predictive_replan(plan::PredictivePlan, local_slots::Integer, n_samples::Integer,
                           constants::Union{Nothing, ParallelCost.MachineConstants})::PredictivePlan
    L = clamp(Int(local_slots), 0, plan.local_slots)
    # The plan's own worker slowdown, not the environment's: a re-plan must
    # price the remainder the same way the first round was priced.
    return _predictive_plan(plan.route, plan.workers, L, Int(n_samples), L == 0, constants,
                            plan.worker_slowdown - 1.0)
end

# Machine constants are read from disk once per process. They are a calibration
# artifact keyed by machine fingerprint: nothing in a running session can
# change them, and re-parsing the TOML at the head of every campaign would put
# file IO on a path whose whole purpose is to cost less than the learner it
# replaces. `nothing` (no calibration on this machine) is cached too, and is
# the common case.
const _PREDICTIVE_CONSTANTS = Ref{Any}(nothing)
const _PREDICTIVE_CONSTANTS_LOADED = Ref(false)
const _PREDICTIVE_CONSTANTS_LOCK = ReentrantLock()

"""
    predictive_machine_constants() -> Union{Nothing, ParallelCost.MachineConstants}

This machine's calibrated cost constants, or `nothing` when it has never been
calibrated (no `output/parallel_policy_state/cost_constants_<fingerprint>.toml`,
or one written by an older schema). Cached for the life of the process; call
[`reset_predictive_machine_constants!`](@ref) in a test that changes the path.
"""
function predictive_machine_constants()
    lock(_PREDICTIVE_CONSTANTS_LOCK) do
        if !_PREDICTIVE_CONSTANTS_LOADED[]
            _PREDICTIVE_CONSTANTS_LOADED[] = true
            _PREDICTIVE_CONSTANTS[] = try
                ParallelCost.load_machine_constants()
            catch err
                @debug "Machine constants could not be loaded; predicting without contention." exception = err
                nothing
            end
        end
        return _PREDICTIVE_CONSTANTS[]
    end
end

"""
    reset_predictive_machine_constants!()

Forget the cached machine constants so the next campaign re-reads them. Exists
for tests; a normal process loads them once.
"""
function reset_predictive_machine_constants!()::Nothing
    lock(_PREDICTIVE_CONSTANTS_LOCK) do
        _PREDICTIVE_CONSTANTS_LOADED[] = false
        _PREDICTIVE_CONSTANTS[] = nothing
    end
    return nothing
end

function _predictive_plan_line(p::PredictivePlan)::String
    return "$(p.route)@w$(p.workers)+l$(p.local_slots) makespan=$(round(p.makespan; digits=3)) " *
           "consumers=$(p.consumers) s_worker=$(round(p.worker_slowdown; digits=3)) " *
           "s_heap=$(round(p.heap_slowdown; digits=3))" * (p.static_equivalent ? " static" : "")
end
