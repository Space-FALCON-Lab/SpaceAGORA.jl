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
#   2. Nothing is fitted to a plan. There is no racing and no per-plan or
#      per-route reward anywhere. What persists between campaigns is a bounded
#      correction to the cost model's shared PARAMETERS (see "Online
#      corrections" below), and the only exploration is the leash: a
#      one-slot step from the previous campaign's plan for the same shape.
#   3. A guard after the first round of samples re-plans the remainder when the
#      observation contradicts the prediction, and it only ever moves TOWARD
#      the static-equivalent plan.
#
# The switch is `SPACEAGORA_CAMPAIGN_PLANNER`; `bandit` (the default) is R6
# unchanged, bit for bit.

using ..SimulationModel: ParallelCost
import ..SimulationEngine
import TOML

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

@inline function _predictive_env_heap_model(name::AbstractString, fallback::Symbol)::Symbol
    raw = lowercase(strip(get(ENV, String(name), "")))
    isempty(raw) && return fallback
    raw == "none" && return :none
    raw == "locals" && return :locals
    raw == "usl" && return :usl
    throw(ArgumentError("$(name) must be \"none\", \"locals\" or \"usl\"; got \"$(raw)\"."))
end

@inline function _predictive_env_bool(name::AbstractString, fallback::Bool)::Bool
    raw = lowercase(strip(get(ENV, String(name), "")))
    isempty(raw) && return fallback
    raw in ("1", "true", "yes", "on") && return true
    raw in ("0", "false", "no", "off") && return false
    throw(ArgumentError("$(name) must be a boolean; got \"$(raw)\"."))
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
- `heap_model` (`SPACEAGORA_PREDICTIVE_HEAP_MODEL`, default `:locals`): whether
  concurrent samples sharing this process's heap are charged a contention term
  at all, and where. `:none` takes `s_heap = 1` at every width, so plans are
  ranked by round count; `:usl` charges `k / usl_speedup(usl_alpha_base,
  usl_beta_alloc, k)` from the machine's calibrated constants to every plan
  that shares this heap; `:locals` charges that term to coordinator local slots
  in a mixed plan only, never to the threads-route static plan. `:locals` is the
  default because two full TRX50 runs measured the term right about local slots
  and wrong about the threads plan -- see the constants table and
  `docs/architecture/predictive_routing_r7.md`. Machine constants are loaded
  and traced under every model; `:none` declines to use them for contention, it
  does not hide them.
- `local_thrash` (`SPACEAGORA_PREDICTIVE_LOCAL_THRASH`, default `3.0`,
  ASSUMED): the observed local-slot slowdown, relative to a pool worker, past
  which the guard trims the local slots even though the pool is unharmed. Well
  above `guard_factor`, because a local slot that is merely slower than a
  worker still adds throughput; only one that is slow enough to be losing more
  than it adds is worth closing.
- `route_switch` (`SPACEAGORA_PREDICTIVE_GUARD_ROUTE_SWITCH`, default `false`):
  whether the guard may move the remainder of a campaign from the pool to the
  threads route. Built, measured, and shipped OFF: on this repo's workstation
  it made `independent_1sat_1hr` roughly twice as slow, because switching
  leaves the pool idle and the next campaign's first round then reads the cost
  of waking it as evidence that it should switch again. The verdict is still
  computed and traced with the switch off, so the evidence is visible without
  being acted on. See `docs/architecture/predictive_routing_r7.md`.
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
    route_switch::Bool
    heap_model::Symbol
    local_thrash::Float64
end

function PredictivePlannerConfig(;
    margin::Real = _predictive_env_float("SPACEAGORA_PREDICTIVE_MARGIN", 0.15),
    guard_factor::Real = _predictive_env_float("SPACEAGORA_PREDICTIVE_GUARD_FACTOR", 1.5),
    local_slots_max::Integer = _predictive_env_int(
        "SPACEAGORA_PREDICTIVE_LOCAL_SLOTS_MAX", max(0, Base.Threads.nthreads() - 1)),
    remote_overhead::Real = _predictive_env_float(
        "SPACEAGORA_PREDICTIVE_REMOTE_OVERHEAD", PREDICTIVE_REMOTE_OVERHEAD),
    route_switch::Bool = _predictive_env_bool(
        "SPACEAGORA_PREDICTIVE_GUARD_ROUTE_SWITCH", false),
    heap_model::Symbol = _predictive_env_heap_model(
        "SPACEAGORA_PREDICTIVE_HEAP_MODEL", :locals),
    local_thrash::Real = _predictive_env_float(
        "SPACEAGORA_PREDICTIVE_LOCAL_THRASH", 3.0),
)
    margin >= 0.0 || throw(ArgumentError("PredictivePlannerConfig margin must be >= 0; got $(margin)."))
    guard_factor >= 1.0 ||
        throw(ArgumentError("PredictivePlannerConfig guard_factor must be >= 1; got $(guard_factor)."))
    local_slots_max >= 0 ||
        throw(ArgumentError("PredictivePlannerConfig local_slots_max must be >= 0; got $(local_slots_max)."))
    remote_overhead >= 0.0 ||
        throw(ArgumentError("PredictivePlannerConfig remote_overhead must be >= 0; got $(remote_overhead)."))
    heap_model in (:none, :locals, :usl) || throw(ArgumentError(
        "PredictivePlannerConfig heap_model must be :none, :locals or :usl; got :$(heap_model)."))
    local_thrash >= 1.0 || throw(ArgumentError(
        "PredictivePlannerConfig local_thrash must be >= 1; got $(local_thrash)."))
    return PredictivePlannerConfig(Float64(margin), Float64(guard_factor), Int(local_slots_max),
                                   Float64(remote_overhead), route_switch, heap_model,
                                   Float64(local_thrash))
end

# The constants the CONTENTION term may use FOR THIS ROUTE, which is not the
# same question as which constants were loaded. The planner can have machine
# constants in hand and decline to charge contention with them; they are still
# reported, so a trace says what the machine knows as well as what the planner
# used.
#
# The route argument is the whole of the `:locals` model: coordinator local
# slots in a mixed plan are charged, threads-route tasks are not, and that
# asymmetry is measured rather than assumed. See the docs.
@inline function _predictive_contention_constants(
    config::PredictivePlannerConfig,
    constants::Union{Nothing, ParallelCost.MachineConstants},
    route::Symbol,
)
    config.heap_model === :usl && return constants
    config.heap_model === :locals && return route === :process ? constants : nothing
    return nothing
end

"""
    PredictiveCostTerms(; heap_scale = 1.0, round_tail = 0.0, startup = 0.0)

The cost-model parameters that are not properties of one plan, applied to every
candidate alike. The default is the model as it was before any of them
existed, and it prices every plan exactly as that model did.

- `heap_scale`: multiplies the heap slowdown the heap model charges a plan,
  wherever it charges one; the result is clamped at `1`. It is the parameter
  the online corrections move when the guard observes the local slots running
  faster or slower than predicted. It never introduces a contention term where
  the heap model charges none.
- `round_tail`: the spread of one sample's time, in sample units, as it shows
  in a FINAL round. Each consumer class (the pool's workers; the consumers on
  this process's heap) finishes its last round of `k` concurrent samples
  `round_tail * (H_k - 1)` after the list schedule says, `H_k` the `k`-th
  harmonic number: the expected excess of the slowest of `k` samples over
  their mean when sample times have an exponential tail. A final round of one
  sample costs nothing extra; a final round of sixteen costs `2.38 x
  round_tail`. Only the final round pays: in earlier rounds a consumer that
  drew a slow sample simply takes fewer of the rest. See "The final round" in
  `docs/architecture/predictive_routing_r7.md` for the measurement.
- `startup`: delay, in sample units, before each pool worker takes its first
  sample, for a pool that has not dispatched in this process yet. The local
  slots start at once, so a cold pool shifts the first samples onto them.
"""
struct PredictiveCostTerms
    heap_scale::Float64
    round_tail::Float64
    startup::Float64
    function PredictiveCostTerms(heap_scale::Real, round_tail::Real, startup::Real)
        (isfinite(heap_scale) && heap_scale > 0.0) || throw(ArgumentError(
            "PredictiveCostTerms heap_scale must be finite and > 0; got $(heap_scale)."))
        (isfinite(round_tail) && round_tail >= 0.0) || throw(ArgumentError(
            "PredictiveCostTerms round_tail must be finite and >= 0; got $(round_tail)."))
        (isfinite(startup) && startup >= 0.0) || throw(ArgumentError(
            "PredictiveCostTerms startup must be finite and >= 0; got $(startup)."))
        return new(Float64(heap_scale), Float64(round_tail), Float64(startup))
    end
end

PredictiveCostTerms(; heap_scale::Real = 1.0, round_tail::Real = 0.0, startup::Real = 0.0) =
    PredictiveCostTerms(heap_scale, round_tail, startup)

"""
    InnerSpeedupCurve(speedup; source = "")

How much faster one sample runs on `b` threads than on one: `speedup[b]` for
`b = 1, ..., length(speedup)`, with `speedup[1] == 1` and never decreasing.
Read from the RHS calibration store's per-candidate timings
(`SimulationEngine.rhs_inner_speedup_curve`); past the widest measured width it
stays at its last value, since nothing measured says more.

It is the RHS evaluation's speedup, applied to the whole sample. That is
ASSUMED: a sample also spends time outside the RHS (the solver's own linear
algebra, callbacks, output), which does not speed up, so the curve is an upper
bound on the sample's speedup and the margin rule carries the difference.
"""
struct InnerSpeedupCurve
    speedup::Vector{Float64}
    source::String
    function InnerSpeedupCurve(speedup::AbstractVector{<:Real}, source::AbstractString)
        isempty(speedup) && throw(ArgumentError("InnerSpeedupCurve needs at least one width."))
        v = Float64.(collect(speedup))
        all(x -> isfinite(x) && x > 0.0, v) || throw(ArgumentError(
            "InnerSpeedupCurve speedups must be finite and positive; got $(v)."))
        # b = 1 is the reference; a wider budget can always run the narrower
        # plan, so the curve never falls.
        v[1] = 1.0
        for b in 2:length(v)
            v[b] = max(v[b], v[b - 1])
        end
        return new(v, String(source))
    end
end

InnerSpeedupCurve(speedup::AbstractVector{<:Real}; source::AbstractString = "") =
    InnerSpeedupCurve(speedup, source)

@inline inner_speedup(c::InnerSpeedupCurve, b::Integer)::Float64 =
    c.speedup[clamp(Int(b), 1, length(c.speedup))]

"""
    PredictivePlan

One candidate campaign plan and its predicted cost.

- `route` -- `:none`, `:threads` or `:process`, as `MonteCarloResult.route`.
- `workers` -- concurrent thread tasks for `:threads`, pool worker processes
  for `:process`, `1` for `:none`.
- `local_slots` -- coordinator samples running beside the pool (mixed
  dispatch); zero for every route but `:process`.
- `consumers` -- `workers + local_slots`, i.e. how many samples are in flight.
- `inner_thread_budget` -- threads one sample may use: `1` for every
  concurrent plan unless an inner-speedup curve is known, in which case the
  planner also offers plans whose samples run on `b > 1` threads (threads-route
  tasks, the local slots of a mixed plan, or the one consumer of `:none`). The
  pool's workers always run one thread each.
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

This function is the `:usl` heap model itself, and it is not on by default.
`PredictivePlannerConfig.heap_model` decides whether the planner feeds it the
machine's constants at all; under `:none` it is called with `nothing` even on a
calibrated machine. The measurement that made `:none` the default is in the
constants table of `docs/architecture/predictive_routing_r7.md`.
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

"""
    predictive_round_excess(k) -> Float64

`H_k - 1`: the expected excess of the slowest of `k` samples over their mean,
in units of the spread, for sample times with an exponential tail. Zero for
`k <= 1`.
"""
function predictive_round_excess(k::Integer)::Float64
    total = 0.0
    for i in 2:Int(k)
        total += 1.0 / i
    end
    return total
end

"""
    predictive_plan_makespan(n_samples, pool_workers, slowdowns, terms) -> Float64

[`predictive_makespan`](@ref) with the plan-independent terms added. The first
`pool_workers` entries of `slowdowns` are pool workers; the rest are consumers
on this process's heap (local slots, or threads-route tasks).

- Each pool worker becomes free for its first sample `terms.startup` after the
  dispatch starts; a heap consumer is free at once.
- Each of the two classes finishes `terms.round_tail * (H_k - 1)` after its
  last consumer does in the list schedule, `k` the number of its consumers that
  ran its final round (those with the class's highest sample count).

With both terms zero this IS `predictive_makespan` on the same arguments, so a
machine without the terms is priced bit for bit as before them.
"""
function predictive_plan_makespan(n_samples::Integer, pool_workers::Integer,
                                  slowdowns::Vector{Float64},
                                  terms::PredictiveCostTerms)::Float64
    W = clamp(Int(pool_workers), 0, length(slowdowns))
    if terms.round_tail == 0.0 && (W == 0 || terms.startup == 0.0)
        return predictive_makespan(n_samples, slowdowns)
    end
    n = Int(n_samples)
    n <= 0 && return 0.0
    isempty(slowdowns) && return Inf
    C = length(slowdowns)
    free = zeros(Float64, C)
    free[1:W] .= terms.startup
    taken = zeros(Int, C)
    for _ in 1:n
        j = argmin(free)
        free[j] += slowdowns[j]
        taken[j] += 1
    end
    finish = 0.0
    for (lo, hi) in ((1, W), (W + 1, C))
        lo > hi && continue
        top = maximum(@view taken[lo:hi])
        top == 0 && continue
        k = count(==(top), @view taken[lo:hi])
        finish = max(finish, maximum(@view free[lo:hi]) + terms.round_tail * predictive_round_excess(k))
    end
    return finish
end

# The heap slowdown a plan is charged: the heap model's value, scaled by the
# correction, never below 1. Where the heap model charges nothing (no
# constants, or a single consumer) the correction has nothing to scale.
@inline function _predictive_scaled_heap(constants::Union{Nothing, ParallelCost.MachineConstants},
                                         k::Integer, heap_scale::Float64)::Float64
    raw = predictive_heap_slowdown(constants, k)
    (constants === nothing || k <= 1) && return raw
    return max(1.0, heap_scale * raw)
end

# `inner_time` is the sample time at the plan's inner budget relative to one
# thread, `1 / speedup(b)`; it multiplies every consumer on this process's heap
# and the `:none` consumer, never a pool worker (one thread each).
function _predictive_slowdowns(route::Symbol, workers::Int, local_slots::Int,
                               constants::Union{Nothing, ParallelCost.MachineConstants},
                               remote_overhead::Float64, heap_scale::Float64 = 1.0,
                               inner_time::Float64 = 1.0)
    worker_s = 1.0 + remote_overhead
    if route === :process
        heap_s = _predictive_scaled_heap(constants, local_slots, heap_scale) * inner_time
        return (vcat(fill(worker_s, workers), fill(heap_s, local_slots)), worker_s, heap_s)
    elseif route === :threads
        # Threads-route tasks and mixed-dispatch local slots are the same
        # thing: samples sharing one process's heap and allocator. That is why
        # both are priced with `predictive_heap_slowdown`.
        heap_s = _predictive_scaled_heap(constants, workers, heap_scale) * inner_time
        return (fill(heap_s, workers), 1.0, heap_s)
    end
    return ([inner_time], 1.0, inner_time)
end

function _predictive_plan(route::Symbol, workers::Int, local_slots::Int, n_samples::Int,
                          static_equivalent::Bool,
                          constants::Union{Nothing, ParallelCost.MachineConstants},
                          remote_overhead::Float64 = PREDICTIVE_REMOTE_OVERHEAD;
                          terms::PredictiveCostTerms = PredictiveCostTerms(),
                          inner_budget::Int = 0,
                          inner_time::Float64 = 1.0)::PredictivePlan
    slowdowns, worker_s, heap_s = _predictive_slowdowns(route, workers, local_slots, constants,
                                                        remote_overhead, terms.heap_scale,
                                                        inner_time)
    consumers = route === :process ? workers + local_slots : workers
    pool = route === :process ? workers : 0
    # `0` is the planner's spelling of "no budget declared" for the serial
    # plan; everything concurrent declares one, `1` unless a curve priced more.
    budget = inner_budget > 0 ? inner_budget : (consumers > 1 ? 1 : 0)
    return PredictivePlan(
        route, workers, local_slots, consumers,
        budget,
        static_equivalent,
        predictive_plan_makespan(n_samples, pool, slowdowns, terms),
        worker_s, heap_s,
    )
end

# `inner_thread_budget = 0` is the planner's spelling of "this plan does not
# declare a budget" for the serial plan, which must leave the sample the whole
# pool. Everything concurrent declares its budget, 1 unless a curve priced more.
@inline function _predictive_declared_budget(plan::PredictivePlan)::Int
    return max(1, plan.inner_thread_budget)
end

"""
    predictive_plan_candidates(; n_samples, threads, process_workers,
                               threads_candidate, local_slots_cap, constants,
                               config, terms, inner_curve) -> Vector{PredictivePlan}

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
    terms::PredictiveCostTerms = PredictiveCostTerms(),
    inner_curve::Union{Nothing, InnerSpeedupCurve} = nothing,
)::Vector{PredictivePlan}
    n = max(0, Int(n_samples))
    T = max(1, Int(threads))
    pool = Int(process_workers) >= 2 ? min(Int(process_workers), n) : 0
    # `constants` says what the machine has measured; these say what this
    # planner will charge for sharing a heap, which under the default heap
    # model is the local slots and nothing else. See
    # `_predictive_contention_constants`.
    contention_process = _predictive_contention_constants(config, constants, :process)
    contention_threads = _predictive_contention_constants(config, constants, :threads)
    plans = PredictivePlan[]
    if n <= 1 || (T <= 1 && pool == 0)
        push!(plans, _predictive_plan(:none, 1, 0, n, true, nothing, config.remote_overhead))
    end
    if n > 1 && T > 1 && threads_candidate
        W = min(n, T)
        W > 1 && push!(plans, _predictive_plan(:threads, W, 0, n, true, contention_threads,
                                              config.remote_overhead; terms = terms))
    end
    if n > 1 && pool >= 2
        lmax = max(0, min(config.local_slots_max, n - pool, Int(local_slots_cap)))
        for L in 0:lmax
            push!(plans, _predictive_plan(:process, pool, L, n, L == 0, contention_process,
                                          config.remote_overhead; terms = terms))
        end
    end
    # A shape with no parallel route at all still has to run.
    isempty(plans) && push!(plans, _predictive_plan(:none, 1, 0, n, true, nothing, config.remote_overhead))
    inner_curve === nothing ||
        _predictive_inner_candidates!(plans, n, T, pool, threads_candidate, local_slots_cap,
                                      contention_process, contention_threads, config, terms,
                                      inner_curve)
    return plans
end

# Plans whose samples run on b > 1 threads, offered only when a measured
# inner-speedup curve says what b buys. All of them are deviations from the
# static plan: whatever a curve predicts, spending it has to clear the margin.
#
#   threads@W, b = fld(T, W), for W = fld(T, b), b = 2, 4, 8, ... and for
#     W = min(n, T) itself (the pinned outer_threads split);
#   none, b = T: one sample at a time on the whole pool;
#   process@W_p + L with the local slots at b threads, L * b <= T - 1 (thread 1
#     stays free for the feeders, as for every mixed plan).
#
# The pool's workers are one-thread processes and are never given b > 1.
function _predictive_inner_candidates!(plans::Vector{PredictivePlan}, n::Int, T::Int, pool::Int,
                                       threads_candidate::Bool, local_slots_cap::Integer,
                                       contention_process, contention_threads,
                                       config::PredictivePlannerConfig,
                                       terms::PredictiveCostTerms,
                                       curve::InnerSpeedupCurve)::Nothing
    (n > 0 && T > 1) || return nothing
    time_at(b) = 1.0 / inner_speedup(curve, b)
    budgets = Int[]
    b = 2
    while b <= T
        push!(budgets, b)
        b *= 2
    end
    if n > 1 && threads_candidate
        widths = Int[]
        for b in budgets
            W = fld(T, b)
            W >= 2 && W <= min(n, T) && push!(widths, W)
        end
        full = min(n, T)
        fld(T, full) >= 2 && push!(widths, full)
        for W in sort!(unique!(widths))
            bW = fld(T, W)
            push!(plans, _predictive_plan(:threads, W, 0, n, false, contention_threads,
                                          config.remote_overhead; terms = terms,
                                          inner_budget = bW, inner_time = time_at(bW)))
        end
    end
    # A static serial plan already runs its one sample on the whole pool.
    any(p -> p.route === :none, plans) ||
        push!(plans, _predictive_plan(:none, 1, 0, n, false, nothing, config.remote_overhead;
                                      terms = terms, inner_budget = T, inner_time = time_at(T)))
    if n > 1 && pool >= 2
        for b in budgets
            lmax = max(0, min(config.local_slots_max, n - pool, Int(local_slots_cap), fld(T - 1, b)))
            for L in 1:lmax
                push!(plans, _predictive_plan(:process, pool, L, n, false, contention_process,
                                              config.remote_overhead; terms = terms,
                                              inner_budget = b, inner_time = time_at(b)))
            end
        end
    end
    return nothing
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
            _predictive_route_rank(p.route), p.local_slots, max(1, p.inner_thread_budget))
end

"""
    predictive_plan(; n_samples, threads, process_workers, threads_candidate,
                    local_slots_cap, constants, config, terms,
                    inner_curve) -> PredictivePlanning

Rank the plan space and choose one.

`terms` carries the parameters shared by every candidate (see
[`PredictiveCostTerms`](@ref)); the default reproduces the model without them.
`inner_curve` (see [`InnerSpeedupCurve`](@ref)) adds plans whose samples run on
more than one thread; without one the plan space is v1's.
The decision rule below is the same whatever `terms` holds: the terms change
the prices, never the rule that turns prices into a plan.

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
    terms::PredictiveCostTerms = PredictiveCostTerms(),
    inner_curve::Union{Nothing, InnerSpeedupCurve} = nothing,
)::PredictivePlanning
    plans = predictive_plan_candidates(
        n_samples = n_samples, threads = threads, process_workers = process_workers,
        threads_candidate = threads_candidate, local_slots_cap = local_slots_cap,
        constants = constants, config = config, terms = terms, inner_curve = inner_curve)
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
    predictive_trim_target(plan, config, constants, local_slowdown) -> Int

How many local slots to keep once the guard has decided to trim.

Not half of them. Halving is a step in a direction, not an answer, and on the
TRX50's P4 at 8 threads it cut seven slots to three where seven was the faster
plan (R6 ran `w8+l7` at 2.563 s against the pinned pool's 3.800 s).

With a calibrated contention curve the answer is the curve's own turnover.
Throughput with `W` workers and `L` local slots is `W/c_w + L/c_l(L)`, and
modeling `c_l(L)` as the observation scaled along the calibrated shape,
`c_l(L) = c_obs * s_heap(L) / s_heap(L_0)`, makes the second term proportional
to `L / s_heap(L)` -- which is `usl_speedup` itself. So the maximizing width is
[`usl_peak_workers`](@ref), `sqrt((1-alpha)/beta)`, and the observation sets
the level while the curve sets the shape; the argmax depends only on the shape.

Two fallbacks, both narrower than the plan by construction. A curve whose peak
is at or beyond the running width says the width is fine, which contradicts the
observation that made the guard fire -- there the observation wins and the
width is halved. A machine with no curve at all has nothing better than halving
either. Both are ASSUMED, and both only ever reduce.
"""
function predictive_trim_target(plan::PredictivePlan, config::PredictivePlannerConfig,
                                constants::Union{Nothing, ParallelCost.MachineConstants},
                                local_slowdown::Real = NaN)::Int
    L0 = plan.local_slots
    L0 > 1 || return 0
    curve = _predictive_contention_constants(config, constants, :process)
    if curve !== nothing
        alpha = max(0.0, curve.usl_alpha_base)
        beta = max(0.0, curve.usl_beta_alloc)
        peak = ParallelCost.usl_peak_workers(alpha, beta)
        if isfinite(peak) && peak >= 1.0
            target = clamp(round(Int, peak), 0, L0)
            target < L0 && return target
        end
    end
    return L0 ÷ 2
end

"""
    predictive_guard_verdict(plan, config; worker_mean_s, local_mean_s,
                             worker_occupancy_s, local_occupancy_s, remaining,
                             threads, threads_candidate, constants,
                             failures) -> NamedTuple

What to do with the rest of the campaign after its first round.

Round one runs exactly one sample per consumer, so it measures both consumer
classes directly and at the same moment. Two different things come out of it,
and they answer two different questions.

`worker_mean_s` / `local_mean_s` are the classes' mean `elapsed_s`: the work
INSIDE the sample, timed where the sample ran. Their ratio tests the plan's
`heap_slowdown / worker_slowdown` -- how much slower a coordinator local slot
should be than a pool worker, which is what the contention model predicted.
`ratio` is observed over predicted; `1.0` is the model being right. Too high
means the local slots cost more than the model said, and the verdict reduces
them toward `:process@L=0`.

`worker_occupancy_s` / `local_occupancy_s` are the classes' mean wall from the
dispatch starting to that sample's result being in hand -- the work PLUS
everything around it, which for a pool worker is the `remotecall_fetch` round
trip and the closure's serialization to that worker, and for a local slot is
little more than task scheduling. This is the measurement the planner's
`remote_overhead` constant stands in for a priori, and it is the one that says
which SIDE the campaign belongs on: a worker that occupies far more than a
local slot for the same work is a worker that is not paying for itself, and
the remainder belongs on the pinned threads route. Nothing inside a sample's
`elapsed_s` can see this, which is why the class means alone were not enough:
measured on this repo's workstation, `independent_1sat_1hr` at 64 samples of
~38 ms, the two classes' `elapsed_s` differ by 21% while the campaign on the
pool takes 1.7x the campaign on threads.

So the guard has two directions, and both of them end on a
static-equivalent plan:

- local slots slower than predicted -> halve them, or drop them entirely past
  twice the factor. The plan that remains is `:process@L=0`.
- pool workers occupying more than `guard_factor` times a local slot, AND the
  pinned threads plan predicted (from the observed occupancies) to finish the
  remainder sooner than continuing -> the remainder runs `:threads` at
  `min(remaining, T)`, budget 1, with the pool left idle.

It still never widens and never invents a plan the planner would not have
enumerated. A failed sample is the strongest signal available and drops the
local slots outright, which is the conservative answer when something about
the deviation is not understood.

When both directions qualify, the route change wins: being on the wrong side
of the machine costs more than holding too many local slots on the right one.

The threads plan is priced from the observed LOCAL class, because a local slot
and a threads-route task are the same thing -- a sample on this process's
heap. Widening from `L` slots to `W` tasks is charged the USL ratio
`s_heap(W) / s_heap(L)` when the machine has constants, and nothing when it
does not, which is ASSUMED and optimistic for the threads plan; it is clamped
at 1 so widening is never predicted to make a sample cheaper.
"""
function predictive_guard_verdict(
    plan::PredictivePlan,
    config::PredictivePlannerConfig;
    worker_mean_s::Float64,
    local_mean_s::Float64,
    worker_occupancy_s::Float64 = NaN,
    local_occupancy_s::Float64 = NaN,
    worker_steady_occupancy_s::Float64 = NaN,
    remaining::Integer = 0,
    threads::Integer = Base.Threads.nthreads(),
    threads_candidate::Bool = true,
    constants::Union{Nothing, ParallelCost.MachineConstants} = nothing,
    failures::Integer = 0,
)
    keep = (; replan = false, route = plan.route, workers = plan.workers,
            local_slots = plan.local_slots, ratio = NaN, occupancy_ratio = NaN,
            local_slowdown = NaN, worker_degradation = NaN,
            continue_s = NaN, threads_s = NaN, reason = :nothing_to_reduce)
    (plan.route === :process && plan.local_slots > 0) || return keep
    if failures > 0
        return (; keep..., replan = true, local_slots = 0, reason = :sample_failure)
    end
    ok = isfinite(worker_mean_s) && worker_mean_s > 0.0 &&
         isfinite(local_mean_s) && local_mean_s > 0.0
    ok || return (; keep..., reason = :not_observed)
    predicted = plan.heap_slowdown / plan.worker_slowdown
    predicted > 0.0 || return (; keep..., reason = :not_observed)
    ratio = (local_mean_s / worker_mean_s) / predicted

    # Direction two: is the campaign on the wrong side of the machine?
    occupancy_ok = isfinite(worker_occupancy_s) && worker_occupancy_s > 0.0 &&
                   isfinite(local_occupancy_s) && local_occupancy_s > 0.0
    occupancy_ratio = occupancy_ok ? worker_occupancy_s / local_occupancy_s : NaN
    n_rest = max(0, Int(remaining))
    width = min(n_rest, max(1, Int(threads)))
    if occupancy_ok && threads_candidate && width > 1 && n_rest > 0 &&
       occupancy_ratio > config.guard_factor
        # The destination is the threads route, so it is priced the way the
        # threads route is priced under this heap model -- which under
        # `:locals` means no contention term, hence no widening charge.
        contention = _predictive_contention_constants(config, constants, :threads)
        scale = if contention === nothing
            1.0
        else
            predictive_heap_slowdown(contention, width) /
                max(1.0e-12, predictive_heap_slowdown(contention, plan.local_slots))
        end
        threads_unit = local_occupancy_s * max(1.0, scale)
        threads_s = predictive_makespan(n_rest, fill(threads_unit, width))
        continue_s = predictive_makespan(
            n_rest, vcat(fill(worker_occupancy_s, plan.workers),
                         fill(local_occupancy_s, plan.local_slots)))
        if threads_s < continue_s
            # Computed and reported either way; acted on only when the switch
            # is enabled. See `PredictivePlannerConfig.route_switch` for the
            # measurement that turned it off.
            config.route_switch || return (; keep..., ratio = ratio,
                    occupancy_ratio = occupancy_ratio, continue_s = continue_s,
                    threads_s = threads_s, reason = :route_switch_disabled)
            return (; keep..., replan = true, route = :threads, workers = width,
                    local_slots = 0, ratio = ratio, occupancy_ratio = occupancy_ratio,
                    continue_s = continue_s, threads_s = threads_s,
                    reason = :workers_occupying_more_than_threads)
        end
        return (; keep..., ratio = ratio, occupancy_ratio = occupancy_ratio,
                continue_s = continue_s, threads_s = threads_s,
                reason = :threads_no_better)
    end

    # Direction one: are the local slots worth keeping?
    #
    # Being slower than a pool worker is not a reason to close one. A local
    # slot at 1.5x a worker still adds two thirds of a worker's throughput, and
    # the old rule -- trim whenever the observed/predicted ratio passed
    # `guard_factor` -- closed slots that were paying their way: on the TRX50's
    # P4 at 8 threads it cut seven to three, where R6's seven ran 2.563 s
    # against the pinned pool's 3.800 s.
    #
    # Two things make a local slot worth closing, and neither is "slower".
    #
    #   (a) The pool is being harmed through the coordinator. The local slots
    #       share this process with the `@async` feeders that keep the workers
    #       supplied, so enough local work stalls the feeders and the workers
    #       go idle waiting for jobs. That shows up as the workers' own STEADY
    #       occupancy -- their cost per sample after their first, which is free
    #       of the one-time per-worker dispatch cost -- inflating against the
    #       work they report doing.
    #
    #   (b) The local slots are thrashing among themselves, past
    #       `local_thrash`, which is far enough above `guard_factor` that a
    #       merely-slower slot is left alone.
    #
    # A steady worker observation does not exist yet when the guard decides on
    # a short campaign -- the decision fires once every consumer has reported
    # ONCE, and for a worker that one sample is its first. Rule (a) then cannot
    # fire, and only (b) can, which is the conservative way round.
    local_slowdown = local_mean_s / worker_mean_s
    worker_degradation = (isfinite(worker_steady_occupancy_s) && worker_steady_occupancy_s > 0.0) ?
        (worker_steady_occupancy_s / worker_mean_s) / max(1.0e-12, plan.worker_slowdown) : NaN
    pool_harmed = isfinite(worker_degradation) && worker_degradation > config.guard_factor
    locals_thrashing = local_slowdown > config.local_thrash
    if pool_harmed || locals_thrashing
        target = predictive_trim_target(plan, config, constants, local_slowdown)
        target < plan.local_slots && return (; keep..., replan = true, local_slots = target,
                ratio = ratio, occupancy_ratio = occupancy_ratio,
                local_slowdown = local_slowdown, worker_degradation = worker_degradation,
                reason = pool_harmed ? :workers_degraded : :local_slots_thrashing)
    end
    return (; keep..., ratio = ratio, occupancy_ratio = occupancy_ratio,
            local_slowdown = local_slowdown, worker_degradation = worker_degradation,
            reason = :observation_matches)
end

"""
    predictive_replan(plan, local_slots, n_samples, constants;
                      route = plan.route, workers = plan.workers) -> PredictivePlan

The plan the remainder of a campaign runs under, re-priced for the samples
that are left so its makespan and the guard's record describe what is about to
run.

Two shapes only, and the function refuses anything else: the same route and
pool with fewer local slots, or the threads plan the guard's second direction
moves to. Both are moves toward a static-equivalent plan, which is the guard's
invariant; enforcing it here rather than at the call site means a future caller
cannot quietly widen through this door.

Since the guard became barrier-free the runtime no longer builds a plan for a
remainder -- it closes consumers inside the one dispatch instead (see
`DispatchAdmission`) -- so nothing in `src/` calls this. It is kept as the
executable statement of what the two legal moves are, and as the constructor
for a post-guard plan that anything reasoning about one should use rather than
assembling a `PredictivePlan` by hand.
"""
function predictive_replan(plan::PredictivePlan, local_slots::Integer, n_samples::Integer,
                           constants::Union{Nothing, ParallelCost.MachineConstants};
                           route::Symbol = plan.route,
                           workers::Integer = plan.workers)::PredictivePlan
    # The plan's own worker slowdown, not the environment's: a re-plan must
    # price the remainder the same way the first round was priced.
    overhead = plan.worker_slowdown - 1.0
    if route === plan.route && Int(workers) == plan.workers
        L = clamp(Int(local_slots), 0, plan.local_slots)
        return _predictive_plan(plan.route, plan.workers, L, Int(n_samples), L == 0,
                                constants, overhead)
    end
    route === :threads || throw(ArgumentError(
        "predictive_replan may only reduce local slots or move to the threads route; " *
        "got route=:$(route) from a :$(plan.route) plan."))
    Int(local_slots) == 0 || throw(ArgumentError(
        "predictive_replan to the threads route must carry no local slots; got $(local_slots)."))
    W = max(1, Int(workers))
    return _predictive_plan(:threads, W, 0, Int(n_samples), true, constants, overhead)
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
        _PREDICTIVE_CAMPAIGN_CONSTANTS[] = nothing
    end
    return nothing
end

# ── Measured campaign-level costs ────────────────────────────────────────────

"""
    PredictiveCampaignConstants

This machine's two measured campaign-level costs, read from the `[campaign]`
table of its fingerprinted machine constants file
(`output/parallel_policy_state/cost_constants_<fingerprint>.toml`):

- `round_tail` -- the spread of a sample's time as it shows in a final round,
  in units of one sample time; the prior for
  [`PredictiveCostTerms`](@ref)`.round_tail`.
- `pool_startup_s` -- the delay, in seconds, before a cold pool's workers
  return their first sample beyond the work it took.
- `source` -- where the two numbers were measured.

Both are measured from campaign dispatch traces, not by
`scripts/calibrate_machine.jl`; `scripts/extract_campaign_cost_terms.py`
extracts them from an archived run and writes this table. A field that is
absent is `nothing`, and the planner then charges no such term: an unmeasured
cost is not modeled.
"""
struct PredictiveCampaignConstants
    round_tail::Union{Nothing, Float64}
    pool_startup_s::Union{Nothing, Float64}
    source::String
end

PredictiveCampaignConstants(; round_tail = nothing, pool_startup_s = nothing,
                            source::AbstractString = "") =
    PredictiveCampaignConstants(round_tail === nothing ? nothing : Float64(round_tail),
                                pool_startup_s === nothing ? nothing : Float64(pool_startup_s),
                                String(source))

"""
    load_predictive_campaign_constants(path = ParallelCost.machine_constants_path())
        -> PredictiveCampaignConstants

Read the `[campaign]` table. A missing file, a file without the table, or a
value that is not a finite non-negative number gives `nothing` for that field.
The table is independent of the rest of the file's schema version: the rest is
rewritten by every machine calibration, and `save_machine_constants` carries
this table across the rewrite.
"""
function load_predictive_campaign_constants(
    path::AbstractString = ParallelCost.machine_constants_path())::PredictiveCampaignConstants
    none = PredictiveCampaignConstants()
    isfile(path) || return none
    parsed = try
        TOML.parsefile(String(path))
    catch
        return none
    end
    table = get(parsed, "campaign", nothing)
    table isa AbstractDict || return none
    field(key) = begin
        v = get(table, key, nothing)
        (v isa Real && isfinite(v) && v >= 0) ? Float64(v) : nothing
    end
    source = get(table, "source", "")
    return PredictiveCampaignConstants(field("round_tail"), field("pool_startup_s"),
                                       source isa AbstractString ? String(source) : "")
end

const _PREDICTIVE_CAMPAIGN_CONSTANTS = Ref{Any}(nothing)

"""
    predictive_campaign_constants() -> PredictiveCampaignConstants

[`load_predictive_campaign_constants`](@ref), cached for the life of the
process alongside the machine constants and forgotten with them by
[`reset_predictive_machine_constants!`](@ref).
"""
function predictive_campaign_constants()::PredictiveCampaignConstants
    lock(_PREDICTIVE_CONSTANTS_LOCK) do
        if _PREDICTIVE_CAMPAIGN_CONSTANTS[] === nothing
            _PREDICTIVE_CAMPAIGN_CONSTANTS[] = try
                load_predictive_campaign_constants()
            catch err
                @debug "Campaign constants could not be loaded; no campaign terms." exception = err
                PredictiveCampaignConstants()
            end
        end
        return _PREDICTIVE_CAMPAIGN_CONSTANTS[]
    end
end

# ── Online corrections to the cost model ─────────────────────────────────────
#
# Every campaign the guard already measures what the model predicted: the
# local slots' slowdown against a pool worker (observed over predicted), the
# pool workers' sample time, and the campaign's wall against its round count.
# The corrections file keeps a BOUNDED running correction to the three model
# parameters those observations bear on, so the next campaign is priced with
# them. What is corrected is always a parameter every candidate shares; nothing
# records how well a plan or a route did, and no campaign is run in order to
# learn. See "Online corrections" in docs/architecture/predictive_routing_r7.md.

const CAMPAIGN_CORRECTIONS_SCHEMA_VERSION = 1

"""
    CampaignCorrectionRules(; step, prior_share, stale_campaigns)

How far and how long a correction may move a parameter away from its prior.
All three are ASSUMED; nothing measures what they should be.

- `step` (`SPACEAGORA_CAMPAIGN_CORRECTION_STEP`, default `0.05`): the largest
  move of a parameter's evidence in one campaign, as a fraction of the prior.
- `prior_share` (`SPACEAGORA_CAMPAIGN_CORRECTION_PRIOR_SHARE`, default `0.25`):
  the weight the prior keeps in the value the planner uses, however much
  evidence has accumulated. At the default a parameter can move at most three
  quarters of the way from its prior to what the campaigns observe.
- `stale_campaigns` (`SPACEAGORA_CAMPAIGN_CORRECTION_STALE_CAMPAIGNS`, default
  `20`): a parameter with no observation in this many campaigns reverts to its
  prior.
"""
struct CampaignCorrectionRules
    step::Float64
    prior_share::Float64
    stale_campaigns::Int
end

function CampaignCorrectionRules(;
    step::Real = _predictive_env_float("SPACEAGORA_CAMPAIGN_CORRECTION_STEP", 0.05),
    prior_share::Real = _predictive_env_float("SPACEAGORA_CAMPAIGN_CORRECTION_PRIOR_SHARE", 0.25),
    stale_campaigns::Integer = _predictive_env_int("SPACEAGORA_CAMPAIGN_CORRECTION_STALE_CAMPAIGNS", 20),
)
    (0.0 < step <= 1.0) || throw(ArgumentError(
        "CampaignCorrectionRules step must be in (0, 1]; got $(step)."))
    (0.0 < prior_share <= 1.0) || throw(ArgumentError(
        "CampaignCorrectionRules prior_share must be in (0, 1]; got $(prior_share)."))
    stale_campaigns >= 1 || throw(ArgumentError(
        "CampaignCorrectionRules stale_campaigns must be >= 1; got $(stale_campaigns)."))
    return CampaignCorrectionRules(Float64(step), Float64(prior_share), Int(stale_campaigns))
end

"""
    CorrectedParameter(evidence, observations, last_observed)

One parameter's running correction: `evidence` is where the observations have
pulled it (it starts at the prior and moves at most `step x prior` per
campaign), `observations` how many campaigns contributed, and `last_observed`
the campaign counter at the latest one.
"""
struct CorrectedParameter
    evidence::Float64
    observations::Int
    last_observed::Int
end

@inline _correction_stale(p::CorrectedParameter, campaign::Integer, rules::CampaignCorrectionRules)::Bool =
    Int(campaign) - p.last_observed > rules.stale_campaigns

"""
    correction_value(p, prior, campaign, rules) -> Float64

The value the planner uses: `prior_share * prior + (1 - prior_share) *
evidence`, or the prior itself when there is no correction or it has gone
stale.
"""
function correction_value(p::Union{Nothing, CorrectedParameter}, prior::Real,
                          campaign::Integer, rules::CampaignCorrectionRules)::Float64
    (p === nothing || _correction_stale(p, campaign, rules)) && return Float64(prior)
    return rules.prior_share * Float64(prior) + (1.0 - rules.prior_share) * p.evidence
end

"""
    correction_observe(p, prior, observed, campaign, rules) -> CorrectedParameter

Fold one campaign's observation of a parameter into its correction. The
evidence moves toward `observed` by at most `rules.step * |prior|`, and never
below zero. A missing or stale correction restarts from the prior.
"""
function correction_observe(p::Union{Nothing, CorrectedParameter}, prior::Real, observed::Real,
                            campaign::Integer, rules::CampaignCorrectionRules)::CorrectedParameter
    base = (p === nothing || _correction_stale(p, campaign, rules)) ?
        CorrectedParameter(Float64(prior), 0, Int(campaign)) : p
    bound = rules.step * abs(Float64(prior))
    evidence = max(0.0, base.evidence + clamp(Float64(observed) - base.evidence, -bound, bound))
    return CorrectedParameter(evidence, base.observations + 1, Int(campaign))
end

# The per-signature sample time has no prior: nothing measures it before the
# campaign runs, and it is not a model verdict -- it only converts the pool's
# start-up latency from seconds into sample units. The first observation sets
# it; later ones move it by at most `step` of its current value.
function _sample_time_observe(p::Union{Nothing, CorrectedParameter}, observed::Real,
                              campaign::Integer, rules::CampaignCorrectionRules)::CorrectedParameter
    (p === nothing || _correction_stale(p, campaign, rules)) &&
        return CorrectedParameter(Float64(observed), 1, Int(campaign))
    bound = rules.step * p.evidence
    evidence = max(0.0, p.evidence + clamp(Float64(observed) - p.evidence, -bound, bound))
    return CorrectedParameter(evidence, p.observations + 1, Int(campaign))
end

@inline _sample_time_value(p::Union{Nothing, CorrectedParameter}, campaign::Integer,
                           rules::CampaignCorrectionRules) =
    (p === nothing || _correction_stale(p, campaign, rules) || !(p.evidence > 0.0)) ?
        nothing : p.evidence

"""
    CampaignCorrections

The corrections file's contents: the machine and code it belongs to, a
campaign counter (the staleness clock), the machine-wide heap-scale and
round-tail corrections, the per-signature sample times, and the plan each
campaign shape last ran (what the leash measures a step from).
"""
mutable struct CampaignCorrections
    fingerprint::String
    code_token::String
    campaigns::Int
    heap_scale::Union{Nothing, CorrectedParameter}
    round_tail::Union{Nothing, CorrectedParameter}
    sample_time_s::Dict{String, CorrectedParameter}
    last_plan::Dict{String, String}
end

CampaignCorrections(fingerprint::AbstractString, code_token::AbstractString) =
    CampaignCorrections(String(fingerprint), String(code_token), 0, nothing, nothing,
                        Dict{String, CorrectedParameter}(), Dict{String, String}())

"""
    campaign_corrections_code_token() -> String

The code version a corrections file is valid for: the RHS calibration store's
code token, so a change to the RHS execution that invalidates the store
invalidates the corrections with it. Sample times and contention measured
against the old code do not describe the new one.
"""
campaign_corrections_code_token()::String = string(SimulationEngine._RHS_CALIB_CODE_TOKEN)

"""
    campaign_corrections_mode() -> Symbol

`SPACEAGORA_CAMPAIGN_CORRECTIONS`: `:on` (the default, and the value of an
unset variable) reads the corrections file and writes it after each campaign
that observed something; `:read` uses it without ever writing it; `:off`
neither reads nor writes, and the planner prices every campaign from the
constants file alone.
"""
function campaign_corrections_mode()::Symbol
    raw = lowercase(strip(get(ENV, "SPACEAGORA_CAMPAIGN_CORRECTIONS", "")))
    isempty(raw) && return :on
    raw in ("1", "on", "true", "yes") && return :on
    raw in ("0", "off", "false", "no") && return :off
    raw == "read" && return :read
    throw(ArgumentError(
        "SPACEAGORA_CAMPAIGN_CORRECTIONS must be \"on\", \"off\" or \"read\"; got \"$(raw)\"."))
end

"""
    campaign_corrections_path() -> String

`SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH` when set; otherwise
`campaign_corrections_<fingerprint>.toml` in the directory of the machine
constants file.
"""
function campaign_corrections_path()::String
    override = strip(get(ENV, "SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH", ""))
    if !isempty(override)
        return normpath(isabspath(override) ? override : joinpath(pwd(), override))
    end
    return joinpath(dirname(ParallelCost.machine_constants_path()),
                    "campaign_corrections_$(ParallelCost.machine_fingerprint()).toml")
end

function _correction_from_toml(d)::Union{Nothing, CorrectedParameter}
    d isa AbstractDict || return nothing
    e = get(d, "evidence", nothing)
    o = get(d, "observations", nothing)
    l = get(d, "last_observed", nothing)
    (e isa Real && isfinite(e) && e >= 0 && o isa Integer && l isa Integer) || return nothing
    return CorrectedParameter(Float64(e), Int(o), Int(l))
end

_correction_to_toml(p::CorrectedParameter) = Dict{String, Any}(
    "evidence" => p.evidence, "observations" => p.observations, "last_observed" => p.last_observed)

"""
    load_campaign_corrections(path; fingerprint, code_token) -> CampaignCorrections

Read a corrections file, or return an empty one when the file is absent,
unreadable, of another schema, or written for another machine fingerprint or
code token. A file that does not describe this machine and this code is
treated as if it did not exist, never reinterpreted.
"""
function load_campaign_corrections(path::AbstractString; fingerprint::AbstractString,
                                   code_token::AbstractString)::CampaignCorrections
    fresh = CampaignCorrections(fingerprint, code_token)
    isfile(path) || return fresh
    parsed = try
        TOML.parsefile(String(path))
    catch
        return fresh
    end
    get(parsed, "schema_version", -1) == CAMPAIGN_CORRECTIONS_SCHEMA_VERSION || return fresh
    string(get(parsed, "fingerprint", "")) == fingerprint || return fresh
    string(get(parsed, "code_token", "")) == code_token || return fresh
    c = fresh
    campaigns = get(parsed, "campaigns", 0)
    c.campaigns = campaigns isa Integer && campaigns >= 0 ? Int(campaigns) : 0
    c.heap_scale = _correction_from_toml(get(parsed, "heap_scale", nothing))
    c.round_tail = _correction_from_toml(get(parsed, "round_tail", nothing))
    times = get(parsed, "sample_time_s", nothing)
    if times isa AbstractDict
        for (k, v) in times
            p = _correction_from_toml(v)
            p === nothing || (c.sample_time_s[String(k)] = p)
        end
    end
    plans = get(parsed, "last_plan", nothing)
    if plans isa AbstractDict
        for (k, v) in plans
            v isa AbstractString && (c.last_plan[String(k)] = String(v))
        end
    end
    return c
end

"""
    save_campaign_corrections(c, path) -> String

Write the corrections file atomically (a temporary file renamed over it).
"""
function save_campaign_corrections(c::CampaignCorrections, path::AbstractString)::String
    payload = Dict{String, Any}(
        "schema_version" => CAMPAIGN_CORRECTIONS_SCHEMA_VERSION,
        "fingerprint" => c.fingerprint,
        "code_token" => c.code_token,
        "campaigns" => c.campaigns,
        "sample_time_s" => Dict{String, Any}(k => _correction_to_toml(v) for (k, v) in c.sample_time_s),
        "last_plan" => Dict{String, Any}(k => v for (k, v) in c.last_plan),
    )
    c.heap_scale === nothing || (payload["heap_scale"] = _correction_to_toml(c.heap_scale))
    c.round_tail === nothing || (payload["round_tail"] = _correction_to_toml(c.round_tail))
    path_s = String(path)
    mkpath(dirname(path_s))
    tmp = path_s * ".tmp"
    open(tmp, "w") do io
        TOML.print(io, payload)
    end
    mv(tmp, path_s; force = true)
    return path_s
end

const _CAMPAIGN_CORRECTIONS = Ref{Any}(nothing)
const _CAMPAIGN_CORRECTIONS_PATH = Ref{String}("")
const _CAMPAIGN_CORRECTIONS_LOCK = ReentrantLock()

"""
    campaign_corrections() -> Union{Nothing, CampaignCorrections}

This process's corrections, loaded from [`campaign_corrections_path`](@ref) on
first use (and again if the path changes); `nothing` when
`SPACEAGORA_CAMPAIGN_CORRECTIONS=off`.
"""
function campaign_corrections()::Union{Nothing, CampaignCorrections}
    campaign_corrections_mode() === :off && return nothing
    path = campaign_corrections_path()
    lock(_CAMPAIGN_CORRECTIONS_LOCK) do
        if _CAMPAIGN_CORRECTIONS[] === nothing || _CAMPAIGN_CORRECTIONS_PATH[] != path
            _CAMPAIGN_CORRECTIONS[] = load_campaign_corrections(
                path; fingerprint = ParallelCost.machine_fingerprint(),
                code_token = campaign_corrections_code_token())
            _CAMPAIGN_CORRECTIONS_PATH[] = path
        end
        return _CAMPAIGN_CORRECTIONS[]
    end
end

"""
    reset_campaign_corrections!()

Forget this process's in-memory corrections so the next campaign re-reads the
file. It does not touch the file; to discard the corrections themselves,
delete the file (see the docs).
"""
function reset_campaign_corrections!()::Nothing
    lock(_CAMPAIGN_CORRECTIONS_LOCK) do
        _CAMPAIGN_CORRECTIONS[] = nothing
        _CAMPAIGN_CORRECTIONS_PATH[] = ""
    end
    return nothing
end

"""
    predictive_cost_terms(corrections, campaign_constants, rules; signature,
                          pool_cold) -> (terms, sample_time_s)

The [`PredictiveCostTerms`](@ref) a campaign is planned with, from the
constants file's priors and this machine's corrections.

- `heap_scale`: prior `1` (the heap model as calibrated), corrected.
- `round_tail`: prior `round_tail` from the constants file,
  corrected; zero, and never corrected, when the file does not have it.
- `startup`: `pool_startup_s / sample_time_s` on a cold pool, when both are
  known; zero otherwise.
"""
function predictive_cost_terms(corrections::Union{Nothing, CampaignCorrections},
                               campaign_constants::PredictiveCampaignConstants,
                               rules::CampaignCorrectionRules;
                               signature::AbstractString, pool_cold::Bool)
    k = corrections === nothing ? 0 : corrections.campaigns
    heap = corrections === nothing ? 1.0 : correction_value(corrections.heap_scale, 1.0, k, rules)
    tail_prior = campaign_constants.round_tail
    tail = if tail_prior === nothing
        0.0
    elseif corrections === nothing
        tail_prior
    else
        correction_value(corrections.round_tail, tail_prior, k, rules)
    end
    t1 = corrections === nothing ? nothing :
        _sample_time_value(get(corrections.sample_time_s, String(signature), nothing), k, rules)
    startup_s = campaign_constants.pool_startup_s
    startup = (pool_cold && t1 !== nothing && startup_s !== nothing) ? startup_s / t1 : 0.0
    return (terms = PredictiveCostTerms(heap, max(0.0, tail), startup), sample_time_s = t1)
end

"""
    predictive_fold_campaign!(corrections, rules, campaign_constants; signature,
                              shape_key, final_plan, worker_sample_s = NaN,
                              heap_scale_observed = NaN,
                              tail_observed = NaN) -> CampaignCorrections

Fold one finished campaign into the corrections. `NaN` means "not observed in
this campaign" and leaves that parameter's correction alone (its staleness
clock keeps running). Parameters that have gone stale are dropped, which is
what reverting to the prior means on disk.

- `worker_sample_s`: mean work of the samples the pool workers ran.
- `heap_scale_observed`: the heap scale that would have predicted the local
  slots' observed slowdown against a pool worker; only from a guarded plan
  whose heap model charged a term.
- `tail_observed`: the round tail implied by a pure process campaign on a
  warm pool; only folded when the constants file has a prior for it.
"""
function predictive_fold_campaign!(c::CampaignCorrections, rules::CampaignCorrectionRules,
                                   campaign_constants::PredictiveCampaignConstants;
                                   signature::AbstractString, shape_key::AbstractString,
                                   final_plan::AbstractString,
                                   worker_sample_s::Real = NaN,
                                   heap_scale_observed::Real = NaN,
                                   tail_observed::Real = NaN)::CampaignCorrections
    c.campaigns += 1
    k = c.campaigns
    sig = String(signature)
    if isfinite(worker_sample_s) && worker_sample_s > 0.0
        c.sample_time_s[sig] = _sample_time_observe(get(c.sample_time_s, sig, nothing),
                                                    worker_sample_s, k, rules)
    end
    if isfinite(heap_scale_observed) && heap_scale_observed > 0.0
        c.heap_scale = correction_observe(c.heap_scale, 1.0, heap_scale_observed, k, rules)
    end
    prior = campaign_constants.round_tail
    if prior !== nothing && isfinite(tail_observed)
        c.round_tail = correction_observe(c.round_tail, prior, tail_observed, k, rules)
    end
    c.last_plan[String(shape_key)] = String(final_plan)
    c.heap_scale !== nothing && _correction_stale(c.heap_scale, k, rules) && (c.heap_scale = nothing)
    c.round_tail !== nothing && _correction_stale(c.round_tail, k, rules) &&
        (c.round_tail = nothing)
    filter!(kv -> !_correction_stale(kv.second, k, rules), c.sample_time_s)
    return c
end

# ── The leash ────────────────────────────────────────────────────────────────

# The leash for a plan without local slots whose samples run on b > 1 threads:
# one step on this shape's budget ladder -- the distinct inner budgets its
# threads-route and serial plans offer, 1 included -- from the previous plan's
# budget toward the chosen one. A previous mixed plan counts as budget 1.
function _predictive_leash_width(planning::PredictivePlanning, chosen::PredictivePlan, prev)
    budget_of(p) = max(1, p.inner_thread_budget)
    cands = filter(p -> p.local_slots == 0 && p.route !== :process, planning.plans)
    ladder = sort!(unique!(vcat(1, [budget_of(p) for p in cands])))
    b_prev = prev[3] == 0 ? prev[4] : 1
    i_prev = something(findfirst(==(b_prev), ladder), 1)
    i_chosen = something(findfirst(==(budget_of(chosen)), ladder), length(ladder))
    abs(i_chosen - i_prev) <= 1 && return (chosen, :within_leash)
    b_step = ladder[i_prev + sign(i_chosen - i_prev)]
    step = findfirst(p -> budget_of(p) == b_step && !p.static_equivalent, cands)
    here = findfirst(p -> budget_of(p) == ladder[i_prev] &&
                          (ladder[i_prev] > 1 || p.static_equivalent), cands)
    from = here === nothing ? planning.best_static : cands[here]
    if step !== nothing && (from === nothing || cands[step].makespan <= from.makespan)
        return (cands[step], :leashed)
    end
    return (from === nothing ? chosen : from, :leash_held)
end

"""
    predictive_plan_key(plan) -> String

`"<route>@w<workers>+l<local_slots>"`, with `"+b<budget>"` appended when the
plan's samples run on more than one thread; the spelling the corrections file
uses for a plan.
"""
predictive_plan_key(p::PredictivePlan)::String =
    "$(p.route)@w$(p.workers)+l$(p.local_slots)" *
    (p.inner_thread_budget > 1 ? "+b$(p.inner_thread_budget)" : "")

function _predictive_parse_plan_key(key::AbstractString)
    m = match(r"^(none|threads|process)@w(\d+)\+l(\d+)(?:\+b(\d+))?$", key)
    m === nothing && return nothing
    budget = m.captures[4] === nothing ? 1 : parse(Int, m.captures[4])
    return (Symbol(m.captures[1]), parse(Int, m.captures[2]), parse(Int, m.captures[3]), budget)
end

"""
    predictive_leash(planning, previous) -> (plan, reason)

The only exploration the planner does. A NON-static plan may differ from the
plan the same campaign shape ran last time (`previous`, a
[`predictive_plan_key`](@ref)) by at most one local slot, and only in the
direction the corrected model favors.

- No previous plan: the chosen plan stands (`:no_previous`). The first
  campaign of a shape is bounded by the margin rule alone, as it always was.
- A static-equivalent chosen plan always stands (`:static`): retreating to the
  static plan is never exploration.
- Within one slot of the previous plan: stands (`:within_leash`).
- Otherwise the plan one slot from the previous one toward the chosen plan is
  taken, if the model prices it no worse than the previous plan
  (`:leashed`); if not, the previous plan is kept (`:leash_held`).

A previous plan on another route or pool width counts as zero local slots on
this one, so every walk away from the static plan starts at one slot. The
cost of a wrong step is therefore bounded by one local slot's share of the
campaign, per campaign.
"""
function predictive_leash(planning::PredictivePlanning,
                          previous::Union{Nothing, AbstractString})
    chosen = planning.chosen
    previous === nothing && return (chosen, :no_previous)
    chosen.static_equivalent && return (chosen, :static)
    prev = _predictive_parse_plan_key(previous)
    prev === nothing && return (chosen, :no_previous)
    chosen_budget = max(1, chosen.inner_thread_budget)
    if chosen.local_slots == 0
        # A threads-route or serial plan with more than one thread per sample:
        # its steps are widths on the budget ladder, not local slots.
        return _predictive_leash_width(planning, chosen, prev)
    end
    same_family = prev[1] === chosen.route && prev[2] == chosen.workers && prev[4] == chosen_budget
    prev_slots = same_family ? prev[3] : 0
    abs(chosen.local_slots - prev_slots) <= 1 && return (chosen, :within_leash)
    step_slots = prev_slots + sign(chosen.local_slots - prev_slots)
    at(L) = L == 0 ?
        findfirst(p -> p.route === chosen.route && p.workers == chosen.workers &&
                       p.local_slots == 0 && p.static_equivalent, planning.plans) :
        findfirst(p -> p.route === chosen.route && p.workers == chosen.workers &&
                       p.local_slots == L && max(1, p.inner_thread_budget) == chosen_budget,
                  planning.plans)
    step_idx = at(step_slots)
    prev_idx = at(prev_slots)
    from = prev_idx === nothing ? planning.best_static : planning.plans[prev_idx]
    if step_idx !== nothing && (from === nothing || planning.plans[step_idx].makespan <= from.makespan)
        return (planning.plans[step_idx], :leashed)
    end
    return (from === nothing ? chosen : from, :leash_held)
end

function _predictive_plan_line(p::PredictivePlan)::String
    return "$(p.route)@w$(p.workers)+l$(p.local_slots)" *
           (p.inner_thread_budget > 1 ? " budget=$(p.inner_thread_budget)" : "") *
           " makespan=$(round(p.makespan; digits=3)) " *
           "consumers=$(p.consumers) s_worker=$(round(p.worker_slowdown; digits=3)) " *
           "s_heap=$(round(p.heap_slowdown; digits=3))" * (p.static_equivalent ? " static" : "")
end
