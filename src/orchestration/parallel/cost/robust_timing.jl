# Measurement primitives for the timings the machine calibration and the RHS
# calibration take.
#
# Wall-clock timing is contaminated by three distinguishable things, and only one
# of them is a real epistemic problem:
#
#   1. Multiplicative, slow-drifting -- thermal throttling, frequency scaling.
#      Scales everything in a measurement window by roughly one factor.
#   2. Additive, one-sided, bursty -- preemption, page faults, GC. Only ever
#      makes a sample larger, never smaller.
#   3. Load-dependent ranking changes -- a co-tenant evicting cache lines can
#      genuinely reorder which candidate is faster.
#
# (1) and (2) are measurement problems with measurement answers, implemented
# here. (3) is not; the paired comparisons below answer it by declining to name
# a winner when the margin is inside the noise.
#
# The estimators below are chosen against those failure modes specifically:
# minimum against (2), paired interleaving against (1).

"""
    timed_min(f; k = 9, warmup = 2, target_window_ns = 20_000) -> Float64

Nanoseconds per call of `f`, estimated as the **minimum** over `k` timing
windows.

The minimum, not the mean. Interference is strictly additive and one-sided --
a preempted sample is longer, never shorter -- so the minimum is the maximum
likelihood estimate of the uncontended cost while the mean estimates "this plus
whatever else the machine was doing". The existing RHS calibration sweep takes a
mean over 10 calls, where a single 2 ms preemption anywhere in that window moves
a candidate by 200 microseconds per call and silently changes which plan wins.

Each window repeats `f` enough times for the window to exceed
`target_window_ns`, so `time_ns` call overhead (tens of ns) stays below a
percent of the measurement. The repeat count is discovered once by doubling and
reused for every window, so all `k` windows are directly comparable.

`f` must return a value; it is accumulated into a sink the compiler cannot see
through, so that a pure kernel is not dead-code eliminated.

!!! warning "Constant folding"
    The sink defeats dead-code elimination, not constant folding. A pure `f`
    whose arguments are compile-time literals -- `() -> kernel(512)` -- can be
    folded to a constant before the loop ever runs, and this then faithfully
    reports the cost of returning that constant: a few tenths of a nanosecond.
    A measurement that reads as nearly free is indistinguishable from a genuinely
    free candidate, and a caller choosing by it would rank it accordingly.

    Callers timing a pure kernel must give it runtime-opaque arguments (a
    `Vector` filled at run time, a value read from a `Ref`). Compare against
    [`timing_noise_floor`](@ref) to detect the failure when in doubt.
"""
function timed_min(f::F; k::Int = 9, warmup::Int = 2, target_window_ns::Int = 20_000)::Float64 where {F}
    k = max(1, k)
    sink = 0.0

    for _ in 1:max(0, warmup)
        sink += _as_sink(f())
    end

    inner = 1
    elapsed = 0
    while inner < 1 << 30
        t0 = time_ns()
        for _ in 1:inner
            sink += _as_sink(f())
        end
        elapsed = Int(time_ns() - t0)
        elapsed >= target_window_ns && break
        inner *= 2
    end

    best = Float64(elapsed) / inner
    for _ in 2:k
        t0 = time_ns()
        for _ in 1:inner
            sink += _as_sink(f())
        end
        sample = Float64(Int(time_ns() - t0)) / inner
        sample < best && (best = sample)
    end

    _consume_sink(sink)
    return best
end

@inline _as_sink(x::Real)::Float64 = Float64(x)
@inline _as_sink(x::Tuple)::Float64 = isempty(x) ? 0.0 : _as_sink(first(x))
@inline _as_sink(x::AbstractArray)::Float64 = isempty(x) ? 0.0 : _as_sink(first(x))
@inline _as_sink(::Any)::Float64 = 0.0

# Opaque sink: keeps the compiler from eliding the work whose cost we are
# measuring, without perturbing the measurement itself.
const _TIMING_SINK = Ref{Float64}(0.0)
@inline function _consume_sink(x::Float64)::Nothing
    if !isfinite(x) || x == 1.2345678e-300
        _TIMING_SINK[] = x
    end
    return nothing
end

"""
    timing_noise_floor(; k = 9) -> Float64

Nanoseconds per call of the measurement harness itself, timing a closure that
does nothing.

The floor a real measurement must clear. A [`timed_min`](@ref) result at or
below this is not a fast kernel, it is an absent one -- constant-folded away, or
a closure that never did the work. Distinguishing "cheap" from "elided" is not
possible from the duration alone, which is why the floor is measured rather than
assumed: it moves with the machine exactly as the measurements it validates do.
"""
function timing_noise_floor(; k::Int = 9)::Float64
    sink = Ref(1.0)
    return timed_min(() -> sink[]; k = k, warmup = 3)
end

"""
    PairedComparison

Outcome of an interleaved A-versus-B comparison: how many pairs each side won,
the median of the per-pair ratios, and whether the win count clears a two-sided
sign test at the requested level.
"""
struct PairedComparison
    wins_a::Int
    wins_b::Int
    median_ratio::Float64
    significant::Bool
end

"""
    paired_compare(fa, fb; pairs = 15, alpha = 0.05) -> PairedComparison

Compare two candidates by measuring them **interleaved** -- A, B, A, B, ... --
and reducing over per-pair differences rather than over separate batches.

Interleaving is what makes this robust to drift. Thermal throttling and
background load change slowly relative to one adjacent pair, so they scale both
halves of a pair by nearly the same factor and cancel in the ratio. Measuring
`AAAA...` then `BBBB...` gives drift a full batch to move between the groups,
which shows up as a difference that is not there.

The decision is a **sign test** on pair wins, not a difference of means: it
needs no distributional assumption, is immune to a single wild pair, and yields
an explicit significance verdict. That last part matters most -- the existing
sweep takes `argmin` over noisy means with no notion of "too close to call", so
it always picks a winner even when the evidence does not support one. A caller
that gets `significant = false` should keep its current choice rather than
switch on noise.
"""
function paired_compare(fa::FA, fb::FB; pairs::Int = 15, alpha::Float64 = 0.05)::PairedComparison where {FA, FB}
    pairs = max(1, pairs)
    ratios = Vector{Float64}(undef, pairs)
    wins_a = 0
    wins_b = 0

    for i in 1:pairs
        # ORDER ALTERNATES WITHIN THE PAIR, and this is not a refinement.
        #
        # A fixed a-then-b order lets any first-position effect register as a
        # difference between the alternatives: a GC landing on the first call,
        # caches or branch predictors warmed by it for the second, a frequency
        # step between them. Measured here on a whole-solve comparison, fixed
        # ordering reported the identification path 59% SLOWER at 15-0 and
        # significant; alternating the order on the same workload in the same
        # session reported it 40% FASTER at 15-0 and significant. The effect was
        # larger than the effect under test and had the opposite sign.
        #
        # `scripts/calibration/paired_profile_probe.jl` and `paired_campaign_probe.jl` both
        # alternate for exactly this reason and say so in their headers; this
        # function did not, so a caller reaching for the library instead of the
        # scripts got the bias silently.
        ta, tb = if isodd(i)
            a = timed_min(fa; k = 1, warmup = i == 1 ? 2 : 0)
            b = timed_min(fb; k = 1, warmup = i == 1 ? 2 : 0)
            (a, b)
        else
            b = timed_min(fb; k = 1, warmup = 0)
            a = timed_min(fa; k = 1, warmup = 0)
            (a, b)
        end
        ratios[i] = tb > 0.0 ? ta / tb : Inf
        if ta < tb
            wins_a += 1
        elseif tb < ta
            wins_b += 1
        end
    end

    sort!(ratios)
    median_ratio = pairs % 2 == 1 ?
        ratios[(pairs + 1) ÷ 2] :
        0.5 * (ratios[pairs ÷ 2] + ratios[pairs ÷ 2 + 1])

    decided = wins_a + wins_b
    winner_wins = max(wins_a, wins_b)
    significant = decided > 0 && _sign_test_two_sided(winner_wins, decided) <= alpha

    return PairedComparison(wins_a, wins_b, median_ratio, significant)
end

# Two-sided sign test: P(as extreme as `wins` out of `n` under Binomial(n, 1/2)).
function _sign_test_two_sided(wins::Int, n::Int)::Float64
    n <= 0 && return 1.0
    wins = max(wins, n - wins)
    tail = 0.0
    for i in wins:n
        tail += binomial(big(n), big(i)) / big(2)^n
    end
    return min(1.0, 2.0 * Float64(tail))
end

# ── Streaming paired identification ───────────────────────────────────────────
#
# `paired_compare` above is a driver: it calls both alternatives itself, in a
# loop it controls. In-run identification needs the inverse. The solver decides
# when the RHS is evaluated; the trial can only say which arm the next call
# should use and be told afterwards what it cost.
#
# Same two properties that make the paired probes trustworthy, kept:
#
#   INTERLEAVING. Arms are visited once per round, and the visiting order
#   rotates between rounds, so anything drifting slower than a round -- thermal
#   state, frequency, a co-tenant, the allocator settling -- scales every arm in
#   a round together and cancels in the within-round comparison. Blocking the
#   arms instead (all of A, then all of B) lets drift masquerade as effect; this
#   repo has one measurement that read 11.1% block-ordered and 1.7% paired.
#
#   A SIGNIFICANCE VERDICT. The reduction is a two-sided sign test on round
#   wins, not a difference of medians, so the trial can return "too close to
#   call" -- which is the right answer for most workloads and the one a median
#   cannot express. An undecided trial keeps the incumbent, which is what makes
#   this safe to run by default: it can only move a decision it can defend.
#
# Rounds, not pairs, because the ladder has more than two arms: a round is one
# observation of every arm, and it is the unit the sign test counts.

"""
    StreamingPairedTrial(n_arms; rounds = 9, alpha = 0.05, warmup_rounds = 1)

An interleaved A/B/.../N trial whose observations arrive one at a time.

`arm 1` is the incumbent: the verdict is expressed as whether some other arm
beats it, so a trial that never reaches significance leaves the caller doing
exactly what it was doing.

`warmup_rounds` observations of each arm are discarded before any are counted.
Switching execution width re-warms caches and can re-enter the dispatch pool, so
the first observation after a switch measures the switch.
"""
mutable struct StreamingPairedTrial
    n_arms::Int
    rounds_target::Int
    warmup_rounds::Int
    alpha::Float64
    obs::Vector{Vector{Float64}}
    order::Vector{Int}
    pos::Int
    round::Int
end

function StreamingPairedTrial(n_arms::Int; rounds::Int = 9, alpha::Float64 = 0.05,
                              warmup_rounds::Int = 1)
    n_arms >= 2 || throw(ArgumentError("A trial needs at least two arms; got $n_arms."))
    return StreamingPairedTrial(
        n_arms, max(1, rounds), max(0, warmup_rounds), alpha,
        [Float64[] for _ in 1:n_arms], collect(1:n_arms), 1, 1)
end

"""
    next_arm(trial) -> Int

Which arm the next observation should use. Cheap and side-effect free, so it can
be called on the hot path before every evaluation.
"""
@inline next_arm(t::StreamingPairedTrial)::Int = @inbounds t.order[t.pos]

"""
    trial_active(trial) -> Bool

Whether the trial still wants observations. False once it has its rounds, at
which point the caller should stop switching arms and use [`trial_verdict`](@ref).
"""
@inline function trial_active(t::StreamingPairedTrial)::Bool
    return t.round <= t.rounds_target + t.warmup_rounds
end

"""
    observe!(trial, elapsed_ns)

Record the cost of the arm [`next_arm`](@ref) named, and advance the schedule.

The arm order is rotated at each round boundary so a systematic first-position
or last-position bias -- a cold cache at the start of a round, a scheduler
preemption at the end -- spreads across arms instead of accumulating on one.
"""
function observe!(t::StreamingPairedTrial, elapsed_ns::Real)::Nothing
    arm = next_arm(t)
    ns = Float64(elapsed_ns)
    if t.round > t.warmup_rounds && isfinite(ns) && ns > 0.0
        push!(t.obs[arm], ns)
    end
    t.pos += 1
    if t.pos > t.n_arms
        t.pos = 1
        t.round += 1
        # Rotate rather than shuffle: deterministic, so a trial is reproducible,
        # and over n_arms rounds every arm has occupied every position once.
        push!(t.order, popfirst!(t.order))
    end
    return nothing
end

"""
    trial_rounds(trial) -> Int

Counted rounds, excluding warm-up.
"""
@inline function trial_rounds(t::StreamingPairedTrial)::Int
    return minimum(length(o) for o in t.obs)
end

"""
    trial_verdict(trial) -> (arm, significant, median_ratio, rounds)

The arm that won, whether the win survives a two-sided sign test against the
incumbent, and the median of the per-round ratios `incumbent / winner`.

`median_ratio > 1` means the winner is that many times faster than arm 1. A
verdict with `significant = false` must be read as "keep the incumbent": it is
returned with the winning arm still named so a caller can log what it nearly
chose, not so it can act on it.
"""
function trial_verdict(t::StreamingPairedTrial)
    n = trial_rounds(t)
    n <= 0 && return (arm = 1, significant = false, median_ratio = 1.0, rounds = 0)

    # Best arm by median, which is insensitive to one wild round in a way a mean
    # is not.
    med(v) = (s = sort(v); length(s) % 2 == 1 ? s[(length(s) + 1) ÷ 2] :
              0.5 * (s[length(s) ÷ 2] + s[length(s) ÷ 2 + 1]))
    medians = [med(t.obs[a][1:n]) for a in 1:t.n_arms]
    best = argmin(medians)
    best == 1 && return (arm = 1, significant = false, median_ratio = 1.0, rounds = n)

    wins = 0
    decided = 0
    ratios = Float64[]
    @inbounds for r in 1:n
        a = t.obs[1][r]
        b = t.obs[best][r]
        (isfinite(a) && isfinite(b) && b > 0.0) || continue
        push!(ratios, a / b)
        if b < a
            wins += 1; decided += 1
        elseif a < b
            decided += 1
        end
    end
    isempty(ratios) && return (arm = 1, significant = false, median_ratio = 1.0, rounds = n)
    significant = decided > 0 && _sign_test_two_sided(max(wins, decided - wins), decided) <= t.alpha
    # Significance alone is not a reason to switch: the winner must be the arm
    # that actually won the rounds, not merely the one with the better median.
    significant &= wins > decided - wins
    return (arm = best, significant = significant, median_ratio = med(ratios), rounds = n)
end

"""
    trial_speedups(trial) -> Vector{Float64}

Per-arm speedup against arm 1, from each arm's minimum observation.

The minimum, not the median: a trial runs inside a live solve, so every
observation carries whatever interruption happened to land on it, and the
minimum is the least contaminated estimate of what the arm costs when nothing
else is happening. Feed these to `_fit_usl` with the arms' widths to identify
the scalability parameters of THIS workload, which is the thing the machine-only
prediction could not do.
"""
function trial_speedups(t::StreamingPairedTrial)::Vector{Float64}
    base = isempty(t.obs[1]) ? NaN : minimum(t.obs[1])
    return [(isempty(t.obs[a]) || !isfinite(base) || base <= 0.0) ? 1.0 :
            max(1e-6, base / minimum(t.obs[a])) for a in 1:t.n_arms]
end
