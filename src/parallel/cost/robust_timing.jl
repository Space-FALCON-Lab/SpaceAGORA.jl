# Measurement primitives for the parts of the cost model that must be measured.
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
# here. (3) is not, and is handled elsewhere by abstaining when the margin is
# inside the model's residual spread.
#
# The estimators below are chosen against those failure modes specifically:
# minimum against (2), paired interleaving against (1), and reference-kernel
# normalisation against both, by turning an absolute duration into a
# dimensionless ratio.

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
    free candidate, and the cost model would rank it accordingly.

    Callers timing a pure kernel must give it runtime-opaque arguments (a
    `Vector` filled at run time, a value read from a `Ref`). Compare against
    [`timing_noise_floor`](@ref) to detect the failure when in doubt; the real
    probe path is safe by construction, since it evaluates effectors against
    live simulation state.
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
        ta = timed_min(fa; k = 1, warmup = i == 1 ? 2 : 0)
        tb = timed_min(fb; k = 1, warmup = i == 1 ? 2 : 0)
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

# Reference kernel: fixed, deterministic, L1-resident FMA throughput.
#
# Serves two purposes, which is why it is worth pinning down precisely.
#
# As a NORMALISER it converts absolute nanoseconds into a dimensionless ratio,
# so a probe measured on a throttled machine is comparable to constants
# calibrated on a cool one -- both numerator and denominator move together.
#
# As a STALENESS CANARY (see MachineConstants.reference_ns) it is re-measured at
# the start of every run and compared against its value at calibration time. If
# it has moved, the cached constants no longer describe this machine -- a
# different cgroup CPU quota, a changed governor, a co-tenant, thermal state --
# and the predictor abstains. It costs microseconds, and it is the only
# freshness check the model needs, because everything else volatile is
# re-derived per run rather than cached.
#
# TWO reference kernels, not one, and the reason is measured rather than assumed.
#
# A single contiguous-FMA reference was tried first. Against a probe kernel
# (strided gather plus a transcendental) under ten memory-bound co-tenants at
# load average 4.6 on this 12-core box:
#
#     raw probe cost        689.8 ns quiet -> 790.6 ns loaded   +14.6%
#     normalised by FMA ref  36812   quiet ->  39220   loaded    +6.5%
#
# Normalisation cut the contention bias by 2.2x but did not remove it, because
# normalisation only cancels the component of a slowdown that both kernels
# share. The probe is memory-sensitive and the FMA reference is not, so they
# degrade at different rates and the ratio drifts.
#
# The fix is to normalise each cost term against a reference with a matching
# resource profile: `ns_per_simd_lane` is arithmetic-bound and normalises
# against the FMA kernel, while `ns_per_coeff_touch` is stride-bound -- the
# coefficient matrices are column-major and the harmonics kernel walks a row --
# and normalises against the memory kernel.
#
# The FMA kernel doubles as the per-run staleness canary because it is the cheap
# one (tens of nanoseconds per call); the memory kernel costs microseconds,
# which is free during a one-shot calibration but not on every run.

# Contiguous SIMD accumulate over L1-resident data. Mirrors the harmonics batch
# inner loop, which is the dominant arithmetic term in the model.
const _REFERENCE_LANES = 1024

function _reference_kernel(buf::Vector{Float64})::Float64
    acc = 0.0
    @inbounds @simd for i in eachindex(buf)
        acc = muladd(buf[i], 1.0000001, acc)
    end
    return acc
end

# One access per 64-byte cache line over a buffer that overflows L1 but stays in
# L2, so the cost is dominated by L1-miss service rather than by arithmetic.
# Mirrors the coefficient-table walk, whose stride guarantees a new line per read.
const _REFERENCE_MEM_DOUBLES = 8192   # 64 KiB
const _REFERENCE_MEM_STRIDE = 8       # 64 B, one cache line
const _REFERENCE_MEM_TOUCHES = _REFERENCE_MEM_DOUBLES ÷ _REFERENCE_MEM_STRIDE

function _reference_mem_kernel(buf::Vector{Float64})::Float64
    acc = 0.0
    @inbounds for i in 1:_REFERENCE_MEM_STRIDE:length(buf)
        acc += buf[i]
    end
    return acc
end

function _make_reference_buffer(n::Int)::Vector{Float64}
    buf = Vector{Float64}(undef, n)
    @inbounds for i in eachindex(buf)
        buf[i] = 1.0 + i * 1e-9
    end
    return buf
end

"""
    reference_kernel_ns(; k = 15) -> Float64

Arithmetic reference: nanoseconds per lane of the L1-resident FMA kernel.

Per lane rather than per call so a stored value stays comparable if
`_REFERENCE_LANES` changes. This is the cheap reading, and the one re-measured
at the start of every run as the staleness canary.
"""
function reference_kernel_ns(; k::Int = 15)::Float64
    buf = _make_reference_buffer(_REFERENCE_LANES)
    return timed_min(() -> _reference_kernel(buf); k = k, warmup = 3) / _REFERENCE_LANES
end

"""
    reference_memory_kernel_ns(; k = 15) -> Float64

Memory reference: nanoseconds per cache-line touch of the L2-resident strided
kernel.

Used to normalise the stride-bound cost terms, which the arithmetic reference
tracks poorly under memory contention -- see the measurement in this file's
header.
"""
function reference_memory_kernel_ns(; k::Int = 15)::Float64
    buf = _make_reference_buffer(_REFERENCE_MEM_DOUBLES)
    return timed_min(() -> _reference_mem_kernel(buf); k = k, warmup = 3) / _REFERENCE_MEM_TOUCHES
end
