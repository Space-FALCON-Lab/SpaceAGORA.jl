# Work-count and machine-constant types for the analytic cost model.
#
# The model is T = sum_k w_k * n_k, split so that the two halves have opposite
# volatility:
#
#   n_k (WorkCounts)       -- countable in closed form from the configuration
#                             before the solve runs. Machine-independent, exact,
#                             zero measurement.
#   w_k (MachineConstants) -- rates in ns per unit of n_k. Measured once per
#                             machine by scripts/calibrate_machine.jl, cached,
#                             and reused by every run on that machine.
#
# That split is what makes a per-solve decision affordable: the part that varies
# per solve needs no measurement, and the part that needs measurement does not
# vary per solve. A routing decision becomes a length-7 dot product instead of a
# timed sweep over candidate plans.

"""
    WorkCounts

Countable work in one full evaluation of a constellation's dynamic effectors,
before any decision about how to split it across threads.

The counts are *per pass over the whole constellation*, not per satellite or per
worker; the predictor applies the split. Fields are `Float64` rather than `Int`
because probe-derived terms are measured, and because the predictor divides them
by worker counts.

Units differ per field and the distinction is the whole model, so they are
stated explicitly:

- `simd_terms`      PER SATELLITE. Vectorized inner-loop iterations, counted as
                    work units per satellite: each `@turbo` loop in the kernel
                    performs one unit for every satellite in its batch, so the
                    loop count and the per-satellite unit count coincide. Scales
                    with the satellites a worker holds, hence shrinks under a
                    split.
- `scalar_items`    PER SATELLITE. Setup/teardown that does not vectorize.
- `coeff_touches`   PER PASS, *not* per satellite. Distinct coefficient-table
                    reads needed to sweep the model once. This is the term that
                    separates the routing candidates: the SIMD batch kernel
                    loads each coefficient once and broadcasts it across its
                    whole satellite batch, so it pays this once per worker,
                    while a per-satellite kernel re-walks the table for every
                    satellite and pays it N times. Everything else in the model
                    is symmetric between the candidates; this is what makes them
                    differ, and what makes the crossover predictable.
- `simd_workspace_bytes_per_sat` SIZE ARGUMENT. Scratch the vectorized kernel
                    holds per satellite in its batch; multiplied by batch width
                    it gives the footprint that decides the SIMD lane rate. For
                    harmonics the A workspace is batch x (L+3) x (M+2), which at
                    batch 1024 and L=20 is 4.1 MB -- nowhere near cache, and the
                    reason a rate indexed by batch width alone under-predicted
                    the serial case by ~2x.
- `queue_nodes`     TOTAL. Flat-queue work items (satellites x effectors),
                    driving per-node dispatch bookkeeping.
- `coeff_table_bytes` SIZE ARGUMENT, not a count. The footprint strided over
                    while accumulating `coeff_touches`, used to look up the
                    right rate on `MachineConstants.coeff_touch` -- a touch into
                    an L1-resident table and one into an L3-resident table
                    differ by 3.4x and the model has to tell them apart. Summed
                    across effectors that contribute touches; zero when none do.
- `probe_ns`        PER SATELLITE. Summed measured cost in ns of effectors the
                    model cannot count analytically. Zero when
                    every effector is either known or has declared its terms.
- `unknown_effectors` how many effectors contributed `probe_ns` rather than
                    analytic counts. Purely informational for telemetry; the
                    probe makes them predictable, so a nonzero value is not by
                    itself a reason to abstain.
- `in_domain`       false when the workload contains a term the linear model
                    provably cannot represent -- today, a density model holding
                    a process-wide native lock, whose cost is superlinear in
                    concurrency. The predictor refuses to answer rather than
                    answering confidently and wrongly; callers fall back to the
                    heuristic route.
"""
Base.@kwdef struct WorkCounts
    simd_terms::Float64 = 0.0
    scalar_items::Float64 = 0.0
    coeff_touches::Float64 = 0.0
    coeff_table_bytes::Float64 = 0.0
    simd_workspace_bytes_per_sat::Float64 = 0.0
    queue_nodes::Float64 = 0.0
    probe_ns::Float64 = 0.0
    unknown_effectors::Int = 0
    in_domain::Bool = true
end

function Base.:+(a::WorkCounts, b::WorkCounts)::WorkCounts
    return WorkCounts(
        simd_terms = a.simd_terms + b.simd_terms,
        scalar_items = a.scalar_items + b.scalar_items,
        coeff_touches = a.coeff_touches + b.coeff_touches,
        coeff_table_bytes = a.coeff_table_bytes + b.coeff_table_bytes,
        simd_workspace_bytes_per_sat = a.simd_workspace_bytes_per_sat + b.simd_workspace_bytes_per_sat,
        queue_nodes = a.queue_nodes + b.queue_nodes,
        probe_ns = a.probe_ns + b.probe_ns,
        unknown_effectors = a.unknown_effectors + b.unknown_effectors,
        in_domain = a.in_domain && b.in_domain,
    )
end

Base.zero(::Type{WorkCounts})::WorkCounts = WorkCounts()

"""
    RateCurve

A machine rate that depends on the size of the thing it operates on, stored as
measured knots and interpolated in log2 of that size.

Scalars were the original design and the measurements rejected them. Both of
the model's discriminating rates move by more than a factor of two across the
sizes real workloads produce, on this 12-core box (L1d 48K, L2 1M, L3 32M):

    coefficient touch, by table footprint
        3.4 KiB  0.103 ns      315 KiB  0.240 ns
       20.3 KiB  0.125 ns     1256 KiB  0.299 ns
       79.7 KiB  0.174 ns    20025 KiB  0.346 ns          -> 3.4x spread

    SIMD lane-term, by batch width
        B=16     0.228 ns      B=4096   0.0227 ns
        B=256    0.0259 ns     B=16384  0.0435 ns
        B=1024   0.0227 ns     B=262144 0.0551 ns

A single number for either would be wrong by 2-3x at the extremes, and a 2x
error in the coefficient-touch rate alone shifts the predicted `satellite_batch`
cost by 49% -- enough to invert the routing decision it exists to make.

The shape is physical, not noise: flat while the working set sits in L1, rising
through the L1-to-L2 and L2-to-L3 transitions, flattening again in L3. The small
batch-width end is a different effect -- per-pass loop overhead that cannot
amortize over a short vector -- and is absorbed into the same curve rather than
modelled as a separate term, because what is measured at each knot is the
*effective* per-lane cost with that overhead already amortized.

Knots are interpolated linearly in log2 of the size and clamped outside the
measured range, so a workload larger or smaller than anything calibrated gets
the nearest measured rate rather than an extrapolation off the end of a curve
whose shape is only known where it was sampled.
"""
struct RateCurve
    log2_size::Vector{Float64}
    ns::Vector{Float64}

    function RateCurve(log2_size::Vector{Float64}, ns::Vector{Float64})
        length(log2_size) == length(ns) ||
            throw(ArgumentError("RateCurve knots and values must have equal length."))
        isempty(log2_size) && throw(ArgumentError("RateCurve needs at least one knot."))
        issorted(log2_size) ||
            throw(ArgumentError("RateCurve knots must be sorted ascending in log2 size."))
        return new(log2_size, ns)
    end
end

"""
    rate_at(curve, size) -> Float64

The rate at `size`, linearly interpolated in log2 and clamped at both ends.
"""
function rate_at(curve::RateCurve, size::Real)::Float64
    n = length(curve.log2_size)
    n == 1 && return curve.ns[1]
    x = log2(max(1.0, Float64(size)))
    x <= curve.log2_size[1] && return curve.ns[1]
    x >= curve.log2_size[n] && return curve.ns[n]
    @inbounds for i in 2:n
        if x <= curve.log2_size[i]
            x0, x1 = curve.log2_size[i - 1], curve.log2_size[i]
            y0, y1 = curve.ns[i - 1], curve.ns[i]
            w = (x - x0) / max(1e-12, x1 - x0)
            return y0 + w * (y1 - y0)
        end
    end
    return curve.ns[n]
end

"""
    MachineConstants

Per-machine rates, in nanoseconds per unit of the corresponding `WorkCounts`
field, plus the reference-kernel reading that validates them.

Produced by `scripts/calibrate_machine.jl` and cached in a fingerprinted TOML.
Never fitted from the benchmark case catalog: the catalog is held out so that a
predicted-versus-actual comparison over it is an out-of-sample test rather than
a restatement of the fit.

- `simd_lane`          [`RateCurve`] cost of one vectorized inner-loop iteration
                       per satellite, as a function of the satellite batch width
                       the loop runs over.
- `coeff_touch`        [`RateCurve`] cost of one coefficient-table read, as a
                       function of the table footprint being strided over.
                       Calibrated rather than derived from element size: the
                       coefficient matrices are column-major and the kernel walks
                       a row, so each read strides by `(L+2)*8` bytes onto its own
                       cache line, and what it costs is set by which cache level
                       holds the table -- not by `sizeof(Float64)`.
- `parallel_speedup`   [`RateCurve`] ACHIEVED speedup against worker count, not
                       assumed linear scaling. Measured here: 0.99 efficiency at
                       one worker falling to 0.26 at twelve, saturating near 3x,
                       and the real harmonics kernel shows the same ~2.3x
                       ceiling. Dividing work by the worker count over-credits
                       every wide plan, and no correction to the other constants
                       can compensate for a term that is simply absent.
- `ns_per_scalar_item` cost of one per-satellite scalar unit.
- `ns_per_queue_node`  per-node flat-queue bookkeeping.
- `dispatch_pool_ns_base` / `dispatch_pool_ns_per_worker`
- `dispatch_batch_ns_base` / `dispatch_batch_ns_per_worker`
                       affine models of one parallel dispatch and join, ONE PER
                       MECHANISM. The routes do not share a dispatcher: the flat
                       routes go through the persistent channel pool
                       (`threaded_foreach_worker_persistent`) while
                       `satellite_batch` goes through `Polyester.@batch`, and
                       they differ by enough to invert a decision. A single
                       constant measured from a third mechanism (`Threads.@spawn`)
                       was the first version of this, and it drove the predictor
                       to 25% decision accuracy against the real kernel: it
                       systematically preferred narrow allotments because it
                       could not see what a wide one really cost on the route it
                       was actually taken on.
- `ns_per_atomic`      one contended atomic read-modify-write, which is what the
                       dynamic scheduler pays per chunk.
- `reference_fma_ns` / `reference_mem_ns`
                       the two reference kernels' costs at calibration time, in
                       ns per lane and ns per cache-line touch. Two rather than
                       one because normalisation only cancels the part of a
                       contention slowdown that the reference and the measured
                       kernel share; an arithmetic reference tracks a
                       stride-bound term poorly (measured: 14.6% raw drift under
                       memory load, still 6.5% after FMA-only normalisation).
                       Arithmetic terms normalise against the first,
                       stride-bound terms against the second.

                       `reference_fma_ns` doubles as the staleness canary:
                       re-measured at the start of every run and compared, and
                       if it has moved beyond tolerance these constants no
                       longer describe this machine (different cgroup quota,
                       changed governor, a co-tenant, thermal state) and the
                       predictor abstains. It is the cheap one, at tens of
                       nanoseconds per call, which is why it and not the memory
                       kernel is on the per-run path. It is the only freshness
                       check the model needs, because everything else volatile
                       is re-derived per run rather than cached.
- `fingerprint`        machine identity these constants were measured on.
- `schema_version`     bumped when term semantics change, so stale files are
                       rejected rather than reinterpreted.
"""
Base.@kwdef struct MachineConstants
    simd_lane::RateCurve
    coeff_touch::RateCurve
    parallel_speedup::RateCurve
    ns_per_scalar_item::Float64
    ns_per_queue_node::Float64
    dispatch_pool_ns_base::Float64
    dispatch_pool_ns_per_worker::Float64
    dispatch_batch_ns_base::Float64
    dispatch_batch_ns_per_worker::Float64
    ns_per_atomic::Float64
    reference_fma_ns::Float64
    reference_mem_ns::Float64
    fingerprint::String = ""
    schema_version::Int = 1
end

"""
    effector_cost_terms(model) -> Union{Nothing, WorkCounts}

Declare the analytic work an effector performs per satellite per RHS pass, so
the cost model can predict its routing without timing it.

Returns `nothing` by default, which means "measure me": the model probes the
effector once at run start and uses the measured per-satellite cost instead.
That fallback is why this is optional -- a user-defined effector works with the
cost model without implementing anything, it just costs one extra call per
effector at setup rather than zero.

Declaring is still worth it for an effector whose cost varies strongly with its
own parameters, because a probe measures only the configuration in front of it
while a declaration extrapolates. A declaration is validated against a single
probe at run start and warns if the two disagree by more than a factor, so a
wrong declaration is caught without being required to be right.

Extend it the same way as [`environment_requirements`](@ref):

```julia
SpaceAGORA.effector_cost_terms(m::MyEffector) =
    WorkCounts(simd_terms = m.n_modes, coeff_touches = m.n_modes)
```

Counts are per satellite per pass; the model multiplies by satellite count and
divides by the worker split.
"""
@inline effector_cost_terms(::Any)::Union{Nothing, WorkCounts} = nothing
