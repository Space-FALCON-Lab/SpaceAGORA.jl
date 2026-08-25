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
- `queue_nodes`     TOTAL. Flat-queue work items (satellites x effectors),
                    driving per-node dispatch bookkeeping.
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
        queue_nodes = a.queue_nodes + b.queue_nodes,
        probe_ns = a.probe_ns + b.probe_ns,
        unknown_effectors = a.unknown_effectors + b.unknown_effectors,
        in_domain = a.in_domain && b.in_domain,
    )
end

Base.zero(::Type{WorkCounts})::WorkCounts = WorkCounts()

"""
    MachineConstants

Per-machine rates, in nanoseconds per unit of the corresponding `WorkCounts`
field, plus the reference-kernel reading that validates them.

Produced by `scripts/calibrate_machine.jl` and cached in a fingerprinted TOML.
Never fitted from the benchmark case catalog: the catalog is held out so that a
predicted-versus-actual comparison over it is an out-of-sample test rather than
a restatement of the fit.

- `ns_per_simd_lane`   cost of one vectorized inner-loop iteration per satellite
                       in the batch.
- `ns_per_scalar_item` cost of one per-satellite scalar unit.
- `ns_per_coeff_touch` cost of one coefficient-table read, including whatever
                       cache-line traffic the access stride implies. Calibrated
                       rather than derived from element size on purpose: the
                       coefficient matrices are column-major and the kernel
                       walks a row, so the bytes actually moved per term depend
                       on stride and table size, not on `sizeof(Float64)`.
- `ns_per_queue_node`  per-node flat-queue bookkeeping.
- `dispatch_ns_base` / `dispatch_ns_per_worker`
                       affine model of one parallel dispatch and join.
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
    ns_per_simd_lane::Float64
    ns_per_scalar_item::Float64
    ns_per_coeff_touch::Float64
    ns_per_queue_node::Float64
    dispatch_ns_base::Float64
    dispatch_ns_per_worker::Float64
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
