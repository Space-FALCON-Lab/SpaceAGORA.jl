# Scoring a routing candidate from work counts and machine constants.
#
# NOT WIRED INTO ANY ROUTING PATH, and not currently fit to be. `select_plan` is
# called by nothing outside tests and the validation script. An earlier version
# of this comment described it as the replacement for the pre-solve sweep, which
# read as shipped intent; it is not, and at its measured accuracy it would be a
# regression against the sweep it was meant to replace.
#
# Measured by scripts/validate_cost_model.jl against the held-out real harmonics
# kernel, over six runs on identical cached constants: 30-45% decision accuracy
# at picking the fastest plan (30/30/35/35/35/45), median regret +2.8% to +11%.
# That is worse than the heuristic it would displace, and the abstention guard
# below cannot rescue it, because where it is wrong the margins are not small --
# it is confidently wrong rather than uncertain.
#
# QUOTE THAT AS A RANGE. Predictions here are deterministic given cached
# constants, so the spread is not the model changing its mind: it is measurement
# noise reordering which candidate counts as the true best. Accuracy is
# therefore unstable exactly where the best is a near-tie, and at those rows it
# measures how often two near-equal plans swap rather than anything about the
# model. Regret is the better headline for the same reason -- continuous, and it
# degrades gracefully at ties where a binary accuracy figure does not.
#
# Scored as the coarser question "serial or threaded?", the same six runs give
# 60-70% and a tighter spread (2/20 against 3/20), because that view discards
# precisely the within-side ties the full-plan metric spends its variance on.
# Serial-vs-threaded plus regret is the more stable pair if a headline is needed.
# One configuration predicts threaded where serial measures faster in every one
# of the six runs -- N=8, L=50, M=50, a 2.4x miss -- which is the shape of error
# no margin can catch, since the model is confident there.
#
# What the same measurement does support is using it as a FILTER rather than a
# chooser. The true best sits in its top-1 only 25% of the time but in its top-2
# 60% and top-3 70%, and a sweep restricted to its top two carries 0.0% median
# regret -- so it can cut a candidate set by most of its width without ever
# being trusted to name the winner. Pruning also fails safely where gating does
# not: a wrongly pruned candidate costs convergence speed, a wrongly opened gate
# costs wall time immediately. If this is used in rhs_calibration.jl at all it
# should prune the sweep's candidates, never replace the sweep.
#
# Mechanically: the sweep evaluates the whole RHS per candidate and takes the
# minimum of repeated timings; this evaluates a handful of multiplies per
# candidate and picks argmin of the predicted cost, with an explicit refusal to
# answer when the candidates are too close to separate.

# Wait/hold above which the model declines to rank plans at all.
#
# Below it, a lock is a serial fraction and `usl_speedup` handles it. Above it
# the workers are spending most of their time queueing, and the cost of that
# queue is set by convoy behaviour the linear-plus-quadratic form has no term
# for. Three is deliberately loose: at wait/hold = 3 the measured ladder is
# already past its turning point, so anything the model says beyond there would
# be extrapolation dressed as prediction.
const LOCK_CONVOY_ABSTAIN_RATIO = 3.0

"""
    PlanCandidate

An RHS execution plan the model can score: which route, how wide, and -- for the
flat route -- which inner scheduler.
"""
struct PlanCandidate
    mode::Symbol
    allotment::Int
    scheduler::Symbol
end

"""
    PlanPrediction

A scored candidate, with the derived quantities that produced the score kept
alongside it. `workers` and `batch_width` are retained because they are what
make a prediction explicable after the fact -- "it chose 4 workers because at 12
the batch width fell to 85 and the lane rate doubled" is a diagnosis, whereas a
bare nanosecond figure is not.
"""
struct PlanPrediction
    candidate::PlanCandidate
    ns::Float64
    workers::Int
    batch_width::Int
end

"""
    predict_plan_ns(counts, mc, candidate; n_active_sats, budget) -> PlanPrediction

Predicted wall time for one pass of the dynamic effectors under `candidate`.

Wall time, not summed work: the per-worker terms are divided by the satellite
split because workers run concurrently, so what the solve waits for is one
worker's share. The terms that do *not* divide are the ones that matter:

  - `coeff_touches` is paid once per worker on the flat route, because the SIMD
    batch kernel loads each coefficient once and broadcasts it across the whole
    satellite batch. On `satellite_batch` it is paid once per satellite, because
    that route re-walks the table for every one. This asymmetry is the entire
    discriminator between the two candidates; everything else is symmetric.
  - the dynamic scheduler's atomic is a single contended cache line, so its
    increments serialise rather than spreading across workers, and the cost
    scales with total items rather than items per worker.

The lane rate is evaluated at the batch width the route actually produces, which
is why `satellite_batch` is expensive in a way a flat-rate model would miss: it
runs at batch width 1, where the measured lane rate is ~90x the plateau because
per-pass loop overhead has nothing to amortise over.
"""
function predict_plan_ns(
    counts::WorkCounts,
    mc::MachineConstants,
    candidate::PlanCandidate;
    n_active_sats::Int,
    budget::Int,
    contention::Union{Nothing, ContentionInputs} = nothing,
)::PlanPrediction
    n_sats = max(1, n_active_sats)
    max_workers = max(1, min(budget, n_sats))

    workers = candidate.mode === :satellite_batch ?
        max_workers :
        max(1, min(candidate.allotment, max_workers))

    batch = max(1, cld(n_sats, workers))
    # Per-mechanism: satellite_batch dispatches through Polyester.@batch, the
    # flat routes through the persistent channel pool. Using one number for both
    # is what took the first predictor to 25% decision accuracy.
    dispatch = if workers <= 1
        0.0
    elseif candidate.mode === :satellite_batch
        mc.dispatch_batch_ns_base + workers * mc.dispatch_batch_ns_per_worker
    else
        mc.dispatch_pool_ns_base + workers * mc.dispatch_pool_ns_per_worker
    end

    scalar = batch * counts.scalar_items * mc.ns_per_scalar_item
    probe = batch * counts.probe_ns

    # Achieved speedup, not assumed linear scaling. `batch` already divides the
    # work by the worker count; this corrects that ideal share back up to what
    # the machine actually delivers.
    #
    # WHICH achieved speedup depends on what is known about the workload, and the
    # difference is the point of the USL path. `parallel_speedup` is a
    # machine-only curve measured from one synthetic kernel, so it says the same
    # thing about a constellation that allocates nothing per pass and one that
    # allocates hundreds of KiB -- and those two measure a plateau and an
    # inversion respectively on the same box. A single curve cannot be right for
    # both, and being wrong here moves the predicted cost of a wide plan by more
    # than the margins the router is trying to resolve.
    #
    # Given a warm probe, the parameters are predicted per workload instead:
    # `alpha` from the dispatch primitive's own serial fraction plus the measured
    # native-lock duty cycle, `beta` from allocation rate and working set. Absent
    # a probe this falls back to the legacy curve, so callers that cannot probe
    # behave exactly as before.
    ideal = Float64(max(1, workers))
    achieved = if workers <= 1
        1.0
    elseif contention === nothing || !contention.valid
        max(1.0, rate_at(mc.parallel_speedup, workers))
    else
        alpha, beta = workload_usl_parameters(mc, contention; n_active_sats = n_sats)
        usl_speedup(alpha, beta, workers)
    end
    contention_factor = ideal / achieved

    total = if candidate.mode === :satellite_batch
        # satellite_batch runs one satellite at a time, so its workspace is one
        # satellite's worth.
        lane1 = rate_at(mc.simd_lane, max(8.0, counts.simd_workspace_bytes_per_sat))
        touch = rate_at(mc.coeff_touch, counts.coeff_table_bytes)
        dispatch + contention_factor * (
            batch * counts.simd_terms * lane1 +
            batch * counts.coeff_touches * touch +
            scalar + probe)
    else
        lane = rate_at(mc.simd_lane, max(8.0, batch * counts.simd_workspace_bytes_per_sat))
        touch = rate_at(mc.coeff_touch, counts.coeff_table_bytes)
        nodes_per_worker = counts.queue_nodes / workers
        atomics = candidate.scheduler === :dynamic && workers > 1 ?
            counts.queue_nodes * mc.ns_per_atomic : 0.0
        dispatch + atomics + contention_factor * (
            batch * counts.simd_terms * lane +
            counts.coeff_touches * touch +
            nodes_per_worker * mc.ns_per_queue_node +
            scalar + probe)
    end

    return PlanPrediction(candidate, total, workers, batch)
end

"""
    plan_candidates(counts; budget, n_active_sats) -> Vector{PlanCandidate}

The candidate set the model scores: `satellite_batch`, plus the flat route over
a geometric allotment ladder crossed with both schedulers.

The ladder is geometric in the thread budget rather than in satellite count.
A ladder built from satellite count collapses under the `min(a, budget)` clamp
-- at 1024 satellites and a 12-thread budget the entries beyond 12 all become
12, so the sweep only ever compares "narrow" against "full budget" and cannot
find an optimum in between. There usually is one, because the arithmetic term
falls as 1/workers while the dispatch and coefficient-walk terms rise with it.
"""
function plan_candidates(
    counts::WorkCounts;
    budget::Int,
    n_active_sats::Int,
)::Vector{PlanCandidate}
    max_workers = max(1, min(budget, n_active_sats))
    out = PlanCandidate[PlanCandidate(:satellite_batch, 1, :static)]

    allotments = Int[1]
    a = 2
    while a < max_workers
        push!(allotments, a)
        a *= 2
    end
    push!(allotments, max_workers)
    sort!(unique!(allotments))

    for allot in allotments
        if allot == 1
            # At width 1 the dispatcher degenerates to a serial loop and no
            # scheduler is observable, so only one entry is meaningful.
            push!(out, PlanCandidate(:flat_constellation_effector_queue, 1, :static))
            continue
        end
        for sched in (:static, :dynamic)
            push!(out, PlanCandidate(:flat_constellation_effector_queue, allot, sched))
        end
    end
    return out
end

"""
    select_plan(counts, mc; budget, n_active_sats, margin = 0.15) -> Union{Nothing, PlanPrediction}

Score every candidate and return the cheapest, or `nothing` when the model
should not be trusted to choose.

Returns `nothing` in three cases, and the distinction between them matters less
than the fact that all three are explicit rather than silent:

  - `counts.in_domain` is false. Something in the workload has no representation
    as a sum of per-unit rates -- native GRAM's process-wide lock, whose cost is
    superlinear in concurrency, or an effector whose probe failed. No
    calibration can fix a missing term, so the model declines.
  - the constants are stale, which the caller checks with
    [`constants_are_current`](@ref) before calling here.
  - the best and second-best predictions are within `margin` of each other. A
    model with residual error cannot rank candidates whose true costs are closer
    together than that error, and picking anyway would be reading noise.

Abstention is what makes "never worse than the static route" structural rather
than hoped-for: the caller keeps its existing heuristic, so the model can only
ever move a decision it is confident about. It is also what makes a
badly-calibrated machine safe -- constants that do not describe the hardware
produce compressed margins, and compressed margins abstain.

`margin` defaults to 15%, chosen from the measured prediction residual against
the held-out real kernel rather than picked for comfort.
"""
function select_plan(
    counts::WorkCounts,
    mc::MachineConstants;
    budget::Int,
    n_active_sats::Int,
    margin::Float64 = 0.15,
    contention::Union{Nothing, ContentionInputs} = nothing,
)::Union{Nothing, PlanPrediction}
    # `in_domain` used to be the whole gate, and it is false for any workload
    # whose density model holds the process-wide native lock -- which refused to
    # answer on exactly the workloads where the split matters most.
    #
    # A measured lock duty cycle makes that refusal too broad. The lock is an
    # Amdahl term: it caps speedup at 1/rho and `workload_usl_parameters` now
    # carries it in `alpha`, which the model can represent. What it still cannot
    # represent is the convoy regime -- once workers spend far longer queueing
    # than holding, the cost stops being linear in occupancy and the fit has no
    # term for it. So the refusal narrows from "this model takes a lock" to
    # "this workload is already past what the lock admits", which is a statement
    # about the regime rather than about the model type.
    if contention !== nothing && contention.valid
        lock_wait_hold_ratio(contention) > LOCK_CONVOY_ABSTAIN_RATIO && return nothing
    else
        counts.in_domain || return nothing
    end

    candidates = plan_candidates(counts; budget=budget, n_active_sats=n_active_sats)
    isempty(candidates) && return nothing

    preds = [
        predict_plan_ns(counts, mc, c; n_active_sats=n_active_sats, budget=budget,
                        contention=contention)
        for c in candidates
    ]
    filter!(p -> isfinite(p.ns) && p.ns > 0.0, preds)
    isempty(preds) && return nothing
    length(preds) == 1 && return preds[1]

    sort!(preds; by = p -> p.ns)
    best, runner_up = preds[1], preds[2]
    (runner_up.ns - best.ns) / best.ns < margin && return nothing
    return best
end
