#!/usr/bin/env julia
#
# Does the native lock explain the multi-effector thread ceiling?
#
#   for t in 1 2 4 8 12; do julia --project=. --threads=$t \
#       scripts/probe_lock_thread_ladder.jl; done
#
# scripts/probe_contention_matrix.jl measures occupancy rho at one thread, which
# gives the Amdahl ceiling 1/rho but cannot show contention: with one worker
# nothing ever waits. This runs the same workload across a thread ladder and
# reads wait/hold, which is the quantity that says whether the width in use is
# already past what the lock admits.
#
# The prediction under test, on the fullstack shape (harmonics + SRP + aero over
# a SPICE-backed planet -- the configuration `heavy_1024sat_fullstack_1hr` uses):
# rho ~ 0.48 at one thread caps speedup near 2.1x, and wait/hold should climb
# with the thread count. If it stays near zero the lock is not the ceiling and
# this line of explanation is dead.
#
# §7.4 of the parallelization findings ruled out GC thread count, false sharing,
# dispatch mechanism, allocation and route choice, and left the root cause open.
# Native-lock occupancy was never on that list because nothing could measure it.

using SpaceAGORA
using StaticArrays
using Printf

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const RS = SpaceAGORA.RuntimeServices

include(joinpath(@__DIR__, "probe_contention_matrix_support.jl"))

const PASSES = parse(Int, get(ENV, "PLT_PASSES", "40"))
const REPEATS = parse(Int, get(ENV, "PLT_REPEATS", "7"))
const CASE_NAME = get(ENV, "PLT_CASE", "fullstack")
# A PLT_CASE the local table does not know is looked up in the benchmark
# harness's catalog, so a published number can be checked against the case
# that produced it rather than against a reconstruction of it.
# Widths swept inside ONE process, so the comparison across widths is paired
# against a single machine state. Sweeping by relaunching Julia per width -- the
# earlier shape of this probe -- lets drift between processes land on the width
# axis, which is the axis under test.
const WIDTHS = [parse(Int, w) for w in split(get(ENV, "PLT_WIDTHS", "1,2,4,8,12"), ",")]

# Widths are INTERLEAVED, with the order rotating between repeats.
#
# The previous shape of this sweep measured width 1 to convergence, then width
# 2, and so on. That is blocked sampling on the axis under test: anything
# drifting over the sweep -- thermal state, frequency, a co-tenant, the
# allocator settling -- lands on the width axis and reads as scaling. It is the
# same defect that made a fresh-process-per-width sweep untrustworthy, moved
# inside one process rather than removed, and this project has now produced
# three retracted results from exactly that family of error.
#
# Rotating the order additionally spreads first-position and last-position bias
# across widths instead of accumulating it on whichever width is measured first.
function sweep(case, du, widths, passes, repeats)
    mission_s = Float64(case.args.mission_configuration.mission_time)
    samples = Dict(w => Float64[] for w in widths)
    rejected = Int[]
    order = collect(widths)
    for rep in 1:repeats
        for w in order
            case.p.shared_buffers.rhs_plan_override[] = SE._make_calib_flat_plan(w, :static)
            # Blocks in which the collector ran are DISCARDED, not minimised
            # over. The fullstack RHS allocates ~790 KiB per pass, so a sweep of
            # this length turns over hundreds of MB and GC lands on blocks at
            # random -- measured spreads of +139% to +1559% across widths, which
            # is an order of magnitude more than the effect under test. Taking a
            # minimum hopes some block escaped; checking the collector's own
            # counter knows which ones did.
            gc0 = Base.gc_num().total_time
            t0 = time_ns()
            for i in 1:passes
                SE.spacecraft_dynamics!(du, case.u0, case.p,
                                        mission_s * (i - 1) / max(1, passes - 1))
            end
            elapsed = Float64(time_ns() - t0)
            if Base.gc_num().total_time == gc0
                push!(samples[w], elapsed / passes)
            else
                push!(rejected, w)
            end
        end
        push!(order, popfirst!(order))
    end
    return samples, rejected
end

function main()
    case = haskey(PROBE_CASES, CASE_NAME) ?
        PROBE_CASES[CASE_NAME]() : build_harness_case(CASE_NAME)
    du = zero(case.u0)
    # Warm every width before timing any of them: the first evaluation at a new
    # width compiles and re-warms the pool, and that cost belongs to neither.
    for w in WIDTHS
        case.p.shared_buffers.rhs_plan_override[] = SE._make_calib_flat_plan(w, :static)
        for _ in 1:3
            SE.spacecraft_dynamics!(du, case.u0, case.p, 0.0)
        end
    end

    RS.reset_native_lock_stats!()
    samples, rejected = sweep(case, du, WIDTHS, PASSES, REPEATS)
    base_samples = samples[first(WIDTHS)]
    base = isempty(base_samples) ? NaN : minimum(base_samples)
    @printf("case=%-30s N=%d L=%s threads=%d passes=%d repeats=%d (interleaved, rotating)\n",
            CASE_NAME, length(case.args.dynamics_model.spacecraft),
            get(ENV, "PCM_L", "20"), Threads.nthreads(), PASSES, REPEATS)
    @printf("  GC-contaminated blocks discarded: %d of %d\n",
            length(rejected), length(WIDTHS) * REPEATS)
    for w in WIDTHS
        v = sort(samples[w])
        if isempty(v)
            @printf("  width %2d   no GC-free block in %d repeats\n", w, REPEATS)
            continue
        end
        lo, hi = v[1], v[end]
        @printf("  width %2d   %9.1f us/call   speedup %5.2fx   spread %+.0f%%\n",
                w, lo / 1e3, base / lo, 100 * (hi - lo) / lo)
    end
end

main()
