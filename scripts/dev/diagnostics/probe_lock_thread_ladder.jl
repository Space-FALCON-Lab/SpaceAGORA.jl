#!/usr/bin/env julia
#
# Compare RHS timing across interleaved execution widths in one Julia process.
#
#   julia --project=. --threads=12 scripts/dev/diagnostics/probe_lock_thread_ladder.jl
#   PLT_WIDTHS=1,2,4 julia --project=. --threads=4 scripts/dev/diagnostics/probe_lock_thread_ladder.jl
#
# Start Julia with enough threads for the requested PLT_WIDTHS. The probe warms
# every width, rotates their measurement order between repeats and discards
# timing blocks in which the collector ran. It reports the minimum retained
# time per RHS call, speedup relative to the first requested width, spread and
# the number of discarded blocks. A width with no retained block is explicit.
#
# This is a timing comparison, not a lock-contention diagnosis: it does not
# report native-lock wait/hold. The separate probe_contention_matrix.jl reports
# allocation pressure and lock statistics by workload family. Both tools use
# probe_contention_matrix_support.jl so their workload construction stays shared.

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
# Measuring each width to completion in sequence can confound width with
# thermal state, frequency, a co-tenant or allocator settling. Interleaving
# spreads those effects across the widths being compared.
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
            # Discard blocks in which the collector ran before taking a minimum.
            # Comparing its counter before and after each block makes that
            # exclusion explicit rather than assuming the fastest block escaped.
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
