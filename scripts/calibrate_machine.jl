#!/usr/bin/env julia
#
# One-shot per-machine calibration for the analytic cost model.
#
#   julia --project=. --threads=auto scripts/calibrate_machine.jl [--verbose] [--repeats=N]
#
# Writes output/parallel_policy_state/cost_constants_<fingerprint>.toml, which
# every later run on this machine reads instead of measuring anything itself.
#
# Run it with the thread count you intend to simulate at: the dispatch and
# atomic constants are properties of the inner dispatch primitive at that width.
#
# --repeats runs the whole pass N times and reports each constant's spread.
#
# Read that spread as diagnostics, NOT as the acceptance criterion. Measured
# here across three passes, the arithmetic rates reproduce to 2-4% while the
# thread-dispatch ones reach 10-29%, because dispatch pays OS wake latency on
# every call and so offers no uncontended window for a minimum to find. That
# looks alarming and mostly is not: dispatch cost is common-mode between the
# routing candidates (both dispatch), so it largely cancels in the comparison,
# while the terms that actually discriminate -- the coefficient-touch curve and
# the SIMD-lane plateau -- are the stable ones.
#
# The criterion that matters is whether calibration noise changes a DECISION.
# Across 100 configurations (satellite count x degree/order x worker count)
# scored under three independent calibration passes, quiet and again under
# eight competing processes, the candidate ranking disagreed zero times; the
# median margin between candidates was ~70%. One configuration sat at a 0.1%
# margin, which is the case the predictor's abstention guard exists for.

using SpaceAGORA

const PC = SpaceAGORA.SimulationModel.ParallelCost

verbose = any(a -> a in ("--verbose", "-v"), ARGS)
function _parse_repeats(args)::Int
    for a in args
        startswith(a, "--repeats=") && return max(1, parse(Int, split(a, "=")[2]))
    end
    return 1
end
repeats = _parse_repeats(ARGS)

println("SpaceAGORA cost-model machine calibration")
println("  fingerprint : ", PC.machine_fingerprint())
println("  threads     : ", Threads.nthreads())
println("  repeats     : ", repeats)
println()

results = PC.MachineConstants[]
for r in 1:repeats
    repeats > 1 && println("--- pass $(r)/$(repeats) ---")
    push!(results, PC.calibrate_machine(; verbose = verbose))
end
mc = results[end]

if repeats > 1
    println()
    println("Reproducibility across $(repeats) passes (max spread, % of min).")
    println("Diagnostics, not an acceptance gate -- see the header comment:")
    spread(v) = isempty(v) || minimum(v) <= 0 ? NaN :
        (maximum(v) - minimum(v)) / minimum(v) * 100
    report(name, f) = println("  ", rpad(name, 24),
        isnan(spread(map(f, results))) ? "n/a" : string(round(spread(map(f, results)), digits=2), "%"))
    report("ns_per_scalar_item", m -> m.ns_per_scalar_item)
    report("ns_per_queue_node", m -> m.ns_per_queue_node)
    report("dispatch_ns_base", m -> m.dispatch_ns_base)
    report("dispatch_ns_per_worker", m -> m.dispatch_ns_per_worker)
    report("ns_per_atomic", m -> m.ns_per_atomic)
    report("reference_fma_ns", m -> m.reference_fma_ns)
    report("reference_mem_ns", m -> m.reference_mem_ns)
    for i in eachindex(mc.coeff_touch.ns)
        report("coeff_touch[$(i)]", m -> m.coeff_touch.ns[i])
    end
    for i in eachindex(mc.simd_lane.ns)
        report("simd_lane[$(i)]", m -> m.simd_lane.ns[i])
    end
end

path = PC.save_machine_constants(mc)
println()
println("Wrote ", path)
