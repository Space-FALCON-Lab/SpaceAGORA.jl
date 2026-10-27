# Diagnostic: is the first RHS calibration sweep's extra cost compilation or
# another first-time effect? One process, 4096-satellite L50 vacuum gravity case
# (the ppc harness case), the sweep driven directly (it never writes the store).
#
#   julia --project=. --threads=16 scripts/verification/first_sweep_cost_probe.jl A   # 3 full sweeps
#   julia --project=. --threads=16 scripts/verification/first_sweep_cost_probe.jl B   # 1-candidate sweep, then full x2
#
# Per sweep it prints wall, compile (total and recompile), GC time, allocated
# bytes, then one line per recorded candidate (ns of its last round). Wall times
# are diagnostic only. Compile time is a process-wide counter, so compilation on
# other threads during the sweep is included.

using SpaceAGORA
using Printf
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PPC = joinpath(ROOT, "benchmarks", "studies", "parallelization_performance")
include(joinpath(PPC, "cli.jl"))
include(joinpath(PPC, "cases.jl"))
include(joinpath(ROOT, "scripts", "dev", "diagnostics", "probe_contention_matrix_support.jl"))

const MODE = isempty(ARGS) ? "A" : ARGS[1]
Base.cumulative_compile_timing(true)
ct() = (x = Base.cumulative_compile_time_ns(); x isa Tuple ? (Float64(x[1]), Float64(x[2])) : (Float64(x), 0.0))

cfg = PPCConfig(profile = "full")
args0 = ppc_single_config("gravity_4096sat_l50_vacuum_5800s", cfg)
(; args, p, u0) = probe_params_from_args(args0)
SM.ParallelPolicy.reset_policy_telemetry!()
SE._spin_barrier_tuner_reset!(p.shared_buffers.spin_barrier_tuner, false)
effs = args.dynamics_model.dynamic_effectors
println("# mode=$MODE threads=$(Threads.nthreads()) budget=$(SM.ParallelPolicy.effective_inner_thread_budget()) n_sats=$(length(args.dynamics_model.spacecraft))")

function one_sweep(label, seed_k)
    GC.gc()
    g0 = Base.gc_num(); c0 = ct(); t0 = time_ns()
    local res
    bytes = @allocated res = SE._run_rhs_sweep_core!(p, u0, effs, false, args, seed_k)
    t1 = time_ns(); c1 = ct(); g1 = Base.gc_num()
    gcns = g1.total_time - g0.total_time
    @printf("SWEEP %s seed_k=%d wall_s=%.3f compile_s=%.3f recompile_s=%.3f gc_s=%.3f alloc_MB=%.1f\n",
            label, seed_k, (t1 - t0) / 1e9, (c1[1] - c0[1]) / 1e9, (c1[2] - c0[2]) / 1e9, gcns / 1e9, bytes / 2^20)
    for t in res[6]
        @printf("CAND %s %s a=%d %s ns=%.0f reps=%d round=%d\n", label, t["mode"], t["allotment"], t["scheduler"], t["ns"], t["reps"], t["round"])
    end
    flush(stdout)
end

if MODE == "A"
    for i in 1:3
        one_sweep("s$i", 0)
    end
else
    one_sweep("k1", 1)
    one_sweep("full1", 0)
    one_sweep("full2", 0)
end
