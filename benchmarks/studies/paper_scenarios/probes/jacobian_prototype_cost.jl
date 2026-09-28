# Cost of the block-diagonal Jacobian prototype at large N, and whether the
# default solver (:tsit5, explicit) ever consumes it.
#
# Prints the active solver mode, then one `JACPROBE ...` line per constellation
# size: per-satellite block width, total state length, build time, allocation,
# nnz and the size of the resulting SparseMatrixCSC.
#
#   julia --project=. --threads=1 <this file>
const PS_REPO_ROOT = get(ENV, "PS_REPO_ROOT", normpath(joinpath(@__DIR__, "..", "..", "..", "..")))
ENV["PS_DENSITY"]="none"; ENV["PS_GRAVITY"]="l20"; ENV["PS_MISSION_S"]="60.0"
ENV["PS_N_SATS"]="4"; ENV["PS_NO_SPICE"]="1"
ENV["PS_WARMUP"]="0"; ENV["PS_REPEATS"]="1"; ENV["PS_WORKLOAD"]="constellation"
include(joinpath(PS_REPO_ROOT, "benchmarks", "studies", "paper_scenarios", "scenario_worker.jl"))

using SparseArrays
const SE = SpaceAGORA.SimulationEngine
println("default solver_mode = ", SpaceAGORA.SimulationModel.SolverConfig().solver_mode)

for n in (1024, 4096, 16384, 32768)
    args = ps_build_config(n_sats=n)
    u0 = SE.build_initial_conditions(args)
    active = trues(n)
    SE._build_block_diagonal_jac_prototype(u0, active)   # warm
    GC.gc()
    s = @timed SE._build_block_diagonal_jac_prototype(u0, active)
    J = s.value
    println("JACPROBE n=$n  block_len=$(length(u0.sc[1]))  n_total=$(length(u0))  " *
            "build_s=$(round(s.time, digits=4))  alloc_mib=$(round(s.bytes/2^20, digits=1))  " *
            "nnz=$(nnz(J))  J_mib=$(round(Base.summarysize(J)/2^20, digits=1))")
end
