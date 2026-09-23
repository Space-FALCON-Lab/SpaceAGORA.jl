# Bit-identity tool: dump the full state history of one
# parallelization_performance case to a raw binary file, for a byte-for-byte
# `cmp` between the base commit and a later tip. Same pattern as
# benchmarks/studies/third_body_cost/variants.jl's `--dump` path (read that
# file first if this one is unclear): full Float64 time vector, then every
# saved state flattened with `getdata`, all raw and untouched by any
# tolerance-bearing comparison.
#
# Usage:
#   julia --project=. --threads=1 benchmarks/studies/heap_contention/dump_state.jl \
#       --case=independent_1sat_1hr --out=/tmp/hc_dump_base --mc-index=1
#
# One case per invocation, one dump file per invocation
# (<out>_<case>.bin). Run once at the base commit and once at the tip, same
# thread count both times, then `cmp` the two files.

const HC_STUDY_DIR = @__DIR__
const HC_REPO_ROOT = normpath(joinpath(HC_STUDY_DIR, "..", "..", ".."))
const HC_PPC_DIR = joinpath(HC_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(HC_PPC_DIR, "cli.jl"))
include(joinpath(HC_PPC_DIR, "modes.jl"))
include(joinpath(HC_PPC_DIR, "cases.jl"))
using ComponentArrays: getdata
using Printf

hc_arg(args, key, default) = begin
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end

const HC_CASE     = hc_arg(ARGS, "case", "independent_1sat_1hr")
const HC_OUT       = hc_arg(ARGS, "out", "/tmp/hc_dump")
const HC_MC_INDEX  = parse(Int, hc_arg(ARGS, "mc-index", "1"))
const HC_SEED       = parse(Int, hc_arg(ARGS, "seed", "20260615"))
const HC_MODE      = hc_arg(ARGS, "mode", "serial")

function main()
    cfg = PPCConfig(profile="full", seed=HC_SEED)
    mode = ppc_mode_specs()[HC_MODE]
    envpairs = ppc_mode_env_pairs(mode, cfg; outer_tasks=1)

    @printf("host=%s threads=%d case=%s mc_index=%d mode=%s\n",
            gethostname(), Threads.nthreads(), HC_CASE, HC_MC_INDEX, HC_MODE)
    println("load-at-start: ", strip(read(`uptime`, String)))
    flush(stdout)

    args = ppc_single_config(HC_CASE, cfg; mc_index=HC_MC_INDEX)

    # Warm-up: JIT the RHS/solver stack on a throwaway short solve of the same
    # case before the timed/dumped solve, same reasoning as the harness's own
    # warm-up requirement (see paper_parallelization_benchmarks/CASES.md).
    withenv(envpairs...) do
        SimulationEngine.run_simulation(args; isolate_state=true, return_solution=true)
    end
    println("warmed"); flush(stdout)

    timed = withenv(envpairs...) do
        @timed SimulationEngine.run_simulation(
            args; isolate_state=false, return_solution=true, return_solver_metadata=true
        )
    end
    r = timed.value
    sol = r.solution
    @printf("wall=%8.3fs gc=%6.3fs alloc=%7.2fMiB steps=%d retcode=%s\n",
            timed.time, timed.gctime, timed.bytes / 2^20, length(sol.t), string(sol.retcode))

    outfile = "$(HC_OUT)_$(HC_CASE).bin"
    open(outfile, "w") do io
        write(io, Float64.(sol.t))
        for u in sol.u
            write(io, Float64.(vec(getdata(u))))
        end
    end
    println("dumped: ", outfile, " (", filesize(outfile), " bytes)")
    println("load-at-end: ", strip(read(`uptime`, String)))
end

main()
