# Collector-settings grid: does varying Julia's `--gcthreads` /
# `--heap-size-hint` move the outer_threads (pinned-threads) route's wall
# time on the P3/P4 shapes? Those two flags are process-launch flags, so this
# script does NOT set them itself -- the driver (collector_grid.sh) launches
# one Julia subprocess per grid point with the flags already on argv, and
# this script just times `repeats` back-to-back solves inside that process
# and prints one machine-parseable line.
#
# Usage (normally invoked by collector_grid.sh, not by hand):
#   julia --project=. --threads=8 [--gcthreads=... --heap-size-hint=...] \
#       benchmarks/studies/heap_contention/collector_grid.jl \
#       --case=independent_1sat_1hr --repeats=5

const HC_STUDY_DIR = @__DIR__
const HC_REPO_ROOT = normpath(joinpath(HC_STUDY_DIR, "..", "..", ".."))
const HC_PPC_DIR = joinpath(HC_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(HC_PPC_DIR, "cli.jl"))
include(joinpath(HC_PPC_DIR, "modes.jl"))
include(joinpath(HC_PPC_DIR, "cases.jl"))
using Printf
using Statistics

hc_arg(args, key, default) = begin
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end

const HC_CASE    = hc_arg(ARGS, "case", "independent_1sat_1hr")
const HC_REPEATS = parse(Int, hc_arg(ARGS, "repeats", "5"))
const HC_MODE    = hc_arg(ARGS, "mode", "outer_threads")

function main()
    cfg = PPCConfig(profile="full")
    mode = ppc_mode_specs()[HC_MODE]
    envpairs = ppc_mode_env_pairs(mode, cfg; outer_tasks=1)
    args = ppc_single_config(HC_CASE, cfg)

    withenv(envpairs...) do
        SimulationEngine.run_simulation(args; isolate_state=true, return_solution=true)
    end

    walls = Float64[]
    for _ in 1:HC_REPEATS
        timed = withenv(envpairs...) do
            @timed SimulationEngine.run_simulation(args; isolate_state=false, return_solution=true)
        end
        push!(walls, timed.time)
    end
    @printf("RESULT case=%s mode=%s threads=%d repeats=%d median_wall_s=%.6f min_wall_s=%.6f max_wall_s=%.6f\n",
            HC_CASE, HC_MODE, Threads.nthreads(), HC_REPEATS, median(walls), minimum(walls), maximum(walls))
end

main()
