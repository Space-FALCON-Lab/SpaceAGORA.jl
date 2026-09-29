# Byte-for-byte dump of a solve's recorded outputs: the results table the run
# writes (every saved field, including drag/lift/cross, heat rate, winds and,
# with the visualization scene on, the per-spacecraft density the viewer
# colors passes by). `variants.jl --dump` covers the integrated state; these
# columns are read from shared buffers at save time, so a change of RHS route
# could move them without moving a state bit. Compare two runs' files with `cmp`.
#
# Usage (one Julia process at a time; check `uptime` first):
#   julia --project=. --threads=1 benchmarks/studies/aero_batch/outputs_dump.jl \
#       --cases=aero_256sat_l50_expatm_100s --dump=/path/prefix

using Printf

const AB_STUDY_DIR = @__DIR__
const AB_REPO_ROOT = normpath(joinpath(AB_STUDY_DIR, "..", "..", ".."))
const AB_PPC_DIR = joinpath(AB_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(AB_PPC_DIR, "cli.jl"))
include(joinpath(AB_PPC_DIR, "modes.jl"))
include(joinpath(AB_PPC_DIR, "cases.jl"))

ab_arg(args, key, default) = begin
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end

const AB_CASES = String.(split(ab_arg(ARGS, "cases", "aero_256sat_l50_expatm_100s"), ","))
const AB_DUMP  = ab_arg(ARGS, "dump", "")
const AB_RHS   = ab_arg(ARGS, "rhs", "auto")

# The catalog config with results writing and the visualization fields on,
# pointed at a scratch directory. Everything else is the catalog's.
function ab_config_with_outputs(case::String, dir::String)
    base = ppc_single_config(case, PPCConfig(profile="full"))
    fields = Dict{Symbol, Any}(f => getfield(base, f) for f in fieldnames(typeof(base)))
    fields[:simulation_settings] = SimulationSettings(
        results=true, verbose=false, generate_plots=false, normalize=false,
        save_csv=true, save_visualization_scene=true, results_directory=dir,
    )
    return SimulationConfiguration(; fields...)
end

function main()
    isempty(AB_DUMP) && error("--dump=<prefix> is required")
    @printf("host=%s threads=%d rhs_mode=%s\n", gethostname(), Threads.nthreads(), AB_RHS)
    for case in AB_CASES
        dir = mktempdir()
        args = ab_config_with_outputs(case, dir)
        withenv("SPACEAGORA_RHS_EXECUTION_MODE" => AB_RHS) do
            SimulationEngine.run_simulation(args; isolate_state=false)
        end
        csvs = String[]
        for (root, _, files) in walkdir(dir), f in files
            endswith(f, ".csv") && push!(csvs, joinpath(root, f))
        end
        length(csvs) == 1 || error("expected one results CSV under $dir, found $(csvs)")
        out = "$(AB_DUMP)_$(case).csv"
        cp(only(csvs), out; force=true)
        header = readline(out)
        has_density = occursin("density", header)
        @printf("%-34s -> %s (%d bytes, %d columns, density column: %s)\n",
                case, out, filesize(out), count(==(','), header) + 1, has_density)
    end
end

main()
