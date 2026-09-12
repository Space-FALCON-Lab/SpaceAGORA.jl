# End-to-end cost of a results-writing run: solves with `results=true`, writes
# the CSV and feather bundle, and reports wall time and allocation.
#
# Two uses: pricing the output path inside a real run, and byte-diffing the
# written bundle between two checkouts (point OUTDIR at a different directory
# per tree, then `cmp` the files).
#
#   OUTDIR=/tmp/e2e N=1024 PS_GRAVITY=l20 PS_MISSION_S=1800 \
#     julia --project=. --threads=12 <this file>
#
# Knobs: OUTDIR (required), N, PS_GRAVITY, PS_MISSION_S, TREE, PS_REPO_ROOT.
const PS_REPO_ROOT = get(ENV, "PS_REPO_ROOT", normpath(joinpath(@__DIR__, "..", "..", "..", "..")))
ENV["PS_DENSITY"]="none"; ENV["PS_GRAVITY"]=get(ENV,"PS_GRAVITY","l20")
ENV["PS_MISSION_S"]=get(ENV,"PS_MISSION_S","600.0")
ENV["PS_N_SATS"]=get(ENV,"N","256"); ENV["PS_NO_SPICE"]="1"
ENV["PS_WARMUP"]="0"; ENV["PS_REPEATS"]="1"; ENV["PS_WORKLOAD"]="constellation"
include(joinpath(PS_REPO_ROOT, "benchmarks", "studies", "paper_scenarios", "scenario_worker.jl"))

n = parse(Int, ENV["PS_N_SATS"])
outdir = ENV["OUTDIR"]; mkpath(outdir)
base = ps_build_config(n_sats=n)
args = SimulationConfiguration(
    file_paths             = base.file_paths,
    simulation_settings    = SimulationSettings(results=true, verbose=false,
                                                results_directory=outdir,
                                                generate_plots=false, normalize=false,
                                                save_csv=true),
    mission_configuration  = base.mission_configuration,
    environment_model      = base.environment_model,
    dynamics_model         = base.dynamics_model,
    guidance_model         = base.guidance_model,
    navigation_model       = base.navigation_model,
    control_model          = base.control_model,
    initial_time           = base.initial_time,
    integration_tolerances = base.integration_tolerances,
    solver_config          = base.solver_config,
)

# Two warm-ups, not one. The results-writing path JITs across more than one
# call, and this probe used to get a second warm solve by accident, from the
# unconditional `main()` that including scenario_worker.jl once ran.
for _ in 1:2
    SpaceAGORA.run_simulation(args; isolate_state=false)
end
ts=Float64[]; al=Float64[]
for _ in 1:2
    GC.gc(); GC.gc()
    s = @timed SpaceAGORA.run_simulation(args; isolate_state=false)
    push!(ts, s.time); push!(al, s.bytes/2^20)
end
using Statistics
println("E2E tree=$(get(ENV, "TREE", "local")) n=$(n) median_s=$(round(median(ts),digits=3)) alloc_mib=$(round(median(al),digits=1))")
