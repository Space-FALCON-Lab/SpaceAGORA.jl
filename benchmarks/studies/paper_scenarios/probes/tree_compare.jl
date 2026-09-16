# Runs one fixed configuration against a given checkout, so two trees can be
# compared on wall time and allocation at the same workload.
#
# Defaults to the shape used for large-N CPU scaling comparisons: inverse-square
# gravity only, no atmosphere, thread count from the caller, one warm-up solve
# and two timed repetitions. Prints one `NCMP ...` line.
#
#   for mode in hist nohist; do
#     PS_REPO_ROOT=<tree> TREE=<label> MODE=$mode N=32768 PS_GRAVITY=invsq \
#     julia --project=<tree> --threads=16 <tree>/<this file>
#   done
#
# Knobs: MODE (hist|nohist), N, PS_GRAVITY, PS_MISSION_S, TREE, PS_REPO_ROOT.
# PS_REPO_ROOT is what lets one invocation measure a checkout other than the one
# the file sits in.
const PS_REPO_ROOT = get(ENV, "PS_REPO_ROOT", normpath(joinpath(@__DIR__, "..", "..", "..", "..")))
ENV["PS_DENSITY"]  = "none"
ENV["PS_GRAVITY"]  = get(ENV, "PS_GRAVITY", "invsq")
ENV["PS_MISSION_S"]= get(ENV, "PS_MISSION_S", "1800.0")
ENV["PS_N_SATS"]   = get(ENV, "N", "32768")
ENV["PS_NO_SPICE"] = "1"
ENV["PS_WARMUP"]   = "0"; ENV["PS_REPEATS"] = "1"; ENV["PS_WORKLOAD"] = "constellation"

include(joinpath(PS_REPO_ROOT, "benchmarks", "studies", "paper_scenarios", "scenario_worker.jl"))

const MODE = get(ENV, "MODE", "hist")          # hist = retain solution+metadata (what S1/main does)
const WANT = MODE == "hist"                    # nohist = ask for neither
n = parse(Int, ENV["PS_N_SATS"])
args = ps_build_config(n_sats=n)
solve() = SpaceAGORA.run_simulation(args; isolate_state=false,
                                    return_solution=WANT, return_solver_metadata=WANT)

solve()                                        # warm-up 1 (JIT + first-touch)
ts = Float64[]; allocs = Float64[]; gcs = Float64[]
for rep in 1:2                                 # 2 timed repetitions
    GC.gc(); GC.gc()
    s = @timed solve()
    push!(ts, s.time); push!(allocs, s.bytes/2^30); push!(gcs, s.gctime)
end
using Statistics
println("NCMP tree=$(get(ENV, "TREE", "local")) mode=$(MODE) n=$(n) gravity=$(ENV["PS_GRAVITY"]) " *
        "threads=$(Threads.nthreads()) median_s=$(round(median(ts),digits=3)) " *
        "times_s=$(join(round.(ts,digits=3),'|')) alloc_gib=$(round(median(allocs),digits=2)) " *
        "gc_s=$(round(median(gcs),digits=3)) maxrss_gib=$(round(Sys.maxrss()/2^30,digits=2))")
