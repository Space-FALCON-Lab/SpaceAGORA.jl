# What a solve's return shape costs, on the S1 constellation workload.
#
#   MODE=hist      return_solution=true,  return_solver_metadata=true
#   MODE=metaonly  return_solution=false, return_solver_metadata=true
#   MODE=none      neither
#
# One mode per process, so Sys.maxrss() is a clean peak for that mode. Prints one
# `M3 ...` line: median wall over REPS timed repeats after one warm-up, plus
# allocation, GC time and peak RSS.
#
#   MODE=metaonly N=32768 REPS=5 PS_GRAVITY=invsq PS_MISSION_S=1800 \
#     julia --project=. --threads=16 <this file>
#
# Knobs: MODE, N, REPS, PS_GRAVITY, PS_MISSION_S, TREE (label only), PS_REPO_ROOT.
const PS_REPO_ROOT = get(ENV, "PS_REPO_ROOT", normpath(joinpath(@__DIR__, "..", "..", "..", "..")))
ENV["PS_DENSITY"]="none"; ENV["PS_GRAVITY"]=get(ENV,"PS_GRAVITY","invsq")
ENV["PS_MISSION_S"]=get(ENV,"PS_MISSION_S","1800.0")
ENV["PS_N_SATS"]=get(ENV,"N","32768"); ENV["PS_NO_SPICE"]="1"
ENV["PS_WARMUP"]="0"; ENV["PS_REPEATS"]="1"; ENV["PS_WORKLOAD"]="constellation"
include(joinpath(PS_REPO_ROOT, "benchmarks", "studies", "paper_scenarios", "scenario_worker.jl"))

const MODE = get(ENV,"MODE","hist")
n = parse(Int, ENV["PS_N_SATS"])
args = ps_build_config(n_sats=n)
solve() = if MODE == "hist"       # solution + metadata: what the harness asks for today
    SpaceAGORA.run_simulation(args; isolate_state=false, return_solution=true,  return_solver_metadata=true)
elseif MODE == "metaonly"          # NEW: metadata without the trajectory
    SpaceAGORA.run_simulation(args; isolate_state=false, return_solution=false, return_solver_metadata=true)
else                               # nothing at all
    SpaceAGORA.run_simulation(args; isolate_state=false, return_solution=false, return_solver_metadata=false)
end

# The check runs inside a function so the returned solution is not retained in a
# global across the timed repeats: holding a full trajectory alive there inflates
# the `hist` row's GC and peak RSS and overstates the metadata-only win.
function check()
    r = solve()
    if MODE == "metaonly"
        println("METAONLY_CHECK retcode=$(r.retcode) solution=$(r.solution === nothing ? "nothing" : "PRESENT") " *
                "solver_mode=$(r.solver_mode) trace_entries=$(length(r.solver_trace))")
    elseif MODE == "hist"
        println("HIST_CHECK retcode=$(r.retcode) solution=$(r.solution === nothing ? "nothing" : "PRESENT") " *
                "saved_steps=$(length(r.solution.t))")
    end
    return nothing
end
check()
GC.gc(); GC.gc()

ts=Float64[]; al=Float64[]; gc=Float64[]
for _ in 1:parse(Int, get(ENV,"REPS","2"))
    GC.gc(); GC.gc()
    s = @timed solve()
    push!(ts,s.time); push!(al,s.bytes/2^30); push!(gc,s.gctime)
end
using Statistics
println("M3 tree=$(get(ENV, "TREE", "local")) mode=$(MODE) n=$(n) gravity=$(ENV["PS_GRAVITY"]) threads=$(Threads.nthreads()) " *
        "median_s=$(round(median(ts),digits=3)) times_s=$(join(round.(ts,digits=3),'|')) " *
        "alloc_gib=$(round(median(al),digits=2)) gc_s=$(round(median(gc),digits=3)) " *
        "maxrss_gib=$(round(Sys.maxrss()/2^30,digits=2))")
