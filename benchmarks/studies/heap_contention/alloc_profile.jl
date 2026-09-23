# WS11c allocation attribution: where does a P3/P4 sample's per-solve
# allocation actually go -- solution storage (the SavedValues/ODESolution
# per-step state history), callbacks, or solver internals?
#
# `Profile.Allocs` samples individual allocations with their originating
# stack. P3 (independent_1sat_1hr, ~15 MiB/sample, ~180k allocations at
# sample_rate=1.0) is small enough to sample exhaustively. P4
# (montecarlo_heavy_aerobraking) is NOT small: a 6 h mission at a 1 s max
# step is ~365 MiB and several million allocations per sample, and running
# `sample_rate=1.0` on it is what took this machine's whole session down via
# the OOM killer on 2026-09-23 (Profile.Allocs' own bookkeeping -- a stack
# capture per sampled allocation -- grows with the *sampled* allocation
# count, not the run's byte total, so an unrestrained rate on a run this size
# multiplies rather than shrinks the footprint being profiled). P4 MUST be
# run at a low `--sample-rate` (default here is 0.005; never above 0.01) and
# every invocation of this script MUST be wrapped in a hard memory cap:
#
#   systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
#       julia --project=. --threads=1 benchmarks/studies/heap_contention/alloc_profile.jl \
#       --case=montecarlo_heavy_aerobraking --sample-rate=0.005 \
#       --out=benchmarks/studies/heap_contention/results/alloc_p4.csv
#
# A sub-1.0 sample rate still gives a valid per-bucket *fraction* of total
# allocation (uniform random sub-sampling preserves the relative split
# between buckets; only the absolute byte/count totals become estimates
# scaled by 1/sample_rate, and this script does not report those as exact).
#
# Every allocation's top few frames are bucketed by the `src/` directory that
# owns the file, using the domain-ownership table in CLAUDE.md, so the report
# reads as "which subsystem", not "which function".
#
# Usage:
#   julia --project=. --threads=1 benchmarks/studies/heap_contention/alloc_profile.jl \
#       --case=independent_1sat_1hr --out=benchmarks/studies/heap_contention/results/alloc_p3.csv

const HC_STUDY_DIR = @__DIR__
const HC_REPO_ROOT = normpath(joinpath(HC_STUDY_DIR, "..", "..", ".."))
const HC_PPC_DIR = joinpath(HC_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(HC_PPC_DIR, "cli.jl"))
include(joinpath(HC_PPC_DIR, "modes.jl"))
include(joinpath(HC_PPC_DIR, "cases.jl"))
using Printf
using Profile
using Profile.Allocs: @profile, fetch, clear

hc_arg(args, key, default) = begin
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end

const HC_CASE = hc_arg(ARGS, "case", "independent_1sat_1hr")
const HC_OUT  = hc_arg(ARGS, "out", joinpath(HC_STUDY_DIR, "results", "alloc_$(HC_CASE).csv"))
const HC_MODE = hc_arg(ARGS, "mode", "serial")
# See the file header: P4-scale cases (tens of MiB or more per sample) must
# not use the default of 1.0. 0.005 is deliberately conservative -- below the
# coordinator's 0.01 ceiling -- since this script has already OOM'd this
# machine once at full sampling.
const HC_SAMPLE_RATE = parse(Float64, hc_arg(ARGS, "sample-rate", "1.0"))
HC_SAMPLE_RATE > 0.01 && HC_CASE != "independent_1sat_1hr" && @warn(
    "sample-rate=$(HC_SAMPLE_RATE) on case=$(HC_CASE): this is only known " *
    "safe for small samples (see file header). Pass --sample-rate=0.005 or " *
    "lower for anything P4-scale."
)

# Bucket a file path into the src/ domain that owns it (see the ownership
# table in CLAUDE.md). Anything outside src/ (DifferentialEquations.jl,
# DiffEqCallbacks.jl, base Julia, ...) is bucketed as "external".
function hc_bucket(path::AbstractString)::String
    isempty(path) && return "unknown"
    # Split the engine/ directory apart: WS11c owns exactly execution.jl and
    # persistence.jl there, and the whole point of this attribution is to say
    # what fraction of allocation those two files (vs. the rest of engine/,
    # owned by other workstreams) are actually responsible for.
    if occursin("src/simulation/engine/execution.jl", path)
        return "engine/execution.jl (WS11c-owned)"
    elseif occursin("src/simulation/engine/persistence.jl", path)
        return "engine/persistence.jl (WS11c-owned)"
    end
    idx = findfirst("src/simulation/engine/", path)
    idx !== nothing && return "engine/other (solver_policy, dynamics_rhs, rhs_calibration, ... -- not WS11c)"
    idx = findfirst("src/simulation/callbacks/", path)
    idx !== nothing && return "callbacks (SavedValues/SaveData/SavingCallback)"
    idx = findfirst("src/io/", path)
    idx !== nothing && return "io (results dataframe/CSV)"
    idx = findfirst("src/dynamics/", path)
    idx !== nothing && return "dynamics (RHS terms)"
    idx = findfirst("src/environment/", path)
    idx !== nothing && return "environment (gravity/atmosphere/ephemerides)"
    idx = findfirst("src/vehicle/", path)
    idx !== nothing && return "vehicle (spacecraft/structure)"
    idx = findfirst("src/core/", path)
    idx !== nothing && return "core"
    idx = findfirst("src/parallel/", path)
    idx !== nothing && return "parallel (routing/policy)"
    occursin("SpaceAGORA.jl/src/", path) && return "src (other)"
    return "external (solver/OrdinaryDiffEq/DiffEqCallbacks/Base)"
end

function hc_first_spaceagora_frame(alloc)
    for frame in alloc.stacktrace
        path = String(frame.file)
        occursin("SpaceAGORA.jl", path) && return path
    end
    isempty(alloc.stacktrace) && return ""
    return String(alloc.stacktrace[1].file)
end

function main()
    cfg = PPCConfig(profile="full")
    mode = ppc_mode_specs()[HC_MODE]
    envpairs = ppc_mode_env_pairs(mode, cfg; outer_tasks=1)

    @printf("host=%s threads=%d case=%s mode=%s\n", gethostname(), Threads.nthreads(), HC_CASE, HC_MODE)
    println("load-at-start: ", strip(read(`uptime`, String)))
    flush(stdout)

    args = ppc_single_config(HC_CASE, cfg)

    # Warm-up (JIT), untimed and unprofiled.
    withenv(envpairs...) do
        SimulationEngine.run_simulation(args; isolate_state=true, return_solution=true)
    end
    println("warmed"); flush(stdout)

    # Byte total from @timed, for the headline number and the
    # bytes-per-accepted-step ratio.
    timed = withenv(envpairs...) do
        @timed SimulationEngine.run_simulation(
            args; isolate_state=false, return_solution=true, return_solver_metadata=true
        )
    end
    sol = timed.value.solution
    naccept = try
        Int(sol.stats.naccept)
    catch
        length(sol.t) - 1
    end
    @printf("measured: wall=%.3fs alloc=%.2fMiB naccept=%d bytes/accepted-step=%.1f\n",
            timed.time, timed.bytes / 2^20, naccept, timed.bytes / max(1, naccept))

    # Allocation-site profile, on a fresh (still-warmed) solve so the JIT
    # compilation itself never enters the sample.
    println("profiling at sample_rate=", HC_SAMPLE_RATE); flush(stdout)
    clear()
    @profile sample_rate=HC_SAMPLE_RATE withenv(envpairs...) do
        SimulationEngine.run_simulation(
            args; isolate_state=false, return_solution=true, return_solver_metadata=true
        )
    end
    results = fetch()
    allocs = results.allocs
    println("profiled allocations: ", length(allocs)); flush(stdout)

    bucket_bytes = Dict{String, Int}()
    bucket_count = Dict{String, Int}()
    for a in allocs
        path = hc_first_spaceagora_frame(a)
        bucket = hc_bucket(path)
        bucket_bytes[bucket] = get(bucket_bytes, bucket, 0) + Int(a.size)
        bucket_count[bucket] = get(bucket_count, bucket, 0) + 1
    end
    total_bytes = sum(values(bucket_bytes); init=0)

    mkpath(dirname(HC_OUT))
    open(HC_OUT, "w") do io
        println(io, "case,mode,sample_rate,bucket,sampled_bytes,sampled_count,fraction_of_sampled_total,estimated_bytes_at_1x")
        for (bucket, bytes) in sort(collect(bucket_bytes); by=kv -> -kv[2])
            frac = total_bytes == 0 ? 0.0 : bytes / total_bytes
            est = round(Int, bytes / HC_SAMPLE_RATE)
            println(io, "$(HC_CASE),$(HC_MODE),$(HC_SAMPLE_RATE),\"$(bucket)\",$(bytes),$(bucket_count[bucket]),$(round(frac; digits=4)),$(est)")
        end
    end
    println("\n-- allocation attribution (sampled) --")
    for (bucket, bytes) in sort(collect(bucket_bytes); by=kv -> -kv[2])
        frac = total_bytes == 0 ? 0.0 : bytes / total_bytes
        @printf("%-55s %10.2f MiB  %6.1f%%  (n=%d)\n", bucket, bytes / 2^20, 100 * frac, bucket_count[bucket])
    end
    println("written: ", HC_OUT)
    println("load-at-end: ", strip(read(`uptime`, String)))
end

main()
