# WS11b step 1: the locked native-GRAM path against the isolated per-worker
# GRAM pool, on one workstation, at one thread count per process.
#
# Native GRAM is single-threaded behind a process-wide lock, so a threaded
# constellation serializes on it. `SPACEAGORA_GRAM_ISOLATED_POOL` replaces the
# one shared model with independent `deepcopy`ed models, one per worker, each
# behind its own lock. This script measures what that buys.
#
# Two density paths, the two the P6/S2 traces use and for the reasons given in
# parallelization_performance/cases.jl (`_ppc_p6_gram_density_env!`):
#
#   lookahead  the vacuum-predicted look-ahead cache (density_callbacks/
#              vacuum_predicted_gram.jl), horizon past the end of the mission and
#              deviation past anything the mission can reach, so only the initial
#              build runs (a rebuild mid-run hangs the job). Not to be confused
#              with the GRAM track cache in callbacks/gram_track_cache/, which is
#              a different mechanism behind a different env var.
#   freeze     direct native GRAM with density frozen per accepted step.
#
# Every (size, density path) group is solved back to back in one process, in one
# process state, locked first and then each pool width. Only ratios within a
# group are meaningful; the absolute seconds are recorded so a group can be
# checked for drift, not quoted as a benchmark. The thread count is the
# process's own `--threads` and so varies across processes, never within one.
#
# Usage (one process per thread count, serialized by the caller):
#   julia --project=. --threads=4 \
#       benchmarks/studies/gram_thread_scaling/run_scaling.jl \
#       --sizes=256 --density=lookahead,freeze --pool=0,2,4,8 \
#       --mission=100 --repeats=1 --out=results/scaling_t4.csv

using Printf
using Statistics

const GTS_DIR = @__DIR__
const GTS_REPO_ROOT = normpath(joinpath(GTS_DIR, "..", "..", ".."))
const GTS_PPC_DIR = joinpath(GTS_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(GTS_PPC_DIR, "cli.jl"))
include(joinpath(GTS_PPC_DIR, "modes.jl"))
include(joinpath(GTS_PPC_DIR, "cases.jl"))
include(joinpath(GTS_DIR, "run_scaling_config.jl"))

ppc_ensure_gramsuite_loaded!()

const RS = SpaceAGORA.RuntimeServices

gts_arg(args, key, default) = begin
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end
gts_ints(s) = [parse(Int, x) for x in split(s, ",") if !isempty(x)]

const GTS_SIZES    = gts_ints(gts_arg(ARGS, "sizes", "256,1024"))
const GTS_DENSITY  = String.(split(gts_arg(ARGS, "density", "lookahead,freeze"), ","))
const GTS_POOL     = gts_ints(gts_arg(ARGS, "pool", "0,2,4,8"))
const GTS_MISSION  = parse(Float64, gts_arg(ARGS, "mission", "100"))
const GTS_REPEATS  = parse(Int, gts_arg(ARGS, "repeats", "1"))
const GTS_WARMUP_S = parse(Float64, gts_arg(ARGS, "warmup-mission", "5"))
const GTS_EI_KM    = parse(Float64, gts_arg(ARGS, "ei-km", "600"))
const GTS_OUT      = gts_arg(ARGS, "out", "")

function gts_solve(n::Int, mission_s::Float64, env)
    args = gts_build_config(n, mission_s, GTS_EI_KM)
    # The counters accumulate for the life of the process, so a run's own
    # occupancy is only readable against a window that starts here.
    RS.reset_native_lock_stats!()
    timed = withenv(env...) do
        @timed SimulationEngine.run_simulation(args; isolate_state=false, return_solution=true)
    end
    sol = timed.value isa NamedTuple ? timed.value.solution : timed.value
    snap = RS.native_lock_stats_snapshot()
    stats = sol.stats
    geti(name) = try Int(getproperty(stats, name)) catch; -1 end
    return (
        wall_s = Float64(timed.time),
        gc_s = Float64(timed.gctime),
        nf = geti(:nf),
        naccept = geti(:naccept),
        gram_density_acq = snap.sites.gram_density.acquisitions,
        gram_density_hold_ns = snap.sites.gram_density.hold_ns,
        gram_density_wait_ns = snap.sites.gram_density.wait_ns,
        gram_cache_acq = snap.sites.gram_cache.acquisitions,
        gram_cache_hold_ns = snap.sites.gram_cache.hold_ns,
        gram_cache_wait_ns = snap.sites.gram_cache.wait_ns,
        total_hold_ns = snap.hold_ns,
        total_wait_ns = snap.wait_ns,
    )
end

const GTS_HEADER = "threads,n_sats,ei_km,density_path,pool_workers,rep,wall_s,gc_s,nf,naccept," *
                   "gram_density_acq,gram_density_hold_ns,gram_density_wait_ns," *
                   "gram_cache_acq,gram_cache_hold_ns,gram_cache_wait_ns," *
                   "total_hold_ns,total_wait_ns"

function main()
    nthreads = Threads.nthreads()
    @printf("threads=%d sizes=%s density=%s pool=%s mission=%.0fs repeats=%d EI=%.0fkm\n",
            nthreads, join(GTS_SIZES, ","), join(GTS_DENSITY, ","),
            join(GTS_POOL, ","), GTS_MISSION, GTS_REPEATS, GTS_EI_KM)

    rows = String[]
    # One warm-up solve per process, at the smallest size and a short mission:
    # the first solve pays compilation and the first native GRAM initialization,
    # and charging that to whichever configuration happens to run first would
    # put the whole difference between locked and pooled inside the noise.
    print("warm-up ... "); flush(stdout)
    gts_solve(minimum(GTS_SIZES), GTS_WARMUP_S,
              vcat(gts_density_env("freeze", GTS_WARMUP_S), gts_pool_env(0), GTS_WIDTH_ENV))
    gts_solve(minimum(GTS_SIZES), GTS_WARMUP_S,
              vcat(gts_density_env("freeze", GTS_WARMUP_S), gts_pool_env(2), GTS_WIDTH_ENV))
    gts_solve(minimum(GTS_SIZES), GTS_WARMUP_S,
              vcat(gts_density_env("lookahead", GTS_WARMUP_S), gts_pool_env(0), GTS_WIDTH_ENV))
    println("done")

    for n in GTS_SIZES, path in GTS_DENSITY
        group = Dict{Int, Vector{Float64}}()
        # Repeat outermost, arm innermost: the arms alternate, so any drift over
        # the group (cache state, thermal, a neighboring process) lands on all
        # of them instead of on whichever ran last. The reported ratio is taken
        # from each arm's minimum, which is the least contaminated repeat.
        for rep in 1:GTS_REPEATS, pool in GTS_POOL
            env = vcat(gts_density_env(path, GTS_MISSION), gts_pool_env(pool), GTS_WIDTH_ENV)
            r = gts_solve(n, GTS_MISSION, env)
            push!(get!(group, pool, Float64[]), r.wall_s)
            push!(rows, join((nthreads, n, GTS_EI_KM, path, pool, rep,
                              @sprintf("%.6f", r.wall_s), @sprintf("%.6f", r.gc_s),
                              r.nf, r.naccept,
                              r.gram_density_acq, r.gram_density_hold_ns, r.gram_density_wait_ns,
                              r.gram_cache_acq, r.gram_cache_hold_ns, r.gram_cache_wait_ns,
                              r.total_hold_ns, r.total_wait_ns), ","))
            @printf("n=%-5d %-9s pool=%-2d rep%d wall=%8.3fs nf=%8d gram_density(acq=%9d hold=%7.3fs wait=%7.3fs) gram_cache(acq=%7d hold=%7.3fs wait=%7.3fs)\n",
                    n, path, pool, rep, r.wall_s, r.nf,
                    r.gram_density_acq, r.gram_density_hold_ns * 1e-9, r.gram_density_wait_ns * 1e-9,
                    r.gram_cache_acq, r.gram_cache_hold_ns * 1e-9, r.gram_cache_wait_ns * 1e-9)
            flush(stdout)
        end
        if haskey(group, 0)
            base = minimum(group[0])
            println("  -- ratios against the locked path, same process, back to back --")
            for pool in GTS_POOL
                haskey(group, pool) || continue
                @printf("  n=%-5d %-9s pool=%-2d  locked/pool x%.3f\n",
                        n, path, pool, base / minimum(group[pool]))
            end
        end
    end

    if !isempty(GTS_OUT)
        mkpath(dirname(GTS_OUT))
        exists = isfile(GTS_OUT)
        open(GTS_OUT, "a") do io
            exists || println(io, GTS_HEADER)
            for row in rows
                println(io, row)
            end
        end
        println("wrote $(length(rows)) rows to $(GTS_OUT)")
    end
    return nothing
end

main()
