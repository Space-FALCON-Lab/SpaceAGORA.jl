# One trial of the process-pool startup measurement. Run in a fresh
# coordinator process per trial, against the tree under test:
#
#   julia --project=<tree> --threads=1 measure.jl <build> <mode> <N> <outdir> [variant] [samples]
#
#   mode = bootstrap  ensure_process_workers!(pool, N), no warm-up
#          warmup     ensure_process_workers!(pool, N; warmup_fn) with one
#                     independent_1sat_1hr sample ("full" profile, 1 h mission)
#          campaign   run_monte_carlo(f, 1:samples; threads=:auto), the shipped
#                     campaign runner, steered to the process route at N workers
#   variant (campaign only) = default | skip | cheap
#          default  no PROCESS_WARMUP override: the runner warms with f(first(seeds))
#          skip     PROCESS_WARMUP => false
#          cheap    PROCESS_WARMUP => the same case at the "test" profile (10 s mission)
#
# The sample closure is anonymous and captures only strings, so Distributed
# ships it whole to a freshly bootstrapped worker; on first call it includes
# the parallelization_performance study files of <tree> (as the harness's own
# pool does) and runs the case through ppc_run_sample_once. That include is
# part of every first call, before and after alike, and is timed separately
# in warmup mode.
#
# Appends one CSV row per trial to <outdir>/<mode>.csv; campaign mode also
# writes every sample's terminal state to <outdir>/campaign_samples.csv.

using Distributed
using Printf
using SpaceAGORA

const BUILD = ARGS[1]
const MODE = ARGS[2]
const N = parse(Int, ARGS[3])
const OUTDIR = abspath(ARGS[4])
const VARIANT = length(ARGS) >= 5 ? ARGS[5] : "default"
const SAMPLES = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 256
const TREE = dirname(Base.active_project())
const STUDY = joinpath(TREE, "benchmarks", "studies", "parallelization_performance")
const CASE = "independent_1sat_1hr"
mkpath(OUTDIR)

Threads.nthreads() == 1 || error("run the coordinator with --threads=1 (got $(Threads.nthreads()))")

# Seed -> NamedTuple of plain values. `profile` selects the mission length.
function sample_closure(study::String, case::String, profile::String)
    return seed -> begin
        t_start = time()
        if !isdefined(Main, :ppc_run_sample_once)
            for f in ("cli.jl", "modes.jl", "cases.jl", "trajectory_parity.jl", "execution.jl")
                Base.include(Main, joinpath(study, f))
            end
        end
        t_loaded = time()
        r = Base.invokelatest() do
            cfg = Main.PPCConfig(; profile=profile)
            Main.ppc_run_sample_once(case, cfg, Int(seed), 20260615 + Int(seed))
        end
        (pid=getpid(), t_start=t_start, t_loaded=t_loaded, t_end=time(),
         success=r.success, solve_s=r.wall_time_s,
         terminal_time_s=Float64(coalesce(r.terminal.terminal_time_s, NaN)),
         pos_norm_m=Float64(coalesce(r.terminal.pos_norm_m, NaN)),
         vel_norm_mps=Float64(coalesce(r.terminal.vel_norm_mps, NaN)))
    end
end

function append_row(name::String, header::String, row::String)
    path = joinpath(OUTDIR, name)
    new = !isfile(path)
    open(path, "a") do io
        new && println(io, header)
        println(io, row)
    end
end

g(x) = @sprintf("%.6f", x)
full = sample_closure(STUDY, CASE, "full")
pool = SpaceAGORA.campaign_process_pool()
has_override = isdefined(SpaceAGORA, :PROCESS_WARMUP)
stamp = string(round(Int, time()))

if MODE == "bootstrap"
    t0 = time()
    ids = SpaceAGORA.ensure_process_workers!(pool, N)
    t1 = time()
    length(ids) == N || error("pool has $(length(ids)) workers, expected $N")
    loaded = all(w -> remotecall_fetch(() -> isdefined(Main, :SpaceAGORA) && isdefined(Main, :GRAMSuite), w), ids)
    append_row("bootstrap.csv", "build,n_workers,stamp,ensure_s,all_loaded_spaceagora_gramsuite",
               join((BUILD, N, stamp, g(t1 - t0), loaded), ","))
    println("RESULT bootstrap build=$BUILD N=$N ensure_s=$(g(t1 - t0)) loaded=$loaded")

elseif MODE == "warmup"
    stampdir = mktempdir()
    warm = let f = full, dir = stampdir
        () -> begin
            v = f(0)
            open(joinpath(dir, "warm_$(getpid()).txt"), "w") do io
                println(io, join((v.t_start, v.t_loaded, v.t_end, v.success, v.solve_s), " "))
            end
            nothing
        end
    end
    t0 = time()
    ids = SpaceAGORA.ensure_process_workers!(pool, N; warmup_fn=warm)
    t1 = time()
    length(ids) == N || error("pool has $(length(ids)) workers, expected $N")
    rows = [split(readline(joinpath(stampdir, f))) for f in readdir(stampdir)]
    length(rows) == N || error("$(length(rows)) of $N warm-ups reported")
    starts = [parse(Float64, r[1]) for r in rows]
    loads = [parse(Float64, r[2]) - parse(Float64, r[1]) for r in rows]
    solves = [parse(Float64, r[5]) for r in rows]
    firsts = [parse(Float64, r[3]) - parse(Float64, r[2]) for r in rows]
    ok = all(r -> r[4] == "true", rows)
    # Bootstrap component: call start to the first warm-up starting (both
    # versions bootstrap every worker before warming any). Warm-up component:
    # from there to the call's return.
    boot = minimum(starts) - t0
    append_row("warmup.csv",
               "build,n_workers,stamp,ensure_s,bootstrap_s,warmup_s,median_include_s,median_first_sample_s,max_first_sample_s,median_first_solve_s,all_success",
               join((BUILD, N, stamp, g(t1 - t0), g(boot), g(t1 - t0 - boot),
                     g(sort(loads)[cld(N, 2)]), g(sort(firsts)[cld(N, 2)]), g(maximum(firsts)),
                     g(sort(solves)[cld(N, 2)]), ok), ","))
    println("RESULT warmup build=$BUILD N=$N ensure_s=$(g(t1 - t0)) bootstrap_s=$(g(boot)) warmup_s=$(g(t1 - t0 - boot)) ok=$ok")

elseif MODE == "campaign"
    SC = SpaceAGORA.SimulationCampaigns
    PP = SpaceAGORA.ParallelProfiles
    VARIANT == "default" || has_override || error("variant $VARIANT needs PROCESS_WARMUP, absent in this build")
    # Route state isolated in a scratch file, so no persisted history is read
    # or written; the route decision is the cold one.
    ENV["SPACEAGORA_OUTER_ROUTE_STATE_PATH"] = joinpath(mktempdir(), "outer_route_state.toml")
    # Route features from the case's own config, as the harness computes them.
    for f in ("cli.jl", "modes.jl", "cases.jl", "trajectory_parity.jl", "execution.jl")
        Base.include(Main, joinpath(STUDY, f))
    end
    features = Base.invokelatest() do
        probe = Main.ppc_single_config(CASE, Main.PPCConfig(; profile="full"); seed=20260615, mc_index=1)
        SC.campaign_route_features(probe; samples=SAMPLES)
    end
    tuning = PP.OuterRouteTuning(process_max_workers=N)
    state = PP.OuterRouteState()
    cheap = sample_closure(STUDY, CASE, "test")
    override = VARIANT == "skip" ? false : VARIANT == "cheap" ? (() -> cheap(0)) : nothing
    run() = SC.run_monte_carlo(full, collect(1:SAMPLES); threads=:auto, fail_fast=true,
                               route_features=features, route_state=state, route_tuning=tuning)
    t0 = time()
    r = has_override ? Base.ScopedValues.with(run, SpaceAGORA.PROCESS_WARMUP => override) : run()
    t1 = time()
    r.route === :process || error("campaign took route $(r.route), not :process")
    vals = [s.value for s in r.samples]
    first_dispatch = minimum(v.t_start for v in vals) - t0
    occupancy = sort([v.t_end - v.t_start for v in vals])
    nworkers = length(unique(v.pid for v in vals))
    ok = all(s -> s.success && s.value.success, r.samples)
    append_row("campaign.csv",
               "build,variant,n_workers,samples,stamp,total_s,first_dispatch_s,dispatch_s,route,threads,workers_used,median_sample_s,max_sample_s,all_success",
               join((BUILD, VARIANT, N, SAMPLES, stamp, g(t1 - t0), g(first_dispatch), g(t1 - t0 - first_dispatch),
                     r.route, r.threads, nworkers, g(occupancy[cld(end, 2)]), g(occupancy[end]), ok), ","))
    for (s, v) in zip(r.samples, vals)
        append_row("campaign_samples.csv",
                   "build,variant,stamp,index,pid,sample_s,solve_s,terminal_time_s,pos_norm_m,vel_norm_mps",
                   join((BUILD, VARIANT, stamp, s.index, v.pid, g(v.t_end - v.t_start), g(v.solve_s),
                         @sprintf("%.17g", v.terminal_time_s), @sprintf("%.17g", v.pos_norm_m),
                         @sprintf("%.17g", v.vel_norm_mps)), ","))
    end
    println("RESULT campaign build=$BUILD variant=$VARIANT N=$N total_s=$(g(t1 - t0)) first_dispatch_s=$(g(first_dispatch)) workers=$nworkers ok=$ok")
else
    error("unknown mode $MODE")
end

nprocs() > 1 && rmprocs(workers())
