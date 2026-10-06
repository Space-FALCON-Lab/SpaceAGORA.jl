# Identity-and-ratio driver for the aerodynamic / density constellation
# cases.
#
# Two jobs, one script, deliberately modeled on
# benchmarks/studies/third_body_cost/variants.jl, whose `--dump` output has
# already been compared with `cmp` byte for byte:
#
#   1. `--dump=<prefix>` writes the FULL state history of each case (every
#      saved time, then every component of every spacecraft at that time) as
#      raw little-endian Float64 to `<prefix>_<case>.bin`, for a byte-for-byte
#      before/after comparison with `cmp`. No tolerance, no summary statistic.
#   2. `--repeats=N` times each case and prints wall seconds, so a before/after
#      pair of runs gives a ratio. Ratios only: these are two runs of the same
#      script on the same machine in one process state, not a benchmark.
#
# Every case is built through the parallelization_performance study's own case
# builders, so the constellation geometry, tolerances, step cap and mission are
# the catalog's and not this script's. The `aero_<N>sat_l50_expatm_<S>s` builder
# branch accepts any N by regex even though the catalog only registers 16 and
# 4096, which is how the 64/256/1024 rungs are reached.
#
# Usage (one Julia process at a time; check `uptime` first):
#   julia --project=. --threads=1 benchmarks/studies/aero_batch/variants.jl \
#       --cases=aero64,atmo256,aero1024 --dump=/path/prefix
#   julia --project=. --threads=8 benchmarks/studies/aero_batch/variants.jl \
#       --cases=atmo256,aero1024 --repeats=3 --mode=inner_only \
#       --csv=benchmarks/studies/aero_batch/results/ratio_8t_after.csv
#
# `--ppc-profile` selects the catalog's mission-length profile (full | smoke |
# test); the default is `full`, i.e. the mission length written into the case
# name. `--no-warmup` skips the untimed 10 s warm-up solves.

using Printf
using Profile

const AB_STUDY_DIR = @__DIR__
const AB_REPO_ROOT = normpath(joinpath(AB_STUDY_DIR, "..", "..", ".."))
const AB_PPC_DIR = joinpath(AB_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(AB_PPC_DIR, "cli.jl"))
include(joinpath(AB_PPC_DIR, "modes.jl"))
include(joinpath(AB_PPC_DIR, "cases.jl"))
using ComponentArrays: getdata

ab_arg(args, key, default) = begin
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end
ab_flag(args, key)::Bool = any(a -> a == "--$key", args)

# Short names for the aero/density reference cases, mapped to the catalog case
# name the builder understands. The mission lengths are the catalog's own
# (100 s for the P6 density ladder, 10 min for the B10 atmosphere ladder);
# nothing here invents a duration.
const AB_CASE_NAMES = Dict(
    "aero64"    => "aero_64sat_l50_expatm_100s",
    "aero256"   => "aero_256sat_l50_expatm_100s",
    "aero1024"  => "aero_1024sat_l50_expatm_100s",
    "atmo256"   => "atmo256_exponential_10min",
    "vacuum256" => "atmo256_vacuum_10min",
    "multi64"   => "multi_64_high_fidelity",
    "mars1"     => "montecarlo_heavy_aerobraking",
)

const AB_CASES   = String.(split(ab_arg(ARGS, "cases", "aero64,atmo256,aero1024"), ","))
const AB_REPEATS = parse(Int, ab_arg(ARGS, "repeats", "1"))
const AB_SOLVER  = ab_arg(ARGS, "solver", "auto_stiff")
const AB_RHS     = ab_arg(ARGS, "rhs", "auto")
const AB_MODE    = ab_arg(ARGS, "mode", "serial")
const AB_PROFILE_NAME = ab_arg(ARGS, "ppc-profile", "full")
const AB_DUMP    = ab_arg(ARGS, "dump", "")
const AB_CSV     = ab_arg(ARGS, "csv", "")
const AB_WARMUP  = !ab_flag(ARGS, "no-warmup")

function ab_case_name(key::String)::String
    haskey(AB_CASE_NAMES, key) && return AB_CASE_NAMES[key]
    # An unrecognized key is passed through verbatim so a catalog case can be
    # named directly without editing this table.
    return key
end

function ab_config(key::String, profile::String)
    cfg = PPCConfig(profile=profile, solver_mode=AB_SOLVER)
    return ppc_single_config(ab_case_name(key), cfg)
end

function ab_solve(key::String, profile::String)
    args = ab_config(key, profile)
    nsats = length(args.dynamics_model.spacecraft)
    timed = @timed SimulationEngine.run_simulation(
        args; isolate_state=false, return_solution=true, return_solver_metadata=true
    )
    r = timed.value
    sol = r.solution
    geti(name) = try Int(getproperty(sol.stats, name)) catch; -1 end
    sc = sol.u[end].sc[1]
    return (
        wall_s=Float64(timed.time), gc_s=Float64(timed.gctime), bytes=Int(timed.bytes),
        nf=geti(:nf), naccept=geti(:naccept), nreject=geti(:nreject),
        steps=length(sol.t), retcode=string(sol.retcode),
        nsats=nsats,
        pos=Float64.(sc.pos), vel=Float64.(sc.vel),
        sol=sol,
    )
end

# Raw Float64 state history: the saved times, then every component of every
# spacecraft at each saved time, in solver order. Identical bytes or the change
# does not ship.
function ab_write_dump(path::String, sol)
    open(path, "w") do io
        write(io, Float64.(sol.t))
        for u in sol.u
            write(io, Float64.(vec(getdata(u))))
        end
    end
    return filesize(path)
end

function main()
    mode = ppc_mode_specs()[AB_MODE]
    cfg = PPCConfig(profile=AB_PROFILE_NAME, solver_mode=AB_SOLVER)
    envpairs = copy(ppc_mode_env_pairs(mode, cfg; outer_tasks=1))
    if AB_RHS != "auto"
        push!(envpairs, "SPACEAGORA_RHS_EXECUTION_MODE" => AB_RHS)
    end

    @printf("host=%s threads=%d cases=%s profile=%s solver=%s rhs_mode=%s parallel_mode=%s\n",
            gethostname(), Threads.nthreads(), join(AB_CASES, ","),
            AB_PROFILE_NAME, AB_SOLVER, AB_RHS, AB_MODE)
    println("load-at-start: ", strip(read(`uptime`, String)))
    flush(stdout)

    # Warm every case at the "test" profile (10 s mission) so the timed solves
    # below are post-JIT. The warm-up solves are never timed and never dumped.
    if AB_WARMUP
        for key in AB_CASES
            withenv(envpairs...) do
                ab_solve(key, "test")
            end
        end
        println("warmed all cases"); flush(stdout)
    end

    rows = Dict{String, Vector{NamedTuple}}()
    for rep in 1:AB_REPEATS, key in AB_CASES
        r = withenv(envpairs...) do
            ab_solve(key, AB_PROFILE_NAME)
        end
        if !isempty(AB_DUMP) && rep == 1
            path = "$(AB_DUMP)_$(key).bin"
            nbytes = ab_write_dump(path, r.sol)
            @printf("dump %-10s -> %s (%d bytes)\n", key, path, nbytes)
        end
        push!(get!(rows, key, NamedTuple[]), r)
        @printf("rep%d %-10s N=%5d wall=%9.3fs gc=%6.3fs alloc=%8.3fGiB nf=%8d acc=%6d rej=%5d steps=%6d %s\n",
                rep, key, r.nsats, r.wall_s, r.gc_s, r.bytes / 2^30,
                r.nf, r.naccept, r.nreject, r.steps, r.retcode)
        flush(stdout)
    end

    println("\n-- terminal state of spacecraft 1 (all 17 digits, for cross-run comparison) --")
    for key in AB_CASES
        haskey(rows, key) || continue
        r = rows[key][1]
        @printf("%-10s pos=[%.17g, %.17g, %.17g] vel=[%.17g, %.17g, %.17g]\n",
                key, r.pos[1], r.pos[2], r.pos[3], r.vel[1], r.vel[2], r.vel[3])
    end

    println("\n-- allocation per spacecraft per derivative evaluation --")
    for key in AB_CASES
        haskey(rows, key) || continue
        r = rows[key][1]
        denom = max(1, r.nf * r.nsats)
        @printf("%-10s %.1f B/sat/eval (total %.3f GiB over nf=%d x N=%d)\n",
                key, r.bytes / denom, r.bytes / 2^30, r.nf, r.nsats)
    end

    if !isempty(AB_CSV)
        mkpath(dirname(AB_CSV))
        open(AB_CSV, "w") do io
            println(io, "host,threads,mode,rhs_mode,profile,case,case_name,rep,nsats,wall_s,gc_s,bytes,nf,naccept,nreject,steps,retcode")
            for key in AB_CASES
                haskey(rows, key) || continue
                for (rep, r) in enumerate(rows[key])
                    @printf(io, "%s,%d,%s,%s,%s,%s,%s,%d,%d,%.6f,%.6f,%d,%d,%d,%d,%d,%s\n",
                            gethostname(), Threads.nthreads(), AB_MODE, AB_RHS, AB_PROFILE_NAME,
                            key, ab_case_name(key), rep, r.nsats, r.wall_s, r.gc_s, r.bytes,
                            r.nf, r.naccept, r.nreject, r.steps, r.retcode)
                end
            end
        end
        println("\ncsv written: ", AB_CSV)
    end

    println("load-at-end: ", strip(read(`uptime`, String)))
end

main()
