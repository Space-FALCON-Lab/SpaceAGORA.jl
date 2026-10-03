# Third-body/SRP cost attribution for constellation RHS evaluation (WS10b).
#
# Reproduces the P6 observation that a degree-50 vacuum constellation solved
# serially costs far less than the same constellation with solar radiation
# pressure and Sun/Moon third-body gravity added, and attributes the gap.
#
# The two cases are the P6 traces themselves, built through the
# parallelization_performance study's own case builders, so the configuration
# (constellation geometry, tolerances, step cap, mission) is the one the
# calibration measured and the only thing that differs is the effector tuple.
#
# Every number this prints is a relative attribution: a ratio between two
# solves timed back to back in the same process, or a share of one profile.
# It is not a machine benchmark and must not be quoted as one.
#
# Usage:
#   julia --project=. --threads=1 benchmarks/studies/third_body_cost/attribution.jl \
#       --n=256 --mission=5800 [--repeats=1] [--mode=serial] [--profile]

using Printf
using Profile

const TBC_STUDY_DIR = @__DIR__
const TBC_REPO_ROOT = normpath(joinpath(TBC_STUDY_DIR, "..", "..", ".."))
const TBC_PPC_DIR = joinpath(TBC_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(TBC_PPC_DIR, "cli.jl"))
include(joinpath(TBC_PPC_DIR, "modes.jl"))
include(joinpath(TBC_PPC_DIR, "cases.jl"))

function tbc_arg(args, key::String, default::String)::String
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    return default
end
tbc_flag(args, key::String)::Bool = any(a -> a == "--$key", args)

const TBC_N        = parse(Int, tbc_arg(ARGS, "n", "256"))
const TBC_MISSION  = parse(Int, tbc_arg(ARGS, "mission", "5800"))
const TBC_REPEATS  = parse(Int, tbc_arg(ARGS, "repeats", "1"))
const TBC_MODE     = tbc_arg(ARGS, "mode", "serial")
const TBC_PROFILE  = tbc_flag(ARGS, "profile")
const TBC_CASES    = split(tbc_arg(ARGS, "cases", "vacuum,srp_nbody"), ",")

tbc_case_name(kind::String) = kind == "vacuum" ?
    "gravity_$(TBC_N)sat_l50_vacuum_$(TBC_MISSION)s" :
    "gravity_$(TBC_N)sat_l50_srp_nbody_vacuum_$(TBC_MISSION)s"

# Counters are per solve; read them from the run's own return value.
function tbc_counter_row(counters)
    counters === nothing && return (nbody=0, srp=0, pxform=0)
    get_i(k) = begin
        v = getproperty(counters, k)
        v isa Base.Threads.Atomic ? v[] : Int(v)
    end
    return (
        nbody  = get_i(:nbody_spkpos_runtime_calls),
        srp    = get_i(:srp_spkpos_runtime_calls),
        pxform = get_i(:planet_pxform_runtime_calls),
    )
end

function tbc_solve(case_name::String, cfg::PPCConfig)
    args = ppc_single_config(case_name, cfg)
    timed = @timed SimulationEngine.run_simulation(
        args; isolate_state=false, return_solution=true, return_solver_metadata=true
    )
    result = timed.value
    sol = result.solution
    nf = try
        Int(sol.stats.nf)
    catch
        -1
    end
    return (
        wall_s = Float64(timed.time),
        gc_s = Float64(timed.gctime),
        bytes = Int(timed.bytes),
        nf = nf,
        steps = length(sol.t),
        retcode = string(sol.retcode),
        counters = tbc_counter_row(get(result, :spice_counters, nothing)),
        terminal = sol.u[end].sc[1],
    )
end

function main()
    mode = ppc_mode_specs()[TBC_MODE]
    warm_cfg = PPCConfig(profile="test")         # 10 s mission: compile everything
    full_cfg = PPCConfig(profile="full")

    @printf("host=%s threads=%d mode=%s N=%d mission=%ds repeats=%d\n",
            gethostname(), Threads.nthreads(), TBC_MODE, TBC_N, TBC_MISSION, TBC_REPEATS)
    println("load-at-start: ", strip(read(`uptime`, String)))

    results = Dict{String, Vector{NamedTuple}}()
    for kind in TBC_CASES
        name = tbc_case_name(String(kind))
        withenv(ppc_mode_env_pairs(mode, warm_cfg; outer_tasks=1)...) do
            tbc_solve(name, warm_cfg)
        end
        println("warmed: ", name)
        flush(stdout)
    end

    for rep in 1:TBC_REPEATS, kind in TBC_CASES
        kind = String(kind)
        name = tbc_case_name(kind)
        r = withenv(ppc_mode_env_pairs(mode, full_cfg; outer_tasks=1)...) do
            tbc_solve(name, full_cfg)
        end
        push!(get!(results, kind, NamedTuple[]), r)
        @printf("rep%d %-12s wall=%.3fs gc=%.3fs alloc=%.2fGiB nf=%d steps=%d retcode=%s spkpos(nbody=%d srp=%d) pxform=%d\n",
                rep, kind, r.wall_s, r.gc_s, r.bytes / 2^30, r.nf, r.steps, r.retcode,
                r.counters.nbody, r.counters.srp, r.counters.pxform)
        flush(stdout)
    end

    if haskey(results, "vacuum") && haskey(results, "srp_nbody")
        v = minimum(r.wall_s for r in results["vacuum"])
        s = minimum(r.wall_s for r in results["srp_nbody"])
        @printf("RATIO srp_nbody/vacuum = %.2f  (best-of-%d each, %.3fs vs %.3fs)\n",
                s / v, TBC_REPEATS, s, v)
        rv = results["vacuum"][1]
        rs = results["srp_nbody"][1]
        @printf("nf: vacuum=%d srp_nbody=%d (ratio %.2f); per-nf wall ratio %.2f\n",
                rv.nf, rs.nf, rs.nf / max(1, rv.nf),
                (s / max(1, rs.nf)) / (v / max(1, rv.nf)))
        @printf("srp_nbody ephemeris calls per RHS evaluation: nbody=%.3f srp=%.3f (N=%d spacecraft)\n",
                rs.counters.nbody / max(1, rs.nf), rs.counters.srp / max(1, rs.nf), TBC_N)
    end

    if TBC_PROFILE
        name = tbc_case_name("srp_nbody")
        Profile.clear()
        Profile.init(n = 20_000_000, delay = 0.001)
        withenv(ppc_mode_env_pairs(mode, full_cfg; outer_tasks=1)...) do
            @profile tbc_solve(name, full_cfg)
        end
        open(joinpath(TBC_STUDY_DIR, "profile_srp_nbody_$(TBC_N)sat_$(TBC_MISSION)s.txt"), "w") do io
            Profile.print(IOContext(io, :displaysize => (24, 2000)); format=:tree, C=true, maxdepth=60, mincount=50)
        end
        open(joinpath(TBC_STUDY_DIR, "profile_flat_srp_nbody_$(TBC_N)sat_$(TBC_MISSION)s.txt"), "w") do io
            Profile.print(IOContext(io, :displaysize => (24, 2000)); format=:flat, C=true, sortedby=:count, mincount=50)
        end
        println("profile written to ", TBC_STUDY_DIR)
    end
    println("load-at-end: ", strip(read(`uptime`, String)))
end

main()
