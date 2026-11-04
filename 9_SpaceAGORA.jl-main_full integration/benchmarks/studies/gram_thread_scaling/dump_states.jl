# WS11b step 2, trajectory half: the full state history of a native-GRAM
# constellation, locked against isolated-pool, for a byte-for-byte `cmp`.
#
# The density-value grid (`density_grid_parity.jl`) answers whether two GRAM
# instances agree on one sample. This answers whether a whole solve agrees, with
# every saved time and every component of every spacecraft written as raw
# Float64 in solver order -- the same dump format
# `benchmarks/studies/third_body_cost/variants.jl --dump` writes, so the two are
# compared the same way.
#
# Two configurations, because they answer different questions:
#
#   reference  the WS11 reference case: `ppc_constellation` at 64 spacecraft,
#              L50 harmonics plus aero, live native GRAM, the look-ahead cache
#              settings the P6 trace uses, entry interface 120 km. This is the
#              shape the other WS11 workstreams dump. Note what it does NOT do:
#              every member starts above the 120 km entry interface, so
#              `in_atmosphere` is false for all of them and the look-ahead cache
#              never builds; and the density callback's automatic minimum thread
#              budget of 16 keeps the pool from engaging on an 8-thread machine.
#              So on this configuration the pool is expected to change nothing
#              because it never runs -- which is worth recording, and is not
#              evidence about the pool.
#
#   engaged    the configuration under which the pool actually evaluates GRAM:
#              the low constellation from `run_scaling.jl`, entry interface
#              above it, and the density callback's minimum thread budget
#              lowered so the pooled batch call clears its `workers > 1` guard.
#              This is the dump that carries the bit-identity claim.
#
# Usage:
#   julia --project=. --threads=4 \
#       benchmarks/studies/gram_thread_scaling/dump_states.jl \
#       --config=engaged --pool=4 --n=64 --mission=100 --out=/path/prefix
#
# Then: cmp <prefix>_locked.bin <prefix>_pool4.bin

using Printf

const GTS_DIR = @__DIR__
const GTS_REPO_ROOT = normpath(joinpath(GTS_DIR, "..", "..", ".."))
const GTS_PPC_DIR = joinpath(GTS_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(GTS_PPC_DIR, "cli.jl"))
include(joinpath(GTS_PPC_DIR, "modes.jl"))
include(joinpath(GTS_PPC_DIR, "cases.jl"))
include(joinpath(GTS_DIR, "run_scaling_config.jl"))

ppc_ensure_gramsuite_loaded!()

using ComponentArrays: getdata

gtsd_arg(args, key, default) = begin
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end

const GTSD_CONFIG  = gtsd_arg(ARGS, "config", "engaged")
const GTSD_N       = parse(Int, gtsd_arg(ARGS, "n", "64"))
const GTSD_MISSION = parse(Float64, gtsd_arg(ARGS, "mission", "100"))
const GTSD_POOL    = parse(Int, gtsd_arg(ARGS, "pool", "4"))
const GTSD_DENSITY = gtsd_arg(ARGS, "density", "lookahead")
const GTSD_OUT     = gtsd_arg(ARGS, "out", joinpath(GTS_DIR, "results", "dump"))

function gtsd_config(n::Int, mission_s::Float64)
    if GTSD_CONFIG == "reference"
        planet = Earth("", PPC_SPICE_PATH)
        return ppc_build_config(
            planet=planet,
            spacecraft=ppc_constellation(planet, n),
            mission_time_s=mission_s,
            orientation_sim=false,
            dynamic_effectors=(ppc_harmonics_model(planet, 50), AerodynamicCoefficientfM()),
            density_model=ppc_gram_atmosphere_model("earth"),
            dt_max_orbit=5.0
        )
    elseif GTSD_CONFIG == "engaged"
        return gts_build_config(n, mission_s, GTS_ENGAGED_EI_KM)
    end
    error("Unknown --config='$(GTSD_CONFIG)'. Use reference or engaged.")
end

function gtsd_env(pool::Int)
    env = gts_density_env(GTSD_DENSITY, GTSD_MISSION)
    append!(env, gts_pool_env(pool))
    # The width override is for the study's own arms only. A "default" dump
    # (pool < 0) must see exactly what a user sees, overrides included.
    GTSD_CONFIG == "engaged" && pool >= 0 && append!(env, GTS_WIDTH_ENV)
    return env
end

function gtsd_run(label::String, pool::Int)
    args = gtsd_config(GTSD_N, GTSD_MISSION)
    r = withenv(gtsd_env(pool)...) do
        SimulationEngine.run_simulation(args; isolate_state=false, return_solution=true)
    end
    sol = r isa NamedTuple ? r.solution : r
    path = "$(GTSD_OUT)_$(label).bin"
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, Float64.(sol.t))
        for u in sol.u
            write(io, Float64.(vec(getdata(u))))
        end
    end
    @printf("%-8s saved %d times x %d state components -> %s (%d bytes)\n",
            label, length(sol.t), length(vec(getdata(sol.u[1]))), path, filesize(path))
    flush(stdout)
    return path
end

function main()
    @printf("config=%s n=%d mission=%.0fs density=%s pool=%d threads=%d\n",
            GTSD_CONFIG, GTSD_N, GTSD_MISSION, GTSD_DENSITY, GTSD_POOL, Threads.nthreads())
    # A warm-up solve first: the dumps must not differ because one of them paid
    # the first native GRAM initialization and the other did not.
    withenv(gtsd_env(0)...) do
        SimulationEngine.run_simulation(gtsd_config(min(GTSD_N, 16), 5.0); isolate_state=false)
    end
    locked = gtsd_run("locked", 0)
    pooled = gtsd_run(GTSD_POOL < 0 ? "default" : "pool$(GTSD_POOL)", GTSD_POOL)
    println("compare with: cmp $(locked) $(pooled)")
    return nothing
end

main()
