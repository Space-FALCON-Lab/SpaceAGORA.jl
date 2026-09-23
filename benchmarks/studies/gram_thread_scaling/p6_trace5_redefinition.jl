# P6 trace 5 before and after its redefinition, back to back in one process.
#
# `aero_<N>sat_l50_gram_lookahead_100s` used to be built on `ppc_constellation`
# with a 120 km entry interface; it is now built on
# `ppc_p6_gram_constellation` with a 600 km one, so that every member calls
# native GRAM and the look-ahead cache actually builds (see
# p6_case_audit.jl). This times the harness's current case against a rebuild of
# the previous definition from the same helpers, under the look-ahead
# environment the harness sets for this case (`_ppc_p6_gram_density_env!`),
# alternating the two across repeats. Only the ratio is meaningful; the thread
# count is the process's own.
#
# Usage (one process, memory-capped; run once per thread count):
#   systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
#     julia --project=. --threads=8 benchmarks/studies/gram_thread_scaling/p6_trace5_redefinition.jl \
#       [--n=256] [--repeats=3] [--out=results/p6_trace5_redefinition.csv]

using Printf

const GTS_DIR = @__DIR__
const GTS_REPO_ROOT = normpath(joinpath(GTS_DIR, "..", "..", ".."))
const GTS_PPC_DIR = joinpath(GTS_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")
include(joinpath(GTS_PPC_DIR, "cli.jl"))
include(joinpath(GTS_PPC_DIR, "modes.jl"))
include(joinpath(GTS_PPC_DIR, "cases.jl"))
ppc_ensure_gramsuite_loaded!()

const RS = SpaceAGORA.RuntimeServices

r5_arg(key, default) = begin
    for a in ARGS
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end
const R5_N = parse(Int, r5_arg("n", "256"))
const R5_REPEATS = parse(Int, r5_arg("repeats", "3"))
const R5_OUT = r5_arg("out", joinpath(GTS_DIR, "results", "p6_trace5_redefinition.csv"))
const R5_CASE = "aero_$(R5_N)sat_l50_gram_lookahead_100s"
const R5_MISSION = 100.0

# The harness's own look-ahead settings for this case, from _ppc_p6_gram_density_env!.
const R5_ENV = [
    "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "0",
    "SPACEAGORA_VACUUM_GRAM_CACHE" => "1",
    "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "20",
    "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => string(R5_MISSION + 500.0),
    "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "1e8",
]

current_config() = ppc_single_config(R5_CASE, PPCConfig(profile="full"))

function previous_config()
    planet = Earth("", PPC_SPICE_PATH)
    return ppc_build_config(
        planet=planet, spacecraft=ppc_constellation(planet, R5_N), mission_time_s=R5_MISSION,
        orientation_sim=false,
        dynamic_effectors=(ppc_harmonics_model(planet, 50), AerodynamicCoefficientfM()),
        density_model=ppc_gram_atmosphere_model("earth"), dt_max_orbit=5.0
    )
end

function solve(builder)
    args = builder()
    RS.reset_native_lock_stats!()
    t = @timed withenv(R5_ENV...) do
        SimulationEngine.run_simulation(args; isolate_state=false, return_solution=true)
    end
    r = t.value
    sol = r isa NamedTuple ? r.solution : r
    snap = RS.native_lock_stats_snapshot()
    return (wall=Float64(t.time), retcode=string(sol.retcode), nf=Int(sol.stats.nf),
            acq=snap.sites.gram_density.acquisitions)
end

function main()
    T = Threads.nthreads()
    @printf("case=%s threads=%d repeats=%d\n", R5_CASE, T, R5_REPEATS)
    # Warm-up: one untimed solve of each definition, so neither timed arm pays
    # compilation or the first native GRAM initialization.
    solve(current_config); solve(previous_config)
    rows = String[]
    best = Dict("current" => Inf, "previous" => Inf)
    for rep in 1:R5_REPEATS, (def, b) in (("previous", previous_config), ("current", current_config))
        r = solve(b)
        best[def] = min(best[def], r.wall)
        push!(rows, join((T, R5_N, def, rep, @sprintf("%.6f", r.wall), r.retcode, r.nf, r.acq), ","))
        @printf("rep%d %-8s wall=%8.3fs retcode=%s nf=%d gram_density acquisitions=%d\n",
                rep, def, r.wall, r.retcode, r.nf, r.acq)
        flush(stdout)
    end
    @printf("current/previous wall ratio (min over repeats): %.2f\n", best["current"] / best["previous"])
    mkpath(dirname(R5_OUT))
    exists = isfile(R5_OUT)
    open(R5_OUT, "a") do io
        exists || println(io, "threads,n_sats,definition,rep,wall_s,retcode,nf,gram_density_acq")
        foreach(r -> println(io, r), rows)
    end
end

main()
