# Steady-state (post-compilation) timing of the real multi_16_gram_live case.
#
# The benchmark harness reports ~50 s wall for this case, but almost all of that
# is Julia compilation: the solve itself is ~1.4 s. Measuring the effect of any
# GRAM-side change on the harness number therefore understates it by ~35x. This
# warms both paths first, then times repeated full-mission solves.

const REPO = "/home/space-falcon-1/Documents/SpaceAGORA.jl"
push!(ARGS, "--case=multi_16_gram_live")
using Printf, Statistics
include(joinpath(REPO, "benchmarks", "studies", "parallelization_performance", "cli.jl"))
include(joinpath(REPO, "benchmarks", "studies", "parallelization_performance", "modes.jl"))
include(joinpath(REPO, "benchmarks", "studies", "parallelization_performance", "cases.jl"))

function build(; density, nsat=16, mission_s=1200.0)
    planet = Earth("", PPC_SPICE_PATH)
    dm = density == :gram ? ppc_gram_atmosphere_model("earth") : ExponentialAtmosphereModel(planet)
    return ppc_build_config(
        planet=planet,
        spacecraft=ppc_constellation(planet, nsat),
        mission_time_s=mission_s,
        orientation_sim=false,
        dynamic_effectors=(ppc_harmonics_model(planet, 20), AerodynamicCoefficientfM()),
        density_model=dm,
        dt_max_orbit=5.0,
    )
end

function solve_time(args)
    t0 = time_ns()
    res = SimulationEngine.run_simulation(args; isolate_state=false,
                                          return_solution=true, return_solver_metadata=true)
    dt = (time_ns() - t0) * 1e-9
    return dt, string(res.solution.retcode), res.solution.u[end].sc[1]
end

function bench(label, envs; density=:gram, reps=5)
    withenv(envs...) do
        args = build(density=density)
        solve_time(args)                       # warm this configuration
        ts = Float64[]
        local rc, sc
        for _ in 1:reps
            t, rc, sc = solve_time(args)
            push!(ts, t)
        end
        pos = sqrt(sum(abs2, sc.pos))
        @printf("  %-34s %6.3f s  (min %.3f, max %.3f)  %s  pos=%.10e\n",
                label, median(ts), minimum(ts), maximum(ts), rc, pos)
        return median(ts)
    end
end

function main()
    @printf("multi_16_gram_live, N=16, full 1200 s mission, threads=%d\n", Threads.nthreads())
    println("steady-state solve time (compilation excluded)\n")

    # global warm-up for both density paths
    withenv("SPACEAGORA_GRAM_EPHEMERIS_CACHE" => "on") do
        solve_time(build(density=:exp)); solve_time(build(density=:gram))
    end

    t_exp  = bench("exponential atmosphere", ("SPACEAGORA_GRAM_EPHEMERIS_CACHE" => "on",); density=:exp)
    t_off  = bench("GRAM, ephemeris cache OFF", ("SPACEAGORA_GRAM_EPHEMERIS_CACHE" => "off",))
    t_on   = bench("GRAM, ephemeris cache ON", ("SPACEAGORA_GRAM_EPHEMERIS_CACHE" => "on",))
    t_pool = bench("GRAM, cache ON + threaded pool",
                   ("SPACEAGORA_GRAM_EPHEMERIS_CACHE" => "on",
                    "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on",
                    "SPACEAGORA_GRAM_ISOLATED_POOL" => "on",
                    "SPACEAGORA_GRAM_LOCK_SCOPE" => "model"))

    println()
    @printf("  GRAM share of solve (cache off): %.1f%%\n", 100 * (t_off - t_exp) / t_off)
    @printf("  ephemeris cache speedup:         %.2fx  (%.3f -> %.3f s)\n", t_off / t_on, t_off, t_on)
    @printf("  + threaded pool:                 %.2fx  (%.3f -> %.3f s)\n", t_off / t_pool, t_off, t_pool)
end

Base.invokelatest(main)
