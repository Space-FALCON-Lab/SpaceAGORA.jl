# Control experiment for multi_16_gram_live.
#
# Same constellation, same harmonics, same tolerances, same dt_max -- the ONLY
# difference is the density model. If the exponential variant is fast and the
# GRAM variant is not, the cost is in the GRAM coupling; if both are slow, it is
# the constellation configuration and GRAM is incidental.
#
# Short missions with solver step counts, so every cell is reachable.

const REPO = "/home/space-falcon-1/Documents/SpaceAGORA.jl"
push!(ARGS, "--case=multi_16_gram_live")
using Printf
include(joinpath(REPO, "benchmarks", "studies", "parallelization_performance", "cli.jl"))
include(joinpath(REPO, "benchmarks", "studies", "parallelization_performance", "modes.jl"))
include(joinpath(REPO, "benchmarks", "studies", "parallelization_performance", "cases.jl"))

function rss_mb()
    for line in eachline("/proc/self/status")
        startswith(line, "VmRSS:") && return parse(Int, split(line)[2]) / 1024
    end
    return NaN
end

function build_case(; density, reltol, mission_s, nsat)
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
        reltol_orbit=reltol,
        abstol_orbit=reltol,
    )
end

function trial(; density, reltol, mission_s, nsat)
    args = build_case(density=density, reltol=reltol, mission_s=mission_s, nsat=nsat)
    GC.gc(true); r0 = rss_mb()
    t0 = time_ns()
    res = try
        SimulationEngine.run_simulation(args; isolate_state=false,
                                        return_solution=true, return_solver_metadata=true)
    catch e
        @printf("  %-5s N=%-3d reltol=%-6.0e mission=%-5.0fs  ERROR after %6.1fs: %s\n",
                density, nsat, reltol, mission_s, (time_ns()-t0)*1e-9,
                first(sprint(showerror, e), 70))
        flush(stdout); return
    end
    dt = (time_ns()-t0)*1e-9
    sol = res.solution
    st = try res.solver_metadata.stats catch; nothing end
    nacc = st === nothing ? missing : try st.naccept catch; missing end
    nrej = st === nothing ? missing : try st.nreject catch; missing end
    nf   = st === nothing ? missing : try st.nf      catch; missing end
    GC.gc(true)
    @printf("  %-5s N=%-3d reltol=%-6.0e mission=%-5.0fs  %7.2fs  %-8s accept=%-7s reject=%-7s rhs=%-8s RSS %+6.0fMB\n",
            density, nsat, reltol, mission_s, dt, string(sol.retcode),
            string(nacc), string(nrej), string(nf), rss_mb()-r0)
    flush(stdout)
end

function main()
    println("Control: identical config, density model swapped. dt_max=5 s, deg/ord-20 harmonics.\n")
    # warm up compilation on both paths so the reported times are steady-state
    trial(density=:exp, reltol=1e-9, mission_s=30.0, nsat=16)
    trial(density=:gram, reltol=1e-9, mission_s=30.0, nsat=16)
    println("--- steady state, N=16, real 1200 s mission ---")
    for _ in 1:2
        trial(density=:exp,  reltol=1e-9, mission_s=1200.0, nsat=16)
        trial(density=:gram, reltol=1e-9, mission_s=1200.0, nsat=16)
    end
end

Base.invokelatest(main)
