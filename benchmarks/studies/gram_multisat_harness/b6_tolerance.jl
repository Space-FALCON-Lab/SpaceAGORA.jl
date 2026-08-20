# Does the multi_16_gram_live failure come from the solver being unable to meet
# reltol against a rough density field, rather than from a leak in the GRAM
# binding? Rebuild the case's exact configuration at several tolerances and
# mission lengths and watch wall time, solver steps, and RSS.
#
# Run: julia --project=<repo> -t 1 b6_tolerance.jl

const REPO = "/home/space-falcon-1/Documents/SpaceAGORA.jl"
push!(ARGS, "--case=multi_16_gram_live")   # makes cases.jl load GRAMSuite up front
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

# The multi_16_gram_live case body, with reltol/abstol and mission time exposed.
function build_case(; reltol, mission_s)
    planet = Earth("", PPC_SPICE_PATH)
    return ppc_build_config(
        planet=planet,
        spacecraft=ppc_constellation(planet, 16),
        mission_time_s=mission_s,
        orientation_sim=false,
        dynamic_effectors=(ppc_harmonics_model(planet, 20), AerodynamicCoefficientfM()),
        density_model=ppc_gram_atmosphere_model("earth"),
        dt_max_orbit=5.0,
        reltol_orbit=reltol,
        abstol_orbit=reltol,
    )
end

function trial(; reltol, mission_s, budget_s)
    args = build_case(reltol=reltol, mission_s=mission_s)
    GC.gc(true); r0 = rss_mb()
    t0 = time_ns()
    local res
    try
        res = SimulationEngine.run_simulation(args; isolate_state=false,
                                              return_solution=true,
                                              return_solver_metadata=true)
    catch e
        dt = (time_ns()-t0)*1e-9
        @printf("  reltol=%-7.0e mission=%-6.0fs  ERROR after %.1f s: %s\n",
                reltol, mission_s, dt, sprint(showerror, e)[1:min(end,90)])
        return nothing
    end
    dt = (time_ns()-t0)*1e-9
    sol = res.solution
    nsteps = length(sol.t)
    nacc = try res.solver_metadata.stats.naccept catch; missing end
    nrej = try res.solver_metadata.stats.nreject catch; missing end
    nf   = try res.solver_metadata.stats.nf      catch; missing end
    GC.gc(true)
    @printf("  reltol=%-7.0e mission=%-6.0fs  %7.1f s  retcode=%-8s saved=%-5d accept=%-8s reject=%-8s rhs=%-9s RSS %+7.0f MB\n",
            reltol, mission_s, dt, string(sol.retcode), nsteps,
            string(nacc), string(nrej), string(nf), rss_mb()-r0)
    return dt
end

# One trial per process, so a hang at tight tolerance cannot starve the rest of
# the sweep and so each line is flushed as soon as it is produced.
function main()
    rt = parse(Float64, get(ENV, "B6_RELTOL", "1e-9"))
    ms = parse(Float64, get(ENV, "B6_MISSION", "60"))
    trial(reltol=rt, mission_s=ms, budget_s=0)
    flush(stdout)
end

Base.invokelatest(main)
