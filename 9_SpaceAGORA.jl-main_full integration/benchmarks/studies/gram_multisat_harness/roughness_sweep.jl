# Characterise the GRAM density field's smoothness in altitude: is it a smooth
# function with a physical gradient, or a staircase? An adaptive ODE solver
# needs the former.

using Printf
const REPO = "/home/space-falcon-1/Documents/SpaceAGORA.jl"
import Pkg; Pkg.activate(REPO; io=devnull)
Base.find_package("GRAMSuite") === nothing && pushfirst!(LOAD_PATH, joinpath(REPO,"data","GRAMSuite.jl"))
using GRAMSuite

build(planet) = GRAMSuite.GRAMAtmosphereModel(;
    planet_name=planet,
    gram_root_directory=joinpath(REPO,"data","GRAMSuite.jl","GRAM Suite 2.0"),
    initial_time=GRAMSuite.InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.0))

function sweep(model, planet)
    lat, lon, t = deg2rad(20.0), deg2rad(40.0), 100.0
    h0 = 300.0e3

    for (span, n, label) in ((1.0, 201, "1 m"), (100.0, 201, "100 m"), (10_000.0, 201, "10 km"))
        hs = range(h0, h0 + span, length=n)
        rho = [GRAMSuite.point_density_state(model, h, lat, lon, t, true)[1] for h in hs]

        # how many DISTINCT values? a staircase has few
        ndistinct = length(unique(rho))
        # successive relative differences
        d = [abs(rho[i+1]-rho[i])/rho[i] for i in 1:n-1]
        nz = count(>(0.0), d)
        # physically expected change across the span (scale height ~7 km)
        expected = span / 7000.0

        @printf("  span %-6s: distinct=%3d/%d   nonzero steps=%3d   ", label, ndistinct, n, nz)
        @printf("total rel change=%.3e (smooth expectation %.3e)  max single jump=%.3e\n",
                abs(rho[end]-rho[1])/rho[1], expected, maximum(d))
    end

    # Same probe in TIME (what the integrator actually advances) at fixed position
    println("  --- in elapsed time, fixed position ---")
    for (span, label) in ((1.0e-3, "1 ms"), (1.0, "1 s"), (60.0, "60 s"))
        ts = range(t, t + span, length=51)
        rho = [GRAMSuite.point_density_state(model, h0, lat, lon, tt, true)[1] for tt in ts]
        d = [abs(rho[i+1]-rho[i])/rho[i] for i in 1:length(rho)-1]
        @printf("  span %-5s: distinct=%2d/51  max single jump=%.3e\n",
                label, length(unique(rho)), maximum(d))
    end
end

let planet = get(ENV, "RHS_PLANET", "earth")
    m = build(planet)
    println("--- $planet density field smoothness ---")
    Base.invokelatest(sweep, m, planet)
end
