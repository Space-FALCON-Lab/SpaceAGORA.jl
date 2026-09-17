# Is flipping the default wind mode a no-op for the published reconstructions?
#
# Odyssey (Mars) and VEX (Venus) run with gram_perturbation_scales = 0. If
# perturbed winds under scales=0 are numerically equal to nominal winds, then
# changing the "auto" default from :perturbed to :nominal cannot move those
# results at all. Dump both modes over a grid and compare.

using Printf
const REPO = "/home/space-falcon-1/Documents/SpaceAGORA.jl"
import Pkg; Pkg.activate(REPO; io=devnull)
Base.find_package("GRAMSuite") === nothing && pushfirst!(LOAD_PATH, joinpath(REPO,"data","GRAMSuite.jl"))
using GRAMSuite

function build(planet, scales)
    kw = scales === nothing ? (;) : (; gram_perturbation_scales=scales)
    GRAMSuite.GRAMAtmosphereModel(;
        planet_name=planet,
        gram_root_directory=joinpath(REPO,"data","GRAMSuite.jl","GRAM Suite 2.0"),
        initial_time=GRAMSuite.InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.0),
        kw...)
end

function sweep(model)
    out = NTuple{5,Float64}[]
    for h in (60.0e3, 100.0e3, 140.0e3, 200.0e3), lat in (-60.0, -20.0, 35.0, 70.0),
        lon in (0.0, 90.0, 187.5, 300.0), t in (0.0, 137.0, 5000.0)
        rho, T, w = GRAMSuite.point_density_state(model, h, deg2rad(lat), deg2rad(lon), t, true)
        push!(out, (rho, T, w[1], w[2], w[3]))
    end
    return out
end

function run(planet, scales)
    m = build(planet, scales)
    pert = withenv("SPACEAGORA_GRAM_WIND_MODE" => "perturbed") do; sweep(m); end
    nom  = withenv("SPACEAGORA_GRAM_WIND_MODE" => "nominal")   do; sweep(m); end
    n = length(pert)
    maxw = 0.0; maxr = 0.0
    for i in 1:n
        maxr = max(maxr, abs(pert[i][1]-nom[i][1]) / max(nom[i][1], 1e-300))
        for k in 3:5
            maxw = max(maxw, abs(pert[i][k] - nom[i][k]))
        end
    end
    @printf("%-6s scales=%-7s n=%d  max |wind_pert - wind_nom| = %.6f m/s   max rel density diff = %.3e  %s\n",
            planet, scales === nothing ? "default" : string(scales), n, maxw, maxr,
            (maxw == 0.0 && maxr == 0.0) ? "IDENTICAL" : "differs")
end

let planet = get(ENV, "EQ_PLANET", "mars"),
    scales = get(ENV, "EQ_SCALES", "zero") == "zero" ? 0.0 : nothing
    Base.invokelatest(run, planet, scales)
end
