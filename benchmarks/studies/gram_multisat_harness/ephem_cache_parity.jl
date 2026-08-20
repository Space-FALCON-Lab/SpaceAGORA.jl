# Epoch-level ephemeris cache: parity + call-count check.
#
# Goes through the real SpaceAGORA path (extension loaded), so the injected
# ephemeris hook is live. Dumps density/temperature/wind over a constellation
# pattern -- many satellites per epoch, several epochs -- with the cache on and
# off. The two must be bit-identical.

using Printf
const REPO = "/home/space-falcon-1/Documents/SpaceAGORA.jl"
import Pkg; Pkg.activate(REPO; io=devnull)
push!(ARGS, "--case=multi_16_gram_live")
include(joinpath(REPO, "benchmarks", "studies", "parallelization_performance", "cli.jl"))
include(joinpath(REPO, "benchmarks", "studies", "parallelization_performance", "modes.jl"))
include(joinpath(REPO, "benchmarks", "studies", "parallelization_performance", "cases.jl"))

const EM = SpaceAGORA.SimulationModel.EnvironmentModels

function dump(model, tag)
    io = IOBuffer()
    # 16 "satellites" per epoch, 12 epochs -- the constellation access pattern.
    for k in 0:11
        el = 5.0 * k
        for s in 1:16
            lat = deg2rad(50.0 * sin(2pi * (k / 11.0) + 2pi * s / 16))
            lon = deg2rad(mod(360.0 * (k / 11.0) + 360.0 * s / 16, 360.0))
            h = 300.0e3 + 20.0e3 * sin(2pi * (k / 6.0) + s)
            rho, T, w = GRAMSuite.point_density_state(model.core, h, lat, lon, el, true)
            @printf(io, "%.17e %.17e %.17e %.17e %.17e\n", rho, T, w[1], w[2], w[3])
        end
    end
    return String(take!(io))
end

function timed(model, n)
    t0 = time_ns()
    for k in 0:(n-1)
        el = 5.0 * k
        for s in 1:16
            lat = deg2rad(50.0 * sin(2pi * (k / 11.0) + 2pi * s / 16))
            lon = deg2rad(mod(360.0 * (k / 11.0) + 360.0 * s / 16, 360.0))
            GRAMSuite.point_density_state(model.core, 300.0e3, lat, lon, el, true)
        end
    end
    return (time_ns() - t0) * 1e-9
end

function main()
    planet = get(ENV, "EC_PLANET", "earth")
    # Constructing a planet furnishes the SPICE kernels (planets.jl `_furnsh_once`);
    # without it utc2et fails for want of a leapseconds kernel.
    planet == "mars" ? Mars("", PPC_SPICE_PATH) : Earth("", PPC_SPICE_PATH)
    model = ppc_gram_atmosphere_model(planet)

    on  = withenv("SPACEAGORA_GRAM_EPHEMERIS_CACHE" => "on")  do; dump(model, "on");  end
    off = withenv("SPACEAGORA_GRAM_EPHEMERIS_CACHE" => "off") do; dump(model, "off"); end
    @printf("%-6s parity: %d chars each, %s\n", planet, length(on),
            on == off ? "BIT-IDENTICAL" : "*** DIFFER ***")

    # warm
    withenv("SPACEAGORA_GRAM_EPHEMERIS_CACHE" => "on") do; timed(model, 5); end
    t_on  = withenv("SPACEAGORA_GRAM_EPHEMERIS_CACHE" => "on")  do; timed(model, 60); end
    t_off = withenv("SPACEAGORA_GRAM_EPHEMERIS_CACHE" => "off") do; timed(model, 60); end
    n = 60 * 16
    @printf("%-6s %d evals: cache off %.3f s (%.0f/s)   cache on %.3f s (%.0f/s)   speedup %.2fx\n",
            planet, n, t_off, n/t_off, t_on, n/t_on, t_off/t_on)
end

Base.invokelatest(main)
