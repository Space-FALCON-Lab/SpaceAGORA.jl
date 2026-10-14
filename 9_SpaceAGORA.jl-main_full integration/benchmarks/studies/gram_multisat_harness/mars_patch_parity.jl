# Dump Mars density/temperature/wind over a grid of states, so the stock and
# patched libGRAM builds can be compared byte-for-byte.
#
#   GRAM_LIB_OVERRIDE=<path to patched libGRAM.so> julia ... mars_patch_parity.jl > patched.txt
#   julia ... mars_patch_parity.jl > stock.txt
#   diff stock.txt patched.txt

using Printf
const REPO = "/home/space-falcon-1/Documents/SpaceAGORA.jl"
import Pkg; Pkg.activate(REPO; io=devnull)
Base.find_package("GRAMSuite") === nothing && pushfirst!(LOAD_PATH, joinpath(REPO,"data","GRAMSuite.jl"))
using GRAMSuite

build(planet) = GRAMSuite.GRAMAtmosphereModel(;
    planet_name=planet,
    gram_root_directory=joinpath(REPO,"data","GRAMSuite.jl","GRAM Suite 2.0"),
    gram_library_path=get(ENV, "GRAM_LIB_OVERRIDE", ""),
    initial_time=GRAMSuite.InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.0))

function dump(model)
    for h in (60.0e3, 100.0e3, 140.0e3, 200.0e3, 300.0e3)
        for lat in (-60.0, -20.0, 0.0, 35.0, 70.0)
            for lon in (0.0, 90.0, 187.5, 300.0)
                for t in (0.0, 137.0, 5000.0, 86400.0)
                    rho, T, w = GRAMSuite.point_density_state(
                        model, h, deg2rad(lat), deg2rad(lon), t, true)
                    @printf("%.1f %.1f %.1f %.1f  %.17e %.17e %.17e %.17e %.17e\n",
                            h, lat, lon, t, rho, T, w[1], w[2], w[3])
                end
            end
        end
    end
end

let
    m = build(get(ENV, "PARITY_PLANET", "mars"))
    Base.invokelatest(dump, m)
end
