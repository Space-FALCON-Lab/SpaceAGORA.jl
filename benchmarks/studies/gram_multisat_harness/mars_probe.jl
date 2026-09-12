# Isolate where the Mars clone path breaks: cloning itself, using a clone
# serially, or only under concurrency.

using Printf, Base.Threads
const REPO = "/home/space-falcon-1/Documents/SpaceAGORA.jl"
import Pkg; Pkg.activate(REPO; io=devnull)
Base.find_package("GRAMSuite") === nothing && pushfirst!(LOAD_PATH, joinpath(REPO,"data","GRAMSuite.jl"))
using GRAMSuite

build(planet) = GRAMSuite.GRAMAtmosphereModel(
    planet_name=planet,
    gram_root_directory=joinpath(REPO,"data","GRAMSuite.jl","GRAM Suite 2.0"),
    gram_library_path=get(ENV, "GRAM_LIB_OVERRIDE", ""),
    initial_time=GRAMSuite.InitialTime(year=2020,month=6,day=1,hour=12,minute=0,second=0.0))

function probe(model, planet)
    G  = getfield(GRAMSuite, :_GRAM_WRAPPER)[]
    gf(n) = Base.invokelatest(getfield, G, n)
    setpos, upd, dyn = gf(:set_position!), gf(:update!), gf(:get_dynamics_state)
    getephem, setephem = gf(:get_ephemeris_state), gf(:set_ephemeris_state!)
    EphT, copyatm = gf(:EphemerisStateC), gf(:copy_atmosphere)
    atm0 = model.gram_atmosphere

    native(atm,h,lat,lon,t) = (setpos(atm;height=h*1e-3,latitude=lat,longitude=lon,elapsed_time=t);
                               upd(atm); dyn(atm).density)

    println("1. native call on original .......... ", native(atm0, 200e3, 10.0, 20.0, 0.0))
    base = getephem(atm0)
    println("   base ephemeris: solarTime=", round(base.solarTime,digits=4),
            " ssLat=", round(base.subsolarLatitude,digits=4),
            " ssLon=", round(base.subsolarLongitude,digits=4))

    println("2. injected call on original ....... ")
    setephem(atm0, EphT(base.solarTime, base.longitudeSun, base.subsolarLatitude,
                        base.subsolarLongitude, base.orbitalRadius, base.oneWayLightTime,
                        base.solarZenithAngle, base.secondsPerSol))
    println("   rho = ", native(atm0, 200e3, 10.0, 20.0, 0.0))

    println("3. clone one atmosphere ............ ")
    c1 = copyatm(atm0); println("   ok, ptr=", c1.ptr, " body=", c1.body)

    println("4. NATIVE call on the clone ........ ")
    flush(stdout)
    println("   rho = ", native(c1, 200e3, 10.0, 20.0, 0.0))

    println("5. INJECTED call on the clone ...... ")
    flush(stdout)
    setephem(c1, EphT(base.solarTime, base.longitudeSun, base.subsolarLatitude,
                      base.subsolarLongitude, base.orbitalRadius, base.oneWayLightTime,
                      base.solarZenithAngle, base.secondsPerSol))
    println("   rho = ", native(c1, 200e3, 10.0, 20.0, 0.0))

    println("6. 500 injected calls on the clone (serial) ")
    flush(stdout)
    for i in 1:500
        setephem(c1, EphT(base.solarTime, base.longitudeSun, base.subsolarLatitude,
                          base.subsolarLongitude, base.orbitalRadius, base.oneWayLightTime,
                          base.solarZenithAngle, base.secondsPerSol))
        native(c1, 200e3 + i, 10.0, 20.0, Float64(i))
    end
    println("   ok")

    println("7. two clones, concurrent injected calls ")
    flush(stdout)
    c2 = copyatm(atm0)
    cl = (c1, c2)
    @sync for w in 1:2
        Threads.@spawn for i in 1:500
            setephem(cl[w], EphT(base.solarTime, base.longitudeSun, base.subsolarLatitude,
                                 base.subsolarLongitude, base.orbitalRadius,
                                 base.oneWayLightTime, base.solarZenithAngle, base.secondsPerSol))
            native(cl[w], 200e3 + i, 10.0, 20.0, Float64(i))
        end
    end
    println("   ok -- concurrent injected calls on separate clones are safe")

    # ---- 8: varying positions (benchmark-like), spawned ONCE ----
    track(sat,k,nsat) = (300.0e3 + 20.0e3*sin(2pi*(k/550.0)+sat),
                         50.0*sin(2pi*(k/1100.0) + 2pi*sat/nsat),
                         mod(360.0*(k/1100.0) + 360.0*sat/nsat, 360.0),
                         5.0*k)
    println("8. varying positions, 2 clones, spawned once ")
    flush(stdout)
    @sync for w in 1:2
        Threads.@spawn for k in 1:200, sat in w:2:16
            h, lat, lon, t = track(sat, k, 16)
            setephem(cl[w], EphT(base.solarTime, base.longitudeSun, base.subsolarLatitude,
                                 base.subsolarLongitude, base.orbitalRadius,
                                 base.oneWayLightTime, base.solarZenithAngle, base.secondsPerSol))
            native(cl[w], h, lat, lon, t)
        end
    end
    println("   ok")

    # ---- 9: varying positions, re-spawned per step (exactly the benchmark) ----
    println("9. varying positions, 2 clones, re-spawned per step ")
    flush(stdout)
    for k in 1:200
        @sync for w in 1:2
            Threads.@spawn for sat in w:2:16
                h, lat, lon, t = track(sat, k, 16)
                setephem(cl[w], EphT(base.solarTime, base.longitudeSun, base.subsolarLatitude,
                                     base.subsolarLongitude, base.orbitalRadius,
                                     base.oneWayLightTime, base.solarZenithAngle, base.secondsPerSol))
                native(cl[w], h, lat, lon, t)
            end
        end
    end
    println("   ok")
end

let planet = get(ENV, "PROBE_PLANET", "mars")
    m = build(planet)
    println("planet=$planet threads=$(nthreads())")
    Base.invokelatest(probe, m, planet)
end
