# Is GRAM's returned density a pure function of (position, time), or does it
# depend on the CALL HISTORY of the instance?
#
# PerturbedAtmosphere::setPosition computes `delta` against previousPosition,
# and updatePerturbations accumulates `perturbationStep` from that delta. If any
# of that reaches the density SpaceAGORA reads, then querying one instance for N
# satellites -- or re-querying it across solver stages -- makes the ODE
# right-hand side non-deterministic, which an adaptive integrator cannot
# distinguish from extreme stiffness.
#
# Run: julia --project=<repo> -t 1 path_dependence.jl

using Printf, Random
const REPO = "/home/space-falcon-1/Documents/SpaceAGORA.jl"
import Pkg; Pkg.activate(REPO; io=devnull)
Base.find_package("GRAMSuite") === nothing && pushfirst!(LOAD_PATH, joinpath(REPO,"data","GRAMSuite.jl"))
using GRAMSuite

build(planet) = GRAMSuite.GRAMAtmosphereModel(;
    planet_name=planet,
    gram_root_directory=joinpath(REPO,"data","GRAMSuite.jl","GRAM Suite 2.0"),
    initial_time=GRAMSuite.InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.0))

rho(model, h, lat, lon, t) = GRAMSuite.point_density_state(model, h, lat, lon, t, true)[1]

reldiff(a, b) = abs(a - b) / max(abs(a), 1e-300)

function run(model, planet)
    LAT, LON, T = deg2rad(20.0), deg2rad(40.0), 100.0
    X = (300.0e3, LAT, LON, T)                     # the probe point
    FAR = (250.0e3, deg2rad(-40.0), deg2rad(200.0), T)   # somewhere else entirely

    println("=== $planet ===")

    # ---- 1. same point, repeated (delta = 0 every time) ----
    a = [rho(model, X...) for _ in 1:5]
    @printf("1. same point x5 ............................ identical=%s\n", all(==(a[1]), a))

    # ---- 2. A, then a FAR point, then A again ----
    r1 = rho(model, X...)
    rho(model, FAR...)
    r2 = rho(model, X...)
    @printf("2. X -> FAR -> X ............................ rel diff %.3e %s\n",
            reldiff(r2, r1), reldiff(r2,r1) == 0 ? "(pure)" : "<-- PATH DEPENDENT")

    # ---- 3. A, then 15 other 'satellites', then A again (the constellation pattern) ----
    r1 = rho(model, X...)
    for s in 2:16
        rho(model, 300.0e3 + 1.0e3*s, LAT + deg2rad(3.0*s), LON + deg2rad(20.0*s), T)
    end
    r2 = rho(model, X...)
    @printf("3. X -> 15 other sats -> X .................. rel diff %.3e %s\n",
            reldiff(r2, r1), reldiff(r2,r1) == 0 ? "(pure)" : "<-- PATH DEPENDENT")

    # ---- 4. the roughness sweep, monotonic vs shuffled order ----
    # If the "staircase" is a property of the FIELD, a given altitude returns the
    # same density regardless of the order in which altitudes are visited.
    hs = collect(range(300.0e3, 300.0e3 + 1.0, length=201))
    mono = [rho(model, h, LAT, LON, T) for h in hs]

    perm = shuffle(MersenneTwister(42), 1:201)
    shuf = Vector{Float64}(undef, 201)
    for i in perm
        shuf[i] = rho(model, hs[i], LAT, LON, T)
    end

    d = [reldiff(shuf[i], mono[i]) for i in 1:201]
    @printf("4. 1 m sweep monotonic vs shuffled .......... max rel diff %.3e  mean %.3e %s\n",
            maximum(d), sum(d)/201, maximum(d) == 0 ? "(field property)" : "<-- ORDER ARTIFACT")

    # jump structure within each ordering
    jm = [reldiff(mono[i+1], mono[i]) for i in 1:200]
    js = [reldiff(shuf[i+1], shuf[i]) for i in 1:200]
    @printf("   largest single step: monotonic %.3e   shuffled-then-reordered %.3e\n",
            maximum(jm), maximum(js))

    # ---- 5. one FRESH instance per point (no shared history at all) ----
    # Approximated by re-visiting each altitude after a long excursion, which
    # resets the accumulated perturbation step the same way a fresh call would.
    isolated = Float64[]
    for h in hs
        rho(model, FAR...)          # flush the history
        push!(isolated, rho(model, h, LAT, LON, T))
    end
    di = [reldiff(isolated[i], mono[i]) for i in 1:201]
    ji = [reldiff(isolated[i+1], isolated[i]) for i in 1:200]
    @printf("5. history-flushed sweep vs monotonic ....... max rel diff %.3e   largest single step %.3e\n",
            maximum(di), maximum(ji))
    println()
end

let planet = get(ENV, "PD_PLANET", "earth")
    m = build(planet)
    Base.invokelatest(run, m, planet)
end
