# Same path-dependence protocol, applied to WINDS.
#
# GRAMSuite._gram_density_state_native returns the PERTURBED winds by default
# (_gram_wind_mode() defaults to :perturbed via "auto"). Those come from
# PerturbedAtmosphere::updatePerturbations, whose random walk advances per call
# and is driven by `delta` = position difference from the PREVIOUS CALL on that
# instance. Winds feed relative velocity and hence drag.
#
# Run: julia --project=<repo> -t 1 wind_path_dependence.jl

using Printf, Random
const REPO = "/home/space-falcon-1/Documents/SpaceAGORA.jl"
import Pkg; Pkg.activate(REPO; io=devnull)
Base.find_package("GRAMSuite") === nothing && pushfirst!(LOAD_PATH, joinpath(REPO,"data","GRAMSuite.jl"))
using GRAMSuite

function build(planet)
    pert = get(ENV, "PD_PERT", "default")
    kw = pert == "zero" ? (; gram_perturbation_scales=0.0) : (;)
    return GRAMSuite.GRAMAtmosphereModel(;
        planet_name=planet,
        gram_root_directory=joinpath(REPO,"data","GRAMSuite.jl","GRAM Suite 2.0"),
        initial_time=GRAMSuite.InitialTime(year=2020,month=1,day=1,hour=0,minute=0,second=0.0),
        kw...)
end

# returns (rho, T, wind::SVector{3})
q(model, h, lat, lon, t) = GRAMSuite.point_density_state(model, h, lat, lon, t, true)

wnorm(w) = sqrt(w[1]^2 + w[2]^2 + w[3]^2)
absdiff(a, b) = maximum(abs.(a .- b))

function run(model, planet)
    LAT, LON, T = deg2rad(20.0), deg2rad(40.0), 100.0
    X   = (300.0e3, LAT, LON, T)
    FAR = (250.0e3, deg2rad(-40.0), deg2rad(200.0), T)

    println("=== $planet : WINDS (m/s)   perturbation_scales=$(get(ENV,"PD_PERT","default")) ===")

    _, _, w0 = q(model, X...)
    @printf("   wind at probe point: (%.3f, %.3f, %.3f)  |w| = %.3f m/s\n",
            w0[1], w0[2], w0[3], wnorm(w0))

    # 1. same point repeated
    ws = [q(model, X...)[3] for _ in 1:5]
    same = all(w -> w == ws[1], ws)
    @printf("1. same point x5 ............................ identical=%s  spread %.4f m/s\n",
            same, maximum(absdiff(w, ws[1]) for w in ws))

    # 2. X -> FAR -> X
    _,_,a = q(model, X...); q(model, FAR...); _,_,b = q(model, X...)
    @printf("2. X -> FAR -> X ............................ max wind diff %.4f m/s %s\n",
            absdiff(b,a), absdiff(b,a) == 0 ? "(pure)" : "<-- PATH DEPENDENT")

    # 3. the constellation pattern: X, then 15 other satellites, then X
    _,_,a = q(model, X...)
    for s in 2:16
        q(model, 300.0e3 + 1.0e3*s, LAT + deg2rad(3.0*s), LON + deg2rad(20.0*s), T)
    end
    _,_,b = q(model, X...)
    @printf("3. X -> 15 other sats -> X .................. max wind diff %.4f m/s %s\n",
            absdiff(b,a), absdiff(b,a) == 0 ? "(pure)" : "<-- PATH DEPENDENT")

    # 4. solver-stage replay: evaluate the SAME constellation state twice in a
    #    row, exactly as a rejected-then-retried RK stage would.
    sats = [(300.0e3 + 1.0e3*s, LAT + deg2rad(3.0*s), LON + deg2rad(20.0*s), T) for s in 1:16]
    pass1 = [q(model, s...)[3] for s in sats]
    pass2 = [q(model, s...)[3] for s in sats]
    dw = [absdiff(pass2[i], pass1[i]) for i in 1:16]
    @printf("4. identical 16-sat state evaluated twice ... max wind diff %.4f m/s  mean %.4f\n",
            maximum(dw), sum(dw)/16)

    # 5. does satellite 1's wind depend on how many OTHER satellites precede it?
    firsts = Float64[]
    for n in (1, 2, 4, 8, 16)
        # flush, then query n satellites, recording satellite 1's wind
        q(model, FAR...)
        w1 = q(model, sats[1]...)[3]
        push!(firsts, wnorm(w1))
        for i in 2:n; q(model, sats[i]...); end
    end
    @printf("5. |wind| of sat 1 with n=1,2,4,8,16 others:  %s\n",
            join((@sprintf("%.3f", f) for f in firsts), "  "))

    # 6. control: nominal (unperturbed) winds, same replay as test 4
    withenv("SPACEAGORA_GRAM_WIND_MODE" => "nominal") do
        p1 = [q(model, s...)[3] for s in sats]
        p2 = [q(model, s...)[3] for s in sats]
        d  = [absdiff(p2[i], p1[i]) for i in 1:16]
        @printf("6. SAME, with WIND_MODE=nominal ............. max wind diff %.4f m/s %s\n",
                maximum(d), maximum(d) == 0 ? "(deterministic)" : "<-- still varying")
    end
    println()
end

let planet = get(ENV, "PD_PLANET", "earth")
    m = build(planet)
    Base.invokelatest(run, m, planet)
end
