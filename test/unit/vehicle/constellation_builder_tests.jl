module ConstellationBuilderTests

# Orbital planes, the Walker patterns and the constellation builder; the mean
# anomaly in initial conditions; spacecraft ids and slots through a run.

using Test
using Logging
using LinearAlgebra
using StaticArrays
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine

# A small bus, built anew on every call.
function bus_factory(plane, slot, index)
    bus = Link(root=true, m=100.0 + index)
    return SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia)
end

# The mean anomaly (degrees) at true anomaly ν (degrees), computed directly.
function mean_anomaly_deg(ν, e)
    E = 2atan(sqrt((1 - e) / (1 + e)) * tan(deg2rad(ν) / 2))
    return mod(rad2deg(E - e * sin(E)), 360)
end

angle_gap(x, y) = abs(mod(x - y + 180, 360) - 180)

# The 50-step Newton loop the conversion used alone before, with whether it
# converged and 1 - e cos E at its end.
function newton_true_anomaly(M, e)
    M = mod(M, 2pi)
    E = e < 0.8 ? M : Float64(pi)
    converged = false
    for _ in 1:50
        step = (E - e * sin(E) - M) / (1.0 - e * cos(E))
        E -= step
        abs(step) <= 1e-14 && (converged = true; break)
    end
    ν = mod(2.0 * atan(sqrt(1.0 + e) * sin(E / 2.0), sqrt(1.0 - e) * cos(E / 2.0)), 2pi)
    return ν, converged, 1.0 - e * cos(E)
end

# The true anomaly (rad) at mean anomaly 0 < M <= π, inverting Kepler's equation
# by bisection in BigFloat.
function bigfloat_true_anomaly(M, e)
    return setprecision(BigFloat, 256) do
        m, ε = BigFloat(M), BigFloat(e)
        lo, hi = big(0.0), big(pi)
        for _ in 1:2000
            mid = (lo + hi) / 2
            mid - ε * sin(mid) < m ? (lo = mid) : (hi = mid)
            hi - lo <= big(2.0)^-240 * hi && break
        end
        E = (lo + hi) / 2
        Float64(2 * atan(sqrt(1 + ε) * sin(E / 2), sqrt(1 - ε) * cos(E / 2)))
    end
end

# Inclination, RAAN and argument of latitude (degrees) of an inertial state.
function plane_angles(r, v)
    h = cross(r, v)
    inclination = acosd(h[3] / norm(h))
    raan = mod(atand(h[1], -h[2]), 360)
    node = SVector(cosd(raan), sind(raan), 0.0)
    latitude = mod(atand(dot(cross(node, r), h) / norm(h), dot(node, r)), 360)
    return inclination, raan, latitude
end

# The warnings `f` gives; under --depwarn=error the deprecation throws instead.
function deprecation_messages(f)
    if Base.JLOptions().depwarn == 2
        try
            f()
            return String[]
        catch err
            return [sprint(showerror, err)]
        end
    end
    logger = Test.TestLogger(min_level=Logging.Warn)
    Logging.with_logger(f, logger)
    return [string(entry.message) for entry in logger.logs]
end

function walker_configuration(spacecraft; duration=60.0,
                              settings=SimulationSettings(results=false, generate_plots=false))
    planet = make_no_gram_planet(:earth)
    return SimulationConfiguration(
        simulation_settings=settings,
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=duration, orientation_sim=false, data_rate=10.0),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]),
        initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-8, dt_max_orbit=10.0),
        solver_config=SM.SolverConfig(solver_mode=:tsit5))
end

@testset "Constellation builder" begin
    @testset "the mean anomaly in initial conditions" begin
        for e in (0.0, 0.01, 0.3, 0.7, 0.95), M in (0.0, 10.0, 90.0, 179.0, 180.0, 270.0, 359.0, 725.0)
            ν = rad2deg(InitialCondition(a=7e6, e=e, M=M).ν)
            @test angle_gap(mean_anomaly_deg(ν, e), M) < 1e-9
        end
        @test InitialCondition(a=7e6, M=33.0) == InitialCondition(a=7e6, ν=33.0)   # equal on a circle
        @test InitialCondition(a=7e6, i=53, Ω=10, ω=5, M=0) == InitialCondition(a=7e6, i=53.0, Ω=10.0, ω=5.0, ν=0.0)
        # The defaults stay: periapsis for a and e, apoapsis for ra and rp.
        @test InitialCondition(a=7e6, e=0.1).ν == 0.0
        @test InitialCondition(ra=8e6, rp=7e6).ν == deg2rad(180.0)
        # With apsides, M is converted on the orbit's own eccentricity.
        apsides = InitialCondition(ra=8e6, rp=7e6, M=40)
        @test angle_gap(mean_anomaly_deg(rad2deg(apsides.ν), apsides.e), 40) < 1e-9
        planet = make_no_gram_planet(:earth)
        oblate = InitialCondition(planet; ra=900e3, hp=400e3, i=30, M=40, L_PI=Matrix(1.0I, 3, 3))
        @test angle_gap(mean_anomaly_deg(rad2deg(oblate.ν), oblate.e), 40) < 1e-9
        @test InitialCondition(planet; ra=900e3, hp=400e3, i=30, L_PI=Matrix(1.0I, 3, 3)).ν == deg2rad(180.0)
        @test_throws ArgumentError InitialCondition(a=7e6, ν=10, M=10)
        @test_throws ArgumentError InitialCondition(a=7e6, e=1.2, M=10)
    end

    @testset "Kepler's equation up to a parabolic orbit" begin
        OG = SM.OrbitalGeometry
        # The apsides are exact at any eccentricity. On a near-parabolic orbit the
        # Newton loop alone ran out of steps at M = 0 and returned about 0.158°.
        @test InitialCondition(a=1e23, e=prevfloat(1.0), M=0.0).ν == 0.0
        @test InitialCondition(a=1e23, e=prevfloat(1.0), M=180.0).ν == deg2rad(180.0)
        for e in (0.0, 0.3, 0.95, 0.999, 1 - 1e-12, prevfloat(1.0))
            @test OG._true_anomaly_from_mean(0.0, e) == 0.0 && OG._true_anomaly_from_mean(2pi, e) == 0.0
            @test OG._true_anomaly_from_mean(Float64(pi), e) == Float64(pi)
        end
        @test_throws ArgumentError OG._true_anomaly_from_mean(Inf, 0.5)
        @test_throws ArgumentError OG._true_anomaly_from_mean(NaN, 0.5)
        # Near-parabolic orbits at small and moderate mean anomalies, and their mirror
        # images below 2π, against an inversion in BigFloat: within a few ulps.
        for e in (1 - 1e-12, prevfloat(1.0)), M in (1e-30, 1e-20, 1e-12, 1e-6, 1e-3, 0.1, 1.0, 3.0)
            ν = OG._true_anomaly_from_mean(M, e)
            reference = bigfloat_true_anomaly(M, e)
            @test abs(ν - reference) <= 8 * eps(reference)
            M >= 1e-12 || continue
            mirrored = 2pi - M
            ν = OG._true_anomaly_from_mean(mirrored, e)
            @test abs(ν - (2 * big(pi) - bigfloat_true_anomaly(2pi - mirrored, e))) <= 2 * eps(2pi)
        end
        # The bracketed solve on its own, at inputs where plain Newton's iterates cycle
        # a few ulps around the root, or where a step lands on 0: the bracket settles them.
        for (M, e) in ((0.37491748451944507, 0.23594472361809046), (0.07641192965098051, 0.6475929648241207),
                       (1e-300, 0.0))
            ν = OG._true_anomaly_from_eccentric(OG._eccentric_anomaly_safeguarded(M, e), e)
            reference = bigfloat_true_anomaly(M, e)
            @test abs(ν - reference) <= 8 * eps(reference)
        end
        # Where the Newton loop converged and 1 - e cos E stays >= 0.05, which holds
        # for every e <= 0.95, the result is the loop's, bit for bit. That includes the
        # cases above, except the apsis M = 0, which the loop missed by a few 1e-37 rad
        # at e >= 0.8 (and returned as 2π at e = 0.999).
        changed = Tuple{Float64, Float64}[]
        for e in (0.01, 0.3, 0.7, 0.8, 0.9, 0.95, 0.99, 1 - 1e-12),
                M in vcat(deg2rad.([10.0, 90.0, 179.0, 180.0, 270.0, 359.0, 725.0]), range(-1.0, 7.0; length=321))
            ν, converged, slope = newton_true_anomaly(M, e)
            (converged && slope >= 0.05 && mod(M, 2pi) != 0.0) || continue
            OG._true_anomaly_from_mean(M, e) === ν || push!(changed, (M, e))
        end
        @test isempty(changed)
        @test newton_true_anomaly(0.0, 0.95)[1] > 0.0 && OG._true_anomaly_from_mean(0.0, 0.95) == 0.0
        @test newton_true_anomaly(0.0, 0.999)[1] == 2pi
    end

    @testset "orbital planes" begin
        plane = OrbitalPlane(a=7e6, i=53, Ω=40, M=[0, 90])
        @test (plane.a, plane.e, plane.i, plane.Ω, plane.ω, plane.M) == (7e6, 0.0, 53.0, 40.0, 0.0, [0.0, 90.0])
        for bad in ((a=-1.0,), (e=1.0,), (M=Float64[],), (M=[NaN],), (i=Inf,))
            @test_throws ArgumentError OrbitalPlane(; merge((a=7e6, i=53, Ω=40, M=[0]), bad)...)
        end
    end

    @testset "Walker patterns" begin
        delta = walker_planes(total=24, planes=3, phasing=1, a=7e6, i=53)
        @test length(delta) == 3 && all(plane -> length(plane.M) == 8 && plane.i == 53.0, delta)
        @test [plane.Ω for plane in delta] == [0.0, 120.0, 240.0]
        @test delta[1].M ≈ 0:45:315
        @test delta[2].M ≈ 15:45:330                       # phasing 1 * 360 / 24 per plane
        @test delta[3].M ≈ 30:45:345
        star = walker_planes(total=6, planes=3, phasing=2, a=7e6, i=86.4, pattern=:star, raan0=10, M0=5)
        @test [plane.Ω for plane in star] == [10.0, 70.0, 130.0]
        @test star[2].M ≈ [125.0, 305.0] && star[3].M ≈ [245.0, 65.0]    # reduced to [0, 360)
        planet = make_no_gram_planet(:earth)
        @test walker_planes(total=2, planes=1, phasing=0, altitude=550e3, planet=planet, i=0)[1].a == planet.Rp_e + 550e3
        for bad in ((total=24, planes=5), (phasing=3,), (phasing=-1,), (pattern=:rosette,), (altitude=550e3,),
                    (a=nothing,), (planes=0,))
            @test_throws ArgumentError walker_planes(; merge((total=24, planes=3, phasing=1, a=7e6, i=53), bad)...)
        end
    end

    @testset "the builder" begin
        calls = Tuple{Int, Int, Int}[]
        factory = (plane, slot, index) -> (push!(calls, (plane, slot, index)); bus_factory(plane, slot, index))
        spacecraft = walker_constellation(factory; total=6, planes=3, phasing=1, a=7e6, i=53, first_id=10)
        @test calls == [(1, 1, 1), (1, 2, 2), (2, 1, 3), (2, 2, 4), (3, 1, 5), (3, 2, 6)]
        @test [sc.id for sc in spacecraft] == 10:15
        @test [sc.constellation_slot for sc in spacecraft] ==
            [ConstellationSlot(p, s) for p in 1:3 for s in 1:2]
        @test [sc.root.m for sc in spacecraft] == 101.0:106.0          # the factory's own spacecraft
        fourth = spacecraft[4].initial_condition                          # plane 2, slot 2
        @test (fourth.a, fourth.e) == (7e6, 0.0)
        @test [rad2deg(fourth.i), rad2deg(fourth.Ω), rad2deg(fourth.ν)] ≈ [53.0, 120.0, 240.0]
        # The factory's attitude and rate are kept; only the orbit is replaced.
        turned = (plane, slot, index) -> begin
            sc = bus_factory(plane, slot, index)
            sc.initial_condition = InitialCondition(a=1.0, q=(0.0, 0.0, sind(15), cosd(15)), ang_vel=(0.0, 0.0, 0.01 * index))
            sc
        end
        kept = build_constellation(turned, [OrbitalPlane(a=7e6, i=10, Ω=0, M=[0, 180])])
        @test kept[2].initial_condition.q == SVector(0.0, 0.0, sind(15), cosd(15))
        @test kept[2].initial_condition.ang_vel == SVector(0.0, 0.0, 0.02) && kept[2].initial_condition.a == 7e6
        cartesian = (plane, slot, index) -> begin
            sc = bus_factory(plane, slot, index)
            sc.initial_condition = CartesianInitialCondition(zeros(3), zeros(3); q=(0.0, 1.0, 0.0, 0.0))
            sc
        end
        @test build_constellation(cartesian, [OrbitalPlane(a=7e6, i=10, Ω=0, M=[0])])[1].initial_condition.q ==
            SVector(0.0, 1.0, 0.0, 0.0)
        # A plain spacecraft has no slot, and keeps it through the GNC flattening.
        @test SpacecraftModel().constellation_slot === nothing
        @test SE._without_gnc(spacecraft[3]).constellation_slot == ConstellationSlot(2, 1)
        @test_throws ArgumentError ConstellationSlot(0, 1)
    end

    @testset "the factory must build fresh objects" begin
        planes = [OrbitalPlane(a=7e6, i=10, Ω=0, M=[0, 120, 240])]
        message(f) = try
            build_constellation(f, planes)
            ""
        catch err
            err isa ArgumentError ? sprint(showerror, err) : rethrow()
        end
        one = bus_factory(1, 1, 1)
        @test occursin("plane 1, slot 2 shares a SpacecraftModel with the one for plane 1, slot 1", message((p, s, i) -> one))
        shared_bus = Link(root=true)
        @test occursin("shares a Link", message((p, s, i) -> SpacecraftModel(root=shared_bus, links=[shared_bus])))
        thruster = Thruster(direction=MVector(1.0, 0.0, 0.0))
        shared_device = (p, s, i) -> begin
            sc = bus_factory(p, s, i)
            push!(sc.root.thrusters, thruster)
            sc
        end
        @test occursin("shares a Thruster", message(shared_device))
        controller = SM.MagneticMomentumManagerModel(commanded_torque=(args...) -> zeros(3), b_field_ii=(args...) -> zeros(3))
        shared_gnc = (p, s, i) -> begin
            bus = Link(root=true)
            SpacecraftModel(root=bus, links=[bus], control=SM.ControlModel(control_effectors=(controller,), control_rates=[1.0]))
        end
        @test occursin("shares a MagneticMomentumManagerModel", message(shared_gnc))
        @test occursin("returned a Int64", message((p, s, i) -> i))
        # Immutable values hold no state and may be shared.
        initial = InitialCondition(a=7e6)
        @test length(build_constellation((p, s, i) -> (b = Link(root=true);
            SpacecraftModel(root=b, links=[b], initial_condition=initial)), planes)) == 3
    end

    @testset "a Walker constellation runs from its elements" begin
        spacecraft = walker_constellation(bus_factory; total=4, planes=2, phasing=1, altitude=700e3,
            planet=make_no_gram_planet(:earth), i=53)
        args = walker_configuration(spacecraft)
        u0 = SE.build_initial_conditions(args)
        for (index, sc) in enumerate(spacecraft)
            inclination, raan, latitude = plane_angles(SVector{3}(u0.sc[index].pos), SVector{3}(u0.sc[index].vel))
            plane, slot = sc.constellation_slot.plane, sc.constellation_slot.slot
            @test inclination ≈ 53.0 atol=1e-9
            @test angle_gap(raan, (plane - 1) * 180.0) < 1e-9
            @test angle_gap(latitude, (slot - 1) * 180.0 + (plane - 1) * 90.0) < 1e-9
        end
        solution = run_simulation(args; return_solution=true)
        ran = solution.prob.p.args.dynamics_model.spacecraft
        @test [sc.id for sc in ran] == 1:4
        @test [sc.constellation_slot for sc in ran] == [sc.constellation_slot for sc in spacecraft]
        @test solution.t[end] ≈ 60.0
    end

    @testset "repeated ids warn" begin
        spacecraft = [bus_factory(1, 1, index) for index in 1:5]
        foreach((sc, id) -> sc.id = id, spacecraft, [1, 1, 1, 4, 4])
        args = walker_configuration(spacecraft)
        messages = deprecation_messages(() -> SE._warn_repeated_spacecraft_ids(args))
        @test length(messages) == 1
        @test occursin("Spacecraft ids repeat (id 1 at positions 1, 2, 3; id 4 at positions 4, 5)", messages[1])
        unique_args = walker_configuration(walker_constellation(bus_factory; total=2, planes=1, phasing=0, a=7e6, i=0))
        @test isempty(deprecation_messages(() -> SE._warn_repeated_spacecraft_ids(unique_args)))
    end

    @testset "ensemble folders sort in member order" begin
        tag = SpaceAGORA.SimulationCampaigns._ensemble_member_tag
        @test tag(3, 9, 7) == "sat_3_id_7"                     # up to nine members: unchanged
        @test tag(3, 12, 7) == "sat_03_id_7"
        @test tag(12, 12, 40) == "sat_12_id_40"
        @test tag(5, 2048, 5) == "sat_0005_id_5"
        @test issorted([tag(index, 12, index) for index in 1:12])
    end

    @testset "a resumed ensemble finds checkpoints in unpadded folders" begin
        mktempdir() do directory
            spacecraft = walker_constellation(bus_factory; total=12, planes=3, phasing=0, a=7e6, i=53)
            settings(resume) = SimulationSettings(results=false, generate_plots=false, checkpoint_enabled=true,
                checkpoint_directory=directory, resume_from_checkpoint=resume)
            member(args, index) = SpaceAGORA.SimulationCampaigns._ensemble_member_configuration(
                args, args.dynamics_model.spacecraft[index], index, 12).simulation_settings.checkpoint_directory
            touch_checkpoint(folder) = (mkpath(joinpath(directory, folder));
                touch(SE._checkpoint_paths((simulation_settings=(checkpoint_directory=joinpath(directory, folder),),)).data))
            resumed = walker_configuration(spacecraft; settings=settings(true))
            @test member(resumed, 3) == joinpath(directory, "sat_03_id_3")       # no checkpoint anywhere
            touch_checkpoint("sat_3_id_3")                                       # written before the padding
            @test member(resumed, 3) == joinpath(directory, "sat_3_id_3")
            @test member(resumed, 4) == joinpath(directory, "sat_04_id_4")
            fresh = walker_configuration(spacecraft; settings=settings(false))
            @test member(fresh, 3) == joinpath(directory, "sat_03_id_3")         # a new run writes the padded name
            touch_checkpoint("sat_03_id_3")
            @test member(resumed, 3) == joinpath(directory, "sat_03_id_3")       # the padded folder wins
        end
    end
end

end # module
