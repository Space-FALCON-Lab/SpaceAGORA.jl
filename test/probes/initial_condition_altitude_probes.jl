using Test
using LinearAlgebra
using StaticArrays
using SPICE: recgeo
using SpaceAGORA
using SpaceAGORA.SimulationModel: InitialCondition, CartesianInitialCondition,
    Earth, Mars, Venus, InitialTime, SimpleEphemeridesModel,
    ephemerides_time_seconds, planet_frame_lpi

const ALT_ENGINE = SpaceAGORA.SimulationEngine
const ALT_GEODESY = SpaceAGORA.SimulationModel.Geodesy
const ALT_GUIDANCE = SpaceAGORA.SimulationModel.GuidanceHooks
const ALT_IDENTITY = SMatrix{3, 3, Float64}(I)

function _altitude_probe_apsis_position(ic, planet, anomaly)
    # Stored angles are radians; the public keyword constructor accepts degrees.
    elements = SVector{7, Float64}(ic.a, ic.e, ic.i, ic.Ω, ic.ω, anomaly, 0.0)
    position, _ = ALT_ENGINE.orbitalelemtorv(elements, planet)
    return SVector{3, Float64}(position)
end

function _altitude_probe_check_apses(ic, planet, l_pi, apo_alt, peri_alt)
    flattening = (planet.Rp_e - planet.Rp_p) / planet.Rp_e
    for (label, anomaly, expected) in (("apoapsis", Float64(pi), apo_alt),
                                      ("periapsis", 0.0, peri_alt))
        @testset "$label" begin
            position_pp = l_pi * _altitude_probe_apsis_position(ic, planet, anomaly)
            # CSPICE's geometric conversion is independent of our Bowring
            # implementation and uses no kernels. Lengths remain in metres.
            _, _, height = recgeo(collect(position_pp), planet.Rp_e, flattening)
            @test height ≈ expected atol=1e-3 rtol=0
            @test ALT_ENGINE.rtolatlong(position_pp, planet)[1] ≈ expected atol=1e-3 rtol=0
        end
    end
end

_altitude_probe_elements(ic) = (ic.a, ic.e, ic.i, ic.ω, ic.Ω, ic.ν, ic.q, ic.ang_vel)

@testset "Planet-aware initial-condition altitudes" begin
    # Plain planet constructors need no external files and have independent
    # frame storage. Do not modify the cached SPICE-backed planet instances.
    for planet in (Earth(), Mars())
        for (label, inclination, argument) in (("north periapsis", 45.0, 90.0),
                                               ("south periapsis", 60.0, 270.0),
                                               ("both poles", 90.0, 90.0),
                                               ("equator", 0.0, 0.0))
            @testset "$(planet.name), $label" begin
                ic = InitialCondition(planet; ra=500e3, hp=150e3,
                    i=inclination, ω=argument, Ω=31.0, L_PI=ALT_IDENTITY)
                _altitude_probe_check_apses(ic, planet, ALT_IDENTITY, 500e3, 150e3)
                @test ic.ν == deg2rad(180.0)
                @test ic.i == deg2rad(inclination)
                @test ic.ω == deg2rad(argument)
                @test ic.Ω == deg2rad(31.0)
            end
        end
        @testset "$(planet.name), surface periapsis" begin
            ic = InitialCondition(planet; ra=500e3, hp=0.0,
                i=60.0, ω=270.0, L_PI=ALT_IDENTITY)
            _altitude_probe_check_apses(ic, planet, ALT_IDENTITY, 500e3, 0.0)
        end
    end

    @testset "Spherical body" begin
        planet = Venus()
        @test planet.Rp_e == planet.Rp_p
        ic = InitialCondition(planet; ra=500e3, hp=150e3,
            i=61.0, ω=240.0, Ω=17.0, L_PI=ALT_IDENTITY)
        @test ic.a * (1 + ic.e) ≈ planet.Rp_e + 500e3 atol=1e-6 rtol=0
        @test ic.a * (1 - ic.e) ≈ planet.Rp_e + 150e3 atol=1e-6 rtol=0
        _altitude_probe_check_apses(ic, planet, ALT_IDENTITY, 500e3, 150e3)
    end

    @testset "Frame selection and precedence" begin
        c, s = cos(0.63), sin(0.63)
        tilted = @SMatrix [1.0 0.0 0.0; 0.0 c -s; 0.0 s c]
        planet = Earth(L_PI=MMatrix{3, 3, Float64}(tilted))
        inputs = (ra=500e3, hp=150e3, i=52.0, ω=71.0, Ω=23.0)
        explicit = InitialCondition(planet; inputs..., L_PI=tilted)
        stored = InitialCondition(planet; inputs...)
        @test isequal(_altitude_probe_elements(stored), _altitude_probe_elements(explicit))
        _altitude_probe_check_apses(stored, planet, tilted, 500e3, 150e3)

        # An uninitialized all-zero planet frame falls back to identity.
        fresh = Earth()
        fallback = InitialCondition(fresh; inputs...)
        identity = InitialCondition(fresh; inputs..., L_PI=ALT_IDENTITY)
        @test isequal(_altitude_probe_elements(fallback), _altitude_probe_elements(identity))
        _altitude_probe_check_apses(fallback, fresh, ALT_IDENTITY, 500e3, 150e3)

        epoch = InitialTime(year=2024, month=3, day=2, hour=4, minute=5, second=6.0)
        ephemerides = SimpleEphemeridesModel()
        time_seconds = ephemerides_time_seconds(epoch, ephemerides)
        epoch_frame = planet_frame_lpi(planet, time_seconds, ephemerides)
        from_epoch = InitialCondition(planet; inputs...,
            initial_time=epoch, ephemerides_model=ephemerides)
        explicit_epoch = InitialCondition(planet; inputs..., L_PI=epoch_frame)
        @test isequal(_altitude_probe_elements(from_epoch), _altitude_probe_elements(explicit_epoch))
        _altitude_probe_check_apses(from_epoch, planet, epoch_frame, 500e3, 150e3)
        override = InitialCondition(planet; inputs..., L_PI=tilted,
            initial_time=epoch, ephemerides_model=ephemerides)
        @test isequal(_altitude_probe_elements(override), _altitude_probe_elements(explicit))
        @test abs(explicit.a - explicit_epoch.a) > 1.0
    end

    @testset "Constructor compatibility" begin
        planet = Earth()
        q = SVector(0.0, 0.0, 0.6, 0.8)
        angular_velocity = SVector(0.01, -0.02, 0.03)
        ic = InitialCondition(planet; ra=500_000, hp=150_000,
            i=45, ω=90, Ω=30, ν=37, q=q, ang_vel=angular_velocity,
            L_PI=ALT_IDENTITY)
        @test ic.ν == deg2rad(37.0)
        @test ic.q === q
        @test ic.ang_vel === angular_velocity
        _altitude_probe_check_apses(ic, planet, ALT_IDENTITY, 500e3, 150e3)
        @test_throws ArgumentError InitialCondition(planet; ra=500e3)
        @test_throws ArgumentError InitialCondition(planet; hp=150e3)
        @test_throws ArgumentError InitialCondition(planet; ra=500e3, hp=-1.0)
        @test_throws ArgumentError InitialCondition(planet; ra=-1.0, hp=150e3)
        @test_throws ArgumentError InitialCondition(planet; ra=150e3, hp=500e3)
        @test_throws ArgumentError InitialCondition(planet; ra=150e3, hp=150e3)

        # Without a planet argument, ra/rp remain radii rather than altitudes.
        radii = InitialCondition(ra=7.0e6, rp=6.8e6, i=45.0, ω=90.0, Ω=30.0)
        @test radii.a == 6.9e6
        @test radii.e ≈ 1 / 69 rtol=1e-15
        @test radii.ν == deg2rad(180.0)
        @test InitialCondition(a=6.9e6, e=0.01).ν == 0.0
        position = SVector(6.9e6, 1.0, -2.0)
        velocity = SVector(0.0, 7.6e3, 3.0)
        cartesian = CartesianInitialCondition(position, velocity; q=q, ang_vel=angular_velocity)
        @test cartesian.pos === position
        @test cartesian.vel === velocity
        @test cartesian.q === q
        @test cartesian.ang_vel === angular_velocity
    end

    @testset "Shared altitude and guidance semantics" begin
        for planet in (Earth(), Mars(), Venus()), latitude in (-90.0, -60.0, 0.0, 45.0, 90.0)
            direction = SVector(cosd(latitude), 0.0, sind(latitude))
            direction /= norm(direction)
            surface = ALT_GEODESY.ellipsoid_surface_radius(direction, planet)
            @test surface == ALT_GUIDANCE._oblate_surface_radius(direction, planet)
            for expected in (0.0, 150e3, 500e3)
                radius = ALT_GEODESY.radius_for_geodetic_altitude(expected, direction, planet)
                guidance_radius = ALT_GUIDANCE._radius_for_oblate_altitude(expected, direction, planet)
                @test reinterpret(UInt64, radius) == reinterpret(UInt64, guidance_radius)
                altitude = ALT_GEODESY.geodetic_altitude(radius, direction, planet)
                guidance_altitude = ALT_GUIDANCE._oblate_altitude_from_radius(radius, direction, planet)
                @test reinterpret(UInt64, altitude) == reinterpret(UInt64, guidance_altitude)
                _, _, independent_altitude = recgeo(collect(radius * direction),
                    planet.Rp_e, (planet.Rp_e - planet.Rp_p) / planet.Rp_e)
                @test independent_altitude ≈ expected atol=1e-3 rtol=0
            end
        end
        # Guidance historically returns NaN for an invalid target; the public
        # initial-condition constructor instead raises ArgumentError above.
        @test isnan(ALT_GUIDANCE._radius_for_oblate_altitude(-1.0, SVector(1.0, 0.0, 0.0), Earth()))
    end
end
