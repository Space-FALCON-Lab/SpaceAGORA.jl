using Test
using LinearAlgebra
using SpaceAGORA
using CSV
using DataFrames

const SM_SFR = SpaceAGORA.SimulationModel
const RVTOOE_SFR = SpaceAGORA.SimulationModel.SimulationCallbacks.rvtoorbitalelement

_sfr_field_names(fields) = Symbol[field.name for field in fields]

function _sfr_config()
    planet = make_no_gram_planet(:earth)
    initial_time = SM_SFR.InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0)
    ephemerides_model = SM_SFR.SimpleEphemeridesModel()
    ic = SM_SFR.CartesianInitialCondition(
        planet;
        lat=0.0,
        lon=0.0,
        alt=200e3,
        speed=sqrt(planet.μ / (planet.Rp_e + 200e3)),
        inclination=89.876,
        initial_time=initial_time,
        ephemerides_model=ephemerides_model
    )
    spacecraft = SpaceAGORA.TelemetryVerification.make_three_body_spacecraft(
        bus_dims=(1.0, 1.0, 1.0),
        panel_dims=(0.01, 1.0, 1.0),
        bus_mass=100.0,
        panel_mass_each=5.0,
        panel_offset_y=1.0,
        ic=ic,
        id=1
    )
    args = SpaceAGORA.TelemetryVerification.make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=60.0,
        initial_time=initial_time,
        density_model=SM_SFR.NoAtmosphereModel(),
        ephemerides_model=ephemerides_model,
        results=false,
        verbose=false
    )
    return planet, initial_time, ephemerides_model, args
end

@testset "SaveFieldRegistry" begin
    _, _, _, args = _sfr_config()
    base = _sfr_field_names(default_save_fields(args))

    # `extra` appends in the order given, after the defaults
    @test _sfr_field_names(default_save_fields(args; extra=(:orbital_elements, :gravity_accel))) ==
        vcat(base, [:orbital_elements, :gravity_accel])

    # a name already in the list is skipped rather than duplicated, so the
    # resolved set never trips run_simulation's uniqueness check
    @test _sfr_field_names(default_save_fields(args; extra=(:gravity_accel, :gravity_accel))) ==
        vcat(base, [:gravity_accel])

    # a SaveField of the caller's own mixes in with the built-in names
    custom = SaveField(:custom_field, (u, t, integrator) -> 1.0)
    @test _sfr_field_names(default_save_fields(args; extra=(:gravity_accel, custom))) ==
        vcat(base, [:gravity_accel, :custom_field])

    # every advertised name builds, and its name is the one asked for
    @test !isempty(available_save_fields())
    for name in available_save_fields()
        field = save_field(name, args)
        @test field isa SaveField
        @test field.name === name
    end

    # an unknown name fails with the list of what is available
    err = try
        save_field(:not_a_save_field, args)
    catch caught
        caught
    end
    @test err isa ArgumentError
    @test occursin("orbital_elements", sprint(showerror, err))
end

@testset "GeodeticInitialCondition" begin
    planet, initial_time, ephemerides_model, _ = _sfr_config()

    geodetic_ic(; kwargs...) = SM_SFR.CartesianInitialCondition(
        planet;
        initial_time=initial_time,
        ephemerides_model=ephemerides_model,
        kwargs...
    )

    # on the equator the geodetic altitude is the radius less the equatorial
    # radius, so the requested point is hit exactly
    equatorial = geodetic_ic(lat=0.0, lon=0.0, alt=200e3, speed=7.8e3, inclination=89.876)
    @test isapprox(norm(equatorial.pos), planet.Rp_e + 200e3; rtol=1e-12)

    # the requested inclination comes out exact away from the equator too, both
    # crossing directions, with and without a flight path angle
    for descending in (false, true), flight_path_angle in (0.0, -5.0)
        ic = geodetic_ic(
            lat=30.0,
            lon=10.0,
            alt=200e3,
            speed=7.6e3,
            inclination=60.0,
            descending=descending,
            flight_path_angle=flight_path_angle
        )
        @test isapprox(rad2deg(RVTOOE_SFR(ic.pos, ic.vel, planet)[3]), 60.0; atol=1e-9)
        # northbound for the ascending crossing, southbound for the descending one
        @test (ic.vel[3] > 0.0) == !descending
    end

    # flight path angle is the orbital one: tan(γ) = r·v / |r × v|
    climbing = geodetic_ic(lat=15.0, lon=-40.0, alt=120e3, speed=7.0e3, azimuth=75.0, flight_path_angle=3.0)
    @test isapprox(
        rad2deg(atan(dot(climbing.pos, climbing.vel), norm(cross(climbing.pos, climbing.vel)))),
        3.0;
        atol=1e-9
    )
    @test isapprox(norm(climbing.vel), 7.0e3; rtol=1e-12)

    # an inclination the point never reaches is refused, as is a missing or a
    # doubled heading input
    @test_throws ArgumentError geodetic_ic(lat=60.0, lon=0.0, alt=200e3, speed=7.6e3, inclination=30.0)
    @test_throws ArgumentError geodetic_ic(lat=0.0, lon=0.0, alt=200e3, speed=7.6e3)
    @test_throws ArgumentError geodetic_ic(lat=0.0, lon=0.0, alt=200e3, speed=7.6e3, azimuth=90.0, inclination=60.0)
    # over a pole there is no local east to measure an azimuth from
    @test_throws ArgumentError geodetic_ic(lat=90.0, lon=0.0, alt=200e3, speed=7.6e3, azimuth=0.0)
end

# The acceleration fields are only worth anything if they actually agree with
# the forces the run integrated, so this one flies a short arc and reads the
# columns back rather than checking that the builders return a SaveField.
@testset "SavedAccelerationColumns" begin
    planet = make_no_gram_planet(:earth)
    initial_time = SM_SFR.InitialTime(year=2014, month=5, day=27, hour=5, minute=0, second=0.0)
    ephemerides_model = SM_SFR.SimpleEphemeridesModel()
    ic = SM_SFR.CartesianInitialCondition(
        planet;
        lat=0.0,
        lon=0.0,
        alt=500e3,
        speed=sqrt(planet.μ / (planet.Rp_e + 500e3)),
        inclination=45.0,
        initial_time=initial_time,
        ephemerides_model=ephemerides_model
    )
    spacecraft = SpaceAGORA.TelemetryVerification.make_three_body_spacecraft(
        bus_dims=(1.0, 1.0, 1.0),
        panel_dims=(0.01, 1.0, 1.0),
        bus_mass=100.0,
        panel_mass_each=5.0,
        panel_offset_y=1.0,
        ic=ic,
        id=1
    )
    results_dir = mktempdir()
    args = SpaceAGORA.TelemetryVerification.make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=120.0,
        initial_time=initial_time,
        dynamic_effectors=(SM_SFR.InverseSquaredJ2GravityModel(),),
        density_model=SM_SFR.NoAtmosphereModel(),
        ephemerides_model=ephemerides_model,
        verbose=false,
        results_directory=results_dir
    )
    run_simulation(args; save_fields=default_save_fields(
        args;
        extra=(:orbital_elements, :gravity_accel, :aero_accel, :total_accel)
    ))

    df = CSV.read(joinpath(results_dir, "simulation_results.csv"), DataFrame)
    @test nrow(df) > 1
    for k in 1:3
        @test Symbol("sc1_gravity_accel_$(k)") in propertynames(df)
        @test Symbol("sc1_aero_accel_$(k)") in propertynames(df)
        @test Symbol("sc1_total_accel_$(k)") in propertynames(df)
    end

    # gravity is the run's only dynamic effector, so the total is exactly it
    for k in 1:3
        @test df[!, Symbol("sc1_total_accel_$(k)")] == df[!, Symbol("sc1_gravity_accel_$(k)")]
        # and with no atmosphere the aerodynamic share is exactly zero
        @test all(iszero, df[!, Symbol("sc1_aero_accel_$(k)")])
    end

    # the gravity column agrees with mu/r^2 to the J2 correction, not by luck
    g = [norm([df[row, Symbol("sc1_gravity_accel_$(k)")] for k in 1:3]) for row in 1:nrow(df)]
    r = [norm([df[row, Symbol("sc1_pos_$(k)")] for k in 1:3]) for row in 1:nrow(df)]
    point_mass = planet.μ ./ r .^ 2
    @test all(@. abs(g - point_mass) / point_mass < 5e-3)
    @test all(>(0.0), g)

    # the first row's elements are the ones the initial condition was built with
    @test isapprox(df[1, :sc1_orbital_elements_3], 45.0; atol=1e-9)
    @test isapprox(df[1, :sc1_orbital_elements_1], planet.Rp_e + 500e3; rtol=1e-6)
end
