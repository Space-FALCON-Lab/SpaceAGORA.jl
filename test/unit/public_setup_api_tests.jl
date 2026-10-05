# Setup API reachable from the root module alone: no SimulationModel import.
using SpaceAGORA, Test

@testset "public setup API" begin
    @test SpaceAGORA.TelemetryVerification.make_example_config === SpaceAGORA.make_example_config
    @test SpaceAGORA.TelemetryVerification.make_three_body_spacecraft === SpaceAGORA.make_three_body_spacecraft

    planet = make_no_gram_planet(:earth)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.0),
        panel_dims=(0.01, 2.0, 1.0),
        bus_mass=500.0,
        panel_mass_each=10.0,
        panel_offset_y=2.0,
        ic=InitialCondition(ra=planet.Rp_e + 500e3, rp=planet.Rp_e + 500e3, i=45.0, ω=0.0, Ω=0.0, ν=0.0),
    )
    config = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=300.0,
        initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(InverseSquaredGravityModel(),),
        density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(),
        verbose=false,
        results=false,
    )
    @test config isa SimulationConfiguration
    @test config.mission_configuration.mission_type == MissionTime
    @test (run_simulation(config); true)
end
