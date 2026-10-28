# Setup API reachable from the root module alone: no SimulationModel import.
using SpaceAGORA, Test

@testset "public setup API" begin
    @test SpaceAGORA.TelemetryVerification.make_example_config === SpaceAGORA.make_example_config
    @test SpaceAGORA.TelemetryVerification.make_three_body_spacecraft === SpaceAGORA.make_three_body_spacecraft

    @test ephemerides_requires_spice(SpiceEphemeridesModel())
    planet = make_no_gram_planet(:earth)
    spacecraft = make_three_body_spacecraft(
        bus_dims=(2.0, 2.0, 2.0),
        panel_dims=(0.01, 2.0, 1.0),
        bus_mass=500.0,
        panel_mass_each=10.0,
        panel_offset_y=2.0,
        ic=InitialCondition(ra=planet.Rp_e + 500e3, rp=planet.Rp_e + 500e3, i=45.0, ω=0.0, Ω=0.0, ν=0.0),
    )
    configure(effectors, ephemerides) = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=300.0,
        initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=effectors,
        density_model=NoAtmosphereModel(),
        ephemerides_model=ephemerides,
        verbose=false,
        results=false,
    )
    config = configure((InverseSquaredGravityModel(),), SimpleEphemeridesModel())
    # Construct the research setup using only root names; no kernels are needed
    # until this SPICE-backed configuration is run.
    research_effectors = (InverseSquaredGravityModel(), SolarRadiationPressureModel(1.2, 4.0))
    research = configure(research_effectors, SpiceEphemeridesModel())
    @test ephemerides_requires_spice(research.environment_model.ephemerides_model)
    # Active sunlight must not silently fall back to the simple backend.
    @test_throws ArgumentError run_simulation(configure(research_effectors, SimpleEphemeridesModel()))
    @test config isa SimulationConfiguration
    @test config.mission_configuration.mission_type == MissionTime
    @test (run_simulation(config); true)
end

# A package built on the root names alone sets an initial attitude and rate with plain
# vectors or tuples; the constructors convert them, and refuse a wrong length.
@testset "initial attitude and rate from any vector" begin
    planet = make_no_gram_planet(:earth)
    q, rate = (0.0, 0.6, 0.0, 0.8), (0.01, -0.02, 0.03)
    elements(; kw...) = InitialCondition(; ra=planet.Rp_e + 700e3, rp=planet.Rp_e + 600e3, i=45.0, ν=10.0, kw...)
    from_tuples, from_vectors = elements(q=q, ang_vel=rate), elements(q=collect(q), ang_vel=collect(rate))
    @test collect(from_tuples.q) == collect(q) && collect(from_tuples.ang_vel) == collect(rate)
    @test eltype(from_tuples.q) === Float64 && eltype(from_tuples.ang_vel) === Float64
    @test isequal(from_tuples.q, from_vectors.q) && isequal(from_tuples.ang_vel, from_vectors.ang_vel)
    @test isequal(elements(q=[0, 0, 0, 1], ang_vel=[0, 0, 0]).q, elements().q)
    @test isequal(elements(q=[0, 0, 0, 1], ang_vel=[0, 0, 0]).ang_vel, elements().ang_vel)
    positional = InitialCondition(7.0e6, 0.01, 45.0, 0.0, 0.0, 10.0, collect(q), rate)
    @test isequal(positional.q, from_tuples.q) && isequal(positional.ang_vel, from_tuples.ang_vel)
    oblate = InitialCondition(planet; ra=700e3, hp=600e3, q=collect(q), ang_vel=rate)
    @test isequal(oblate.q, from_tuples.q) && isequal(oblate.ang_vel, from_tuples.ang_vel)
    cartesian = CartesianInitialCondition([planet.Rp_e + 600e3, 0.0, 0.0], [0.0, 7600.0, 0.0]; q=q, ang_vel=collect(rate))
    @test isequal(cartesian.q, from_tuples.q) && isequal(cartesian.ang_vel, from_tuples.ang_vel)
    @test collect(CartesianInitialCondition([7.0e6, 0.0, 0.0], [0.0, 7600.0, 0.0]).q) == [0.0, 0.0, 0.0, 1.0]
    @test_throws ArgumentError elements(q=[0.0, 0.0, 1.0])
    @test_throws ArgumentError elements(ang_vel=(0.0, 0.0))
    @test_throws ArgumentError InitialCondition(planet; ra=700e3, hp=600e3, ang_vel=zeros(4))
    @test_throws ArgumentError CartesianInitialCondition([7.0e6, 0.0, 0.0], [0.0, 7600.0, 0.0]; q=zeros(5))
end
