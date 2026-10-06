module ConfigDefaultsTests
# GuidanceModel/NavigationModel/ControlModel and the matching SimulationConfiguration
# fields default to the empty model, and omitting them changes nothing.
using Test
using SpaceAGORA
using SpaceAGORA.SimulationModel

function cd_config(; gnc=NamedTuple())
    planet = Earth()
    root = Link(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(ra=planet.Rp_e + 700.0e3, rp=planet.Rp_e + 600.0e3,
        i=53.0, ω=0.0, Ω=10.0, ν=0.0)
    sc = SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, 1)
    return SimulationConfiguration(;
        simulation_settings=SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=MissionConfiguration(MissionTime, true, 1, 300.0, false, 20, 2.0),
        environment_model=EnvironmentModel(planet=planet, EI=600.0,
            density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false, ephemerides_model=SimpleEphemeridesModel()),
        dynamics_model=DynamicsModel([sc], (InverseSquaredGravityModel(),)),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-10,
            dt_max_orbit=1.0),
        gnc...)
end

function cd_final_state(args)
    recorder = TrajectoryRecorder(args; capacity=1)
    result = withenv("SPACEAGORA_RHS_CALIBRATE" => "off", "SPACEAGORA_RHS_IDENTIFY" => "0") do
        run_simulation(args; isolate_state=false, return_solver_metadata=true,
            visualization=false, extra_callbacks=(get_trajectory_recorder_callback(recorder),))
    end
    @test result.retcode == "Success"
    return (copy(trajectory_positions(recorder)), copy(trajectory_velocities(recorder)))
end

@testset "ConfigDefaults" begin
    @testset "zero-argument GNC models are empty and unshared" begin
        for (m, eff, rates) in ((GuidanceModel(), :guidance_effectors, :guidance_rates),
                (NavigationModel(), :navigation_effectors, :navigation_rates),
                (ControlModel(), :control_effectors, :control_rates))
            @test getfield(m, eff) === ()
            @test isempty(getfield(m, rates))
            @test getfield(Base.typename(typeof(m)).wrapper(), rates) !== getfield(m, rates)
        end
        @test GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]) isa GuidanceModel
        @test_throws ArgumentError GuidanceModel(guidance_effectors=(), guidance_rates=[1.0])
    end

    @testset "omitting GNC models matches explicit empty models bit-for-bit" begin
        explicit = cd_config(gnc=(
            guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
            navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
            control_model=ControlModel(control_effectors=(), control_rates=Float64[])))
        omitted = cd_config()
        @test typeof(omitted) == typeof(explicit)
        @test cd_final_state(omitted) == cd_final_state(explicit)
    end
end
end
