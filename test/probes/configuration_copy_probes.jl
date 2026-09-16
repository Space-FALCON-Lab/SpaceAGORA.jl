using Test
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const SC = SM.SimConfig
const TV = SpaceAGORA.TelemetryVerification

function copy_probe_config(; solver_config=SM.SolverConfig(solver_mode=:rodas5p, maxiters=12345))
    planet = SM.make_no_gram_planet(:earth)
    root = SM.Link(root=true, m=500.0, ref_area=12.0)
    ic = SM.InitialCondition(ra=planet.Rp_e + 550_000.0, rp=planet.Rp_e + 550_000.0,
                             i=53.0, ω=0.0, Ω=10.0, ν=15.0)
    spacecraft = SM.SpacecraftModel(SM.Joint[], [root], root, true, 500.0, 0.0,
                                    root.inertia, 0, 0, ic, 1)
    return TV.make_example_config(
        planet=planet, spacecraft=spacecraft, mission_time=600.0,
        initial_time=SM.InitialTime(year=2020, month=1, day=1),
        dynamic_effectors=(SM.InverseSquaredGravityModel(),),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        results=false, verbose=false, results_directory="custom-results",
        solver_config=solver_config,
    )
end

@testset "Configuration updates preserve fields and model ownership" begin
    args = copy_probe_config()
    unchanged = SC._with_configuration(args)
    for name in fieldnames(typeof(args))
        @test getfield(unchanged, name) === getfield(args, name)
    end
    @test SC._with_configuration(args; solver_config=nothing).solver_config === nothing
    @test_throws MethodError SC._with_configuration(args; solver_confg=nothing)

    # Both parameterized model types must be inferred anew by the constructor.
    env = args.environment_model
    replacement_env = SM.EnvironmentModel(
        planet=env.planet, EI=env.EI,
        density_model=SM.ExponentialAtmosphereModel(1e-12, 550_000.0, 60_000.0),
        ephemerides_model=env.ephemerides_model, thermal_model=env.thermal_model,
    )
    replacement_dynamics = SM.DynamicsModel(args.dynamics_model.spacecraft,
                                            (SM.InverseSquaredJ2GravityModel(),))
    updated = SC._with_configuration(args; environment_model=replacement_env,
                                     dynamics_model=replacement_dynamics)
    @test typeof(updated) !== typeof(args)
    @test updated.environment_model === replacement_env
    @test updated.dynamics_model === replacement_dynamics
    @test updated.dynamics_model.spacecraft === args.dynamics_model.spacecraft
    for name in fieldnames(typeof(args))
        name in (:environment_model, :dynamics_model) && continue
        @test getfield(updated, name) === getfield(args, name)
    end
    @test args.environment_model === env
    @test args.environment_model.density_model isa SM.NoAtmosphereModel
end

@testset "Smoke updates retain the existing explicit overrides" begin
    args = SC._with_configuration(copy_probe_config();
        mission_configuration=SM.MissionConfiguration(
            mission_type=SM.MissionOrbits, keplerian=false, number_of_orbits=7,
            mission_time=600.0, orientation_sim=true, num_steps_to_save=999, data_rate=3.0),
        simulation_settings=SM.SimulationSettings(
            results=true, verbose=true, results_directory="custom-results",
            generate_plots=true, generate_filenames=true, normalize=true, save_csv=true,
            checkpoint_enabled=true, checkpoint_interval_s=17.0,
            checkpoint_directory="custom-checkpoints", resume_from_checkpoint=true))
    withenv("SPACEAGORA_EXAMPLE_SMOKE" => "0") do
        @test TV._example_smoke_args(args) === args
    end
    for keep_results in (false, true)
        withenv("SPACEAGORA_EXAMPLE_SMOKE" => "1",
                "SPACEAGORA_EXAMPLE_SMOKE_RESULTS" => (keep_results ? "1" : "0"),
                "SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME" => "2.0") do
            smoke = TV._example_smoke_args(args)
            for name in fieldnames(typeof(args))
                name in (:mission_configuration, :simulation_settings) && continue
                @test getfield(smoke, name) === getfield(args, name)
            end
            mc = smoke.mission_configuration
            @test mc.mission_type === SM.MissionOrbits
            @test !mc.keplerian
            @test mc.number_of_orbits == 1
            @test mc.mission_time == 2.0
            @test mc.orientation_sim
            @test mc.num_steps_to_save == 200
            @test mc.data_rate == 3.0
            ss = smoke.simulation_settings
            @test ss.results == keep_results
            @test !ss.verbose
            @test ss.results_directory == joinpath(pwd(), "output")
            @test !ss.generate_plots
            @test ss.generate_filenames
            @test !ss.normalize
            @test ss.save_csv == keep_results
            # Preserve smoke's existing restart defaults, not the full-run state.
            @test !ss.checkpoint_enabled
            @test ss.checkpoint_interval_s == 300.0
            @test ss.checkpoint_directory == ""
            @test !ss.resume_from_checkpoint
        end
    end
    @test args.simulation_settings.checkpoint_enabled
    @test args.simulation_settings.results_directory == "custom-results"
    @test args.mission_configuration.number_of_orbits == 7
end

@testset "Smoke solver choice survives a conflicting environment" begin
    withenv("SPACEAGORA_EXAMPLE_SMOKE" => "1",
            "SPACEAGORA_EXAMPLE_SMOKE_RESULTS" => "0",
            "SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME" => "2.0",
            "SPACEAGORA_SOLVER_MODE" => "tsit5",
            "SPACEAGORA_SOLVER_MAXITERS" => "99999") do
        args = copy_probe_config()
        smoke = TV._example_smoke_args(args)
        @test smoke.solver_config === args.solver_config
        @test smoke.solver_config.maxiters == 12345
        @test smoke.environment_model.planet === args.environment_model.planet
        @test smoke.dynamics_model.spacecraft === args.dynamics_model.spacecraft
        frame_before = copy(args.environment_model.planet.L_PI)
        typed = SpaceAGORA.run_simulation(smoke; return_solver_metadata=true)
        @test typed.retcode == "Success"
        @test typed.solver_mode == "rodas5p"
        @test !isempty(typed.solver_trace)
        @test all(step -> step.solver == "Rodas5P", typed.solver_trace)
        @test args.environment_model.planet.L_PI == frame_before

        # An explicit nothing retains the existing environment-selected behavior.
        env_args = SC._with_configuration(args; solver_config=nothing)
        env_smoke = TV._example_smoke_args(env_args)
        @test env_smoke.solver_config === nothing
        from_env = SpaceAGORA.run_simulation(env_smoke; return_solver_metadata=true)
        @test from_env.retcode == "Success"
        @test from_env.solver_mode == "tsit5"
        @test !isempty(from_env.solver_trace)
        @test all(step -> step.solver == "Tsit5", from_env.solver_trace)
        @test args.environment_model.planet.L_PI == frame_before
    end
end

println("configuration_copy_probes_ok")
