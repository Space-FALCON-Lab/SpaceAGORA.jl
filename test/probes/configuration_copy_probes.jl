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

function copy_probe_preserved_fields(updated, original, changed::Tuple)
    for name in fieldnames(typeof(original))
        name in changed && continue
        @test getfield(updated, name) === getfield(original, name)
    end
end

function copy_probe_maneuver_scenario(; enabled::Bool)
    return TV.OrbitEventsScenarioConfig(
        name="configuration_copy_probe", planet_name="earth",
        telemetry_peri_path="", telemetry_apo_path="",
        target_orbits_quick=2, target_orbits_full=3,
        compare_points_quick=2, compare_points_full=3, min_eval_points=1,
        units_x="orbit", units_y=Dict("peri" => "km", "apo" => "km"),
        tolerances_quick=Dict{String, TV.EventTolerance}(),
        tolerances_full=Dict{String, TV.EventTolerance}(),
        initial_time=SM.InitialTime(year=2020, month=1, day=1),
        ra_m=7.1e6, rp_altitude_m=550_000.0, i_deg=53.0,
        aop_deg=0.0, raan_deg=10.0, ta_deg=15.0,
        spacecraft=TV.SpacecraftConfig(
            bus_dims=(1.0, 1.0, 1.0), panel_dims=(0.1, 0.2, 0.3),
            bus_mass_kg=500.0, panel_mass_each_kg=0.0,
            panel_offset_y_m=0.5, prop_mass_kg=10.0, id=1),
        gravity_model=:inverse_squared, EI_km=300.0,
        maneuver_orbit_numbers=enabled ? Int64[2, 4] : Int64[],
        maneuver_delta_v_mps=enabled ? [0.1, 0.2] : Float64[],
        maneuver_thrust_n=1.5, maneuver_isp_s=230.0,
        maneuver_guidance_rate_s=13.0, maneuver_control_rate_s=7.0,
    )
end

@testset "Telemetry configuration adapters preserve unrelated fields" begin
    args = SC._with_configuration(copy_probe_config();
        file_paths=SM.FilePaths(results="configuration-copy-results"),
        mission_configuration=SM.MissionConfiguration(
            mission_type=SM.MissionTime, keplerian=false, number_of_orbits=7,
            mission_time=17.0, orientation_sim=true, num_steps_to_save=29, data_rate=0.25),
        simulation_settings=SM.SimulationSettings(
            results=false, verbose=true, results_directory="adapter-results",
            generate_plots=true, generate_filenames=true, normalize=true, save_csv=false,
            checkpoint_enabled=true, checkpoint_interval_s=17.0,
            checkpoint_directory="adapter-checkpoints", resume_from_checkpoint=true),
        integration_tolerances=SM.IntegrationTolerances(
            reltol=2e-6, abstol=3e-8, reltol_mass=4e-7, dt_max=0.75),
    )

    @testset "wind" begin
        for wind in (true, false)
            updated = TV._with_environment_wind(args, wind)
            copy_probe_preserved_fields(updated, args, (:environment_model,))
            copy_probe_preserved_fields(updated.environment_model, args.environment_model, (:wind,))
            @test updated.environment_model.wind == wind
        end
    end

    @testset "campaign maneuvers" begin
        @test TV._with_campaign_maneuvers(args, copy_probe_maneuver_scenario(enabled=false)) === args
        scenario = copy_probe_maneuver_scenario(enabled=true)
        updated = TV._with_campaign_maneuvers(args, scenario)
        copy_probe_preserved_fields(updated, args, (:guidance_model, :control_model))
        guidance = only(updated.guidance_model.guidance_effectors)
        thruster = only(updated.control_model.control_effectors)
        @test guidance.maneuver_orbit_number == [2, 4]
        @test guidance.maneuver_Δv == [0.1, 0.2]
        @test isempty(guidance.maneuver_flight_apoapsis_radius_m)
        @test updated.guidance_model.guidance_rates == [13.0]
        @test updated.control_model.control_rates == [7.0]
        @test thruster.thrust == [1.5]
        @test thruster.Isp == [230.0]
        @test thruster.Δv == [0.0]
        @test thruster.start_burn_time == [-1.0]
        @test thruster.stop_burn_time == [-1.0]
        @test isempty(args.guidance_model.guidance_effectors)
        @test isempty(args.control_model.control_effectors)
    end

    @testset "orbit mission" begin
        updated = TV._with_orbit_mission(args, 3, 42.0)
        copy_probe_preserved_fields(updated, args, (:mission_configuration,))
        copy_probe_preserved_fields(updated.mission_configuration, args.mission_configuration,
                                   (:mission_type, :number_of_orbits, :mission_time))
        @test updated.mission_configuration.mission_type === SM.MissionOrbits
        @test updated.mission_configuration.number_of_orbits == 3
        @test updated.mission_configuration.mission_time == 42.0
    end

    @testset "study settings" begin
        # Pin study-specific controls so ambient preferences cannot change the
        # expected override values; all component-tolerance defaults stay intact.
        withenv("SPACEAGORA_TELEMETRY_RELTOL_ORBIT" => "2e-8",
                "SPACEAGORA_TELEMETRY_ABSTOL_ORBIT" => "3e-10",
                "SPACEAGORA_TELEMETRY_RELTOL_ATM" => "4e-8",
                "SPACEAGORA_TELEMETRY_ABSTOL_ATM" => "5e-10",
                "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => "7.0",
                "SPACEAGORA_TELEMETRY_DT_MAX_ATM" => "0.1") do
            for quick in (false, true)
                updated = TV._with_study_settings(args; quick=quick)
                copy_probe_preserved_fields(updated, args,
                    (:simulation_settings, :mission_configuration, :integration_tolerances))
                copy_probe_preserved_fields(updated.mission_configuration, args.mission_configuration,
                                           (:num_steps_to_save,))
                @test updated.mission_configuration.num_steps_to_save == 2000
                expected_settings = SM.SimulationSettings(
                    results=true, verbose=false, results_directory="adapter-results",
                    generate_plots=false, generate_filenames=false, normalize=false, save_csv=true)
                for name in fieldnames(typeof(expected_settings))
                    @test getfield(updated.simulation_settings, name) == getfield(expected_settings, name)
                end
                expected_tolerances = SM.IntegrationTolerances(
                    reltol_orbit=2e-8, abstol_orbit=3e-10, dt_max_orbit=7.0,
                    reltol_atmosphere=4e-8, abstol_atmosphere=5e-10, dt_max_atmosphere=0.1)
                for name in fieldnames(typeof(expected_tolerances))
                    @test getfield(updated.integration_tolerances, name) == getfield(expected_tolerances, name)
                end
            end
        end
    end
end

@testset "Ensemble configuration copies retain sharing and member paths" begin
    campaigns = SpaceAGORA.SimulationCampaigns
    args = copy_probe_config()
    spacecraft = only(args.dynamics_model.spacecraft)
    member_tag = "sat_1_id_1"
    # With no outputs/checkpoints, preserve the settings object as before.
    unchanged_paths = campaigns._ensemble_member_configuration(args, spacecraft, member_tag)
    copy_probe_preserved_fields(unchanged_paths, args, (:dynamics_model,))
    @test unchanged_paths.dynamics_model.spacecraft !== args.dynamics_model.spacecraft
    @test only(unchanged_paths.dynamics_model.spacecraft) === spacecraft
    @test unchanged_paths.dynamics_model.dynamic_effectors === args.dynamics_model.dynamic_effectors

    for checkpoint_dir in ("", "campaign-checkpoints")
        settings = SM.SimulationSettings(
            results=false, verbose=true, results_directory="campaign-results",
            generate_plots=false, generate_filenames=true, normalize=false, save_csv=false,
            checkpoint_enabled=true, checkpoint_interval_s=17.0,
            checkpoint_directory=checkpoint_dir, resume_from_checkpoint=true)
        configured = SC._with_configuration(args; simulation_settings=settings)
        member = campaigns._ensemble_member_configuration(configured, spacecraft, member_tag)
        copy_probe_preserved_fields(member, configured, (:simulation_settings, :dynamics_model))
        copy_probe_preserved_fields(member.simulation_settings, settings,
                                   (:results_directory, :checkpoint_directory))
        @test member.simulation_settings.results_directory ==
            (isempty(checkpoint_dir) ? joinpath("campaign-results", member_tag) : "campaign-results")
        @test member.simulation_settings.checkpoint_directory ==
            (isempty(checkpoint_dir) ? "" : joinpath(checkpoint_dir, member_tag))
        @test only(member.dynamics_model.spacecraft) === spacecraft
        @test member.dynamics_model.dynamic_effectors === configured.dynamics_model.dynamic_effectors
        @test configured.simulation_settings.results_directory == "campaign-results"
        @test configured.simulation_settings.checkpoint_directory == checkpoint_dir
    end
end

@testset "Telemetry runner retains its solver and retry policy" begin
    # The harness deliberately selects its own solver, even when its input
    # configuration carries a typed solver. A one-iteration first attempt plus
    # a small maximum step guarantees that the real MaxIters retry is exercised.
    withenv("SPACEAGORA_SOLVER_MODE" => "rodas5p",
            "SPACEAGORA_SOLVER_MAXITERS" => "12345",
            "SPACEAGORA_TELEMETRY_SOLVER_MODE" => "tsit5",
            "SPACEAGORA_TELEMETRY_SOLVER_MAXITERS" => "1",
            "SPACEAGORA_TELEMETRY_SOLVER_MAXITERS_RETRY" => "10000",
            "SPACEAGORA_SOLVER_SAVE_EVERYSTEP" => nothing,
            "SPACEAGORA_SOLVER_SAVE_ON" => nothing,
            "SPACEAGORA_RHS_CALIBRATE" => "off") do
        mktempdir() do tmp
            checkpoint_dir = joinpath(tmp, "unused-input-checkpoints")
            args = SC._with_configuration(copy_probe_config();
                mission_configuration=SM.MissionConfiguration(
                    mission_time=2.0, keplerian=true, data_rate=0.1, num_steps_to_save=50),
                integration_tolerances=SM.IntegrationTolerances(
                    reltol_orbit=1e-8, abstol_orbit=1e-8, dt_max_orbit=0.05,
                    reltol_atmosphere=1e-8, abstol_atmosphere=1e-8, dt_max_atmosphere=0.05),
                simulation_settings=SM.SimulationSettings(
                    results=false, verbose=false, generate_plots=false, save_csv=false,
                    checkpoint_enabled=true, checkpoint_directory=checkpoint_dir,
                    resume_from_checkpoint=true),
            )
            # Capture the entire environment, including absent keys, to detect
            # leaks from the runner's nested withenv scopes after its retry.
            environment_before = Dict(ENV)
            directory_before = pwd()
            result = TV._run_simulation_dataframe(args, "configuration_copy_retry",
                                                  TV.AtmosphereTruthConfig(), :quick)
            @test result.solver_info.solver_retcode == "Success"
            @test result.solver_info.solver_mode == "tsit5"
            @test result.solver_info.solver_sequence == "Tsit5"
            @test result.solver_info.solver_maxiters_retry_used
            @test result.solver_info.solver_maxiters == 10000
            @test size(result.results_df, 1) > 1
            @test args.solver_config.solver_mode === :rodas5p
            @test args.solver_config.maxiters == 12345
            @test args.simulation_settings.checkpoint_enabled
            @test args.simulation_settings.resume_from_checkpoint
            @test !ispath(checkpoint_dir)
            @test pwd() == directory_before
            @test Dict(ENV) == environment_before
        end
    end
end

println("configuration_copy_probes_ok")
