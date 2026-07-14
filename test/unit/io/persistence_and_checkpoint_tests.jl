@testset "Deterministic Smoke + No-Drag Energy Invariant" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=400e3, ν_deg=175.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=1200.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df = run_case(args)
    @test nrow(df) > 100
    @test all(isfinite, df.sc1_pos_1)
    @test all(isfinite, df.sc1_vel_1)

    eps = specific_energy(df, EARTH.μ)
    energy_drift = maximum(abs.(eps .- first(eps))) / abs(first(eps))
    @test energy_drift < 1e-8
end

@testset "Deterministic Replay (No-Drag)" begin
    sc = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df1 = run_case_silent(args)
    df2 = run_case_silent(args)
    @test nrow(df1) == nrow(df2)
    @test nrow(df1) > 10

    sample_idxs = round.(Int, range(1, nrow(df1), length=8))
    for idx in sample_idxs
        t = Float64(df1.time[idx])
        p1 = SVector{3, Float64}(Float64(df1.sc1_pos_1[idx]), Float64(df1.sc1_pos_2[idx]), Float64(df1.sc1_pos_3[idx]))
        v1 = SVector{3, Float64}(Float64(df1.sc1_vel_1[idx]), Float64(df1.sc1_vel_2[idx]), Float64(df1.sc1_vel_3[idx]))
        p2 = SVector{3, Float64}(
            interp_linear(df2.time, df2.sc1_pos_1, t),
            interp_linear(df2.time, df2.sc1_pos_2, t),
            interp_linear(df2.time, df2.sc1_pos_3, t)
        )
        v2 = SVector{3, Float64}(
            interp_linear(df2.time, df2.sc1_vel_1, t),
            interp_linear(df2.time, df2.sc1_vel_2, t),
            interp_linear(df2.time, df2.sc1_vel_3, t)
        )
        @test norm(p1 - p2) < 0.1
        @test norm(v1 - v2) < 1e-4
    end
end


@testset "Typed Results Bundle Persistence" begin
    sc = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    settings = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    mktempdir() do tmp
        cd(tmp) do
            run_simulation(args)
            @test isfile(joinpath("output", "simulation_results.csv"))

            bundle_prefix = joinpath("output", "simulation_results")
            feather_path = bundle_prefix * ".feather"
            bundle_csv_path = bundle_prefix * ".csv"
            manifest_path = bundle_prefix * ".manifest.toml"

            @test isfile(feather_path)
            @test isfile(bundle_csv_path)
            @test isfile(manifest_path)

            manifest = TOML.parsefile(manifest_path)
            @test get(manifest, "schema_version", "") == "1"
            @test get(manifest, "steps", 0) > 10
            @test haskey(manifest, "files")

            files = manifest["files"]
            @test haskey(files, "feather")
            @test haskey(files, "csv")
            @test get(files["feather"], "size_bytes", 0) > 0
            @test length(get(files["feather"], "sha256", "")) == 64
        end
    end
end

@testset "IOConfig Owner Coverage" begin
    mktempdir() do tmp
        results_dir = joinpath(tmp, "results")
        mkpath(results_dir)

        args_default_ckpt = (; simulation_settings=(; results_directory=results_dir, checkpoint_directory=""))
        args_custom_ckpt = (; simulation_settings=(; results_directory=results_dir, checkpoint_directory=joinpath(tmp, "custom_ckpt")))

        bundle_prefix = SimulationModel.IOConfig._results_bundle_prefix(args_default_ckpt)
        @test bundle_prefix == joinpath(results_dir, "simulation_results")

        results_csv = SimulationModel.IOConfig._results_csv_path(args_default_ckpt)
        @test results_csv == joinpath(results_dir, "simulation_results.csv")

        collision_csv = SimulationModel.IOConfig._collision_results_csv_path(args_default_ckpt)
        @test startswith(collision_csv, joinpath(results_dir, "simulation_results."))
        @test endswith(collision_csv, ".csv")

        @test SimulationModel.IOConfig._checkpoint_directory(args_default_ckpt) == joinpath(results_dir, "checkpoints")
        @test SimulationModel.IOConfig._checkpoint_directory(args_custom_ckpt) == joinpath(tmp, "custom_ckpt")

        paths_default = SimulationModel.IOConfig._checkpoint_paths(args_default_ckpt)
        @test paths_default.data == joinpath(results_dir, "checkpoints", "simulation_checkpoint.bin")
        @test paths_default.manifest == joinpath(results_dir, "checkpoints", "simulation_checkpoint.manifest.toml")

        paths_custom = SimulationModel.IOConfig._checkpoint_paths(args_custom_ckpt)
        @test paths_custom.data == joinpath(tmp, "custom_ckpt", "simulation_checkpoint.bin")
        @test paths_custom.manifest == joinpath(tmp, "custom_ckpt", "simulation_checkpoint.manifest.toml")
    end
end

@testset "Checkpoint Resume Parity" begin
    sc = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    baseline_settings = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )
    baseline_args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=240.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=baseline_settings,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    df_full = run_case_silent(baseline_args)

    checkpoint_dir = joinpath("output", "checkpoints")
    checkpoint_data_path = joinpath(checkpoint_dir, "simulation_checkpoint.bin")
    checkpoint_manifest_path = joinpath(checkpoint_dir, "simulation_checkpoint.manifest.toml")

    phase1_settings = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false,
        checkpoint_enabled=true,
        checkpoint_interval_s=40.0,
        checkpoint_directory=checkpoint_dir,
        resume_from_checkpoint=false
    )
    phase1_args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=phase1_settings,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    resume_settings = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false,
        checkpoint_enabled=true,
        checkpoint_interval_s=40.0,
        checkpoint_directory=checkpoint_dir,
        resume_from_checkpoint=true
    )
    resume_args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=240.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=resume_settings,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    mktempdir() do tmp
        cd(tmp) do
            run_simulation(phase1_args)
            @test isfile(checkpoint_data_path)
            @test isfile(checkpoint_manifest_path)

            run_simulation(resume_args)
            @test isfile(checkpoint_data_path)
            @test isfile(checkpoint_manifest_path)
            @test isfile(joinpath("output", "simulation_results.csv"))

            df_resume = CSV.read(joinpath("output", "simulation_results.csv"), DataFrame)
            @test nrow(df_resume) >= 2
            @test issorted(df_resume.time)
            @test Float64(df_resume.time[1]) > 0.0
            @test Float64(df_resume.time[end]) > Float64(df_resume.time[1])
            @test abs(Float64(df_resume.time[end]) - Float64(df_full.time[end])) < 1e-8

            p_full = SVector{3, Float64}(Float64(df_full.sc1_pos_1[end]), Float64(df_full.sc1_pos_2[end]), Float64(df_full.sc1_pos_3[end]))
            v_full = SVector{3, Float64}(Float64(df_full.sc1_vel_1[end]), Float64(df_full.sc1_vel_2[end]), Float64(df_full.sc1_vel_3[end]))
            p_resume = SVector{3, Float64}(Float64(df_resume.sc1_pos_1[end]), Float64(df_resume.sc1_pos_2[end]), Float64(df_resume.sc1_pos_3[end]))
            v_resume = SVector{3, Float64}(Float64(df_resume.sc1_vel_1[end]), Float64(df_resume.sc1_vel_2[end]), Float64(df_resume.sc1_vel_3[end]))

            @test norm(p_resume - p_full) < 1.0
            @test norm(v_resume - v_full) < 1e-3
            @test abs(Float64(df_resume.sc1_mass[end]) - Float64(df_full.sc1_mass[end])) < 1e-8
        end
    end
end

@testset "Gravity Backbone Checkpoint Resume + Save Fields" begin
    sc = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    baseline_settings = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )
    baseline_args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=180.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel(),
        simulation_settings=baseline_settings,
        tolerances=IntegrationTolerances(reltol_orbit=1e-8, abstol_orbit=1e-8, dt_max_orbit=20.0)
    )
    df_full = withenv(
        "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
        "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "2.0"
    ) do
        run_case_silent(baseline_args)
    end
    @test all(abs.(Float64.(df_full.sc1_heat_rate)) .< 1e-12)
    @test all(abs.(Float64.(df_full.sc1_heat_load)) .< 1e-12)

    checkpoint_dir = joinpath("output", "checkpoints")
    checkpoint_data_path = joinpath(checkpoint_dir, "simulation_checkpoint.bin")
    checkpoint_manifest_path = joinpath(checkpoint_dir, "simulation_checkpoint.manifest.toml")

    phase1_args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=90.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel(),
        simulation_settings=SimulationSettings(
            results=true,
            verbose=false,
            generate_plots=false,
            results_directory="output",
            save_csv=true,
            normalize=false,
            checkpoint_enabled=true,
            checkpoint_interval_s=30.0,
            checkpoint_directory=checkpoint_dir,
            resume_from_checkpoint=false
        ),
        tolerances=IntegrationTolerances(reltol_orbit=1e-8, abstol_orbit=1e-8, dt_max_orbit=20.0)
    )

    resume_args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=180.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel(),
        simulation_settings=SimulationSettings(
            results=true,
            verbose=false,
            generate_plots=false,
            results_directory="output",
            save_csv=true,
            normalize=false,
            checkpoint_enabled=true,
            checkpoint_interval_s=30.0,
            checkpoint_directory=checkpoint_dir,
            resume_from_checkpoint=true
        ),
        tolerances=IntegrationTolerances(reltol_orbit=1e-8, abstol_orbit=1e-8, dt_max_orbit=20.0)
    )

    tsit_resume_args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=180.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel(),
        simulation_settings=SimulationSettings(
            results=true,
            verbose=false,
            generate_plots=false,
            results_directory="output",
            save_csv=true,
            normalize=false,
            checkpoint_enabled=true,
            checkpoint_interval_s=30.0,
            checkpoint_directory=checkpoint_dir,
            resume_from_checkpoint=true
        ),
        tolerances=IntegrationTolerances(reltol_orbit=1e-8, abstol_orbit=1e-8, dt_max_orbit=20.0)
    )

    mktempdir() do tmp
        cd(tmp) do
            withenv(
                "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
                "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "2.0"
            ) do
                run_simulation(phase1_args)
            end
            @test isfile(checkpoint_data_path)
            @test isfile(checkpoint_manifest_path)

            ckpt_manifest = TOML.parsefile(checkpoint_manifest_path)
            @test get(ckpt_manifest, "solver_mode", "") == "gravity_backbone_split"

            ckpt = _load_checkpoint(resume_args)
            @test ckpt !== nothing
            @test ckpt.solver_mode == "gravity_backbone_split"
            @test _is_gravity_backbone_state(ckpt.u)

            withenv(
                "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
                "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "2.0"
            ) do
                run_simulation(resume_args)
            end
            df_resume = CSV.read(joinpath("output", "simulation_results.csv"), DataFrame)
            @test nrow(df_resume) >= 2

            p_full = SVector{3, Float64}(Float64(df_full.sc1_pos_1[end]), Float64(df_full.sc1_pos_2[end]), Float64(df_full.sc1_pos_3[end]))
            v_full = SVector{3, Float64}(Float64(df_full.sc1_vel_1[end]), Float64(df_full.sc1_vel_2[end]), Float64(df_full.sc1_vel_3[end]))
            p_resume = SVector{3, Float64}(Float64(df_resume.sc1_pos_1[end]), Float64(df_resume.sc1_pos_2[end]), Float64(df_resume.sc1_pos_3[end]))
            v_resume = SVector{3, Float64}(Float64(df_resume.sc1_vel_1[end]), Float64(df_resume.sc1_vel_2[end]), Float64(df_resume.sc1_vel_3[end]))
            @test norm(p_resume - p_full) < 5.0
            @test norm(v_resume - v_full) < 5e-3

            withenv("SPACEAGORA_SOLVER_MODE" => "tsit5") do
                @test_throws ArgumentError run_simulation(tsit_resume_args)
            end
        end
    end
end

@testset "Gravity Backbone Kick Checkpoint Resume Parity" begin
    sc = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    baseline_args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=90.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            NBodyGravityModel(["moon"], "Earth", SPICE_PATH),
            SolarRadiationPressureModel(1.2, 12.0; direct=true, albedo=true, ir=true),
        ),
        keplerian=true,
        ephemerides_model=SpiceEphemeridesModel(),
        simulation_settings=SimulationSettings(
            results=true,
            verbose=false,
            generate_plots=false,
            results_directory="output",
            save_csv=true,
            normalize=false
        ),
        tolerances=IntegrationTolerances(reltol_orbit=1e-8, abstol_orbit=1e-8, dt_max_orbit=20.0)
    )

    df_full = withenv(
        "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
        "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "2.0"
    ) do
        run_case_silent(baseline_args)
    end
    @test all(abs.(Float64.(df_full.sc1_heat_rate)) .< 1e-12)
    @test all(abs.(Float64.(df_full.sc1_heat_load)) .< 1e-12)

    checkpoint_dir = joinpath("output", "checkpoints")
    phase1_args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=45.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            NBodyGravityModel(["moon"], "Earth", SPICE_PATH),
            SolarRadiationPressureModel(1.2, 12.0; direct=true, albedo=true, ir=true),
        ),
        keplerian=true,
        ephemerides_model=SpiceEphemeridesModel(),
        simulation_settings=SimulationSettings(
            results=true,
            verbose=false,
            generate_plots=false,
            results_directory="output",
            save_csv=true,
            normalize=false,
            checkpoint_enabled=true,
            checkpoint_interval_s=15.0,
            checkpoint_directory=checkpoint_dir,
            resume_from_checkpoint=false
        ),
        tolerances=IntegrationTolerances(reltol_orbit=1e-8, abstol_orbit=1e-8, dt_max_orbit=20.0)
    )

    resume_args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=90.0,
        EI_km=120.0,
        dynamic_effectors=(
            InverseSquaredGravityModel(),
            NBodyGravityModel(["moon"], "Earth", SPICE_PATH),
            SolarRadiationPressureModel(1.2, 12.0; direct=true, albedo=true, ir=true),
        ),
        keplerian=true,
        ephemerides_model=SpiceEphemeridesModel(),
        simulation_settings=SimulationSettings(
            results=true,
            verbose=false,
            generate_plots=false,
            results_directory="output",
            save_csv=true,
            normalize=false,
            checkpoint_enabled=true,
            checkpoint_interval_s=15.0,
            checkpoint_directory=checkpoint_dir,
            resume_from_checkpoint=true
        ),
        tolerances=IntegrationTolerances(reltol_orbit=1e-8, abstol_orbit=1e-8, dt_max_orbit=20.0)
    )

    mktempdir() do tmp
        cd(tmp) do
            withenv(
                "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
                "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "2.0"
            ) do
                run_simulation(phase1_args)
                run_simulation(resume_args)
            end

            df_resume = CSV.read(joinpath("output", "simulation_results.csv"), DataFrame)
            @test nrow(df_resume) >= 2
            @test all(abs.(Float64.(df_resume.sc1_heat_rate)) .< 1e-12)
            @test all(abs.(Float64.(df_resume.sc1_heat_load)) .< 1e-12)

            p_full = SVector{3, Float64}(Float64(df_full.sc1_pos_1[end]), Float64(df_full.sc1_pos_2[end]), Float64(df_full.sc1_pos_3[end]))
            v_full = SVector{3, Float64}(Float64(df_full.sc1_vel_1[end]), Float64(df_full.sc1_vel_2[end]), Float64(df_full.sc1_vel_3[end]))
            p_resume = SVector{3, Float64}(Float64(df_resume.sc1_pos_1[end]), Float64(df_resume.sc1_pos_2[end]), Float64(df_resume.sc1_pos_3[end]))
            v_resume = SVector{3, Float64}(Float64(df_resume.sc1_vel_1[end]), Float64(df_resume.sc1_vel_2[end]), Float64(df_resume.sc1_vel_3[end]))

            @test norm(p_resume - p_full) < 25.0
            @test norm(v_resume - v_full) < 2.5e-2
        end
    end
end

@testset "Checkpoint Guards + Missing Resume + Bundle Toggle" begin
    base_spacecraft() = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    base_settings(; kwargs...) = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false;
        kwargs...
    )

    args_bad_interval = build_config(
        spacecraft=base_spacecraft(),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=base_settings(
            checkpoint_enabled=true,
            checkpoint_interval_s=0.0,
            checkpoint_directory=joinpath("output", "checkpoints")
        ),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    @test_throws ArgumentError run_simulation(args_bad_interval)

    args_resume_missing = build_config(
        spacecraft=base_spacecraft(),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=base_settings(
            checkpoint_enabled=true,
            checkpoint_interval_s=20.0,
            checkpoint_directory=joinpath("output", "checkpoints"),
            resume_from_checkpoint=true
        ),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    mktempdir() do tmp
        cd(tmp) do
            @test_logs (:warn, r"resume_from_checkpoint=true but no checkpoint file was found") run_simulation(args_resume_missing)
            @test isfile(joinpath("output", "simulation_results.csv"))

            df = CSV.read(joinpath("output", "simulation_results.csv"), DataFrame)
            expected_rows = Int(floor(
                args_resume_missing.mission_configuration.mission_time /
                args_resume_missing.mission_configuration.data_rate
            )) + 1
            @test nrow(df) >= expected_rows
            @test abs(Float64(df.time[end]) - 60.0) < 1e-8

            ckpt_manifest_path = joinpath("output", "checkpoints", "simulation_checkpoint.manifest.toml")
            @test isfile(ckpt_manifest_path)
            ckpt_manifest = TOML.parsefile(ckpt_manifest_path)
            @test get(ckpt_manifest, "schema_version", "") == "1"
            @test get(ckpt_manifest, "time_s", 0.0) > 0.0
            @test get(ckpt_manifest, "data_size_bytes", 0) > 0
            @test length(get(ckpt_manifest, "data_sha256", "")) == 64
        end
    end

    args_bundle_disabled = build_config(
        spacecraft=base_spacecraft(),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=base_settings(),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    mktempdir() do tmp
        cd(tmp) do
            withenv("SPACEAGORA_SAVE_BUNDLE" => "0") do
                run_simulation(args_bundle_disabled)
            end
            @test isfile(joinpath("output", "simulation_results.csv"))
            @test !isfile(joinpath("output", "simulation_results.feather"))
            @test !isfile(joinpath("output", "simulation_results.manifest.toml"))
        end
    end
end


@testset "Units/Normalization Consistency Audit" begin
    function strip_comments(src::String)
        no_block = replace(src, r"#=.*?=#"s => "")
        no_line = map(line -> first(split(line, '#'; limit=2)), split(no_block, '\n'; keepempty=true))
        return join(no_line, "\n")
    end

    run_src = strip_comments(read(joinpath(REPO_ROOT, "src", "simulation", "engine", "execution.jl"), String))
    legacy_elements_path = joinpath(REPO_ROOT, "src", "simulation", "execution", "simulation_elements.jl")

    # Canonical engine path should stay SI-native.
    @test !occursin("cnf.DU", run_src)
    @test !occursin("cnf.TU", run_src)
    @test !occursin("cnf.MU", run_src)
    @test !isfile(legacy_elements_path)

    sc = make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0)
    settings_norm_true = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        normalize=true
    )
    settings_norm_false = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        normalize=false
    )

    args_norm_true = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=700.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_norm_true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    args_norm_false = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=700.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_norm_false,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    # Build-state sanity: typed initial conditions are SI-scale, not O(1) normalized values.
    u0 = build_initial_conditions(args_norm_true)
    @test norm(SVector{3, Float64}(u0.sc[1].pos)) > 1e6
    @test norm(SVector{3, Float64}(u0.sc[1].vel)) > 1e3
    @test u0.sc[1].mass > 1.0

    @test_throws ArgumentError run_case_silent(args_norm_true)

    df_norm_true = withenv(
        "SPACEAGORA_ALLOW_TYPED_NORMALIZE" => "1",
        "SPACEAGORA_WARN_NORMALIZE" => "0"
    ) do
        run_case_silent(args_norm_true)
    end
    df_norm_false = run_case_silent(args_norm_false)

    @test nrow(df_norm_true) == nrow(df_norm_false)
    @test nrow(df_norm_true) > 10
    sample_idxs = round.(Int, range(1, nrow(df_norm_true), length=8))
    for idx in sample_idxs
        t = Float64(df_norm_true.time[idx])
        p_true = SVector{3, Float64}(Float64(df_norm_true.sc1_pos_1[idx]), Float64(df_norm_true.sc1_pos_2[idx]), Float64(df_norm_true.sc1_pos_3[idx]))
        v_true = SVector{3, Float64}(Float64(df_norm_true.sc1_vel_1[idx]), Float64(df_norm_true.sc1_vel_2[idx]), Float64(df_norm_true.sc1_vel_3[idx]))
        p_false = SVector{3, Float64}(
            interp_linear(df_norm_false.time, df_norm_false.sc1_pos_1, t),
            interp_linear(df_norm_false.time, df_norm_false.sc1_pos_2, t),
            interp_linear(df_norm_false.time, df_norm_false.sc1_pos_3, t)
        )
        v_false = SVector{3, Float64}(
            interp_linear(df_norm_false.time, df_norm_false.sc1_vel_1, t),
            interp_linear(df_norm_false.time, df_norm_false.sc1_vel_2, t),
            interp_linear(df_norm_false.time, df_norm_false.sc1_vel_3, t)
        )
        @test norm(p_true - p_false) < 0.1
        @test norm(v_true - v_false) < 1e-4
    end
end

@testset "Normalize Flag Runtime Policy" begin
    args_warn = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=true),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    withenv("SPACEAGORA_ALLOW_TYPED_NORMALIZE" => "0") do
        @test_throws ArgumentError run_simulation(args_warn)
    end

    withenv(
        "SPACEAGORA_ALLOW_TYPED_NORMALIZE" => "1",
        "SPACEAGORA_WARN_NORMALIZE" => "1"
    ) do
        _normalize_warning_emitted[] = false
        @test_logs (:warn, r"normalize=true is legacy-only") run_simulation(args_warn)
        @test _normalize_warning_emitted[] == true
        @test_logs run_simulation(args_warn)
    end

    withenv(
        "SPACEAGORA_ALLOW_TYPED_NORMALIZE" => "1",
        "SPACEAGORA_WARN_NORMALIZE" => "0"
    ) do
        _normalize_warning_emitted[] = false
        @test_logs run_simulation(args_warn)
        @test _normalize_warning_emitted[] == false
    end
end

@testset "Run Simulation Debug Branches" begin
    debug_thruster = TimedTangentialThrusterModel(5.0, 1.0, 0.0, 120.0)
    args_debug_control = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(debug_thruster,),
        control_rates=[1.0],
        keplerian=true,
        simulation_settings=SimulationSettings(results=true, verbose=false, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    withenv("SPACEAGORA_DEBUG_CONTROL" => "1") do
        _, output = run_case_capture_stdout(args_debug_control)
        @test occursin("Applying control effect for spacecraft 1", output)
        @test occursin("Control force:", output)
    end

    args_debug_throw = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=20.0,
        EI_km=120.0,
        dynamic_effectors=(ThrowingForceModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    withenv("SPACEAGORA_DEBUG_INITIAL_DERIVATIVE" => "1") do
        @test_logs (:error, r"The derivative function itself crashed!") begin
            @test_throws ErrorException run_simulation(args_debug_throw)
        end
    end

    args_debug_nan = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=20.0,
        EI_km=120.0,
        dynamic_effectors=(NaNForceModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    nan_output = ""
    withenv("SPACEAGORA_DEBUG_INITIAL_DERIVATIVE" => "1") do
        mktemp() do path, io
            redirect_stdout(io) do
                @test_throws ErrorException run_simulation(args_debug_nan)
            end
            flush(io)
            seekstart(io)
            nan_output = read(io, String)
        end
    end
    @test occursin("--- INITIAL NaN DETECTED ---", nan_output)
    @test occursin("NaN found in Satellite 1 derivative!", nan_output)
    @test occursin("Pos:", nan_output)
    @test occursin("Vel:", nan_output)

    args_debug_nan_param = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=420e3, ν_deg=165.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=20.0,
        EI_km=120.0,
        dynamic_effectors=(NaNParamForceModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    nan_param_output = ""
    withenv("SPACEAGORA_DEBUG_INITIAL_DERIVATIVE" => "1") do
        mktemp() do path, io
            redirect_stdout(io) do
                @test_throws ErrorException run_simulation(args_debug_nan_param)
            end
            flush(io)
            seekstart(io)
            nan_param_output = read(io, String)
        end
    end
    @test occursin("NaN found in parameter: p.shared_buffers.current_time[]", nan_param_output)

    nan_scan_output = ""
    mktemp() do path, io
        redirect_stdout(io) do
            _debug_print_nan_parameter_paths!(NaN, "p.scalar_test")
            _debug_print_nan_parameter_paths!([1.0, NaN], "p.array_test")
        end
        flush(io)
        seekstart(io)
        nan_scan_output = read(io, String)
    end
    @test occursin("NaN found in parameter: p.scalar_test", nan_scan_output)
    @test occursin("NaN found in parameter: p.array_test[2]", nan_scan_output)
end

@testset "Verbose Gating for Callback/Runtime Logs" begin
    settings_quiet = SimulationSettings(
        results=true,
        verbose=false,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )
    settings_verbose = SimulationSettings(
        results=true,
        verbose=true,
        generate_plots=false,
        results_directory="output",
        save_csv=true,
        normalize=false
    )

    args_orbit_quiet = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=400e3, ν_deg=175.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=400.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_quiet,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    _, out_orbit_quiet = run_case_capture_stdout(args_orbit_quiet)
    @test !occursin("Initial conditions:", out_orbit_quiet)
    @test !occursin("Orbit ", out_orbit_quiet)

    args_orbit_verbose = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=400e3, ν_deg=175.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=400.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=settings_verbose,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    _, out_orbit_verbose = run_case_capture_stdout(args_orbit_verbose)
    @test occursin("Initial conditions:", out_orbit_verbose)

    args_drag_quiet = build_config(
        spacecraft=make_spacecraft(ra_alt_m=220e3, rp_alt_m=100e3, ν_deg=180.0),
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=700.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        simulation_settings=settings_quiet,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )
    _, out_drag_quiet = run_case_capture_stdout(args_drag_quiet)
    @test !occursin("Switching to atmosphere integration", out_drag_quiet)
    @test !occursin("Impact detected", out_drag_quiet)

    args_drag_verbose = build_config(
        spacecraft=make_spacecraft(ra_alt_m=220e3, rp_alt_m=100e3, ν_deg=180.0),
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=700.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        simulation_settings=settings_verbose,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )
    _, out_drag_verbose = run_case_capture_stdout(args_drag_verbose)
    @test occursin("Switching to atmosphere integration", out_drag_verbose)
    @test occursin("Impact detected", out_drag_verbose)
end

