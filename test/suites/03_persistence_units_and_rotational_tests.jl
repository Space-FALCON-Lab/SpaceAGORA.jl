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
            @test nrow(df) > 10
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

@testset "Orbital Elements Round-Trip Invariant" begin
    ic = InitialCondition(
        ra=EARTH.Rp_e + 550e3,
        rp=EARTH.Rp_e + 300e3,
        i=37.0,
        ω=55.0,
        Ω=25.0,
        ν=130.0
    )
    r, v = orbitalelemtorv(ic, EARTH)
    oe = rvtoorbitalelement(SVector{3, Float64}(r), SVector{3, Float64}(v), EARTH)

    @test isapprox(oe[1], ic.a; rtol=1e-10, atol=1e-6)
    @test isapprox(oe[2], ic.e; rtol=1e-10, atol=1e-10)
    @test angle_distance(oe[3], ic.i) < 1e-8
    @test angle_distance(oe[4], ic.Ω) < 1e-8
    @test angle_distance(oe[5], ic.ω) < 1e-8
    @test angle_distance(oe[6], ic.ν) < 1e-8
end

@testset "Reference System Branch Coverage" begin
    ic_arg_wrap = InitialCondition(
        ra=EARTH.Rp_e + 700e3,
        rp=EARTH.Rp_e + 300e3,
        i=55.0,
        ω=300.0,
        Ω=40.0,
        ν=25.0
    )
    r_wrap, v_wrap = orbitalelemtorv(ic_arg_wrap, EARTH)
    oe_wrap = rvtoorbitalelement(SVector{3, Float64}(r_wrap), SVector{3, Float64}(v_wrap), EARTH)
    @test oe_wrap[5] > pi

    ic_circ_incl = InitialCondition(
        ra=EARTH.Rp_e + 500e3,
        rp=EARTH.Rp_e + 500e3,
        i=40.0,
        ω=0.0,
        Ω=15.0,
        ν=250.0
    )
    r_circ, v_circ = orbitalelemtorv(ic_circ_incl, EARTH)
    oe_circ = rvtoorbitalelement(SVector{3, Float64}(r_circ), SVector{3, Float64}(v_circ), EARTH)
    @test abs(oe_circ[2]) < 1e-10
    @test oe_circ[6] > pi

    planet_topo = (
        Rp_e=10.0,
        Rp_p=9.0,
        topography_function=(a, Clm, Slm, lat, lon, A) -> 9.5,
        Clm_topo=zeros(1, 1),
        Slm_topo=zeros(1, 1),
        A_topo=0.0
    )
    rp_topo = SVector{3, Float64}(0.1, 0.5, -8.944723618090451)
    lla = rtolatlong(rp_topo, planet_topo)
    @test all(isfinite, lla)

    topo_module_name = gensym(:ReferenceSystemTopoSandbox)
    Core.eval(Main, :(module $topo_module_name
        using LinearAlgebra
        using StaticArrays
        args = :legacy_topography_args
        include(joinpath(Main.REPO_ROOT, "src", "core", "interfaces", "reference_system.jl"))
    end))
    topo_sandbox = getfield(Main, topo_module_name)
    lla_topo = topo_sandbox.rtolatlong(rp_topo, planet_topo, true)
    @test isapprox(lla_topo[1], norm(rp_topo) - 9.5; atol=1e-12, rtol=0.0)

    q_branch2 = orbital_elements_to_lvlh_quaternion(
        1.186823891356144,
        0.13305160284932352,
        0.0,
        2.0245819323134224
    )
    q_branch3 = orbital_elements_to_lvlh_quaternion(
        0.0,
        0.1,
        0.0,
        0.10471975511965977
    )
    @test isapprox(norm(q_branch2), 1.0; atol=1e-12, rtol=0.0)
    @test isapprox(norm(q_branch3), 1.0; atol=1e-12, rtol=0.0)
end

@testset "Circular Orbit Robustness" begin
    ic_eq = InitialCondition(
        ra=EARTH.Rp_e + 500e3,
        rp=EARTH.Rp_e + 500e3,
        i=0.0,
        ω=0.0,
        Ω=0.0,
        ν=30.0
    )
    r, v = orbitalelemtorv(ic_eq, EARTH)
    oe = rvtoorbitalelement(SVector{3, Float64}(r), SVector{3, Float64}(v), EARTH)
    @test all(isfinite, oe)
    @test abs(oe[2]) < 1e-10

    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=500e3,
        i_deg=0.0,
        ω_deg=0.0,
        Ω_deg=0.0,
        ν_deg=30.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=300.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    df = run_case(args)
    @test nrow(df) > 10
end

@testset "Quaternion Norm Invariant" begin
    q0 = 2.3 * normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    w0 = SVector{3, Float64}(0.01, -0.015, 0.02)
    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=420e3,
        i_deg=40.0,
        ω_deg=10.0,
        Ω_deg=20.0,
        ν_deg=175.0,
        orientation_state=(q0, w0)
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-9,
            abstol_orbit=1e-9,
            reltol_quaternion=1e-9,
            abstol_quaternion=1e-9,
            dt_max_orbit=5.0
        )
    )

    u0 = build_initial_conditions(args)
    @test isapprox(norm(SVector{4, Float64}(u0.sc[1].q)), 1.0; atol=1e-12, rtol=0.0)

    df = run_case(args)
    qnorm = sqrt.(df.sc1_q_1.^2 .+ df.sc1_q_2.^2 .+ df.sc1_q_3.^2 .+ df.sc1_q_4.^2)
    @test maximum(abs.(qnorm .- 1.0)) < 1e-6
end

@testset "DynamicsRotational Core Contracts" begin
    rotational = SimulationModel.DynamicsRotational

    q = normalize(SVector{4, Float64}(0.15, -0.25, 0.35, 0.85))
    ω = SVector{3, Float64}(0.02, -0.03, 0.01)

    qdot_new = rotational.quaternion_derivative(ω, q)
    qdot_legacy = 0.5 * SVector{4, Float64}(SimulationModel.quat_mult(SVector{4, Float64}(ω..., 0.0), q))
    @test isapprox(qdot_new, qdot_legacy; atol=1e-12, rtol=1e-12)
    @test abs(dot(q, qdot_new)) < 1e-12

    J = SMatrix{3, 3, Float64}(2.0, 0.1, 0.0, 0.1, 3.0, 0.0, 0.0, 0.0, 4.0)
    τ = SVector{3, Float64}(0.12, -0.08, 0.05)
    ωdot_new = rotational.angular_acceleration(ω, J, τ; include_gyroscopic=true)
    ωdot_legacy = J \ (τ - cross(ω, J * ω))
    @test isapprox(ωdot_new, ωdot_legacy; atol=1e-12, rtol=1e-12)

    ωdot_no_gyro = rotational.angular_acceleration(ω, J, τ; include_gyroscopic=false)
    @test isapprox(ωdot_no_gyro, J \ τ; atol=1e-12, rtol=1e-12)

    @test rotational.body_torque(MVector{3, Float64}(τ)) == τ
    @test rotational.body_angular_velocity(MVector{3, Float64}(ω)) == ω
    @test rotational.combine_torques(τ, -τ) == SVector{3, Float64}(0.0, 0.0, 0.0)

    J_spherical = SMatrix{3, 3, Float64}(2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0)
    ω_any = SVector{3, Float64}(0.3, -0.25, 0.11)
    ωdot_torque_free = rotational.angular_acceleration(
        ω_any,
        J_spherical,
        SVector{3, Float64}(0.0, 0.0, 0.0);
        include_gyroscopic=true
    )
    @test isapprox(ωdot_torque_free, SVector{3, Float64}(0.0, 0.0, 0.0); atol=1e-12, rtol=0.0)
end

@testset "Torque-Free Rigid-Body Regression Cases" begin
    function make_torque_free_spacecraft(
        inertia_tensor::SMatrix{3, 3, Float64},
        q0::SVector{4, Float64},
        w0::SVector{3, Float64};
        id::Int64=1
    )
        root = Link{0}(root=true, m=500.0, ref_area=12.0, inertia=inertia_tensor)
        ic = InitialCondition(
            ra=EARTH.Rp_e + 500e3,
            rp=EARTH.Rp_e + 500e3,
            i=0.0,
            ω=0.0,
            Ω=0.0,
            ν=0.0,
            q=q0,
            ang_vel=w0
        )
        return SpacecraftModel(
            joints=Joint[],
            links=[root],
            root=root,
            prop_mass=0.0,
            inertia_tensor=inertia_tensor,
            initial_condition=ic,
            id=id
        )
    end

    function run_torque_free_case(
        spacecraft::SpacecraftModel;
        mission_time::Float64,
        dt_max_orbit::Float64=0.5
    )
        args = build_config(
            spacecraft=spacecraft,
            density_model=NoAtmosphereModel(),
            orientation_sim=true,
            mission_time=mission_time,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            keplerian=true,
            simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
            tolerances=IntegrationTolerances(
                reltol_orbit=1e-10,
                abstol_orbit=1e-10,
                reltol_quaternion=1e-10,
                abstol_quaternion=1e-10,
                reltol_angular_rate=1e-10,
                abstol_angular_rate=1e-10,
                dt_max_orbit=dt_max_orbit
            ),
            ephemerides_model=SimpleEphemeridesModel()
        )
        return run_simulation(args; return_solution=true)
    end

    withenv("SPACEAGORA_SOLVER_MODE" => "tsit5") do
        J_axis = SMatrix{3, 3, Float64}(8.0, 0.0, 0.0, 0.0, 12.0, 0.0, 0.0, 0.0, 20.0)
        q_axis0 = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
        ω_axis0 = SVector{3, Float64}(0.0, 0.0, 0.2)
        sol_axis = run_torque_free_case(
            make_torque_free_spacecraft(J_axis, q_axis0, ω_axis0);
            mission_time=40.0,
            dt_max_orbit=0.2
        )

        @test length(sol_axis.t) > 10
        for (t, state) in zip(sol_axis.t, sol_axis.u)
            q_state = SimulationModel.project_unit_quaternion(state.sc[1].q)
            ω_state = SVector{3, Float64}(state.sc[1].ω)
            q_expected = SVector{4, Float64}(0.0, 0.0, sin(0.5 * ω_axis0[3] * Float64(t)), cos(0.5 * ω_axis0[3] * Float64(t)))

            @test isapprox(ω_state, ω_axis0; atol=1e-9, rtol=1e-9)
            @test abs(dot(q_state, q_expected)) > 1.0 - 5e-8
        end

        J_free = SMatrix{3, 3, Float64}(8.0, 0.0, 0.0, 0.0, 11.0, 0.0, 0.0, 0.0, 15.0)
        q_free0 = normalize(SVector{4, Float64}(0.2, -0.3, 0.1, 0.92))
        ω_free0 = SVector{3, Float64}(0.19, -0.27, 0.31)
        sol_free = run_torque_free_case(
            make_torque_free_spacecraft(J_free, q_free0, ω_free0);
            mission_time=180.0
        )

        energies = Float64[]
        angular_momentum_norms = Float64[]
        quaternion_norm_errors = Float64[]
        for state in sol_free.u
            ω_state = SVector{3, Float64}(state.sc[1].ω)
            push!(energies, 0.5 * dot(ω_state, J_free * ω_state))
            push!(angular_momentum_norms, norm(J_free * ω_state))
            push!(quaternion_norm_errors, abs(norm(SVector{4, Float64}(state.sc[1].q)) - 1.0))
        end

        energy_drift = (maximum(energies) - minimum(energies)) / abs(first(energies))
        angular_momentum_drift = (maximum(angular_momentum_norms) - minimum(angular_momentum_norms)) / first(angular_momentum_norms)

        @test energy_drift < 1e-7
        @test angular_momentum_drift < 1e-7
        @test maximum(quaternion_norm_errors) < 1e-10
    end
end

@testset "Solver Tolerances Apply Quaternion Overrides" begin
    q0 = normalize(SVector{4, Float64}(0.1, -0.2, 0.3, 0.9))
    w0 = SVector{3, Float64}(0.01, -0.015, 0.02)
    sc = make_spacecraft(
        ra_alt_m=500e3,
        rp_alt_m=420e3,
        i_deg=40.0,
        ω_deg=10.0,
        Ω_deg=20.0,
        ν_deg=175.0,
        orientation_state=(q0, w0)
    )
    tols = IntegrationTolerances(
        reltol_orbit=1e-5,
        abstol_orbit=2e-6,
        reltol_quaternion=3e-7,
        abstol_quaternion=4e-8,
        reltol_mass=5e-7,
        abstol_mass=6e-8,
        reltol_heat_load=7e-7,
        abstol_heat_load=8e-8,
        reltol_angular_rate=9e-7,
        abstol_angular_rate=1e-7,
        dt_max_orbit=5.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=tols
    )

    u0 = build_initial_conditions(args)
    reltol_vec, abstol_vec = _build_solver_tolerances(u0, args)
    @test all(isapprox.(reltol_vec.sc[1].pos, tols.reltol_orbit; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec.sc[1].pos, tols.abstol_orbit; atol=0.0, rtol=0.0))
    @test all(isapprox.(reltol_vec.sc[1].q, tols.reltol_quaternion; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec.sc[1].q, tols.abstol_quaternion; atol=0.0, rtol=0.0))
    @test all(isapprox.(reltol_vec.sc[1].ω, tols.reltol_angular_rate; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec.sc[1].ω, tols.abstol_angular_rate; atol=0.0, rtol=0.0))
    @test isapprox(reltol_vec.sc[1].mass, tols.reltol_mass; atol=0.0, rtol=0.0)
    @test isapprox(abstol_vec.sc[1].mass, tols.abstol_mass; atol=0.0, rtol=0.0)
    @test all(isapprox.(reltol_vec.sc[1].heat_loads, tols.reltol_heat_load; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec.sc[1].heat_loads, tols.abstol_heat_load; atol=0.0, rtol=0.0))

    tols_no_orient = IntegrationTolerances(
        reltol_orbit=1e-5,
        abstol_orbit=2e-6,
        reltol_mass=0.0,
        abstol_mass=0.0,
        reltol_heat_load=0.0,
        abstol_heat_load=0.0,
        reltol_angular_rate=0.0,
        abstol_angular_rate=0.0,
        dt_max_orbit=5.0
    )
    args_no_orient = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=tols_no_orient
    )
    u0_no_orient = build_initial_conditions(args_no_orient)
    reltol_scalar, abstol_scalar = _build_solver_tolerances(u0_no_orient, args_no_orient)
    @test reltol_scalar == tols_no_orient.reltol_orbit
    @test abstol_scalar == tols_no_orient.abstol_orbit

    tols_no_orient_component = IntegrationTolerances(
        reltol_orbit=1e-5,
        abstol_orbit=2e-6,
        reltol_mass=3e-7,
        abstol_mass=4e-8,
        reltol_heat_load=5e-7,
        abstol_heat_load=6e-8,
        dt_max_orbit=5.0
    )
    args_no_orient_component = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=tols_no_orient_component
    )
    u0_no_orient_component = build_initial_conditions(args_no_orient_component)
    reltol_vec_no_orient, abstol_vec_no_orient = _build_solver_tolerances(u0_no_orient_component, args_no_orient_component)
    @test isapprox(reltol_vec_no_orient.sc[1].mass, tols_no_orient_component.reltol_mass; atol=0.0, rtol=0.0)
    @test isapprox(abstol_vec_no_orient.sc[1].mass, tols_no_orient_component.abstol_mass; atol=0.0, rtol=0.0)
    @test all(isapprox.(reltol_vec_no_orient.sc[1].heat_loads, tols_no_orient_component.reltol_heat_load; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_vec_no_orient.sc[1].heat_loads, tols_no_orient_component.abstol_heat_load; atol=0.0, rtol=0.0))

    tols_ω_no_orient = IntegrationTolerances(
        reltol_orbit=1e-5,
        abstol_orbit=2e-6,
        reltol_mass=0.0,
        abstol_mass=0.0,
        reltol_heat_load=0.0,
        abstol_heat_load=0.0,
        reltol_angular_rate=3e-7,
        dt_max_orbit=5.0
    )
    args_ω_no_orient = build_config(
        spacecraft=make_single_link_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=tols_ω_no_orient
    )
    u0_ω_no_orient = build_initial_conditions(args_ω_no_orient)
    reltol_ω_no_orient, abstol_ω_no_orient = _build_solver_tolerances(u0_ω_no_orient, args_ω_no_orient)
    @test reltol_ω_no_orient == tols_ω_no_orient.reltol_orbit
    @test abstol_ω_no_orient == tols_ω_no_orient.abstol_orbit
end
