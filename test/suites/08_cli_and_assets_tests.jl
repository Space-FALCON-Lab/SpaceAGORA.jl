import GRAMSuite

@testset "CLI and Asset Surface" begin
    @testset "Asset report covers baseline and optional roots" begin
        report = SpaceAGORA.check_assets(repo_root=REPO_ROOT)
        names = Set(item.name for item in report.items)
        @test "no_gram_mode" in names
        @test "gram_root" in names
        @test isfile(joinpath(REPO_ROOT, "data", "assets_manifest.toml"))
        @test isfile(joinpath(REPO_ROOT, "scripts", "assets", "check_assets.jl"))
        @test isfile(joinpath(REPO_ROOT, "scripts", "assets", "show_asset_manifest.jl"))
        @test isfile(joinpath(REPO_ROOT, "scripts", "assets", "setup_open_assets.jl"))
        text = sprint(io -> SpaceAGORA.render_asset_report(report; io=io))
        @test occursin("no_gram_mode", text)
        @test occursin("gram_root", text)
    end

    @testset "CLI help and print-only dispatch" begin
        help_text = sprint(io -> @test SpaceAGORA.run_cli(["help"]; io=io, errio=io) == 0)
        @test occursin("spaceagora run", help_text)

        run_text = sprint(io -> @test SpaceAGORA.run_cli([
            "run",
            "--example=AGORA_Earth_NoGRAM.jl",
            "--output-dir=$(mktempdir())",
            "--smoke",
            "--print-only",
        ]; io=io, errio=io) == 0)
        @test occursin("SPACEAGORA_CLI_OUTPUT_DIR", run_text)
        @test occursin("AGORA_Earth_NoGRAM.jl", run_text)

        telemetry_text = sprint(io -> @test SpaceAGORA.run_cli([
            "telemetry",
            "quick",
            "--output-dir=$(mktempdir())",
            "--enforce=1",
            "--print-only",
        ]; io=io, errio=io) == 0)
        @test occursin("SPACEAGORA_TELEMETRY_PLOTS=0", telemetry_text)
        @test occursin("telemetry_orbit_accuracy_study.jl", telemetry_text)

        perf_text = sprint(io -> @test SpaceAGORA.run_cli([
            "benchmark",
            "runtime-analysis",
            "smoke",
            "--output-dir=$(mktempdir())",
            "--print-only",
        ]; io=io, errio=io) == 0)
        @test occursin("SPACEAGORA_PERF_OUTDIR", perf_text)
        @test occursin("performance_runtime_analysis.jl", perf_text)

        ladder_text = sprint(io -> @test SpaceAGORA.run_cli([
            "benchmark",
            "smart-parallel-ladder",
            "smoke",
            "--output-dir=$(mktempdir())",
            "--print-only",
        ]; io=io, errio=io) == 0)
        @test occursin("SPACEAGORA_SMART_LADDER_OUTDIR", ladder_text)
        @test occursin("performance_smart_parallel_ladder.jl", ladder_text)

        manifest_text = sprint(io -> @test SpaceAGORA.run_cli(["assets", "manifest"]; io=io, errio=io) == 0)
        @test occursin("SpaceAGORA asset manifest", manifest_text)
        @test occursin("data/assets_manifest.toml", read(joinpath(REPO_ROOT, "docs", "src", "assets.md"), String))

        setup_text = sprint(io -> @test SpaceAGORA.run_cli(["assets", "setup-open"]; io=io, errio=io) == 0)
        @test occursin("No downloads are required for baseline no-GRAM mode", setup_text)
    end

    @testset "Earth GRAM density channels" begin
        withenv(
            "SPACEAGORA_GRAM_STATIC_GRID" => "0",
            "SPACEAGORA_GRAM_OFFLINE_SURROGATE" => "off",
        ) do
            it = InitialTime(year=2026, month=3, day=20, hour=12, minute=0, second=0.0)
            kwargs = (
                planet_name="earth",
                initial_time=it,
                earth_daily_f10=120.0,
                earth_mean_f10=120.0,
                earth_ap=10.0,
            )
            low_model = GRAMAtmosphereModel(; kwargs..., gram_density_channel=:low)
            nominal_model = GRAMAtmosphereModel(; kwargs..., gram_density_channel=:nominal)
            high_model = GRAMAtmosphereModel(; kwargs..., gram_density_channel=:high)

            h_m = 200.0e3
            lat_rad = deg2rad(0.0)
            lon_rad = deg2rad(0.0)
            rho_low, _, _ = GRAMSuite.point_density_state(low_model.core, h_m, lat_rad, lon_rad, 0.0, false)
            rho_nominal, _, _ = GRAMSuite.point_density_state(nominal_model.core, h_m, lat_rad, lon_rad, 0.0, false)
            rho_high, _, _ = GRAMSuite.point_density_state(high_model.core, h_m, lat_rad, lon_rad, 0.0, false)

            @test rho_low <= rho_nominal <= rho_high
        end
    end

    @testset "VLEO drag trade smoke" begin
        study_path = joinpath(REPO_ROOT, "benchmarks", "studies", "vleo_drag_trade.jl")
        if !isdefined(@__MODULE__, :run_vleo_drag_trade)
            include(study_path)
        end

        route_state = SpaceAGORA.OuterRouteState()
        features = _case_outer_features(150.0, 350.0)
        @test _select_case_outer_route(:auto, route_state, features) == :process

        mktempdir() do tmp
            result = run_vleo_drag_trade(smoke=true, out_dir=tmp, generate_plots=true, outer_route=:none)
            @test isfile(result.cases_csv_path)
            @test isfile(result.summary_csv_path)
            @test isfile(result.report_path)
            @test isfile(result.outer_route_state_path)
            @test length(result.plot_paths) == 4
            @test all(isfile, result.plot_paths)

            cases_df = CSV.read(result.cases_csv_path, DataFrame)
            @test nrow(cases_df) == 3
            @test Set(String.(cases_df.channel)) == Set(["low", "nominal", "high"])
            @test Set(String.(cases_df.outer_route)) == Set(["none"])
            @test all(Float64.(cases_df.drag_impulse_ns) .> 0.0)
            @test all(Float64.(cases_df.required_reboost_dv_mps) .> 0.0)

            report_text = read(result.report_path, String)
            @test occursin("AerodynamicCoefficientfM()", report_text)
            @test occursin("2026-03-20 12:00:00 UTC", report_text)
            @test occursin("Starlink-like", report_text)
        end
    end
end
