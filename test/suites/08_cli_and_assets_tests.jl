@testset "CLI and Asset Surface" begin
    @testset "Asset report covers baseline and optional roots" begin
        report = SpaceAGORA.check_assets(repo_root=REPO_ROOT)
        names = Set(item.name for item in report.items)
        @test "no_gram_mode" in names
        @test "gram_root" in names
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
    end
end
