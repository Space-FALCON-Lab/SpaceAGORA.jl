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

        telemetry_smoke_text = sprint(io -> @test SpaceAGORA.run_cli([
            "telemetry",
            "smoke",
            "--output-dir=$(mktempdir())",
            "--print-only",
        ]; io=io, errio=io) == 0)
        @test occursin("telemetry_orbit_accuracy_study.jl", telemetry_smoke_text)
        @test occursin(" quick", telemetry_smoke_text)

        telemetry_smoke_profile_text = sprint(io -> @test SpaceAGORA.run_cli([
            "telemetry",
            "--profile=smoke",
            "--output-dir=$(mktempdir())",
            "--print-only",
        ]; io=io, errio=io) == 0)
        @test occursin("telemetry_orbit_accuracy_study.jl", telemetry_smoke_profile_text)
        @test occursin(" quick", telemetry_smoke_profile_text)

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

        ladder_pass_through_text = sprint(io -> @test SpaceAGORA.run_cli([
            "benchmark",
            "smart-parallel-ladder",
            "full",
            "--passes=3",
            "--output-dir=$(mktempdir())",
            "--print-only",
        ]; io=io, errio=io) == 0)
        @test occursin(" full --passes=3", ladder_pass_through_text)
        @test !occursin("full--passes=3", ladder_pass_through_text)

        manifest_text = sprint(io -> @test SpaceAGORA.run_cli(["assets", "manifest"]; io=io, errio=io) == 0)
        @test occursin("SpaceAGORA asset manifest", manifest_text)
        @test occursin("data/assets_manifest.toml", read(joinpath(REPO_ROOT, "docs", "src", "assets.md"), String))

        setup_text = sprint(io -> @test SpaceAGORA.run_cli(["assets", "setup-open"]; io=io, errio=io) == 0)
        @test occursin("No downloads are required for baseline no-GRAM mode", setup_text)
    end

    @testset "CLI children run under the repository project" begin
        # The child launched by `run` must find SpaceAGORA through the project the
        # CLI hands it and nothing else. The probe loads the package directly (no
        # examples/common.jl re-activation), and the child's load path is reduced
        # to its own project so a globally installed SpaceAGORA cannot satisfy the
        # import; both would have hidden the original defect (a child launched with
        # a `.AGORA` project directory that does not exist).
        probe_dir = mktempdir()
        probe = joinpath(probe_dir, "cli_project_probe.jl")
        write(probe, """
            using SpaceAGORA
            println("PROBE_PROJECT=", Base.active_project())
            println("PROBE_PACKAGE=", pathof(SpaceAGORA))
            """)
        expected_project = joinpath(REPO_ROOT, "Project.toml")
        expected_package = joinpath(REPO_ROOT, "src", "SpaceAGORA.jl")

        log_path = joinpath(probe_dir, "child.log")
        code = withenv("JULIA_LOAD_PATH" => "@", "JULIA_PROJECT" => nothing) do
            open(log_path, "w") do log
                SpaceAGORA.run_cli(["run", "--example=$(probe)"]; io=log, errio=log)
            end
        end
        child_output = read(log_path, String)
        @test code == 0
        @test occursin("PROBE_PROJECT=$(expected_project)", child_output)
        @test occursin("PROBE_PACKAGE=$(expected_package)", child_output)

        # The command the CLI builds: compare its actual arguments, not a parse of
        # the printed line. The `--project=` argument must be the repository
        # root, and a script path containing spaces must survive as one argument.
        CLI = SpaceAGORA.SpaceAGORACLI
        same_dir(a, b) = rstrip(normpath(a), '/') == rstrip(normpath(b), '/')
        spaced_dir = joinpath(mktempdir(), "path with spaces")
        mkpath(spaced_dir)
        spaced_script = joinpath(spaced_dir, "spaced probe.jl")
        cp(probe, spaced_script)
        for (script, script_args) in (
            (probe, String[]),
            (spaced_script, String[]),
            (CLI.TELEMETRY_LAUNCHER, ["quick"]),
            (CLI.PERF_RUNTIME_LAUNCHER, ["smoke"]),
            (CLI.SMART_LADDER_LAUNCHER, ["full", "--passes=3"]),
        )
            argv = collect(CLI._child_command(script, script_args).exec)
            project_args = filter(a -> startswith(a, "--project="), argv)
            @test length(project_args) == 1
            @test same_dir(chopprefix(only(project_args), "--project="), REPO_ROOT)
            @test count(==(script), argv) == 1            # one intact argument, spaces included
            @test argv[end - length(script_args):end] == vcat([script], script_args)
            @test !any(a -> occursin(".AGORA", a), argv)
        end

        # The printed launcher must agree with what launches: the `cmd=` line is
        # Julia's rendering of the very command object the CLI would run, and the
        # `project=` line names the same directory.
        for (args, script, script_args) in (
            (["run", "--example=$(spaced_script)", "--print-only"], spaced_script, String[]),
            (["telemetry", "quick", "--output-dir=$(mktempdir())", "--print-only"], CLI.TELEMETRY_LAUNCHER, ["quick"]),
            (["benchmark", "runtime-analysis", "smoke", "--output-dir=$(mktempdir())", "--print-only"], CLI.PERF_RUNTIME_LAUNCHER, ["smoke"]),
            (["benchmark", "smart-parallel-ladder", "smoke", "--output-dir=$(mktempdir())", "--print-only"], CLI.SMART_LADDER_LAUNCHER, ["smoke"]),
        )
            printed = sprint(io -> @test SpaceAGORA.run_cli(args; io=io, errio=io) == 0)
            expected_cmd = CLI._child_command(script, script_args)
            @test occursin("cmd=$(expected_cmd)\n", printed)
            @test occursin("project=$(CLI.CHILD_PROJECT)\n", printed)
            @test same_dir(CLI.CHILD_PROJECT, REPO_ROOT)
            @test !occursin(".AGORA", printed)
        end

        # And the spaced path really launches: the same probe, from a directory
        # whose name contains a space, through the real child process.
        spaced_log = joinpath(spaced_dir, "child.log")
        spaced_code = withenv("JULIA_LOAD_PATH" => "@", "JULIA_PROJECT" => nothing) do
            open(spaced_log, "w") do log
                SpaceAGORA.run_cli(["run", "--example=$(spaced_script)"]; io=log, errio=log)
            end
        end
        spaced_output = read(spaced_log, String)
        @test spaced_code == 0
        @test occursin("PROBE_PROJECT=$(expected_project)", spaced_output)
        @test occursin("PROBE_PACKAGE=$(expected_package)", spaced_output)
    end
end
