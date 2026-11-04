include("production_pipeline_tests.jl")

function helper_extension_fixtures(campaign="helpers400500")
    if campaign == "helpers400500"
        counts, altitudes = (400, 500), (850, 950, 1000, 1050, 1150)
        previous = ("main",)
        results = joinpath(@__DIR__, "results", "helper_extension")
    elseif campaign == "helpers600800"
        counts, altitudes = (600, 800), (850, 900, 950, 1000, 1050, 1100, 1150)
        previous = ("main", "helpers400500", "targets9001100km")
        results = joinpath(@__DIR__, "results", "helper_extension600800")
    elseif campaign == "helpers700"
        counts, altitudes = (700,), (850, 900, 950, 1000, 1050, 1100, 1150)
        previous = ("main", "helpers400500", "targets9001100km", "helpers600800")
        results = joinpath(@__DIR__, "results", "helper_extension700")
    else
        throw(ArgumentError("Use helpers400500, helpers600800 or helpers700"))
    end
    expected_cases = length(altitudes) * 3 * length(counts) * 2
    original = vcat([grid_cases(name) for name in previous]...)
    extension = grid_cases(campaign)
    fingerprints = Dict{String,String}()

    @testset "$campaign grid and unchanged previous settings" begin
        @test length(grid_cases()) == 210
        @test length(grid_cases("helpers400500")) == 60
        @test length(grid_cases("targets9001100km")) == 108
        @test length(grid_cases("helpers600800")) == 84
        @test length(grid_cases("helpers700")) == 42
        @test length(extension) == expected_cases
        @test campaign_settings(campaign).helper_counts == counts
        @test campaign_settings(campaign).target_altitudes_km == altitudes
        @test all(c.orbits == 100 && c.power_W == 10_000 && c.helper_altitude_km == 1000 for c in extension)
        @test length(unique(case_signature.(extension))) == expected_cases
        @test length(unique(case_directory.(Ref(DEFAULT_OUTPUT), extension))) == expected_cases
        @test isempty(intersect(case_signature.(original), case_signature.(extension)))
        @test isempty(intersect(case_directory.(Ref(DEFAULT_OUTPUT), original),
            case_directory.(Ref(DEFAULT_OUTPUT), extension)))
        for altitude in altitudes, helpers in counts
            @test count(c -> c.altitude_km == altitude && c.helpers == helpers, extension) == 6
        end
        @test_throws ArgumentError grid_cases("unknown")
        @test_throws ArgumentError StudyCase(1050, 1, 401, :gve_sma)
        @test_throws ArgumentError StudyCase(1050, 1, 701, :gve_sma)
        for case in original
            directory = case_directory(DEFAULT_OUTPUT, case)
            info = TOML.parsefile(joinpath(directory, "run_info.toml"))
            @test info["status"] == "verified"
            @test info["run_id"] == case_signature(case)
            @test info["settings"] == case_settings(case)
            for name in ("run_info.toml", "verification.toml")
                path = joinpath(directory, name)
                fingerprints[relpath(path, DEFAULT_OUTPUT)] = filehash(path)
            end
        end
        for case in extension
            command = collect(BatchCompatibility.simulation_command(DEFAULT_OUTPUT, case))
            @test command[end-7:end] == string.([DEFAULT_OUTPUT, case.altitude_km,
                case.inclination_deg, case.helpers, case.policy, 100.0, 10_000.0, 1000.0])
        end
        for case in extension[1:2:end]
            args = study_configuration(case)
            state = SE.build_initial_conditions(args)
            @test length(state.sc) == case.helpers + 1
            @test args.scheduling_policy_model.target_idx == case.helpers + 1
            @test args.interlink_diagnostics
            @test args.integration_tolerances.dt_max_orbit == 10.0
            @test length(trajectory_times(case)) == 1001
            @test isapprox(args.dynamics_model.spacecraft[end].initial_condition.ν, deg2rad(1.0))
            @test all(isapprox(norm(state.sc[k].pos), EARTH_RADIUS_M + 1000e3; atol=1e-8, rtol=0)
                for k in 1:case.helpers)
            @test all(isapprox(args.dynamics_model.spacecraft[k].initial_condition.ν,
                deg2rad(360.0 * (k - 1) / case.helpers)) for k in 1:case.helpers)
            @test isapprox(norm(state.sc[end].pos), EARTH_RADIUS_M + case.altitude_km * 1000; atol=1e-8, rtol=0)
            @test OracleStudy.check_initial_geometry(state) > 1e-3
        end
    end

    @testset "$campaign logs, status, locks and stop controls are separate" begin
        mktempdir() do root
            saved = Dict{String,String}()
            for name in (previous..., "helpers9001100km")
                selection = campaign_settings(name)
                batch = joinpath(root, selection.batch_directory)
                mkpath(joinpath(batch, "lock"))
                write(joinpath(root, selection.progress_file), "$name log\n")
                atomic_toml(joinpath(batch, "status.toml"), Dict("status" => "complete"))
                atomic_toml(joinpath(batch, "lock", "owner.toml"), Dict("pid" => getpid()))
                write(joinpath(batch, "STOP_AFTER_CASE"), "$name stop flag\n")
                for path in (joinpath(root, selection.progress_file), joinpath(batch, "status.toml"),
                    joinpath(batch, "lock", "owner.toml"), joinpath(batch, "STOP_AFTER_CASE"))
                    saved[path] = filehash(path)
                end
            end
            BatchCompatibility.run_grid(root; campaign, available_space=_ -> 0)
            selection = campaign_settings(campaign)
            batch = joinpath(root, selection.batch_directory)
            progress = joinpath(root, selection.progress_file)
            status = TOML.parsefile(joinpath(batch, "status.toml"))
            @test status["total_cases"] == expected_cases
            @test status["helper_counts"] == collect(counts)
            @test status["target_altitudes_km"] == collect(altitudes)
            @test status["helper_altitudes_km"] == [1000]
            @test status["status"] == "paused_low_disk"
            @test occursin("/N$(first(counts))/", status["next_case"])
            @test !haskey(status, "child_pid")
            @test !ispath(joinpath(root, "paper_plot_mode"))
            @test !ispath(joinpath(batch, "lock"))
            @test occursin("BATCH_PAUSED_LOW_DISK", read(progress, String))
            scenario = "N$(last(counts)) fixture"
            BatchCompatibility.tracked_child(
                `$(Base.julia_cmd()) --startup-file=no -e 'println("extension fixture")'`,
                joinpath(batch, "fixture.log"), joinpath(batch, "fixture.toml"),
                Dict{String,Any}(), root, scenario, "SIMULATION"; campaign)
            @test occursin("SIMULATION_STARTED | $scenario", read(progress, String))
            @test occursin("SIMULATION_FINISHED | $scenario | elapsed=", read(progress, String))
            script = joinpath(@__DIR__, "..", "Simulation", "Code", "launch_study.sh")
            run(`bash $script stop-after-case $root $campaign`)
            @test isfile(joinpath(batch, "STOP_AFTER_CASE"))
            BatchCompatibility.run_grid(root; campaign)
            @test TOML.parsefile(joinpath(batch, "status.toml"))["status"] == "stopped_after_case"
            @test occursin("stopped_after_case", read(`bash $script status $root $campaign`, String))
            @test occursin("complete", read(`bash $script status $root`, String))
            @test !success(pipeline(`bash $script status $root invalid`; stdout=devnull, stderr=devnull))
            @test all(filehash(path) == hash for (path, hash) in saved)
        end
    end

    for case in grid_cases("helpers9001100km")
        for name in ("run_info.toml", "verification.toml", "progress.toml")
            path = joinpath(case_directory(DEFAULT_OUTPUT, case), name)
            isfile(path) && (fingerprints[relpath(path, DEFAULT_OUTPUT)] = filehash(path))
        end
    end
    for name in (previous..., "helpers9001100km")
        selection = campaign_settings(name)
        for relative in (selection.progress_file, joinpath(selection.batch_directory, "status.toml"),
            joinpath(selection.batch_directory, "STOP_AFTER_CASE"))
            path = joinpath(DEFAULT_OUTPUT, relative)
            isfile(path) && (fingerprints[relative] = filehash(path))
        end
    end
    mkpath(results)
    atomic_toml(joinpath(results, "fixtures.toml"), Dict(
        "passed" => true, "new_cases" => expected_cases, "helpers" => collect(counts),
        "target_altitudes_km" => collect(altitudes), "preserved_files" => fingerprints,
        "maximum_case_disk_requirement_bytes" =>
            maximum(BatchCompatibility.case_disk_requirement, BatchCompatibility.grid_cases(campaign))))
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) <= 1 || error("Usage: julia extended_helper_tests.jl [helpers400500|helpers600800|helpers700]")
    helper_extension_fixtures(isempty(ARGS) ? "helpers400500" : only(ARGS))
end
