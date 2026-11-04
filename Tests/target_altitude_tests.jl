include("production_pipeline_tests.jl")

const ALTITUDE_RESULTS = joinpath(@__DIR__, "results", "target_altitudes")
const ALTITUDE_CAMPAIGN = "targets9001100km"
const ALTITUDE_COUNTS = (1, 50, 100, 150, 200, 250, 300, 400, 500)
const ORIGINAL_FINGERPRINTS = Dict{String,String}()

@testset "900/1100 km targets with 1000 km helpers and unchanged original cases" begin
    original = vcat(grid_cases(), grid_cases("helpers400500"))
    cases = grid_cases(ALTITUDE_CAMPAIGN)
    @test length(grid_cases()) == 210
    @test length(grid_cases("helpers400500")) == 60
    @test length(cases) == 108
    @test length(unique(case_signature.(cases))) == 108
    @test length(unique(case_directory.(Ref(DEFAULT_OUTPUT), cases))) == 108
    @test isempty(intersect(case_signature.(original), case_signature.(cases)))
    @test isempty(intersect(case_directory.(Ref(DEFAULT_OUTPUT), original),
        case_directory.(Ref(DEFAULT_OUTPUT), cases)))
    @test all(c.orbits == 100 && c.power_W == 10_000 for c in cases)
    @test all(c.helper_altitude_km == 1000 for c in original)
    @test all(c.helper_altitude_km == 1000 for c in cases)
    @test length(grid_cases("helpers9001100km")) == 540
    @test isempty(intersect(case_directory.(Ref(DEFAULT_OUTPUT), grid_cases("helpers9001100km")),
        case_directory.(Ref(DEFAULT_OUTPUT), cases)))
    for altitude in (900, 1100), helpers in ALTITUDE_COUNTS
        @test count(c -> c.altitude_km == altitude && c.helpers == helpers, cases) == 6
    end
    for altitude in (800, 1001, 1200, NaN, Inf)
        @test_throws ArgumentError StudyCase(altitude, 1, 50, :gve_sma)
    end
    for case in original
        directory = case_directory(DEFAULT_OUTPUT, case)
        info = TOML.parsefile(joinpath(directory, "run_info.toml"))
        @test info["status"] == "verified"
        @test info["run_id"] == case_signature(case)
        @test info["settings"] == case_settings(case)
        for name in ("run_info.toml", "verification.toml")
            path = joinpath(directory, name)
            ORIGINAL_FINGERPRINTS[relpath(path, DEFAULT_OUTPUT)] = filehash(path)
        end
    end
    for case in cases
        settings = case_settings(case)
        @test settings["target_altitude_km"] == case.altitude_km
        @test settings["helper_altitude_km"] == 1000
        command = collect(BatchCompatibility.simulation_command(DEFAULT_OUTPUT, case))
        @test command[end-7:end] == string.([DEFAULT_OUTPUT, case.altitude_km,
            case.inclination_deg, case.helpers, case.policy, 100.0, 10_000.0, 1000.0])
    end
    for case in cases[1:2:end]
        args = study_configuration(case)
        state = SE.build_initial_conditions(args)
        helper_radius = EARTH_RADIUS_M + 1000e3
        target_radius = EARTH_RADIUS_M + case.altitude_km * 1000
        @test length(state.sc) == case.helpers + 1
        @test all(isapprox(norm(state.sc[k].pos), helper_radius; atol=1e-8, rtol=0)
            for k in 1:case.helpers)
        @test all(isapprox(norm(state.sc[k].vel), sqrt(MU / helper_radius); atol=1e-10, rtol=0)
            for k in 1:case.helpers)
        @test isapprox(norm(state.sc[end].pos), target_radius; atol=1e-8, rtol=0)
        @test isapprox(norm(state.sc[end].vel), sqrt(MU / target_radius); atol=1e-10, rtol=0)
        @test mission_duration(case) ≈ 100 * 2pi * sqrt(norm(state.sc[end].pos)^3 / MU)
        @test args.scheduling_policy_model.target_idx == case.helpers + 1
        @test args.interlink_diagnostics
        @test args.integration_tolerances.dt_max_orbit == 10.0
        @test length(trajectory_times(case)) == 1001
        @test isapprox(args.dynamics_model.spacecraft[end].initial_condition.ν, deg2rad(1.0))
        @test OracleStudy.check_initial_geometry(state) > 1e-3
    end
end

@testset "Altitude campaign has isolated logs, locks and stop controls" begin
    mktempdir() do root
        saved = Dict{String,String}()
        for campaign in ("main", "helpers400500", "helpers9001100km")
            selection = campaign_settings(campaign)
            batch = joinpath(root, selection.batch_directory)
            mkpath(joinpath(batch, "lock"))
            write(joinpath(root, selection.progress_file), "$campaign log\n")
            atomic_toml(joinpath(batch, "status.toml"), Dict("status" => "complete"))
            atomic_toml(joinpath(batch, "lock", "owner.toml"), Dict("pid" => getpid()))
            write(joinpath(batch, "STOP_AFTER_CASE"), "$campaign stop flag\n")
            for path in (joinpath(root, selection.progress_file), joinpath(batch, "status.toml"),
                joinpath(batch, "lock", "owner.toml"), joinpath(batch, "STOP_AFTER_CASE"))
                saved[path] = filehash(path)
            end
        end
        BatchCompatibility.run_grid(root; campaign=ALTITUDE_CAMPAIGN, available_space=_ -> 0)
        selection = campaign_settings(ALTITUDE_CAMPAIGN)
        batch = joinpath(root, selection.batch_directory)
        progress = joinpath(root, selection.progress_file)
        status = TOML.parsefile(joinpath(batch, "status.toml"))
        @test status["total_cases"] == 108
        @test status["helper_counts"] == collect(ALTITUDE_COUNTS)
        @test status["helper_altitudes_km"] == [1000]
        @test status["target_altitudes_km"] == [900, 1100]
        @test status["status"] == "paused_low_disk"
        @test occursin("/h1000km_t900km/", status["next_case"])
        @test !haskey(status, "child_pid")
        @test !ispath(joinpath(root, "paper_plot_mode"))
        @test !ispath(joinpath(batch, "lock"))
        @test occursin("BATCH_PAUSED_LOW_DISK", read(progress, String))
        BatchCompatibility.tracked_child(
            `$(Base.julia_cmd()) --startup-file=no -e 'println("altitude fixture")'`,
            joinpath(batch, "fixture.log"), joinpath(batch, "fixture.toml"),
            Dict{String,Any}(), root, "h1000km_t1100km N500 fixture", "SIMULATION"; campaign=ALTITUDE_CAMPAIGN)
        @test occursin("SIMULATION_STARTED | h1000km_t1100km N500 fixture", read(progress, String))
        @test occursin("SIMULATION_FINISHED | h1000km_t1100km N500 fixture | elapsed=", read(progress, String))
        script = joinpath(@__DIR__, "..", "Simulation", "Code", "launch_study.sh")
        run(`bash $script stop-after-case $root $ALTITUDE_CAMPAIGN`)
        @test isfile(joinpath(batch, "STOP_AFTER_CASE"))
        BatchCompatibility.run_grid(root; campaign=ALTITUDE_CAMPAIGN)
        @test TOML.parsefile(joinpath(batch, "status.toml"))["status"] == "stopped_after_case"
        @test occursin("stopped_after_case", read(`bash $script status $root $ALTITUDE_CAMPAIGN`, String))
        @test all(filehash(path) == hash for (path, hash) in saved)
    end
end

for case in grid_cases("helpers9001100km")
    for name in ("run_info.toml", "verification.toml", "progress.toml")
        path = joinpath(case_directory(DEFAULT_OUTPUT, case), name)
        isfile(path) && (ORIGINAL_FINGERPRINTS[relpath(path, DEFAULT_OUTPUT)] = filehash(path))
    end
end
for campaign in ("main", "helpers400500", "helpers9001100km")
    selection = campaign_settings(campaign)
    for relative in (selection.progress_file, joinpath(selection.batch_directory, "status.toml"))
        ORIGINAL_FINGERPRINTS[relative] = filehash(joinpath(DEFAULT_OUTPUT, relative))
    end
end
stopped_batch = joinpath(DEFAULT_OUTPUT, campaign_settings("helpers9001100km").batch_directory)
@test TOML.parsefile(joinpath(stopped_batch, "status.toml"))["status"] == "stopped_by_user"
@test isfile(joinpath(stopped_batch, "STOP_AFTER_CASE"))
ORIGINAL_FINGERPRINTS[relpath(joinpath(stopped_batch, "STOP_AFTER_CASE"), DEFAULT_OUTPUT)] =
    filehash(joinpath(stopped_batch, "STOP_AFTER_CASE"))
mkpath(ALTITUDE_RESULTS)
atomic_toml(joinpath(ALTITUDE_RESULTS, "fixtures.toml"), Dict(
    "passed" => true, "new_cases" => 108, "helper_altitudes_km" => [1000],
    "target_altitudes_km" => [900, 1100],
    "helper_counts" => collect(ALTITUDE_COUNTS), "preserved_files" => ORIGINAL_FINGERPRINTS,
    "maximum_case_disk_requirement_bytes" =>
        maximum(BatchCompatibility.case_disk_requirement, BatchCompatibility.grid_cases(ALTITUDE_CAMPAIGN))))
