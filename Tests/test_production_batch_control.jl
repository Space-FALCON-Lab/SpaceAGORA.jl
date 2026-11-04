module ProductionBatchControlTests

include(joinpath(@__DIR__, "..", "Simulation", "Code", "run_grid.jl"))
using Test

root = joinpath(@__DIR__, "results", "production_pipeline", "batch_control_fixture")
batch = joinpath(root, "_batch")
mkpath(batch)
open(joinpath(batch, "STOP_AFTER_CASE"), "w") do io
    println(io, "Test fixture: do not launch a simulation")
end

@testset "Batch stop, lock and source compatibility" begin
    run_grid(root)
    status = TOML.parsefile(joinpath(batch, "status.toml"))
    @test status["status"] == "stopped_after_case"
    @test status["verified_cases"] == 0
    @test !ispath(joinpath(root, "paper_plot_mode"))
    @test !ispath(joinpath(batch, "lock"))
    rm(joinpath(batch, "STOP_AFTER_CASE"))
    run_grid(root; available_space=_ -> 0)
    paused = TOML.parsefile(joinpath(batch, "status.toml"))
    @test paused["status"] == "paused_low_disk"
    @test paused["verified_cases"] == 0
    @test !haskey(paused, "child_pid")
    @test !ispath(joinpath(root, "paper_plot_mode"))
    @test paused["required_free_bytes"] > 5 * 1024^3
    progress = read(joinpath(root, "simulation_progress.log"), String)
    @test occursin("BATCH_STARTED", progress)
    @test occursin("BATCH_PAUSED_LOW_DISK", progress)
    child_status = Dict{String,Any}()
    child_log = joinpath(batch, "child_test.log")
    child_status_path = joinpath(batch, "child_test.toml")
    elapsed = tracked_child(`$(Base.julia_cmd()) --startup-file=no -e 'println("fixture")'`,
        child_log, child_status_path, child_status, root, "fixture-success", "SIMULATION")
    @test elapsed > 0
    @test_throws ErrorException tracked_child(
        `$(Base.julia_cmd()) --startup-file=no -e 'exit(1)'`,
        child_log, child_status_path, child_status, root, "fixture-failure", "SIMULATION")
    progress = read(joinpath(root, "simulation_progress.log"), String)
    @test occursin("SIMULATION_STARTED | fixture-success", progress)
    @test occursin("SIMULATION_FINISHED | fixture-success | elapsed=", progress)
    @test occursin("SIMULATION_FAILED | fixture-failure | elapsed=", progress)
    @test !occursin("SIMULATION_FINISHED | fixture-failure", progress)
    write(joinpath(batch, "STOP_AFTER_CASE"), "Test fixture\n")
    lock = joinpath(batch, "lock")
    mkdir(lock)
    atomic_toml(joinpath(lock, "owner.toml"), Dict("pid" => getpid()))
    @test_throws ErrorException run_grid(root)
    @test isfile(joinpath(lock, "owner.toml"))
    rm(joinpath(lock, "owner.toml"))
    rm(lock)
    hashes = source_hashes()
    @test same_simulation_sources(hashes, copy(hashes))
    for name in ("Code/plot_orbital_elements.jl", "Code/plot_combined_orbital_elements.jl",
        "Code/plot_pooled_orbital_elements.jl", "Code/plot_new_fixture.jl")
        changed = copy(hashes)
        changed[name] = "different-postprocessing"
        @test batch_source_hashes(hashes) == batch_source_hashes(changed)
        @test same_simulation_sources(hashes, changed)
    end
    removed_plot = copy(hashes)
    delete!(removed_plot, "Code/plot_orbital_elements.jl")
    @test batch_source_hashes(hashes) == batch_source_hashes(removed_plot)
    for name in ("Code/StudyConfig.jl", "Code/StudyRunner.jl", "Code/VerifyOutput.jl",
        "Code/run_grid.jl", "Code/launch_study.sh")
        changed = copy(hashes)
        changed[name] = "different-production"
        @test batch_source_hashes(hashes) != batch_source_hashes(changed)
    end
    changed_config = copy(hashes)
    changed_config["Code/StudyConfig.jl"] = "different-config"
    @test !same_simulation_sources(hashes, changed_config)
    changed = copy(hashes)
    changed["Code/VerifyOutput.jl"] = "different-verifier"
    @test same_simulation_sources(hashes, changed)
    changed["Code/StudyRunner.jl"] = "different-simulation"
    @test !same_simulation_sources(hashes, changed)
    for case in (StudyCase(1050, 1, 300, :gve_sma),
        StudyCase(1050, 1, 50, :gve_inc), StudyCase(850, 1, 50, :gve_sma))
        info = TOML.parsefile(joinpath(case_directory(DEFAULT_OUTPUT, case), "run_info.toml"))
        @test info["settings"] == case_settings(case)
        @test info["run_id"] == case_signature(case)
    end
end
atomic_toml(joinpath(batch, "test_result.toml"), Dict("passed" => true))

end
