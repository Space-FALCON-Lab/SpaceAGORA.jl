include(joinpath(@__DIR__, "..", "Simulation", "Code", "StudyConfig.jl"))
using .StudyConfig
import Pkg
Pkg.activate(LIBRARY_ROOT; io=devnull)
using Arrow, Test, TOML

function compare_formats(case)
    compact_dir = case_directory(DEFAULT_OUTPUT, case)
    relative = relpath(compact_dir, joinpath(DEFAULT_OUTPUT, "paper_plot_mode"))
    dense_dir = joinpath(@__DIR__, "results", "production_pipeline", "dense_v1", relative)
    compact_info = TOML.parsefile(joinpath(compact_dir, "run_info.toml"))
    dense_info = TOML.parsefile(joinpath(dense_dir, "run_info.toml"))
    compact, dense = Arrow.Table(joinpath(compact_dir, "trajectory.feather")),
        Arrow.Table(joinpath(dense_dir, "trajectory.feather"))
    @testset "Compact saving preserves the simulated trajectory and firing" begin
        @test length(compact.time) == 1001
        @test compact_info["schedule_sha256"] == dense_info["schedule_sha256"]
        @test compact_info["solver_accepted_steps"] == dense_info["solver_accepted_steps"]
        @test compact_info["solver_rejected_steps"] == dense_info["solver_rejected_steps"]
        @test propertynames(compact) == propertynames(dense)
        @test all(first(compact[name]) == first(dense[name]) &&
            last(compact[name]) == last(dense[name]) for name in propertynames(compact))
        for file in ("encounters.feather", "laser_firing.feather", "firing_range_audit.feather")
            old = Arrow.Table(joinpath(dense_dir, file))
            new = Arrow.Table(joinpath(compact_dir, file))
            @test all(isequal(old[name], new[name]) for name in propertynames(old) if name != :run_id)
        end
    end
end

for case in (StudyCase(850, 1, 50, :gve_sma), StudyCase(1050, 1, 50, :gve_inc),
    StudyCase(1050, 1, 300, :gve_sma))
    compare_formats(case)
end
