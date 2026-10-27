module BenchmarkSupportTests
using Test
const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(REPO, "benchmarks", "support", "artifacts.jl"))
include(joinpath(REPO, "benchmarks", "studies", "performance_runtime_analysis", "measurement", "routing.jl"))
using .BenchmarkArtifacts
using .BenchmarkRouting

@testset "artifact selection and missing-output behavior" begin
    mktempdir() do dir
        @test _latest_artifact_path_optional(dir, "summary", "quick", ".csv") === nothing
        @test_throws ErrorException _latest_artifact_path(dir, "summary", "quick", ".csv")
        # Names deliberately disagree with creation order; selection is by mtime.
        old = joinpath(dir, "summary_quick_z.csv")
        latest = joinpath(dir, "summary_quick_a.csv")
        write(old, "old")
        old_mtime = stat(old).mtime
        # Wait only for the filesystem timestamp to advance, with a bounded deadline.
        deadline = time() + 5.0
        while true
            write(latest, "latest")
            stat(latest).mtime > old_mtime && break
            time() >= deadline && error("filesystem timestamp did not advance")
            sleep(0.01)
        end
        write(joinpath(dir, "summary_other_new.csv"), "wrong profile")
        write(joinpath(dir, "report_quick_new.csv"), "wrong prefix")
        write(joinpath(dir, "summary_quick_new.txt"), "wrong suffix")
        @test _latest_artifact_path(dir, "summary", "quick", ".csv") == latest
        @test _latest_artifact_path_optional(dir, "summary", "quick", ".csv") == latest
        rm(latest)
        @test _latest_artifact_path(dir, "summary", "quick", ".csv") == old
        # Optional means no matching file, not an unreadable/missing directory.
        missing_dir = joinpath(dir, "missing")
        @test_throws Base.IOError _latest_artifact_path_optional(missing_dir, "summary", "quick", ".csv")
        @test_throws Base.IOError _latest_artifact_path(missing_dir, "summary", "quick", ".csv")
    end
end

@testset "route partition preserves order, safety and caller state" begin
    selected = [(id=1, route=:process, safe=false),
                (id=2, route=:threads, safe=false),
                (id=3, route=:threads, safe=true),
                (id=4, route=:none, safe=true),
                (id=5, route=:process, safe=true),
                (id=6, route=:custom, safe=true)]
    original = copy(selected)
    environment = copy(ENV)
    decisions, safety_checks = Int[], Int[]
    result = partition_sweep_routes(selected;
        choose_route=case -> (push!(decisions, case.id); case.route),
        thread_safe=case -> (push!(safety_checks, case.id); case.safe))
    @test result.chosen_routes == [:process, :none, :threads, :none, :process, :custom]
    @test result.process_tasks == [(1, selected[1]), (5, selected[5])]
    @test eltype(result.process_tasks) == Tuple{Int, eltype(selected)}
    @test result.thread_indices == [3]
    @test result.serial_indices == [2, 4, 6]
    @test decisions == collect(1:6)
    @test safety_checks == [2, 3]
    @test selected == original
    @test ENV == environment
    empty = partition_sweep_routes(selected[1:0];
        choose_route=_ -> error("empty selection routed"),
        thread_safe=_ -> error("empty selection checked"))
    @test isempty(empty.chosen_routes)
    @test isempty(empty.process_tasks)
    @test isempty(empty.thread_indices)
    @test isempty(empty.serial_indices)
    @test_throws ErrorException partition_sweep_routes(selected;
        choose_route=_ -> error("route failure"), thread_safe=_ -> true)
    @test_throws ErrorException partition_sweep_routes(selected;
        choose_route=_ -> :threads, thread_safe=_ -> error("safety failure"))
    @test ENV == environment
end
end
