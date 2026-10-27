module BenchmarkLauncherTests

using Test

struct PipelineIncludeBoundary <: Exception
    path::String
end

@testset "Static-versus-parallel launcher routing" begin
    repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
    launcher = joinpath(repo, "benchmarks", "studies", "performance", "performance_static_vs_parallel.jl")
    owner = joinpath(repo, "benchmarks", "scripts", "performance_paper_pipeline.jl")
    probe = Module(gensym(:BenchmarkLauncherProbe))
    support = joinpath(repo, "benchmarks", "support", "artifacts.jl")
    # Load the pure shared helper, then stop before the pipeline loads GRAM.
    # Any unexpected include must still fail at the boundary assertion below.
    Core.eval(probe, :(function include(path::AbstractString)
        normalized = normpath(String(path))
        normalized == $support && return Base.include(@__MODULE__, normalized)
        throw($PipelineIncludeBoundary(normalized))
    end))
    original_dir = pwd()
    original_project = Base.active_project()

    mktempdir() do unrelated_dir
        cd(unrelated_dir) do
            failure = try
                Base.include(probe, launcher)
                nothing
            catch err
                err
            end
            @test failure isa LoadError
            if failure isa LoadError
                @test failure.error isa PipelineIncludeBoundary
                if failure.error isa PipelineIncludeBoundary
                    @test failure.error.path == owner
                    @test isfile(failure.error.path)
                end
            end
            @test Base.invokelatest(getfield, probe, :_WRAPPER_DIR) == dirname(launcher)
            @test isdefined(probe, :BenchmarkArtifacts)
            @test !isdefined(probe, :REPO_ROOT)
            @test !isdefined(probe, :Random)
            @test !isdefined(probe, :StaticVsParallelConfig)
            @test pwd() == realpath(unrelated_dir)
            @test Base.active_project() == original_project
        end
    end

    @test pwd() == original_dir
    @test Base.active_project() == original_project
end

end
