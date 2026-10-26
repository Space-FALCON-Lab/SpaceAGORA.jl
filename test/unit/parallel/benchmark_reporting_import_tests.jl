module BenchmarkReportingImportTests
using Test

@testset "report module imports without simulation or process setup" begin
    source = normpath(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies",
        "performance_runtime_analysis", "reporting", "RuntimeReports.jl"))
    # A fresh process detects transitive native imports even when the enclosing
    # unit shard has already loaded SpaceAGORA or plotting packages.
    script = raw"""
        using Test, Distributed
        environment = Dict(ENV)
        load_path = copy(LOAD_PATH)
        project = Base.active_project()
        directory = pwd()
        worker_ids = workers()
        include(ARGS[1])
        @test Dict(ENV) == environment
        @test LOAD_PATH == load_path
        @test Base.active_project() == project
        @test pwd() == directory
        @test workers() == worker_ids
        forbidden = Set(["SpaceAGORA", "GRAMSuite", "SPICE", "Plots"])
        @test isempty(intersect(forbidden, Set(id.name for id in keys(Base.loaded_modules))))
        @test !isdefined(RuntimeReports, :SimulationModel)
        @test !isdefined(RuntimeReports, :SimulationEngine)
        @test parentmodule(RuntimeReports.ProfileSpec) === RuntimeReports
    """
    command = `$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) -e $script $source`
    @test success(command)
end
end
