using Test, TOML, SpaceAGORA, SpaceAGORAEDG
const CLI = SpaceAGORA.SpaceAGORACLI
const environment = dirname(Base.active_project())
@testset "EDG CLI retains the explicit companion environment" begin
    withenv("SPACEAGORA_EDG_ENV" => nothing) do
        @test CLI._edg_project() == environment
        printed = sprint(io -> @test run_cli(["run", "--example=AGORA_Odyssey_Energy_Depletion_Control_Test", "--print-only"]; io, errio=io) == 0)
        @test occursin("project=$environment\n", printed)
    end
    mktempdir() do temporary
        invalid = joinpath(temporary, "invalid")
        mkpath(invalid)
        open(joinpath(invalid, "Project.toml"), "w") do io
            TOML.print(io, Dict("deps" => Dict(
                "SpaceAGORA" => "afbfb69f-5c0b-4832-b760-43725dff8540",
                "SpaceAGORAEDG" => "00000000-0000-0000-0000-000000000000")))
        end
        @test_throws ArgumentError CLI._edg_project(invalid)
        @test_throws ArgumentError CLI._edg_project(joinpath(temporary, "missing"))
        withenv("SPACEAGORA_EDG_ENV" => environment) do
            @test CLI._edg_project() == environment
        end
        probe = joinpath(temporary, "EDG child with spaces.jl")
        write(probe, """
            using SpaceAGORAEDG
            include($(repr(joinpath(CLI.REPO_ROOT, "examples", "common.jl"))))
            @assert SpaceAGORA.EDGServices.available()
            # The EDG environment lacks the example packages: common.jl continues in the examples environment.
            @assert Base.active_project() == $(repr(joinpath(CLI.REPO_ROOT, "examples", "Project.toml")))
            @assert abspath(PROGRAM_FILE) == @__FILE__
            @assert ARGS == ["argument with spaces"]
            println("EDG_CHILD_OK")
            """)
        command = CLI._child_command(probe, ["argument with spaces"]; project=environment)
        @test count(==("--project=$environment"), command.exec) == 1
        @test command.exec[end-1:end] == [probe, "argument with spaces"]
        output = read(command, String)
        @test occursin("EDG_CHILD_OK", output)
    end
end

