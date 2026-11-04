using Test, TOML, SpaceAGORA, SpaceAGORAEDG

# Include the complete smoke driver without selecting any native examples.
# The tests below then call its actual run_example path with owned probe scripts.
module ExampleSmokeHarness
    using Test
    withenv("SPACEAGORA_EXAMPLE_FILTER" => "__edg_smoke_environment_regression_59fb7e3b__") do
        include(joinpath(@__DIR__, "..", "smoke", "ci_examples_suite_smoke.jl"))
        @test Base.invokelatest(() -> isempty(list_examples()))
    end
end

@testset "Example smoke preserves the paired EDG project" begin
    harness = ExampleSmokeHarness
    environment = dirname(Base.active_project())
    native_example = joinpath(harness.EXAMPLES_DIR, "AGORA_Odyssey_Energy_Depletion_Control_Test.jl")
    withenv("SPACEAGORA_EXAMPLE_FILTER" => "") do
        @test native_example in harness.list_examples()
        # The separate surrogate workflow remains responsible for this example.
        @test !(joinpath(harness.EXAMPLES_DIR, "odyssey_surrogate.jl") in harness.list_examples())
    end
    withenv("SPACEAGORA_EDG_ENV" => nothing) do
        @test harness.example_project(native_example) == environment
    end
    mktempdir() do temporary
        directory = joinpath(temporary, "example probes with spaces")
        mkpath(directory)
        probe = joinpath(directory, "AGORA_Odyssey_Energy_Depletion_Control_Test.jl")
        write(probe, """
            using SpaceAGORAEDG
            @assert Base.active_project() == $(repr(joinpath(environment, "Project.toml")))
            @assert LOAD_PATH == ["@", "@stdlib"]
            @assert abspath(PROGRAM_FILE) == @__FILE__
            @assert ENV["SPACEAGORA_EXAMPLE_SMOKE"] == "1"
            @assert ENV["SPACEAGORA_EXAMPLE_SMOKE_RESULTS"] == "1"
            @assert ENV["SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME"] == "120.0"
            include($(repr(joinpath(harness.EXAMPLES_DIR, "common.jl"))))
            @assert SpaceAGORA.EDGServices.available()
            @assert Base.active_project() == $(repr(joinpath(harness.REPO_ROOT, "Project.toml")))
            @assert LOAD_PATH == ["@", "@stdlib"]
            println("EDG_EXAMPLE_SMOKE_CHILD_OK")
            """)
        withenv("SPACEAGORA_EDG_ENV" => environment) do
            command = harness.example_command(probe)
            @test count(==("--project=$environment"), command.exec) == 1
            @test command.exec[end] == probe
            ok, has_nan, output = harness.run_example(probe)
            ok || println(output)
            @test ok
            @test !has_nan
            @test occursin("EDG_EXAMPLE_SMOKE_CHILD_OK", output)

            # Selecting an EDG environment must not reroute ordinary or HYPR examples.
            ordinary = joinpath(directory, "Ordinary_Project_Probe.jl")
            write(ordinary, """
                @assert Base.active_project() == $(repr(joinpath(harness.REPO_ROOT, "Project.toml")))
                @assert LOAD_PATH == ["@", "@stdlib"]
                println("ORDINARY_EXAMPLE_SMOKE_CHILD_OK")
                """)
            ok, has_nan, output = harness.run_example(ordinary)
            @test ok
            @test !has_nan
            @test occursin("ORDINARY_EXAMPLE_SMOKE_CHILD_OK", output)
            @test harness.example_project(joinpath(directory, "Earth_RPO_CubeSat_MPC.jl")) ==
                joinpath(harness.REPO_ROOT, "examples", "rpo_planner_env")
        end

        withenv("SPACEAGORA_EDG_ENV" => joinpath(temporary, "missing")) do
            @test_throws ArgumentError harness.run_example(probe)
        end
        invalid = joinpath(temporary, "invalid")
        mkpath(invalid)
        for dependencies in (
            Dict("SpaceAGORA" => "afbfb69f-5c0b-4832-b760-43725dff8540"),
            Dict("SpaceAGORA" => "afbfb69f-5c0b-4832-b760-43725dff8540",
                 "SpaceAGORAEDG" => "00000000-0000-0000-0000-000000000000"),
        )
            open(joinpath(invalid, "Project.toml"), "w") do io
                TOML.print(io, Dict("deps" => dependencies))
            end
            withenv("SPACEAGORA_EDG_ENV" => invalid) do
                @test_throws ArgumentError harness.run_example(probe)
            end
        end
    end
end
