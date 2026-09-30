using Test

@testset "Planning owners retain architecture enforcement" begin
    repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
    gates = (
        "ci_no_legacy_include_chains_gate.jl" => "__legacy_probe = nothing\n",
        "ci_no_guidance_control_cross_include_gate.jl" => "include(\"control/probe.jl\")\n",
        "ci_gnc_aerobraking_boundary_gate.jl" => "using DynamicEffectors\n",
        "ci_gnc_typed_command_boundary_gate.jl" => "using ThrusterModels\n",
    )
    # Only copy the small source files required by these gates. Each fixture is
    # isolated, so a forbidden source file never enters the working checkout.
    required = (
        "src/gnc/command_types.jl",
        "src/gnc/guidance/guidance_hooks.jl",
        "src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl",
        "src/gnc/control/control_hooks.jl",
        "src/gnc/control/propulsive_maneuvers.jl",
        "src/gnc/guidance/aerobraking/interfaces.jl",
        "src/gnc/guidance/aerobraking/dispatcher.jl",
        "src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl",
        "src/gnc/guidance/aerobraking/t_edg/t_edg_strategy.jl",
        "src/gnc/control/aerobraking/tracking_executor.jl",
        "src/gnc/control/aerobraking/control_commands.jl",
        "src/gnc/control/aerobraking/constraint_tracking.jl",
        "src/mission/operations/aerobraking_policy/policy_types.jl",
    )
    mktempdir() do fixture
        for rel in required
            target = joinpath(fixture, rel)
            mkpath(dirname(target))
            cp(joinpath(repo, rel), target)
        end
        for owner in ("shared", "hypr", "rrt")
            mkpath(joinpath(fixture, "src", "gnc", owner))
        end
        for (gate, forbidden) in gates
            path = joinpath(repo, "test", "gates", gate)
            source = read(path, String)
            root_line = "const REPO_ROOT = normpath(joinpath(@__DIR__, \"..\", \"..\"))"
            @test occursin(root_line, source)
            isolated_source = replace(source, root_line => "const REPO_ROOT = $(repr(fixture))"; count=1)
            run_gate() = redirect_stdout(devnull) do
                Base.include_string(Module(gensym(:PlanningGate)), isolated_source, path)
            end
            # Establish that unrelated required-file checks do not cause the failure.
            @test isnothing(run_gate())
            for owner in ("shared", "hypr", "rrt")
                probe = joinpath(fixture, "src", "gnc", owner, "boundary_probe.jl")
                write(probe, forbidden)
                @test_throws LoadError run_gate()
                rm(probe)
                @test isnothing(run_gate())
            end
        end
    end
end
