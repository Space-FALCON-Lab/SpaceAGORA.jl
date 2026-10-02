using Test

@testset "Planning owners retain architecture enforcement" begin
    repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
    # Include a nested owner so the fixtures also enforce recursive scanning.
    owner_paths = ("shared", "hypr", "rrt", joinpath("shared", "rpo"))
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
        for owner in owner_paths
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
            for owner in owner_paths
                probe = joinpath(fixture, "src", "gnc", owner, "boundary_probe.jl")
                write(probe, forbidden)
                @test_throws LoadError run_gate()
                rm(probe)
                @test isnothing(run_gate())
            end
        end
    end
end

@testset "Shared metric boundary rejects an algorithm configuration dependency" begin
    repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
    path = joinpath(repo, "src", "gnc", "shared", "rpo", "path_metrics.jl")
    source = read(path, String)
    function load_metric_source(text)
        isolated = Module(gensym(:MetricBoundary))
        Core.eval(isolated, :(using LinearAlgebra, StaticArrays))
        Base.include_string(isolated, text, path)
        return isolated
    end
    valid = load_metric_source(source)
    @test !isdefined(valid, :RPOPSOConfig)
    @test isdefined(valid, :rpo_path_cost_normalization_refs)
    @test isdefined(valid, :rpo_fuel_proxy_from_samples)
    # The forbidden annotation exists only in this isolated fixture, never in src.
    forbidden = source * "\nmetric_boundary_probe(points, cfg::RPOPSOConfig) = nothing\n"
    rejected = try load_metric_source(forbidden); nothing catch error; error end
    @test rejected isa LoadError && rejected.error isa UndefVarError
    @test rejected isa LoadError && rejected.error isa UndefVarError && rejected.error.var === :RPOPSOConfig
    @test !isdefined(load_metric_source(source), :RPOPSOConfig)
end
