using Test

@testset "Planning owners retain architecture enforcement" begin
    repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
    # Include a nested owner so the fixtures also enforce recursive scanning.
    owner_paths = (
        joinpath("src", "models", "gnc", "planning", "common"),
        joinpath("src", "models", "gnc", "companions", "hypr"),
        joinpath("src", "models", "gnc", "companions", "edg"),
        joinpath("src", "models", "gnc", "planning", "rpo"),
        joinpath("src", "models", "gnc", "planning", "robot_arm"),
        joinpath("src", "models", "gnc", "planning", "rpo", "nested_probe"),
        joinpath("src", "interfaces", "analysis", "planning"),
    )
    gates = (
        ("ci_src_completeness_contract_gate.jl", "boundary_stub() = nothing\n",
            "Source completeness contract gate failed", "contains forbidden silent stub pattern"),
        ("ci_no_legacy_include_chains_gate.jl", "__legacy_probe = nothing\n",
            "Legacy include-chain gate failed", "contains forbidden legacy token '__legacy_'"),
        ("ci_no_guidance_control_cross_include_gate.jl", "include(\"control/probe.jl\")\n",
            "Guidance/control cross-include gate failed", "include(\"control/probe.jl\")"),
        ("ci_gnc_aerobraking_boundary_gate.jl", "using DynamicEffectors\n",
            "GNC aerobraking boundary gate failed", "GNC source still depends on DynamicEffectors directly"),
        ("ci_gnc_typed_command_boundary_gate.jl", "using ThrusterModels\n",
            "GNC typed-command boundary gate failed", "guidance must not depend on ThrusterModels directly"),
    )
    # Only copy the small source files required by these gates. Each fixture is
    # isolated, so a forbidden source file never enters the working checkout.
    required = (
        "src/models/gnc/contracts/command_types.jl",
        "src/models/gnc/guidance/hooks.jl",
        "src/models/gnc/guidance/maneuvers/thruster_guidance_functions.jl",
        "src/models/gnc/control/hooks.jl",
        "src/models/gnc/control/maneuvers/propulsive_maneuvers.jl",
        # The aerobraking gate also verifies the complete typed EDG owner.
        "src/models/gnc/guidance/aerobraking/typed_edg/algorithms.jl",
        "src/models/gnc/guidance/aerobraking/typed_edg/services.jl",
        "src/models/gnc/guidance/aerobraking/typed_edg/adapter_queries.jl",
        "src/models/gnc/companions/edg/services.jl",
    )
    mktempdir() do fixture
        for rel in required
            target = joinpath(fixture, rel)
            mkpath(dirname(target))
            cp(joinpath(repo, rel), target)
        end
        # The gates require every scan root to exist so a source move cannot
        # silently disable a scan. Unrelated owners can be empty in this fixture;
        # the copied source and injected probes remain the only files inspected.
        scan_roots = (
            joinpath("src", "foundation", "core"),
            joinpath("src", "orchestration", "simulation", "composition"),
            joinpath("src", "orchestration", "simulation", "state"),
            joinpath("src", "orchestration", "simulation", "config", "presets"),
            joinpath("src", "models", "environment"),
            joinpath("src", "models", "gnc", "control"),
            joinpath("src", "models", "gnc", "guidance"),
            joinpath("src", "models", "gnc", "contracts"),
            joinpath("src", "orchestration", "simulation", "engine"),
            joinpath("src", "orchestration", "simulation", "callbacks"),
            joinpath("src", "orchestration", "parallel"),
        )
        for owner in (owner_paths..., scan_roots...)
            mkpath(joinpath(fixture, owner))
        end
        for (gate, forbidden, failure_message, violation_message) in gates
            path = joinpath(repo, "test", "gates", gate)
            source = read(path, String)
            root_line = "const REPO_ROOT = normpath(joinpath(@__DIR__, \"..\", \"..\"))"
            @test occursin(root_line, source)
            isolated_source = replace(source, root_line => "const REPO_ROOT = $(repr(fixture))"; count=1)
            function run_gate()
                isolated = Module(gensym(:PlanningGate))
                # Module(name) has no include binding. Gate helpers must load into
                # the same isolated module, with paths resolved from the gate file.
                Core.eval(isolated, :(include(path) = Base.include(@__MODULE__, path)))
                failure = try
                    redirect_stdout(devnull) do
                        Base.include_string(isolated, isolated_source, path)
                    end
                    nothing
                catch err
                    err
                end
                # Read newly evaluated bindings in the gate's own module.
                violations = Core.eval(isolated,
                    :(isdefined(@__MODULE__, :violations) ? copy(violations) : String[]))
                return (; failure, violations)
            end
            # Establish that unrelated required-file checks do not cause the failure.
            @test isnothing(run_gate().failure)
            for owner in owner_paths
                probe = joinpath(fixture, owner, "boundary_probe.jl")
                write(probe, forbidden)
                try
                    rejected = run_gate()
                    @test rejected.failure isa LoadError
                    @test rejected.failure isa LoadError &&
                        rejected.failure.error isa ErrorException &&
                        rejected.failure.error.msg == failure_message
                    # Reject this injected dependency specifically, rather than an
                    # unrelated include error or another incomplete source fixture.
                    prefix = relpath(probe, fixture) * ":"
                    @test any(v -> startswith(v, prefix) && occursin(violation_message, v),
                        rejected.violations)
                finally
                    rm(probe)
                end
                @test isnothing(run_gate().failure)
            end
            if gate == "ci_gnc_typed_command_boundary_gate.jl"
                @testset "Comparison tracking permits only its six-axis type import" begin
                    rel = joinpath("src", "interfaces", "analysis", "planning", "rpo_comparison_tracking.jl")
                    tracking = joinpath(fixture, rel)
                    original = read(joinpath(repo, rel), String)
                    allowed_import = "using ..ThrusterModels: SixAxisThrusterModel"
                    @test count(allowed_import, original) == 1
                    write(tracking, original)
                    try
                        @test isnothing(run_gate().failure)
                        thruster_dependency = "guidance must not depend on ThrusterModels directly"
                        base_dependency = "guidance must not depend on BaseThrusterModel directly"
                        abstract_dependency = "guidance must not depend on AbstractThrusterModel directly"
                        # Exercise the actual gate on actual source, not a second
                        # implementation of its allowlist or forbidden patterns.
                        cases = [
                            ("duplicate import", original * "\n" * allowed_import * "\n", [thruster_dependency]),
                            ("semicolon suffix", replace(original, allowed_import => allowed_import * "; nothing"; count=1), [thruster_dependency]),
                            ("combined import", replace(original, allowed_import => allowed_import * ", BaseThrusterModel"; count=1), [thruster_dependency, base_dependency]),
                            ("other imported type", replace(original, allowed_import => "using ..ThrusterModels: OtherThrusterModel"; count=1), [thruster_dependency]),
                            ("abstract imported type", replace(original, allowed_import => "using ..ThrusterModels: AbstractThrusterModel"; count=1), [thruster_dependency, abstract_dependency]),
                            ("qualified access", original * "\nprobe = ThrusterModels.SixAxisThrusterModel\n", [thruster_dependency]),
                            ("base type access", original * "\nprobe = BaseThrusterModel\n", [base_dependency]),
                            ("legacy effector access", original * "\nusing DynamicEffectors\n", ["GNC source still depends on DynamicEffectors directly"]),
                            ("command buffer access", original * "\nprobe = p.shared_buffers.maneuver_commands\n", ["maneuver_commands ownership is restricted to typed guidance writer/control consumer files"]),
                        ]
                        for (field, message) in (
                            ("Δv", "guidance must not mutate thruster delta-v state directly"),
                            ("direction", "guidance must not mutate thruster direction state directly"),
                            ("start_burn_time", "guidance must not mutate thruster burn scheduling directly"),
                            ("stop_burn_time", "guidance must not mutate thruster burn scheduling directly"),
                            ("thrust", "guidance must not mutate thruster thrust state directly"),
                        )
                            push!(cases, ("direct $field mutation", original * "\ncontrol.$field[1] = 0.0\n", [message]))
                        end
                        for (label, changed, messages) in cases
                            @testset "$label" begin
                                @test changed != original
                                write(tracking, changed)
                                rejected = run_gate()
                                @test rejected.failure isa LoadError &&
                                    rejected.failure.error isa ErrorException &&
                                    rejected.failure.error.msg == failure_message
                                @test sort(unique(rejected.violations)) == sort(["$rel: $message" for message in messages])
                                write(tracking, original)
                                @test isnothing(run_gate().failure)
                            end
                        end
                        # Identical source in another file is not allowlisted.
                        rm(tracking)
                        wrong_rel = joinpath("src", "interfaces", "analysis", "planning", "tracking_import_probe.jl")
                        wrong_path = joinpath(fixture, wrong_rel)
                        write(wrong_path, original)
                        try
                            rejected = run_gate()
                            @test rejected.failure isa LoadError &&
                                rejected.failure.error isa ErrorException &&
                                rejected.failure.error.msg == failure_message
                            @test rejected.violations == ["$wrong_rel: $thruster_dependency"]
                        finally
                            rm(wrong_path)
                        end
                    finally
                        isfile(tracking) && rm(tracking)
                    end
                    @test isnothing(run_gate().failure)
                end
            end
            if gate == "ci_src_completeness_contract_gate.jl"
                dispatcher = joinpath(fixture, "src", "models", "gnc", "guidance",
                    "aerobraking", "typed_edg", "algorithms.jl")
                original = read(dispatcher, String)
                try
                    write(dispatcher, replace(original,
                        "_invoke_kernel(:guidance_decision!" => "_invoke_kernel(:_edg_wrong_kernel"))
                    rejected = run_gate()
                    @test rejected.failure isa LoadError
                    @test any(v -> occursin("must be an exact package forwarder", v), rejected.violations)
                finally
                    write(dispatcher, original)
                end
                @test isnothing(run_gate().failure)
            end
        end
    end
end

@testset "Shared metric boundary rejects an algorithm configuration dependency" begin
    repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
    path = joinpath(repo, "src", "models", "gnc", "planning", "rpo", "path_metrics.jl")
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
