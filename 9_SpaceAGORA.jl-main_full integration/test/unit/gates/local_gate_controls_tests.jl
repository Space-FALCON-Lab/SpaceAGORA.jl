using Test
using TOML

# Negative controls for the four contract gates that the full local runner
# (test/contracts/runtests.jl) adds on top of test/contracts/pr_runtests.jl.
# Each gate runs against an isolated copy of the source tree, so a forbidden
# construct never enters the working checkout. The copy is the whole src/ tree
# plus the small inputs the gates read, and injection targets are located by
# file name, so the controls follow the gates across source moves.

const _GATE_ROOT_LINE = "const REPO_ROOT = normpath(joinpath(@__DIR__, \"..\", \"..\"))"

function _copy_into!(fixture::String, repo::String, rel::String)
    source = joinpath(repo, rel)
    target = joinpath(fixture, rel)
    mkpath(dirname(target))
    cp(source, target)
    return target
end

function _unique_source(fixture::String, name::String)
    hits = String[]
    for (root, _, files) in walkdir(joinpath(fixture, "src"))
        name in files && push!(hits, joinpath(root, name))
    end
    length(hits) == 1 || error("Expected exactly one $(name) under the fixture source tree; found $(length(hits)).")
    return only(hits)
end

function _gate_runner(repo::String, fixture::String, gate::String)
    path = joinpath(repo, "test", "gates", gate)
    source = read(path, String)
    occursin(_GATE_ROOT_LINE, source) || error("$(gate) does not define the expected REPO_ROOT line.")
    isolated_source = replace(source, _GATE_ROOT_LINE => "const REPO_ROOT = $(repr(fixture))"; count=1)
    return function ()
        isolated = Module(gensym(Symbol(gate)))
        # Module(name) has no include binding. Gate helpers must load into the
        # same isolated module, with paths resolved from the gate file.
        Core.eval(isolated, :(include(path) = Base.include(@__MODULE__, path)))
        failure = try
            redirect_stdout(devnull) do
                Base.include_string(isolated, isolated_source, path)
            end
            nothing
        catch err
            err
        end
        violations = Core.eval(isolated,
            :(isdefined(@__MODULE__, :violations) ? copy(violations) : String[]))
        return (; failure, violations)
    end
end

_failed_with(result, message::String) =
    result.failure isa LoadError &&
    result.failure.error isa ErrorException &&
    result.failure.error.msg == message

_failed_starting(result, prefix::String) =
    result.failure isa LoadError &&
    result.failure.error isa ErrorException &&
    startswith(result.failure.error.msg, prefix)

_reports(result, probe::String, fixture::String, fragment::String) =
    any(v -> startswith(v, relpath(probe, fixture) * ":") && occursin(fragment, v), result.violations)

function _with_text(f, path::String, text::String)
    original = read(path, String)
    write(path, text)
    try
        f()
    finally
        write(path, original)
    end
end

function _with_probe(f, path::String, text::String)
    mkpath(dirname(path))
    write(path, text)
    try
        f()
    finally
        rm(path)
    end
end

@testset "Local contract gates reject deliberate violations" begin
    repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
    mktempdir() do fixture
        _copy_into!(fixture, repo, "src")
        _copy_into!(fixture, repo, "Project.toml")
        _copy_into!(fixture, repo, joinpath("packages", "SpaceAGORAHYPR", "src"))

        @testset "Project compat gate" begin
            run_gate = _gate_runner(repo, fixture, "ci_project_compat_gate.jl")
            @test isnothing(run_gate().failure)

            project_path = joinpath(fixture, "Project.toml")
            project = TOML.parsefile(project_path)
            compat = project["compat"]
            # A real package dependency whose compat entry the gate must require.
            @test haskey(compat, "StaticArrays")
            # A standard library with a compat entry; the gate exempts standard libraries.
            @test haskey(compat, "Base64")
            original = read(project_path, String)
            # Remove the entry from the [compat] table only; [deps] lists the same names first.
            compat_start = first(findfirst("[compat]", original))
            function without(name)
                head = original[1:compat_start - 1]
                tail = original[compat_start:end]
                return head * replace(tail, Regex("^" * name * " = .*\\n", "m") => ""; count=1)
            end
            @test without("StaticArrays") != original
            @test without("Base64") != original

            _with_text(project_path, without("StaticArrays")) do
                @test _failed_with(run_gate(), "Missing compat entries in Project.toml: StaticArrays")
            end
            _with_text(project_path, without("Base64")) do
                @test isnothing(run_gate().failure)
            end
            @test isnothing(run_gate().failure)
        end

        @testset "Translational ownership gate" begin
            run_gate = _gate_runner(repo, fixture, "ci_translational_ownership_gate.jl")
            @test isnothing(run_gate().failure)

            owner_dir = dirname(_unique_source(fixture, "translational_models.jl"))
            outsider_dir = dirname(_unique_source(fixture, "guidance_hooks.jl"))
            forbidden = "du_view.pos .= sc_view.vel\n"
            outsider_probe = joinpath(outsider_dir, "translational_probe.jl")

            _with_probe(outsider_probe, forbidden) do
                rejected = run_gate()
                @test _failed_with(rejected, "Translational ownership gate failed")
                @test _reports(rejected, outsider_probe, fixture, "forbidden inline translational dynamics pattern")
            end
            # Comments are not active source.
            _with_probe(outsider_probe, "# " * forbidden) do
                @test isnothing(run_gate().failure)
            end
            # The owner directory may contain the equations.
            _with_probe(joinpath(owner_dir, "translational_probe.jl"), forbidden) do
                @test isnothing(run_gate().failure)
            end
            @test isnothing(run_gate().failure)
        end

        @testset "Runtime legacy surface gate" begin
            run_gate = _gate_runner(repo, fixture, "ci_runtime_legacy_surface_gate.jl")
            @test isnothing(run_gate().failure)

            # A shared forbidden token in a listed file.
            stub = _unique_source(fixture, "selector_stub.jl")
            _with_text(stub, read(stub, String) * "\ngravity_const(1.0)\n") do
                rejected = run_gate()
                @test _failed_with(rejected, "Runtime legacy surface gate failed")
                @test _reports(rejected, stub, fixture, "gravity_const")
            end
            # A file-specific forbidden token.
            setup = _unique_source(fixture, "setup.jl")
            _with_text(setup, read(setup, String) * "\n_compat_probe = nothing\n") do
                rejected = run_gate()
                @test _failed_with(rejected, "Runtime legacy surface gate failed")
                @test _reports(rejected, setup, fixture, "_compat_")
            end
            # Directory targets are scanned recursively.
            control_dir = dirname(_unique_source(fixture, "tracking_executor.jl"))
            nested_probe = joinpath(control_dir, "nested", "legacy_probe.jl")
            _with_probe(nested_probe, "gravity_GRAM(x)\n") do
                rejected = run_gate()
                @test _failed_with(rejected, "Runtime legacy surface gate failed")
                @test _reports(rejected, nested_probe, fixture, "gravity_GRAM")
            end
            rm(dirname(nested_probe))
            # A missing listed target is an explicit failure, not a silent skip.
            helpers = _unique_source(fixture, "bridge_helpers.jl")
            helpers_source = read(helpers, String)
            rm(helpers)
            try
                @test _failed_with(run_gate(),
                    "Missing runtime legacy-surface target: " * relpath(helpers, fixture))
            finally
                write(helpers, helpers_source)
            end
            @test isnothing(run_gate().failure)
        end

        @testset "Ownership gates contain by path, not by name prefix" begin
            for (gate, owner_file, forbidden, failure, message) in (
                ("ci_translational_ownership_gate.jl", "translational_models.jl",
                    "du_view.pos .= sc_view.vel\n", "Translational ownership gate failed",
                    "forbidden inline translational dynamics pattern"),
                ("ci_rotational_ownership_gate.jl", "rotational_models.jl",
                    "dq = 0.5 * quat_mult(q, w)\n", "Rotational ownership gate failed",
                    "forbidden inline rotational dynamics pattern"),
            )
                run_gate = _gate_runner(repo, fixture, gate)
                @test isnothing(run_gate().failure)
                owner_dir = dirname(_unique_source(fixture, owner_file))
                # A sibling whose name merely starts with the owner's name is not the owner.
                sibling = joinpath(dirname(owner_dir), basename(owner_dir) * "_outside", "probe.jl")
                _with_probe(sibling, forbidden) do
                    rejected = run_gate()
                    @test _failed_with(rejected, failure)
                    @test _reports(rejected, sibling, fixture, message)
                end
                rm(dirname(sibling))
                # Nested directories inside the owner remain exempt.
                nested = joinpath(owner_dir, "nested", "probe.jl")
                _with_probe(nested, forbidden) do
                    @test isnothing(run_gate().failure)
                end
                rm(dirname(nested))
                # A missing owner root is an explicit failure, not an empty scan.
                parked = owner_dir * "_parked"
                mv(owner_dir, parked)
                try
                    @test _failed_with(run_gate(), "Missing " * split(basename(owner_dir), "_")[1] * " ownership root")
                finally
                    mv(parked, owner_dir)
                end
                @test isnothing(run_gate().failure)
            end
        end

        @testset "Benchmark mode parity gate" begin
            # The gate loads SpaceAGORA through the fixture root on LOAD_PATH; keep
            # the process load path unchanged afterwards.
            modes_rel = joinpath("benchmarks", "studies", "parallelization_performance", "modes.jl")
            _copy_into!(fixture, repo, modes_rel)
            modes = joinpath(fixture, modes_rel)
            original = read(modes, String)
            gate = _gate_runner(repo, fixture, "ci_ppc_mode_profile_parity_gate.jl")
            run_gate = () -> (result = gate(); filter!(!=(fixture), LOAD_PATH); result)
            @test isnothing(run_gate().failure)
            disagree = "Benchmark mode table disagrees with the shipped profiles"
            # An inline mode block that drifts from its profile.
            inline = replace(original, "\"outer_threads\" => PPCModeSpec(\n            name=\"outer_threads\",\n            profile=\"R1a\",\n            backend=\"threads\",\n            outer_active=true,\n            policy_adaptive=false,\n            rhs_batch=\"auto\",\n            density=\"off\"" =>
                "\"outer_threads\" => PPCModeSpec(\n            name=\"outer_threads\",\n            profile=\"R1a\",\n            backend=\"threads\",\n            outer_active=true,\n            policy_adaptive=false,\n            rhs_batch=\"auto\",\n            density=\"bogus\""; count=1)
            @test inline != original
            _with_text(modes, inline) do
                @test _failed_starting(run_gate(), disagree)
            end
            # The serial mode is bound to a local name before it enters the table;
            # the gate must follow that binding rather than skip the mode.
            serial = replace(original, "    serial = PPCModeSpec(\n        name=\"serial\",\n        profile=\"R0\",\n        backend=\"none\",\n        outer_active=false,\n        policy_adaptive=false,\n        rhs_batch=\"off\",\n        density=\"off\"" =>
                "    serial = PPCModeSpec(\n        name=\"serial\",\n        profile=\"R0\",\n        backend=\"none\",\n        outer_active=false,\n        policy_adaptive=false,\n        rhs_batch=\"off\",\n        density=\"bogus\""; count=1)
            @test serial != original
            _with_text(modes, serial) do
                @test _failed_starting(run_gate(), disagree)
            end
            # A checked mode that disappears from the table is an explicit failure.
            _with_text(modes, replace(original, "\"serial\" => serial," => ""; count=1)) do
                @test _failed_starting(run_gate(), "Benchmark mode table has no PPCModeSpec for checked mode \"serial\"")
            end
            # A missing table is an explicit failure.
            rm(modes)
            try
                @test _failed_starting(run_gate(), "modes.jl not found")
            finally
                write(modes, original)
            end
            @test isnothing(run_gate().failure)
        end

        @testset "Typed command boundary gate" begin
            run_gate = _gate_runner(repo, fixture, "ci_gnc_typed_command_boundary_gate.jl")
            @test isnothing(run_gate().failure)

            # The control owner must import the typed command from CommandTypes.
            hooks = _unique_source(fixture, "control_hooks.jl")
            hooks_source = read(hooks, String)
            needle = "using ..CommandTypes: PropulsiveManeuverCommand"
            @test occursin(needle, hooks_source)
            _with_text(hooks, replace(hooks_source, needle => "using ..CommandTypes: AerobrakingControlCommand"; count=1)) do
                @test _failed_with(run_gate(),
                    "Typed-command boundary violation: ControlHooks is not importing PropulsiveManeuverCommand from CommandTypes.")
            end
            # maneuver_commands may only be touched by the typed writer and consumer.
            guidance_dir = dirname(_unique_source(fixture, "guidance_hooks.jl"))
            probe = joinpath(guidance_dir, "ownership_probe.jl")
            _with_probe(probe, "x = p.shared_buffers.maneuver_commands\n") do
                rejected = run_gate()
                @test _failed_with(rejected, "GNC typed-command boundary gate failed")
                @test _reports(rejected, probe, fixture, "maneuver_commands ownership is restricted")
            end
            @test isnothing(run_gate().failure)
        end
    end
end
