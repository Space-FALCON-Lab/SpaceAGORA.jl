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

# The one source file whose path ends with `name`: a file name, or a folder and a
# file name ("guidance/hooks.jl") where the file name alone is not unique.
function _unique_source(fixture::String, name::String)
    tail = splitpath(name)
    hits = String[]
    for (root, _, files) in walkdir(joinpath(fixture, "src")), file in files
        parts = splitpath(joinpath(root, file))
        length(parts) >= length(tail) && parts[end-length(tail)+1:end] == tail && push!(hits, joinpath(root, file))
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
        _copy_into!(fixture, repo, joinpath("docs", "Project.toml"))
        _copy_into!(fixture, repo, joinpath("test", "Project.toml"))
        _copy_into!(fixture, repo, joinpath("examples", "Project.toml"))
        _copy_into!(fixture, repo, joinpath("examples", "odyssey_surrogate_env", "Project.toml"))
        _copy_into!(fixture, repo, joinpath("test", "grid_atmosphere", "Project.toml"))
        _copy_into!(fixture, repo, joinpath("build", "build_sysimage.jl"))

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

            # The docs environment's own dependencies need compat entries too.
            docs_project_path = joinpath(fixture, "docs", "Project.toml")
            docs_original = read(docs_project_path, String)
            @test haskey(TOML.parsefile(docs_project_path)["compat"], "Documenter")
            docs_without = replace(docs_original, r"^Documenter = \"1\"\n"m => ""; count=1)
            @test docs_without != docs_original
            _with_text(docs_project_path, docs_without) do
                @test _failed_with(run_gate(), "Missing compat entries in docs/Project.toml: Documenter")
            end
            @test isnothing(run_gate().failure)

            # So do the test environment's.
            test_project_path = joinpath(fixture, "test", "Project.toml")
            test_original = read(test_project_path, String)
            @test haskey(TOML.parsefile(test_project_path)["compat"], "Aqua")
            # Remove the entry from the [compat] table only; [deps] lists the same name first.
            test_compat_start = first(findfirst("[compat]", test_original))
            test_without = test_original[1:test_compat_start - 1] *
                replace(test_original[test_compat_start:end], r"^Aqua = .*\n"m => ""; count=1)
            @test !haskey(TOML.parse(test_without)["compat"], "Aqua")
            @test haskey(TOML.parse(test_without)["deps"], "Aqua")
            _with_text(test_project_path, test_without) do
                @test _failed_with(run_gate(), "Missing compat entries in test/Project.toml: Aqua")
            end
            @test isnothing(run_gate().failure)

            # And the examples environment's.
            examples_project_path = joinpath(fixture, "examples", "Project.toml")
            examples_original = read(examples_project_path, String)
            @test haskey(TOML.parsefile(examples_project_path)["compat"], "Plots")
            examples_compat_start = first(findfirst("[compat]", examples_original))
            examples_without = examples_original[1:examples_compat_start - 1] *
                replace(examples_original[examples_compat_start:end], r"^Plots = .*\n"m => ""; count=1)
            @test !haskey(TOML.parse(examples_without)["compat"], "Plots")
            @test haskey(TOML.parse(examples_without)["deps"], "Plots")
            _with_text(examples_project_path, examples_without) do
                @test _failed_with(run_gate(), "Missing compat entries in examples/Project.toml: Plots")
            end
            @test isnothing(run_gate().failure)

            # A present bound can still exclude the checkout, including in the
            # nested smoke environments that the ordinary dependency scan omits.
            core_version = VersionNumber(project["version"])
            excluded_bound = "=$(core_version.major + 1).0.0"
            for relative in ("docs", "test", "examples", "examples/odyssey_surrogate_env", "test/grid_atmosphere")
                bundled_path = joinpath(fixture, relative, "Project.toml")
                bundled = TOML.parsefile(bundled_path)
                bundled["compat"]["SpaceAGORA"] = excluded_bound
                _with_text(bundled_path, sprint(io -> TOML.print(io, bundled))) do
                    @test _failed_with(run_gate(),
                        "$relative/Project.toml excludes checked-out SpaceAGORA $core_version (compat: $excluded_bound)")
                end
                @test isnothing(run_gate().failure)
            end
        end

        @testset "Sysimage package lists gate" begin
            run_gate = _gate_runner(repo, fixture, "ci_sysimage_package_lists_gate.jl")
            @test isnothing(run_gate().failure)

            builder = joinpath(fixture, "build", "build_sysimage.jl")
            original = read(builder, String)
            prefix = "Sysimage package lists in build/build_sysimage.jl name packages that are not " *
                     ":SpaceAGORA or a Project.toml [deps] entry: "
            full_head = "const FULL_PACKAGES = [\n"
            ci_head = "const CI_PACKAGES = [\n"
            @test occursin(full_head, original)
            @test occursin(ci_head, original)

            # The two names both lists carried after their dependencies were removed.
            stale = "    :Quaternions,\n    :SatelliteToolbox,\n"
            both = replace(replace(original, full_head => full_head * stale; count=1), ci_head => ci_head * stale; count=1)
            _with_text(builder, both) do
                @test _failed_with(run_gate(), prefix *
                    "FULL_PACKAGES: Quaternions, SatelliteToolbox; CI_PACKAGES: Quaternions, SatelliteToolbox")
            end
            # Each list is checked on its own.
            _with_text(builder, replace(original, ci_head => ci_head * "    :Quaternions,\n"; count=1)) do
                @test _failed_with(run_gate(), prefix * "CI_PACKAGES: Quaternions")
            end
            _with_text(builder, replace(original, full_head => full_head * "    :SatelliteToolbox,\n"; count=1)) do
                @test _failed_with(run_gate(), prefix * "FULL_PACKAGES: SatelliteToolbox")
            end
            # A list that is not a literal vector of symbols, or is missing, cannot be checked: fail.
            _with_text(builder, replace(original, ci_head => "const CI_PACKAGES = Symbol[\n"; count=1)) do
                @test _failed_with(run_gate(),
                    "build/build_sysimage.jl: `const CI_PACKAGES` must be a literal vector of symbols")
            end
            _with_text(builder, replace(original, "const FULL_PACKAGES = [" => "const FULL_PACKAGES_RENAMED = ["; count=1)) do
                @test _failed_with(run_gate(),
                    "build/build_sysimage.jl must define `const FULL_PACKAGES` exactly once; found 0")
            end

            # Dropping a dependency from the project makes its builder entries stale.
            project_path = joinpath(fixture, "Project.toml")
            project = read(project_path, String)
            without_dep = replace(project, r"^SpecialFunctions = .*\n"m => ""; count=1)
            @test without_dep != project
            @test !haskey(TOML.parse(without_dep)["deps"], "SpecialFunctions")
            _with_text(project_path, without_dep) do
                @test _failed_with(run_gate(), prefix * "FULL_PACKAGES: SpecialFunctions; CI_PACKAGES: SpecialFunctions")
            end
            @test isnothing(run_gate().failure)
        end

        @testset "Translational ownership gate" begin
            run_gate = _gate_runner(repo, fixture, "ci_translational_ownership_gate.jl")
            @test isnothing(run_gate().failure)

            owner_dir = dirname(_unique_source(fixture, "translational_models.jl"))
            outsider_dir = dirname(_unique_source(fixture, "guidance/hooks.jl"))
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
            listed = _unique_source(fixture, "io_paths.jl")
            _with_text(listed, read(listed, String) * "\ngravity_const(1.0)\n") do
                rejected = run_gate()
                @test _failed_with(rejected, "Runtime legacy surface gate failed")
                @test _reports(rejected, listed, fixture, "gravity_const")
            end
            # A file-specific forbidden token.
            setup = _unique_source(fixture, "setup.jl")
            _with_text(setup, read(setup, String) * "\n_compat_probe = nothing\n") do
                rejected = run_gate()
                @test _failed_with(rejected, "Runtime legacy surface gate failed")
                @test _reports(rejected, setup, fixture, "_compat_")
            end
            # Directory targets are scanned recursively.
            guidance_dir = dirname(_unique_source(fixture, "adapter_queries.jl"))
            nested_probe = joinpath(guidance_dir, "nested", "legacy_probe.jl")
            _with_probe(nested_probe, "gravity_GRAM(x)\n") do
                rejected = run_gate()
                @test _failed_with(rejected, "Runtime legacy surface gate failed")
                @test _reports(rejected, nested_probe, fixture, "gravity_GRAM")
            end
            rm(dirname(nested_probe))
            # A missing listed target is an explicit failure, not a silent skip.
            persistence = _unique_source(fixture, "persistence.jl")
            persistence_source = read(persistence, String)
            rm(persistence)
            try
                @test _failed_with(run_gate(),
                    "Missing runtime legacy-surface target: " * relpath(persistence, fixture))
            finally
                write(persistence, persistence_source)
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

        @testset "Extension hook reference gate" begin
            for rel in (joinpath("docs", "extension_hooks.jl"), joinpath("docs", "public_api_symbols.jl"),
                        joinpath("docs", "src"))
                _copy_into!(fixture, repo, rel)
            end
            run_gate = _gate_runner(repo, fixture, "ci_extension_hook_reference_gate.jl")
            @test isnothing(run_gate().failure)
            registry = joinpath(fixture, "docs", "extension_hooks.jl")
            original = read(registry, String)
            rejects(result, fragment) = _failed_with(result, "Extension hook reference gate failed") &&
                any(v -> occursin(fragment, v), result.violations)
            function with_registry(f, pairs::Pair...)
                edited = original
                for (needle, replacement) in pairs
                    @test occursin(needle, edited)
                    edited = replace(edited, needle => replacement; count=1)
                end
                _with_text(f, registry, edited)
            end
            # Faults go into rows and list heads that later pull requests leave alone,
            # never into the known gaps, which those pull requests close.
            row(symbol) = match(Regex("        hook\\(:" * symbol * ";.*?\\),\\n", "s"), original).match
            in_row(symbol, change::Pair) = row(symbol) => replace(row(symbol), change; count=1)
            not_hooks_head = "const NOT_HOOKS = [\n"
            gaps_head = "const KNOWN_GAPS = [\n"
            # Rule 1: a row for a name SpaceAGORA does not export.
            with_registry(in_row("touchdown_spec", "hook(:touchdown_spec;" => "hook(:touchdown_specs;")) do
                @test rejects(run_gate(), "touchdown_specs: not exported from SpaceAGORA.")
            end
            # Rule 2: an exported hook with no row, NOT_HOOKS entries that are rows or
            # outside the hook sections, and a renamed hook section.
            with_registry(row("touchdown_spec") => "") do
                @test rejects(run_gate(), "touchdown_spec (\"Model Extension Hooks\"): exported function with no registry row")
            end
            with_registry(not_hooks_head => not_hooks_head * "    (symbol = :wrench, reason = \"probe\"),\n    (symbol = :run_simulation, reason = \"probe\"),\n") do
                rejected = run_gate()
                @test rejects(rejected, "wrench: both a registry row and in NOT_HOOKS.")
                @test rejects(rejected, "NOT_HOOKS lists run_simulation, which is not an exported function in a hook section.")
            end
            with_registry("    \"Sensors\",\n" => "    \"Sensor\",\n") do
                rejected = run_gate()
                @test rejects(rejected, "HOOK_SECTIONS names \"Sensor\"")
                @test rejects(rejected, "sensor_measure: listed in \"Sensors\", which is not in HOOK_SECTIONS.")
            end
            # Rule 3: a signature that is not the docstring's first line.
            wrench_signature = "wrench(model, x::StateSample, env::EnvironmentSample, t::Float64) -> (force_ii, torque_body)"
            with_registry(in_row("wrench", wrench_signature => "wrench(model, x, env, t) -> (force_ii, torque_body)")) do
                @test rejects(run_gate(), "wrench: signature cell \"wrench(model, x, env, t) -> (force_ii, torque_body)\" differs from the docstring's first line")
            end
            # Rule 4: arity, dispatch argument and both directions of the fallback.
            with_registry(in_row("wrench", wrench_signature => "wrench(model, x, env) -> (force_ii, torque_body)")) do
                @test rejects(run_gate(), "wrench: no method takes 3 positional arguments")
            end
            with_registry(in_row("wrench", wrench_signature => "wrench(m, x, env, t) -> (force_ii, torque_body)")) do
                @test rejects(run_gate(), "wrench: the signature names no dispatch argument")
            end
            with_registry(in_row("environment_requirements", "fallback = \"returns no sampled fields\"" => "fallback = \"none\"")) do
                @test rejects(run_gate(), "environment_requirements: fallback cell says none, but")
            end
            with_registry(in_row("wrench", "fallback = \"none\"" => "fallback = \"returns `0.0`\"")) do
                @test rejects(run_gate(), "wrench: fallback cell says \"returns `0.0`\", but no method with `model` untyped or typed with AbstractForceTorqueModel or a supertype takes 4 positional arguments.")
            end
            with_registry(in_row("environment_requirements", "fallback = \"returns no sampled fields\"" => "fallback = \"defaults\"")) do
                @test rejects(run_gate(), "environment_requirements: fallback cell must be none")
            end
            # Rule 5: deprecated rows whose docstring does not say so, with a stray index-word gap,
            # and with a missing replacement and the wrong removal release.
            with_registry(in_row("wrench", "hook(:wrench; shape = \"R\"," =>
                              "hook(:wrench; shape = \"R\", status = \"deprecated\", replacement = :environment_requirements, removed_in = \"0.5.0\","),
                          gaps_head => gaps_head * "    (symbol = :wrench, rule = :index_word, reason = \"probe\"),\n") do
                rejected = run_gate()
                @test rejects(rejected, "wrench: the docstring's first paragraph does not open with \"Deprecated\".")
                @test rejects(rejected, "wrench: the docstring's first paragraph does not name the replacement environment_requirements.")
                @test rejects(rejected, "wrench: the docstring does not give the removal release 0.5.0.")
                @test rejects(rejected, "wrench: deprecated rows need no :index_word entry in KNOWN_GAPS.")
            end
            with_registry(in_row("wrench", "hook(:wrench; shape = \"R\"," =>
                              "hook(:wrench; shape = \"R\", status = \"deprecated\", replacement = :wrenches, removed_in = \"0.4.0\",")) do
                rejected = run_gate()
                @test rejects(rejected, "wrench: a deprecated row must name a replacement that is another registry row.")
                @test rejects(rejected, "wrench: removed_in must be 0.5.0.")
            end
            # Rule 6: a missing page, a page that does not name the hook, no page without a gap,
            # an unknown gap rule, and a deleted page.
            with_registry(in_row("touchdown_spec", "page = \"user/stop_conditions.md\"" => "page = \"user/missing_page.md\"")) do
                @test rejects(run_gate(), "touchdown_spec: page docs/src/user/missing_page.md does not exist")
            end
            with_registry(in_row("environment_requirements", "page = \"user/custom_effector.md\"" => "page = \"publications.md\"")) do
                @test rejects(run_gate(), "environment_requirements: page docs/src/publications.md does not name the hook.")
            end
            with_registry(in_row("wrench", "page = \"user/custom_effector.md\"" => "page = nothing"),
                          gaps_head => gaps_head * "    (symbol = :wrench, rule = :unknown, reason = \"probe\"),\n") do
                rejected = run_gate()
                @test rejects(rejected, "wrench: no page; name the page that explains it.")
                @test rejects(rejected, "KNOWN_GAPS entry wrench: unknown rule :unknown.")
            end
            actuator_page = joinpath(fixture, "docs", "src", "user", "custom_actuator.md")
            actuator_text = read(actuator_page, String)
            rm(actuator_page)
            try
                @test rejects(run_gate(), "actuator_wrench: page docs/src/user/custom_actuator.md does not exist")
            finally
                write(actuator_page, actuator_text)
            end
            # Rule 7: `sat_idx` in a current row, and an index-word gap that no longer applies.
            with_registry(in_row("touchdown_spec", "touchdown_spec(effector, spacecraft_index::Int)" => "touchdown_spec(effector, sat_idx::Int)")) do
                @test rejects(run_gate(), "touchdown_spec: a current row names the spacecraft index sat_idx; use spacecraft_index.")
            end
            with_registry(gaps_head => gaps_head * "    (symbol = :touchdown_spec, rule = :index_word, reason = \"probe\"),\n") do
                @test rejects(run_gate(), "touchdown_spec: the signature now uses spacecraft_index; remove its :index_word entry from KNOWN_GAPS.")
            end
            # A facade page that does not exist.
            with_registry("page = \"reference/hypr_services.md\"" => "page = \"reference/missing_services.md\"") do
                @test rejects(run_gate(), "HYPRServices: page docs/src/reference/missing_services.md does not exist or does not name the facade.")
            end
            # Package faults that no registry edit can express. A child process changes
            # SpaceAGORA in memory and runs the gate on this fixture, so this process's
            # package is never touched:
            # - status: environment_requirements's docstring opens with a deprecation while
            #   its row still says current;
            # - arity: environment_requirements loses its one-argument fallback and gains an
            #   untyped two-argument one, so a type with no method of its own, called with
            #   the documented one argument, has no fallback;
            # - foreign: wrench gains a method typed with AbstractRPOPlanner, which no force
            #   model reaches, so it is no fallback: a `none` row passes and a `returns` row
            #   fails;
            # - other_arity / documented_arity: an untyped wrench at an arity the row does
            #   not document leaves a `none` row passing; one at the documented arity fails it.
            # The hooks are ones the hook design keeps unchanged, so later PRs leave them alone.
            mutations = joinpath(fixture, "hook_gate_mutations.jl")
            write(mutations, raw"""
                using SpaceAGORA
                fixture, gate, root_line = ARGS
                source = replace(read(gate, String), root_line => "const REPO_ROOT = " * repr(fixture); count=1)
                function report(name)
                    m = Module(:MutatedHookGate)
                    failed = try
                        redirect_stdout(devnull) do
                            Base.include_string(m, source, gate)
                        end
                        false
                    catch
                        true
                    end
                    println(name, " failed=", failed)
                    foreach(v -> println(name, " violation: ", v), Core.eval(m, :violations))
                end
                hook = SpaceAGORA.environment_requirements
                owner = parentmodule(hook)
                saved = []
                for mod in Base.Docs.modules, (binding, multidoc) in Base.Docs.meta(mod)
                    binding.var === :environment_requirements || continue
                    for d in values(multidoc.docs)
                        text = join(string.(d.text))
                        startswith(lstrip(text), "environment_requirements(") || continue
                        push!(saved, (d, d.text, d.object))
                        d.text = Core.svec(replace(text, "\n\n" => "\n\nDeprecated: use `wrench` instead; `environment_requirements` will be removed in SpaceAGORA 0.5.0.\n\n"; count=1))
                        d.object = nothing
                    end
                end
                isempty(saved) && error("environment_requirements's docstring was not found")
                report("status")
                for (d, text, object) in saved
                    d.text, d.object = text, object
                end
                registry = joinpath(fixture, "docs", "extension_hooks.jl")
                registry_text = read(registry, String)
                wrench_row = match(r"        hook\(:wrench;.*?\),\n"s, registry_text).match
                function with_wrench_methods(f, definition)
                    before = Set(methods(SpaceAGORA.wrench))
                    Core.eval(parentmodule(SpaceAGORA.wrench), definition)
                    try
                        f()
                    finally
                        foreach(Base.delete_method, setdiff(Set(methods(SpaceAGORA.wrench)), before))
                    end
                end
                with_wrench_methods(:(wrench(model::$(SpaceAGORA.AbstractRPOPlanner), x, env, t) = (nothing, nothing))) do
                    report("foreign_none")
                    write(registry, replace(registry_text, wrench_row =>
                        replace(wrench_row, "fallback = " * repr("none") => "fallback = " * repr("returns a zero wrench"); count=1); count=1))
                    try
                        report("foreign_returns")
                    finally
                        write(registry, registry_text)
                    end
                end
                with_wrench_methods(:(wrench(model, x) = (nothing, nothing))) do
                    report("other_arity")
                end
                with_wrench_methods(:(wrench(model, x, env, t) = (nothing, nothing))) do
                    report("documented_arity")
                end
                generic = only(m for m in methods(hook) if m.module === owner && m.nargs == 2 &&
                    all(T -> T === Any, Base.unwrap_unionall(m.sig).parameters[2:end]))
                Base.delete_method(generic)
                Core.eval(owner, :(environment_requirements(model, extra) = nothing))
                report("arity")
                """)
            gate_path = joinpath(repo, "test", "gates", "ci_extension_hook_reference_gate.jl")
            child = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $mutations $fixture $gate_path $_GATE_ROOT_LINE`
            lines = split(read(child, String), '\n')
            @test "status failed=true" in lines
            @test "status violation: environment_requirements: the docstring's first paragraph opens with \"Deprecated\", but the row says current." in lines
            @test count(startswith("status violation: "), lines) == 1
            @test "arity failed=true" in lines
            @test "arity violation: environment_requirements: fallback cell says \"returns no sampled fields\", but no method with `model` untyped or typed with AbstractForceTorqueModel or a supertype takes 1 positional arguments." in lines
            @test count(startswith("arity violation: "), lines) == 1
            @test "foreign_none failed=false" in lines
            @test "foreign_returns failed=true" in lines
            @test "foreign_returns violation: wrench: fallback cell says \"returns a zero wrench\", but no method with `model` untyped or typed with AbstractForceTorqueModel or a supertype takes 4 positional arguments." in lines
            @test count(startswith("foreign_returns violation: "), lines) == 1
            @test "other_arity failed=false" in lines
            @test "documented_arity failed=true" in lines
            @test any(startswith("documented_arity violation: wrench: fallback cell says none, but wrench(model, x, env, t)"), lines)
            @test count(startswith("documented_arity violation: "), lines) == 1
            rm(mutations)
            @test isnothing(run_gate().failure)
        end

        @testset "Typed command boundary gate" begin
            run_gate = _gate_runner(repo, fixture, "ci_gnc_typed_command_boundary_gate.jl")
            @test isnothing(run_gate().failure)

            # The control owner must import the typed command from CommandTypes.
            hooks = _unique_source(fixture, "control/hooks.jl")
            hooks_source = read(hooks, String)
            needle = "using ..CommandTypes: PropulsiveManeuverCommand"
            @test occursin(needle, hooks_source)
            _with_text(hooks, replace(hooks_source, needle => "using ..CommandTypes: AerobrakingControlCommand"; count=1)) do
                @test _failed_with(run_gate(),
                    "Typed-command boundary violation: ControlHooks is not importing PropulsiveManeuverCommand from CommandTypes.")
            end
            # maneuver_commands may only be touched by the typed writer and consumer.
            guidance_dir = dirname(_unique_source(fixture, "guidance/hooks.jl"))
            probe = joinpath(guidance_dir, "ownership_probe.jl")
            _with_probe(probe, "x = p.shared_buffers.maneuver_commands\n") do
                rejected = run_gate()
                @test _failed_with(rejected, "GNC typed-command boundary gate failed")
                @test _reports(rejected, probe, fixture, "maneuver_commands ownership is restricted")
            end
            @test isnothing(run_gate().failure)
        end

        @testset "Removed hook names gate" begin
            run_gate = _gate_runner(repo, fixture, "ci_removed_hook_names_gate.jl")
            @test isnothing(run_gate().failure)
            message = "Removed hook names found in src/, ext/, templates/ or examples/."
            probe = joinpath(dirname(_unique_source(fixture, "gravity_gradient_effector.jl")), "removed_hook_probe.jl")
            # A method of a removed hook is refused in every spelling.
            for text in ("calcForceTorque(model::Probe, x, p, i) = (0.0, 0.0)\n",
                         "function SimulationModel.calcForceTorque(model::Probe, x::AbstractVector, p, i::Int)\n    return 0.0, 0.0\nend\n",
                         "@inline SpaceAGORA.calcForceTorque(::Probe, x, p, i)::Tuple = (0.0, 0.0)\n",
                         "(f::typeof(calcForceTorque))(model::Probe, x, p, i) = (0.0, 0.0)\n",
                         "(::typeof(SimulationModel.calcForceTorque))(model::Probe, x, p, i) = (0.0, 0.0)\n",
                         "calcForceTorque(model, x, p, i) where {} = (0.0, 0.0)\n")
                _with_probe(probe, text) do
                    rejected = run_gate()
                    @test _failed_with(rejected, message)
                    @test _reports(rejected, probe, fixture, "defines a method of the removed hook calcForceTorque")
                end
            end
            # Every removed name, in a method, a declaration and the other uses.
            for (hook, text, kind) in (
                    ("calcControlForceTorque", "SpaceAGORA.calcControlForceTorque(model::Probe, u, p, i, t) = (0.0, 0.0)\n", "defines a method of"),
                    ("calcControlMassFlowRate", "function ControlHooks.calcControlMassFlowRate(model, u, p, i, t)\n    return -1.0\nend\n", "defines a method of"),
                    ("calcReactionWheelTorque", "(::typeof(calcReactionWheelTorque))(model::Probe, u, p, i, t) = nothing\n", "defines a method of"),
                    ("calcControlEffect!", "function calcControlEffect! end\n", "declares"),
                    ("calcGuidanceEffect!", "\"\"\"\n    calcGuidanceEffect!(model, u, p, t, i)\n\"\"\"\nfunction calcGuidanceEffect! end\n", "declares"),
                    ("calcNavigationEffect!", "tick!(m, u, p, t, i) = SpaceAGORA.calcNavigationEffect!(m, u, p, t, i)\n", "uses"),
                    ("planet_frame_lpi", "using ..EphemeridesModels: planet_frame_lpi\n", "uses"),
                    ("ephemerides_sun_direction_ii", "export ephemerides_sun_direction_ii\n", "uses"),
                    ("ephemerides_time_seconds", "hook = :ephemerides_time_seconds\n", "uses"))
                _with_probe(probe, text) do
                    rejected = run_gate()
                    @test _failed_with(rejected, message)
                    @test _reports(rejected, probe, fixture, "$(kind) the removed hook $(hook)")
                end
            end
            # A five-argument wrench is the removed plume form, as a method or a call;
            # the four-argument hook is not. In a signature a defaulted parameter is
            # positional and makes one method per arity, as a trailing `args...` does.
            for text in ("SpaceAGORA.wrench(model::Probe, x::StateSample, env::EnvironmentSample, t::Float64, i::Int) = (0.0, 0.0)\n",
                         "function wrench(model::Probe, x, env, t, i; extra=nothing)\n    return 0.0, 0.0\nend\n",
                         "wrench(model::Probe, x, env, t, i=1) = (0.0, 0.0)\n",
                         "function wrench(model::Probe, x, env, t, i, j=2)\n    return 0.0, 0.0\nend\n",
                         "wrench(model::Probe, x, env, t, rest...) = (0.0, 0.0)\n")
                _with_probe(probe, text) do
                    rejected = run_gate()
                    @test _failed_with(rejected, message)
                    @test _reports(rejected, probe, fixture, "defines a five-argument method of wrench")
                    @test !_reports(rejected, probe, fixture, "calls the five-argument form")
                end
            end
            # In a call a keyword is not positional, after a semicolon or a comma.
            for text in ("f, τ = wrench(model, x, env, t, 2)\n",
                         "SimulationModel.wrench(model, x, env, t, i; extra=1)\n",
                         "SimulationModel.wrench(model, x, env, t, i, extra=1)\n")
                _with_probe(probe, text) do
                    rejected = run_gate()
                    @test _failed_with(rejected, message)
                    @test _reports(rejected, probe, fixture, "calls the five-argument form of wrench")
                end
            end
            _with_probe(probe, "f, τ = wrench(model, x, env, t)\nf, τ = wrench(model, x, env, t; i=2)\n" *
                               "f, τ = wrench(model, x, env, t, option=true)\nargs = (1, 2, 3, 4, 5)\nwrench(args...)\n" *
                               "wrench(model::Probe, x, env, t; i=1) = (0.0, 0.0)\n" *
                               "wrench(model::Probe, x, env, t=0.0) = (0.0, 0.0)\n" *
                               "function wrench(model::Probe, x, env, t, i, j, k=3)\n    return 0.0, 0.0\nend\n") do
                @test isnothing(run_gate().failure)
            end
            # Shipped extension code under ext/ is scanned too.
            extension = joinpath(fixture, "ext", "SpaceAGORAProbeExt", "probe.jl")
            _with_probe(extension, "SpaceAGORA.calcGuidanceEffect!(m::Probe, u, p, t, i) = nothing\n") do
                rejected = run_gate()
                @test _failed_with(rejected, message)
                @test _reports(rejected, extension, fixture, "defines a method of the removed hook calcGuidanceEffect!")
            end
            rm(joinpath(fixture, "ext"); recursive=true, force=true)
            # Comments, strings and longer names are not uses, and the current hooks are fine.
            _with_probe(probe, "# calcForceTorque was removed in 0.5.0.\nnote = \"calcControlEffect! is now update_control!\"\n" *
                               "my_calcForceTorque(x) = x\nSpaceAGORA.wrench(model::Probe, x, env, t::Float64) = (0.0, 0.0)\n") do
                @test isnothing(run_gate().failure)
            end
            # A file that does not parse is reported.
            _with_probe(probe, "function broken(\n") do
                rejected = run_gate()
                @test _failed_with(rejected, message)
                @test _reports(rejected, probe, fixture, "does not parse")
            end
            # The starters name no removed hook anywhere, in any file type; a longer
            # name that only contains one is not a mention.
            starter = joinpath(fixture, "templates", "starters", "package")
            for (file, text) in (("README.md", "Implement `calcForceTorque` here.\n"),
                                 (joinpath("test", "runtests.jl"), "# x = SpaceAGORA.calcControlMassFlowRate\n"),
                                 (joinpath(".github", "workflows", "ci.yml"), "# calcReactionWheelTorque\n"),
                                 (joinpath("src", "tick_notes.jl"), "# The tick used to be calcGuidanceEffect!.\n"))
                probe = joinpath(starter, file)
                _with_probe(probe, text) do
                    rejected = run_gate()
                    @test _failed_with(rejected, message)
                    @test _reports(rejected, probe, fixture, "names the removed hook")
                end
            end
            _with_probe(joinpath(starter, "README.md"), "my_calcForceTorque_notes and calcForceTorques\n") do
                @test isnothing(run_gate().failure)
            end
            rm(joinpath(fixture, "templates"); recursive=true)
            @test isnothing(run_gate().failure)
        end
    end
end

@testset "Random-stream gate distinguishes in-place destinations" begin
    repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
    mktempdir() do fixture
        _copy_into!(fixture, repo, joinpath("src", "models", "gnc"))
        run_gate = _gate_runner(repo, fixture, "ci_gnc_random_stream_gate.jl")
        @test isnothing(run_gate().failure)
        probe = joinpath(fixture, "src", "models", "gnc", "random_stream_probe.jl")
        message = "draw without an explicit stream as its first argument"
        failure = "GNC random-stream gate failed: take an explicit rng argument and pass it to every draw"

        for draw in ("rand!", "randn!", "randexp!", "shuffle!", "randperm!", "randcycle!")
            @testset "$draw" begin
                for source in (
                    "rng_buffer = zeros(3); $draw(rng_buffer)\n",
                    "Random.$draw(\n    rng_buffer,\n)\n",
                    "$draw(view(rng_buffer, :, 1))\n",
                    "$draw((rng_buffer,))\n",
                    "$draw(rng_buffer; keyword=rng)\n",
                    "$draw(rng_buffer, keyword=rng)\n",
                    "$draw(rng_buffer, rest...)\n",
                    "$draw(rng..., buffer)\n",
                )
                    _with_probe(probe, source) do
                        rejected = run_gate()
                        @test _failed_with(rejected, failure)
                        @test _reports(rejected, probe, fixture, message)
                    end
                end
                for source in (
                    "$draw(rng, buffer)\n",
                    "Random.$draw(\n    rng, # caller's stream\n    view(buffer, :, 1),\n)\n",
                    "$draw(state.master_rng, buffer)\n",
                )
                    _with_probe(probe, source) do
                        @test isnothing(run_gate().failure)
                    end
                end
            end
        end
        # The permutation APIs have exactly two explicit-stream arguments.
        for draw in ("randperm!", "randcycle!"), source in (
            "$draw(rng, buffer, extra)\n",
            "$draw(rng, buffer, rest...)\n",
        )
            _with_probe(probe, source) do
                rejected = run_gate()
                @test _failed_with(rejected, failure)
                @test _reports(rejected, probe, fixture, message)
            end
        end
        @testset "randsubseq!" begin
            # A destination named rng_output is not the first of four arguments.
            for source in (
                "randsubseq!(rng_output, values, p)\n",
                "Random.randsubseq!(\n    rng_output,\n    view(values, :, 1),\n    p,\n)\n",
                "randsubseq!(rng, output)\n",
                "randsubseq!(rng_output, values, p; keyword=rng)\n",
                "randsubseq!(rng, output, values, p=p)\n",
                "randsubseq!(rng..., output, values, p)\n",
                "randsubseq!(rng, outputs..., values, p)\n",
                "randsubseq!(rng, output, values..., p)\n",
                "randsubseq!(rng, output, values, probabilities...)\n",
                "randsubseq!(rng, output, values, p, extra)\n",
            )
                _with_probe(probe, source) do
                    rejected = run_gate()
                    @test _failed_with(rejected, failure)
                    @test _reports(rejected, probe, fixture, message)
                end
            end
            for source in (
                "randsubseq!(rng, output, values, p)\n",
                "Random.randsubseq!(\n    rng, # caller's stream\n    output,\n    view(values, :, 1),\n    p,\n)\n",
                "randsubseq!(state.master_rng, output, values, p)\n",
            )
                _with_probe(probe, source) do
                    @test isnothing(run_gate().failure)
                end
            end
        end
        # Existing variable-arity forms keep optional arguments after rng and A.
        for draw in ("rand!", "randn!", "randexp!", "shuffle!"), source in (
            "$draw(rng, buffer, extra)\n",
            "$draw(rng, buffer, rest...)\n",
        )
            _with_probe(probe, source) do
                @test isnothing(run_gate().failure)
            end
        end
        # Non-mutating calls still require the named first stream argument.
        for source in ("rand()\n", "randn(3)\n", "randexp(Float64)\n", "rand!(buffer, rng)\n")
            _with_probe(probe, source) do
                rejected = run_gate()
                @test _failed_with(rejected, failure)
                @test _reports(rejected, probe, fixture, message)
            end
        end
        _with_probe(probe, "rand(rng)\nrandn(master_rng, 3)\nrandexp(state.rng)\n") do
            @test isnothing(run_gate().failure)
        end
        # Comments and string examples are not draw call sites.
        _with_probe(probe, "# rand!(rng_buffer)\nexample = \"rand!(rng_buffer)\"\n") do
            @test isnothing(run_gate().failure)
        end
        # Both roots remain covered, with the draw's source line in the report.
        hypr_probe = joinpath(fixture, "src", "models", "gnc", "companions", "hypr", "random_stream_probe.jl")
        _with_probe(hypr_probe, "rng_buffer = zeros(3)\nrand!(\n    rng_buffer\n)\n") do
            rejected = run_gate()
            @test _failed_with(rejected, failure)
            @test _reports(rejected, hypr_probe, fixture, "random_stream_probe.jl:2: " * message)
        end
        @test isnothing(run_gate().failure)
    end
end
