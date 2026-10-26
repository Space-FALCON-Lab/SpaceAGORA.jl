module PPCModeProfileParityGateTests
using Test

const ROOT_LINE = "const REPO_ROOT = normpath(joinpath(@__DIR__, \"..\", \"..\"))"
const MODES_REL = joinpath("benchmarks", "studies", "parallelization_performance", "modes.jl")
# The gate names its canonical owner explicitly. The fixture locates that one
# source file so it remains usable when an independently reviewed owner move
# updates the gate path; duplicate owners are an error, not an arbitrary pick.
function profile_definitions_rel(repo::String)
    paths = String[]
    for (root, _, files) in walkdir(joinpath(repo, "src"))
        "profile_definitions.jl" in files && push!(paths, joinpath(root, "profile_definitions.jl"))
    end
    length(paths) == 1 || error("Expected one profile_definitions.jl source owner; found $(length(paths))")
    return relpath(only(paths), repo)
end

function replace_once(source::AbstractString, old::AbstractString, new::AbstractString)
    length(findall(old, source)) == 1 || error("Expected one mutation target: $(repr(old))")
    return replace(source, old => new; count=1)
end

function with_text(f, path::String, text::String)
    original = read(path, String)
    write(path, text)
    try
        f()
    finally
        write(path, original)
    end
end

function run_gate(repo::String, fixture::String)
    gate_path = joinpath(repo, "test", "gates", "ci_ppc_mode_profile_parity_gate.jl")
    source = replace_once(read(gate_path, String), ROOT_LINE,
                          "const REPO_ROOT = $(repr(fixture))")
    isolated = Module(gensym(:PPCParityGateControl))
    failure = try
        redirect_stdout(devnull) do
            Base.include_string(isolated, source, gate_path)
        end
        nothing
    catch err
        err
    end
    violations = Core.eval(isolated,
        :(isdefined(@__MODULE__, :violations) ? copy(violations) : String[]))
    return (; failure, violations)
end

function rejects(result, diagnostic::String)
    return result.failure isa LoadError &&
           result.failure.error isa ErrorException &&
           startswith(result.failure.error.msg,
                      "Benchmark mode table disagrees with the shipped profiles:\n") &&
           any(v -> occursin(diagnostic, v), result.violations)
end

function fails_starting(result, diagnostic::String)
    return result.failure isa LoadError &&
           result.failure.error isa ErrorException &&
           startswith(result.failure.error.msg, diagnostic)
end

@testset "PPC mode/profile gate checks actual values" begin
    repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
    mktempdir() do fixture
        profiles_rel = profile_definitions_rel(repo)
        for rel in (MODES_REL, profiles_rel)
            target = joinpath(fixture, rel)
            mkpath(dirname(target))
            cp(joinpath(repo, rel), target)
        end
        modes_path = joinpath(fixture, MODES_REL)
        profiles_path = joinpath(fixture, profiles_rel)
        modes = read(modes_path, String)
        profiles = read(profiles_path, String)

        # Covers R1a/R1b aliases, implicit scheduler defaults, attribution arms,
        # and the deliberate R3 backend difference without changing them.
        @test isnothing(run_gate(repo, fixture).failure)

        @testset "Separately constructed serial mode" begin
            changed = replace_once(modes,
                "        thermal=\"off\",\n        multibody=\"off\",",
                "        thermal=\"on\",\n        multibody=\"off\",")
            with_text(modes_path, changed) do
                @test rejects(run_gate(repo, fixture), "mode \"serial\" sets thermal=\"on\"")
            end
            changed = replace_once(modes,
                "        density=\"off\",\n        control=\"off\",",
                "        density=string(:on),\n        control=\"off\",")
            with_text(modes_path, changed) do
                @test rejects(run_gate(repo, fixture), "mode \"serial\" sets density=\"on\"")
            end
            # Equal expression-valued fields are valid, rather than forbidden syntax.
            changed = replace_once(modes,
                "        density=\"off\",\n        control=\"off\",",
                "        density=string(:off),\n        control=\"off\",")
            with_text(modes_path, changed) do
                @test isnothing(run_gate(repo, fixture).failure)
            end
        end

        @testset "Inline mode and scheduler drift" begin
            # An expression-valued field must be evaluated, not silently
            # omitted because it does not match a quoted-literal pattern.
            changed = replace_once(modes,
                "            name=\"outer_threads\",",
                "            name=\"outer_threads\",\n            scheduler=string(:dynamic),")
            with_text(modes_path, changed) do
                @test rejects(run_gate(repo, fixture), "mode \"outer_threads\" sets scheduler=\"dynamic\"")
            end
            changed = replace_once(modes,
                "            name=\"full_smart\",",
                "            name=\"full_smart\",\n            # Preserved constructor layout.\n")
            with_text(modes_path, changed) do
                @test isnothing(run_gate(repo, fixture).failure)
            end
            full_smart = match(r"        \"full_smart\" => PPCModeSpec\(.*?\n        \),"s, modes).match
            changed = replace_once(modes, full_smart,
                replace_once(full_smart, "thermal=\"auto\"", "thermal=string(:on)"))
            with_text(modes_path, changed) do
                @test rejects(run_gate(repo, fixture), "mode \"full_smart\" sets thermal=\"on\"")
            end
            changed = replace_once(modes, full_smart,
                replace_once(full_smart, "scheduler=\"static\"", "scheduler=\"dynamic\""))
            with_text(modes_path, changed) do
                @test rejects(run_gate(repo, fixture), "mode \"full_smart\" sets scheduler=\"dynamic\"")
            end
            # The serial constructor omits scheduler and must inherit its actual default.
            changed = replace_once(modes, "scheduler::String = \"static\"",
                "scheduler::String = \"dynamic\"")
            with_text(modes_path, changed) do
                @test rejects(run_gate(repo, fixture), "mode \"serial\" sets scheduler=\"dynamic\"")
            end
        end

        @testset "Missing expected mode and profile identity" begin
            changed = replace_once(modes, "\"serial\" => serial", "\"renamed_serial\" => serial")
            with_text(modes_path, changed) do
                @test fails_starting(run_gate(repo, fixture),
                    "Benchmark mode table has no PPCModeSpec for checked mode \"serial\"")
            end
            changed = replace_once(modes, "        profile=\"R0\",", "        profile=\"unknown\",")
            with_text(modes_path, changed) do
                @test rejects(run_gate(repo, fixture), "mode \"serial\" has invalid profile \"unknown\"")
            end
            changed = replace_once(modes, "        profile=\"R0\",", "        profile=\"R5\",")
            with_text(modes_path, changed) do
                @test rejects(run_gate(repo, fixture), "mode \"serial\" names profile \"R5\" but requires R0")
            end
            # R5 and R6 share the checked values but retain distinct identities.
            changed = replace_once(modes,
                "            name=\"full_smart\",\n            profile=\"R5\",",
                "            name=\"full_smart\",\n            profile=\"R6\",")
            with_text(modes_path, changed) do
                @test rejects(run_gate(repo, fixture), "mode \"full_smart\" names profile \"R6\" but requires R5")
            end
        end

        @testset "Missing fields fail closed" begin
            # Rename declarations and constructor keywords together: actual
            # construction succeeds, but the required field is absent.
            for field in ("thermal", "density", "control", "multibody", "effector", "scheduler")
                with_text(modes_path, replace(modes, Regex("\\b" * field * "\\b") => "removed_" * field)) do
                    @test rejects(run_gate(repo, fixture), "mode \"serial\" is missing " * field)
                end
            end
            with_text(modes_path, replace(modes, r"\bprofile\b" => "removed_profile")) do
                @test rejects(run_gate(repo, fixture), "mode \"serial\" is missing profile")
            end
            with_text(profiles_path, replace(profiles, r"\bthermal_mode\b" => "removed_thermal_mode")) do
                @test rejects(run_gate(repo, fixture), "profile R0 is missing thermal_mode")
            end
        end

        @testset "Missing source files fail closed" begin
            for (path, original, diagnostic) in (
                (modes_path, modes, "modes.jl not found"),
                (profiles_path, profiles, "profile_definitions.jl not found"),
            )
                rm(path)
                try
                    @test fails_starting(run_gate(repo, fixture), diagnostic)
                finally
                    write(path, original)
                end
            end
        end

        @test isnothing(run_gate(repo, fixture).failure)
    end
end
end # module
