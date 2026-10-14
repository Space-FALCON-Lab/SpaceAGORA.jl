using Test

@testset "Runtime Any gate covers moved and nested owners" begin
    repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
    gate_path = joinpath(repo, "test", "gates", "ci_runtime_any_hotpath_gate.jl")
    source = read(gate_path, String)
    root_line = "const REPO_ROOT = normpath(joinpath(@__DIR__, \"..\", \"..\"))"
    @test occursin(root_line, source)
    mktempdir() do fixture
        paths = (
            "src/orchestration/simulation/state/runtime_types.jl",
            "src/orchestration/simulation/engine/setup.jl",
            "src/models/dynamics/coupled/aerodynamic_wrench_models.jl",
            "src/models/dynamics/coupled/perturbations.jl",
            "src/orchestration/parallel/policy/multibody_routing.jl",
            "src/models/dynamics/coupled/perturbations/gravity/nested/probe.jl",
            "src/orchestration/simulation/callbacks/event_callbacks.jl",
            "src/orchestration/simulation/callbacks/nested/probe.jl",
        )
        baseline = Dict(relative => "# valid fixture\n" for relative in paths)
        baseline[paths[1]] = "walk_models::Vector{Any}\nexternal_runtimes::Vector{Any} = Any[]\n"
        baseline[paths[7]] = join((
            "storage::Dict{Symbol, Any}", "fallback::Dict{Symbol, Vector{Any}}",
            "storage = Dict{Symbol, Any}()", "fallback = Dict{Symbol, Vector{Any}}()",
            "fallback[field.name] = Vector{Any}(undef, cap)", "grown = Vector{Any}(undef, new_cap)",
        ), "\n") * "\n"
        for relative in paths
            path = joinpath(fixture, relative)
            mkpath(dirname(path))
            write(path, baseline[relative])
        end
        isolated_source = replace(source, root_line => "const REPO_ROOT = $(repr(fixture))"; count=1)
        function run_gate()
            target = Module(gensym(:RuntimeAnyGate))
            failure = try
                redirect_stdout(devnull) do
                    Base.include_string(target, isolated_source, gate_path)
                end
                nothing
            catch error
                error
            end
            violations = Core.eval(target,
                :(isdefined(@__MODULE__, :violations) ? copy(violations) : String[]))
            return (; failure, violations)
        end

        @test run_gate().failure === nothing
        # A closing brace followed by a delimiter, whitespace or newline is
        # not a regex word boundary. Test constructors and typed declarations.
        for relative in paths, forbidden in (
            "probe = Vector{Any}()\n",
            "probe::Vector{Any}\n",
            "probe = Dict{Symbol, Any}()\n",
            "probe::Dict{Symbol,Any}\n",
        )
            path = joinpath(fixture, relative)
            write(path, baseline[relative] * forbidden)
            result = run_gate()
            @test result.failure isa LoadError
            @test any(v -> startswith(v, relative * ":"), result.violations)
            write(path, baseline[relative])
            @test run_gate().failure === nothing
        end

        # Exact baseline entries must neither disappear nor admit duplicates.
        for relative in (paths[1], paths[7]), statement in split(chomp(baseline[relative]), '\n')
            path = joinpath(fixture, relative)
            write(path, baseline[relative] * statement * "\n")
            result = run_gate()
            @test result.failure isa LoadError
            @test any(v -> occursin("expected 1 occurrence(s), observed 2", v), result.violations)
            write(path, baseline[relative])
        end
        state = joinpath(fixture, first(paths))
        write(state, "external_runtimes::Vector{Any} = Any[]\n")
        missing = run_gate()
        @test missing.failure isa LoadError
        @test any(v -> occursin("expected 1 occurrence(s), observed 0", v), missing.violations)
        state_baseline = baseline[first(paths)]
        write(state, state_baseline * "const SaveData = Dict{Symbol, Any}\n")
        @test run_gate().failure === nothing
        # The one documented SaveData alias must not exempt another Any field.
        write(state, state_baseline * "const SaveData = Dict{Symbol, Any}\nprobe::Dict{Symbol,Any}\n")
        @test run_gate().failure isa LoadError
        write(state, state_baseline * "probe = Vector{Float64}()\nother = Dict{Symbol, Float64}()\n# ignored = Vector{Any}()\n")
        @test run_gate().failure === nothing
    end
end
