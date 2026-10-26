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
            "src/models/dynamics/effectors/aerodynamic_wrench_models.jl",
            "src/models/dynamics/effectors/perturbations.jl",
            "src/orchestration/parallel/policy/multibody_routing.jl",
            "src/models/dynamics/effectors/perturbations/gravity/nested/probe.jl",
            "src/orchestration/simulation/callbacks/event_callbacks.jl",
            "src/orchestration/simulation/callbacks/nested/probe.jl",
        )
        baseline = Dict(relative => "# valid fixture\n" for relative in paths)
        baseline[paths[1]] = """
            module ConfigTypes
            mutable struct GramDensityPerturbationState
                walk_models::Vector{Any}
            end
            @kwdef struct SharedBuffers
                external_runtimes::Vector{Any} = Any[]
            end
            const SaveData = Dict{Symbol, Any}
            end
            """
        baseline[paths[7]] = """
            mutable struct TrajectoryRecorder
                storage::Dict{Symbol, Any}
                fallback::Dict{Symbol, Vector{Any}}
            end
            function TrajectoryRecorder(args)
                storage = Dict{Symbol, Any}()
                fallback = Dict{Symbol, Vector{Any}}()
                for field in args
                    fallback[field.name] = Vector{Any}(undef, cap)
                end
            end
            function _grow_trajectory_recorder!(rec::TrajectoryRecorder)::Nothing
                for (name, data) in rec.fallback
                    grown = Vector{Any}(undef, new_cap)
                end
            end
            """
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
        # Duplicate within the original owner; flattening a field into top-level
        # text would test an owner violation instead of the occurrence bound.
        boundary_statements = (
            (paths[1], "walk_models::Vector{Any}"),
            (paths[1], "external_runtimes::Vector{Any} = Any[]"),
            (paths[1], "const SaveData = Dict{Symbol, Any}"),
            (paths[7], "storage::Dict{Symbol, Any}"),
            (paths[7], "fallback::Dict{Symbol, Vector{Any}}"),
            (paths[7], "storage = Dict{Symbol, Any}()"),
            (paths[7], "fallback = Dict{Symbol, Vector{Any}}()"),
            (paths[7], "fallback[field.name] = Vector{Any}(undef, cap)"),
            (paths[7], "grown = Vector{Any}(undef, new_cap)"),
        )
        for (relative, statement) in boundary_statements
            path = joinpath(fixture, relative)
            write(path, replace(baseline[relative], statement => statement * "; " * statement; count=1))
            result = run_gate()
            @test result.failure isa LoadError
            @test any(v -> occursin("expected 1 occurrence(s), observed 2", v), result.violations)
            write(path, baseline[relative])
        end
        state = joinpath(fixture, first(paths))
        state_baseline = baseline[first(paths)]
        write(state, replace(state_baseline, "walk_models::Vector{Any}" => ""; count=1))
        missing = run_gate()
        @test missing.failure isa LoadError
        @test any(v -> occursin("expected 1 occurrence(s), observed 0", v), missing.violations)
        write(state, state_baseline)
        @test run_gate().failure === nothing
        # SaveData belongs to ConfigTypes, not a second top-level alias in the
        # same file, and never exempts an unrelated heterogeneous field.
        write(state, state_baseline * "const SaveData = Dict{Symbol, Any}\n")
        @test run_gate().failure isa LoadError
        write(state, state_baseline * "probe::Dict{Symbol,Any}\n")
        @test run_gate().failure isa LoadError
        write(state, state_baseline * "probe = Vector{Float64}()\nother = Dict{Symbol, Float64}()\n# ignored = Vector{Any}()\n")
        @test run_gate().failure === nothing
    end
end
