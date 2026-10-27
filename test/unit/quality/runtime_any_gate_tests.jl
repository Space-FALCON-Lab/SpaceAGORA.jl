using Test
if !isdefined(@__MODULE__, :RuntimeAnyGateFixtures)
    include(joinpath(@__DIR__, "..", "..", "support", "runtime_any_gate_fixtures.jl"))
end
using .RuntimeAnyGateFixtures

@testset "Runtime Any gate covers moved and nested owners" begin
    repo = normpath(joinpath(@__DIR__, "..", "..", ".."))
    gate_path = joinpath(repo, "test", "gates", "ci_runtime_any_hotpath_gate.jl")
    source = read(gate_path, String)
    root_line = "const REPO_ROOT = normpath(joinpath(@__DIR__, \"..\", \"..\"))"
    @test occursin(root_line, source)
    mktempdir() do fixture
        p = _rah_fixture(fixture)
        nested_gravity = joinpath(p.split, "gravity", "nested", "probe.jl")
        nested_callback = joinpath(p.callbacks, "nested", "probe.jl")
        for path in (nested_gravity, nested_callback)
            mkpath(dirname(path))
            write(path, "# valid fixture\n")
        end
        paths = [_rah_relative_path(path, fixture) for path in
            (p.files..., nested_gravity, p.event, nested_callback)]
        baseline = Dict(relative => read(joinpath(fixture, relative), String) for relative in paths)
        run_gate() = _rah_gate_result(fixture)

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
