using Test

const _RAH_REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(_RAH_REPO, "test", "gates", "runtime_any_hotpath_checks.jl"))

# Follow the actual scan owners across reviewed source moves; ambiguous names
# fail rather than silently selecting a different owner.
const _RAH_SOURCE_NAMES = ("runtime_types.jl", "setup.jl", "aerodynamic_wrench_models.jl",
    "perturbations.jl", "event_callbacks.jl", "multibody_routing.jl")
const _RAH_RELATIVE_PATHS = Dict(name => only([
    relpath(joinpath(dir, name), _RAH_REPO)
    for (dir, _, files) in walkdir(joinpath(_RAH_REPO, "src")) if name in files
]) for name in _RAH_SOURCE_NAMES)

_rah_relative_path(path, root) = replace(relpath(path, root), '\\' => '/')

function _rah_paths(root)
    files = [joinpath(root, _RAH_RELATIVE_PATHS[name]) for name in _RAH_SOURCE_NAMES[1:4]]
    push!(files, joinpath(root, _RAH_RELATIVE_PATHS["multibody_routing.jl"]))
    event = joinpath(root, _RAH_RELATIVE_PATHS["event_callbacks.jl"])
    split = joinpath(dirname(files[4]), "perturbations")
    return (; runtime=files[1], setup=files[2], files, callbacks=dirname(event), event, split)
end
function _rah_scan(root)
    p = _rah_paths(root)
    isdir(p.split) || error("Missing runtime hot-path scan root: $(relpath(p.split, root))")
    files = copy(p.files)
    for (dir, _, names) in walkdir(p.split), name in names
        endswith(name, ".jl") && push!(files, joinpath(dir, name))
    end
    return RuntimeAnyHotpathChecks.scan(root, files, p.callbacks)
end

function _rah_fixture(root)
    p = _rah_paths(root)
    mkpath(p.split)
    for file in (p.files..., p.event)
        mkpath(dirname(file))
        write(file, "# parse-only fixture\n")
    end
    write(p.runtime, """
        module ConfigTypes
        mutable struct GramDensityPerturbationState
            walk_models::Vector{Any}
        end
        @kwdef struct SharedBuffers
            external_runtimes::Vector{Any} = Any[]
        end
        const SaveData = Dict{Symbol, Any}
        end
        """)
    write(p.event, """
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
        """)
    return p
end

function _rah_with_source(f, path, source)
    original = isfile(path) ? read(path, String) : nothing
    mkpath(dirname(path))
    write(path, source)
    try
        f()
    finally
        original === nothing ? rm(path) : write(path, original)
    end
end

function _rah_actual_gate(root)
    gate = joinpath(_RAH_REPO, "test", "gates", "ci_runtime_any_hotpath_gate.jl")
    needle = "const REPO_ROOT = normpath(joinpath(@__DIR__, \"..\", \"..\"))"
    source = read(gate, String)
    @assert occursin(needle, source)
    source = replace(source, needle => "const REPO_ROOT = $(repr(root))"; count=1)
    mod = Module(gensym(:RuntimeAnyGate))
    Core.eval(mod, :(include(path) = Base.include(@__MODULE__, path)))
    try
        redirect_stdout(devnull) do
            Base.include_string(mod, source, gate)
        end
        return nothing
    catch err
        return err
    end
end

@testset "Runtime Any explicit-type syntax fence" begin
    # Exercise real source and the actual gate entrypoint without loading packages.
    @test isempty(_rah_scan(_RAH_REPO))
    @test isnothing(_rah_actual_gate(_RAH_REPO))
    mktempdir() do root
        p = _rah_fixture(root)
        @test isempty(_rah_scan(root))
        @test isnothing(_rah_actual_gate(root))

        @testset "New explicit types are rejected" begin
            for source in (
                "x = Vector{Any}()",
                "x::Vector{Any}",
                "x = Dict{Symbol, Any}()",
                "x::Dict{Symbol,Any}",
                "x = Dict{Symbol, Vector{Any}}()",
                "x = Vector{\n Any\n}()",
                "x = Vector{ Any }()",
                "x = Dict{\n Symbol,\n Any\n}()",
                "x = Base.Vector{Core.Any}()",
                "x = Base.Dict{Core.Symbol,Core.Any}()",
                "f(x::Vector{Any}) = x",
                "function f(x::Dict{Symbol, Any}); x; end",
                "struct Unsafe <: AbstractVector{Vector{Any}}; end",
            )
                _rah_with_source(p.setup, source) do
                    found = _rah_scan(root)
                    @test !isempty(found)
                    @test all(v -> startswith(v, _rah_relative_path(p.setup, root) * ":"), found)
                    @test _rah_actual_gate(root) isa LoadError
                end
            end
            _rah_with_source(joinpath(p.callbacks, "nested", "new_callback.jl"), "x = Vector{Any}()") do
                @test any(v -> occursin("new_callback.jl", v), _rah_scan(root))
            end
            for path in (joinpath(p.split, "gravity", "nested", "probe.jl"), last(p.files))
                _rah_with_source(path, "x = Base.Vector{Core.Any}()") do
                    @test any(v -> startswith(v, _rah_relative_path(path, root) * ":"), _rah_scan(root))
                    @test _rah_actual_gate(root) isa LoadError
                end
            end
            _rah_with_source(p.setup, "x = Vector{Any}()") do
                failure = _rah_actual_gate(root)
                @test failure isa LoadError
                @test failure.error isa ErrorException
                @test failure.error.msg == "Runtime Any hot-path gate failed"
            end
        end

        @testset "Exceptions have path, owner, statement and occurrence bounds" begin
            original = read(p.event, String)
            # Formatting is not a different use.
            _rah_with_source(p.event, replace(original,
                "storage::Dict{Symbol, Any}" => "storage :: Dict{Symbol,\n Any}")) do
                @test isempty(_rah_scan(root))
            end
            # The same complete declaration in another file or enclosing owner fails.
            _rah_with_source(p.setup, read(p.runtime, String)) do
                @test !isempty(_rah_scan(root))
            end
            _rah_with_source(p.event, replace(original,
                "mutable struct TrajectoryRecorder" => "mutable struct OtherRecorder")) do
                @test !isempty(_rah_scan(root))
            end
            _rah_with_source(p.event, replace(original,
                "function TrajectoryRecorder(args)" => "function OtherRecorder(args)")) do
                @test !isempty(_rah_scan(root))
            end
            _rah_with_source(p.event, replace(original,
                "storage = Dict{Symbol, Any}()" =>
                "storage = Dict{Symbol, Any}(); storage = Dict{Symbol, Any}()")) do
                @test !isempty(_rah_scan(root))
            end
            # Keep the allowed prefix, but add another Any use on the same line.
            _rah_with_source(p.event, replace(original,
                "storage = Dict{Symbol, Any}()" =>
                "storage = Dict{Symbol, Any}(); extra = Vector{Any}()")) do
                @test !isempty(_rah_scan(root))
            end
            # Exception matching applies to complete statements, never operands.
            for wrapper in (
                "extra = (storage = Dict{Symbol, Any}())",
                "consume((storage = Dict{Symbol, Any}()))",
                "consume(1, (storage = Dict{Symbol, Any}()))",
                "identity(begin; storage = Dict{Symbol, Any}(); end)",
                "extra = begin; storage = Dict{Symbol, Any}(); end",
                "extra = if true; storage = Dict{Symbol, Any}(); end",
                "extra = let; storage = Dict{Symbol, Any}(); end",
                "return (storage = Dict{Symbol, Any}())",
                "quote; storage = Dict{Symbol, Any}(); end",
                ":(begin; storage = Dict{Symbol, Any}(); end)",
                "quote; begin; storage = Dict{Symbol, Any}(); end; end",
            )
                _rah_with_source(p.event, replace(original,
                    "storage = Dict{Symbol, Any}()" => wrapper)) do
                    found = _rah_scan(root)
                    @test any(v -> occursin("forbidden Dict{Symbol, Any}", v), found)
                    @test all(v -> startswith(v, _rah_relative_path(p.event, root) * ":"), found)
                end
            end
            # Ordinary statement blocks retain eligibility, including existing
            # loop bodies and @kwdef fields in the baseline fixture above.
            for statement in (
                "begin; storage = Dict{Symbol, Any}(); end",
                "if true; storage = Dict{Symbol, Any}(); end",
                "let; storage = Dict{Symbol, Any}(); end",
                "try; storage = Dict{Symbol, Any}(); finally; nothing; end",
            )
                _rah_with_source(p.event, replace(original,
                    "storage = Dict{Symbol, Any}()" => statement)) do
                    @test isempty(_rah_scan(root))
                end
            end
            # A larger expression containing the former allowed statement is not exempt.
            _rah_with_source(p.event, replace(original,
                "storage = Dict{Symbol, Any}()" =>
                "storage = Dict{Symbol, Any}(:extra => Vector{Any}())")) do
                @test length(_rah_scan(root)) >= 2
            end
            _rah_with_source(p.event, original * "\nmodule OtherOwner\n" * original * "\nend\n") do
                @test !isempty(_rah_scan(root))
            end
            # An anonymous function is not the permitted constructor owner.
            _rah_with_source(p.event, replace(original,
                "storage = Dict{Symbol, Any}()" =>
                "f = () -> begin; storage = Dict{Symbol, Any}(); end")) do
                @test !isempty(_rah_scan(root))
            end
            # Preserve the current presence/occurrence policy. A legitimate
            # boundary change needs an explicit inventory update, not a bypass.
            _rah_with_source(p.event, "# recorder removed\n") do
                @test any(v -> occursin("expected 1 occurrence(s), observed 0", v), _rah_scan(root))
            end
            _rah_with_source(p.runtime, "module ConfigTypes; const SaveData = Dict{Symbol,Float64}; end") do
                @test any(v -> occursin("expected 1 occurrence(s), observed 0", v), _rah_scan(root))
            end
            runtime = read(p.runtime, String)
            savedata = "const SaveData = Dict{Symbol, Any}"
            for path in (p.setup, p.event, joinpath(p.callbacks, "nested", "misplaced.jl"))
                _rah_with_source(path, (isfile(path) ? read(path, String) : "") * "\n" * savedata) do
                    @test any(v -> occursin("forbidden Dict{Symbol, Any}", v), _rah_scan(root))
                    @test _rah_actual_gate(root) isa LoadError
                end
            end
            for replacement in (
                savedata * "; " * savedata,
                "module Other; " * savedata * "; end",
                "function other(); " * savedata * "; end",
            )
                _rah_with_source(p.runtime, replace(runtime, savedata => replacement)) do
                    @test !isempty(_rah_scan(root))
                    @test _rah_actual_gate(root) isa LoadError
                end
            end
            _rah_with_source(p.runtime, replace(runtime, savedata => "")) do
                @test any(v -> occursin("expected 1 occurrence(s), observed 0", v), _rah_scan(root))
            end
        end

        @testset "Inert text, typed code and missing inputs" begin
            _rah_with_source(p.setup, """
                # Vector{Any}() and Dict{Symbol, Any}()
                #= nested #= Vector{Any}() =# Dict{Symbol,Any}() =#
                text = "Vector{Any}() # Dict{Symbol, Any}()"
                x = Vector{Float64}()
                y = Dict{Symbol, Float64}()
                z = Dict{String, Any}()
                """) do
                @test isempty(_rah_scan(root))
            end
            # A '#' inside a string must not hide subsequent active code.
            _rah_with_source(p.setup, "text = \"#\"; x = Vector{Any}()") do
                @test !isempty(_rah_scan(root))
            end
            _rah_with_source(p.setup, "x = Vector{") do
                @test any(v -> occursin("cannot parse source", v), _rah_scan(root))
            end
            setup = read(p.setup, String)
            rm(p.setup)
            try
                @test_throws r"Missing runtime Any hot-path target" _rah_scan(root)
            finally
                write(p.setup, setup)
            end
            moved = p.callbacks * "_removed"
            mv(p.callbacks, moved)
            try
                @test_throws r"Missing runtime callback scan root" _rah_scan(root)
            finally
                mv(moved, p.callbacks)
            end
            split_moved = p.split * "_removed"
            mv(p.split, split_moved)
            try
                @test_throws r"Missing runtime hot-path scan root" _rah_scan(root)
                @test _rah_actual_gate(root) isa LoadError
            finally
                mv(split_moved, p.split)
            end
        end
        @test isempty(_rah_scan(root))
    end
end
