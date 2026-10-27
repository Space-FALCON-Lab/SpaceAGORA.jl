module RuntimeAnyGateFixtures
export _RAH_REPO, _rah_relative_path, _rah_paths, _rah_scan, _rah_fixture,
    _rah_with_source, _rah_actual_gate, _rah_gate_result
const _RAH_REPO = normpath(joinpath(@__DIR__, "..", ".."))
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

function _rah_gate_result(root)
    gate = joinpath(_RAH_REPO, "test", "gates", "ci_runtime_any_hotpath_gate.jl")
    needle = "const REPO_ROOT = normpath(joinpath(@__DIR__, \"..\", \"..\"))"
    source = read(gate, String)
    @assert occursin(needle, source)
    source = replace(source, needle => "const REPO_ROOT = $(repr(root))"; count=1)
    mod = Module(gensym(:RuntimeAnyGate))
    Core.eval(mod, :(include(path) = Base.include(@__MODULE__, path)))
    failure = try
        redirect_stdout(devnull) do
            Base.include_string(mod, source, gate)
        end
        nothing
    catch err
        err
    end
    # Read dynamically included globals in the module's current world.
    violations = Core.eval(mod,
        :(isdefined(@__MODULE__, :violations) ? copy(violations) : String[]))
    return (; failure, violations)
end
_rah_actual_gate(root) = _rah_gate_result(root).failure

end
