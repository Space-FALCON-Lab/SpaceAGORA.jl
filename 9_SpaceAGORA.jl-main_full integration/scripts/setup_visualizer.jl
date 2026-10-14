"""Install this core checkout and the pinned optional renderer in a separate project."""
module VisualizerInstallation
using Pkg, TOML, UUIDs
const ROOT = normpath(joinpath(@__DIR__, ".."))
const VISUALIZER_UUID = UUID("a4677e92-0b7f-4d10-9ebd-8cde6a49f03a")
const SUPPORTED_VERSION = "0.1.0"

function source_spec(root=ROOT; local_path=get(ENV, "SPACEAGORA_VISUALIZER_PATH", ""))
    if !isempty(local_path)
        path = abspath(local_path)
        project = TOML.parsefile(joinpath(path, "Project.toml"))
        get(project, "name", nothing) == "SpaceAGORAVisualizer" &&
            get(project, "uuid", nothing) == string(VISUALIZER_UUID) &&
            get(project, "version", nothing) == SUPPORTED_VERSION ||
            error("The explicit visualizer override must be SpaceAGORAVisualizer $SUPPORTED_VERSION with its declared UUID.")
        return PackageSpec(path=path)
    end
    pin = TOML.parsefile(joinpath(root, "packages", "VisualizerSource.toml"))
    get(pin, "uuid", nothing) == string(VISUALIZER_UUID) &&
        get(pin, "version", nothing) == SUPPORTED_VERSION || error("Unsupported visualizer source pin")
    occursin(r"^[0-9a-f]{40}$", get(pin, "rev", "")) ||
        error("Visualizer source must pin a full immutable commit")
    get(pin, "url", "") == "https://github.com/Space-FALCON-Lab/SpaceAGORAVisualizer.jl.git" ||
        error("Unexpected visualizer source repository")
    return PackageSpec(url=pin["url"], rev=pin["rev"])
end

# Pure preparation allows source/protection tests without installing packages.
function installation_plan(environment; root=ROOT, with_visualizer=true,
                           local_path=get(ENV, "SPACEAGORA_VISUALIZER_PATH", ""))
    root = abspath(root)
    env = abspath(environment)
    ispath(env) && samefile(env, root) &&
        error("Use a separate visualizer project, not the SpaceAGORA root project.")
    spec = with_visualizer ? source_spec(root; local_path) : nothing
    if with_visualizer && !isempty(local_path) && ispath(env) && samefile(env, spec.path)
        error("Use a separate visualizer project, not the companion source project.")
    end
    project_path = joinpath(env, "Project.toml")
    project = isfile(project_path) ? TOML.parsefile(project_path) : Dict{String,Any}()
    dependencies = get!(project, "deps", Dict{String,Any}())
    isfile(project_path) || merge!(dependencies, TOML.parsefile(joinpath(root, "Project.toml"))["deps"])
    !with_visualizer && haskey(dependencies, "SpaceAGORAVisualizer") &&
        error("The core-only project must not already contain SpaceAGORAVisualizer.")
    dependencies["SpaceAGORA"] = "afbfb69f-5c0b-4832-b760-43725dff8540"
    sources = get!(project, "sources", Dict{String,Any}())
    sources["SpaceAGORA"] = Dict("path" => root)
    if with_visualizer
        dependencies["SpaceAGORAVisualizer"] = string(VISUALIZER_UUID)
        sources["SpaceAGORAVisualizer"] = isempty(local_path) ?
            Dict("url" => spec.url, "rev" => spec.rev) : Dict("path" => spec.path)
    end
    return (; environment=env, root, project, spec, with_visualizer)
end

function setup(environment; root=ROOT, with_visualizer=true,
               local_path=get(ENV, "SPACEAGORA_VISUALIZER_PATH", ""))
    plan = installation_plan(environment; root, with_visualizer, local_path)
    mkpath(plan.environment)
    open(joinpath(plan.environment, "Project.toml"), "w") do io
        TOML.print(io, plan.project; sorted=true)
    end
    manifest = joinpath(plan.environment, "Manifest.toml")
    isfile(manifest) || cp(joinpath(plan.root, "Manifest.toml"), manifest)
    previous = Base.active_project()
    try
        Pkg.activate(plan.environment)
        specs = [PackageSpec(path=plan.root)]
        if with_visualizer
            if isempty(local_path)
                Pkg.add(plan.spec; preserve=Pkg.PRESERVE_ALL)
            else
                push!(specs, plan.spec)
            end
        end
        Pkg.develop(specs; preserve=Pkg.PRESERVE_ALL)
        Pkg.instantiate()
    finally
        previous === nothing ? Pkg.activate() : Pkg.activate(dirname(previous))
    end
    return plan.environment
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: julia scripts/setup_visualizer.jl <separate-environment>")
    VisualizerInstallation.setup(only(ARGS))
end
