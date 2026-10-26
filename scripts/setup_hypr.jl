"""Install the checkout and its pinned optional HYPR dependency in an isolated project."""
module HYPRInstallation
using Pkg, TOML, UUIDs
const SUPPORTED_HYPR_VERSIONS = ("0.5.0",)
const ROOT = normpath(joinpath(@__DIR__, ".."))
const HYPR_UUID = UUID("342c8629-b435-4bda-9d24-004c735fb297")

function source_spec(root=ROOT; local_path=get(ENV, "SPACEAGORA_HYPR_PATH", ""))
    if !isempty(local_path)
        path = abspath(local_path)
        project = TOML.parsefile(joinpath(path, "Project.toml"))
        project["uuid"] == string(HYPR_UUID) && project["version"] in SUPPORTED_HYPR_VERSIONS ||
            error("The explicit HYPR development override must be a supported HYPR version: $SUPPORTED_HYPR_VERSIONS.")
        return PackageSpec(path=path)
    end
    pin = TOML.parsefile(joinpath(root, "packages", "HYPRSource.toml"))
    pin["uuid"] == string(HYPR_UUID) && pin["version"] in SUPPORTED_HYPR_VERSIONS || error("Unsupported HYPR source pin")
    occursin(r"^[0-9a-f]{40}$", pin["rev"]) || error("HYPR source must pin a full immutable commit")
    return PackageSpec(url=pin["url"], rev=pin["rev"])
end

# Pkg preserves path-developed manifest entries during resolution, even after
# their direct dependency is removed. Retire this one entry before Pkg tries
# to read the deleted package directory; keep optional weak bindings by UUID.
function retire_loader_manifest!(path)
    data = TOML.parsefile(path)
    entries = get(data, "deps", data)
    haskey(entries, "SpaceAGORAHYPR") || return nothing
    for (name, records) in entries, record in records
        dependencies = get(record, "deps", String[])
        depends_on_loader = dependencies isa AbstractDict ?
            haskey(dependencies, "SpaceAGORAHYPR") : "SpaceAGORAHYPR" in dependencies
        depends_on_loader && error("$name still depends on the retired SpaceAGORAHYPR loader. Migrate that package to direct SpaceAGORA and HYPR imports before rerunning setup.")
        weak = get(record, "weakdeps", nothing)
        if weak isa AbstractVector && "SpaceAGORAHYPR" in weak
            record["weakdeps"] = Dict(dependency => only(entries[dependency])["uuid"] for dependency in weak)
        end
    end
    delete!(entries, "SpaceAGORAHYPR")
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return nothing
end

function setup(environment; root=ROOT, with_hypr=true, local_path=get(ENV, "SPACEAGORA_HYPR_PATH", ""))
    env = abspath(environment)
    ispath(env) && samefile(env, root) &&
        error("Use a separate HYPR project, not the SpaceAGORA root project.")
    spec = with_hypr ? source_spec(root; local_path) : nothing
    mkpath(env)
    project_path = joinpath(env, "Project.toml")
    project = isfile(project_path) ? TOML.parsefile(project_path) : Dict{String,Any}()
    dependencies = get!(project, "deps", Dict{String,Any}())
    isfile(project_path) || merge!(dependencies, TOML.parsefile(joinpath(root, "Project.toml"))["deps"])
    # Reusing an environment must not retain the removed loader or its source path.
    pop!(dependencies, "SpaceAGORAHYPR", nothing)
    haskey(project, "compat") && pop!(project["compat"], "SpaceAGORAHYPR", nothing)
    dependencies["SpaceAGORA"] = "afbfb69f-5c0b-4832-b760-43725dff8540"
    if with_hypr
        dependencies["HYPR"] = string(HYPR_UUID)
    end
    sources = get!(project, "sources", Dict{String,Any}())
    pop!(sources, "SpaceAGORAHYPR", nothing)
    sources["SpaceAGORA"] = Dict("path" => abspath(root))
    if with_hypr
        sources["HYPR"] = isempty(local_path) ? Dict("url" => spec.url, "rev" => spec.rev) : Dict("path" => spec.path)
    end
    manifest = joinpath(env, "Manifest.toml")
    isfile(manifest) || cp(joinpath(root, "Manifest.toml"), manifest)
    retire_loader_manifest!(manifest)
    open(project_path, "w") do io
        TOML.print(io, project; sorted=true)
    end
    previous = Base.active_project()
    try
        Pkg.activate(env)
        specs = [PackageSpec(path=root)]
        if with_hypr
            if isempty(local_path)
                Pkg.add(spec; preserve=Pkg.PRESERVE_ALL)
            else
                push!(specs, spec)
            end
        end
        Pkg.develop(specs; preserve=Pkg.PRESERVE_ALL)
        Pkg.instantiate()
    finally
        previous === nothing ? Pkg.activate() : Pkg.activate(dirname(previous))
    end
    return env
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: julia scripts/setup_hypr.jl <project-directory>")
    HYPRInstallation.setup(only(ARGS))
end
