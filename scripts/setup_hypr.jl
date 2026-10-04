"""Install the checkout and its pinned optional HYPR dependency in an isolated project."""
module HYPRInstallation
using Pkg, TOML, UUIDs
const ROOT = normpath(joinpath(@__DIR__, ".."))
const HYPR_UUID = UUID("342c8629-b435-4bda-9d24-004c735fb297")

function source_spec(root=ROOT; local_path=get(ENV, "SPACEAGORA_HYPR_PATH", ""))
    if !isempty(local_path)
        path = abspath(local_path)
        project = TOML.parsefile(joinpath(path, "Project.toml"))
        project["uuid"] == string(HYPR_UUID) && project["version"] == "0.1.0" ||
            error("The explicit HYPR development override must be HYPR 0.1.0.")
        return PackageSpec(path=path)
    end
    pin = TOML.parsefile(joinpath(root, "packages", "SpaceAGORAHYPR", "HYPRSource.toml"))
    pin["uuid"] == string(HYPR_UUID) && pin["version"] == "0.1.0" || error("Unsupported HYPR source pin")
    occursin(r"^[0-9a-f]{40}$", pin["rev"]) || error("HYPR source must pin a full immutable commit")
    return PackageSpec(url=pin["url"], rev=pin["rev"])
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
    dependencies["SpaceAGORA"] = "afbfb69f-5c0b-4832-b760-43725dff8540"
    if with_hypr
        dependencies["HYPR"] = string(HYPR_UUID)
        dependencies["SpaceAGORAHYPR"] = "4096c904-79d0-4d9a-9f5d-f1d6c7a2f1b4"
    end
    sources = get!(project, "sources", Dict{String,Any}())
    sources["SpaceAGORA"] = Dict("path" => abspath(root))
    if with_hypr
        sources["SpaceAGORAHYPR"] = Dict("path" => joinpath(abspath(root), "packages", "SpaceAGORAHYPR"))
        sources["HYPR"] = isempty(local_path) ? Dict("url" => spec.url, "rev" => spec.rev) : Dict("path" => spec.path)
    end
    open(project_path, "w") do io
        TOML.print(io, project; sorted=true)
    end
    manifest = joinpath(env, "Manifest.toml")
    isfile(manifest) || cp(joinpath(root, "Manifest.toml"), manifest)
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
            push!(specs, PackageSpec(path=joinpath(root, "packages", "SpaceAGORAHYPR")))
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
