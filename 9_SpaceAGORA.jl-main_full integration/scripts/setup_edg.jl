"""Prepare an isolated core test/development environment with the EDG provider."""
module EDGTestInstallation
using Pkg, TOML
include(joinpath(@__DIR__, "optional_packages.jl"))
const ROOT = normpath(joinpath(@__DIR__, ".."))
const UUID = "a782e8ee-96d7-4dd6-8231-bb64649e6ad4"

function provider_path(root; local_path=get(ENV, "SPACEAGORA_EDG_PATH", ""))
    path = isempty(local_path) ? OptionalPackages.check("edg"; root) : realpath(local_path)
    project = TOML.parsefile(joinpath(path, "Project.toml"))
    get(project, "name", "") == "SpaceAGORAEDG" &&
        get(project, "uuid", "") == UUID && get(project, "version", "") == "0.1.0" ||
        throw(ArgumentError("Explicit EDG source must be SpaceAGORAEDG 0.1.0 with its declared UUID."))
    return path
end

function setup(environment; root=ROOT, with_edg=true,
               local_path=get(ENV, "SPACEAGORA_EDG_PATH", ""))
    root = realpath(root)
    destination = OptionalPackages.destination(environment; root)
    if !isempty(local_path)
        OptionalPackages.destination(destination; root=realpath(local_path))
    end
    marker = joinpath(destination, ".spaceagora-edg-test-environment.toml")
    if isdir(destination) && !isempty(readdir(destination))
        isfile(marker) || throw(ArgumentError("Refusing an existing unmanaged environment: $destination"))
        old = TOML.parsefile(marker)
        get(old, "root", "") == root && get(old, "with_edg", nothing) == with_edg &&
            get(old, "local_path", "") == local_path ||
            throw(ArgumentError("Existing environment belongs to a different EDG source pair."))
    end
    # Validate the destination before any network or source acquisition.
    with_edg && isempty(local_path) && OptionalPackages.checkout("edg"; root)
    provider = with_edg ? provider_path(root; local_path) : nothing
    project = Dict{String,Any}("deps" => copy(TOML.parsefile(joinpath(root, "Project.toml"))["deps"]))
    project["deps"]["SpaceAGORA"] = "afbfb69f-5c0b-4832-b760-43725dff8540"
    with_edg && (project["deps"]["SpaceAGORAEDG"] = UUID)
    mkpath(destination)
    open(marker, "w") do io
        TOML.print(io, Dict("root" => root, "with_edg" => with_edg, "local_path" => local_path))
    end
    open(joinpath(destination, "Project.toml"), "w") do io
        TOML.print(io, project; sorted=true)
    end
    manifest = joinpath(destination, "Manifest.toml")
    isfile(manifest) || cp(joinpath(root, "Manifest.toml"), manifest)
    previous = Base.active_project()
    try
        Pkg.activate(destination)
        specs = [PackageSpec(path=root)]
        with_edg && push!(specs, PackageSpec(path=provider))
        Pkg.develop(specs; preserve=Pkg.PRESERVE_ALL)
        Pkg.instantiate()
    finally
        previous === nothing ? Pkg.activate() : Pkg.activate(dirname(previous))
    end
    return destination
end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 1 || error("Usage: julia scripts/setup_edg.jl <separate-environment>")
    EDGTestInstallation.setup(only(ARGS))
end

