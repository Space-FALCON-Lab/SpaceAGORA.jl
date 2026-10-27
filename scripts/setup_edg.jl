"""Compatibility entrypoint for the EDG companion's shared standalone installer."""
module EDGTestInstallation
using TOML
include(joinpath(@__DIR__, "setup_core.jl"))
include(joinpath(@__DIR__, "optional_packages.jl"))
const ROOT = normpath(joinpath(@__DIR__, ".."))
const UUID = "a782e8ee-96d7-4dd6-8231-bb64649e6ad4"

function provider_path(root; local_path=get(ENV, "SPACEAGORA_EDG_PATH", ""))
    path = isempty(local_path) ? OptionalPackages.check("edg"; root) : realpath(local_path)
    project = TOML.parsefile(joinpath(path, "Project.toml"))
    get(project, "name", "") == "SpaceAGORAEDG" &&
        get(project, "uuid", "") == UUID && get(project, "version", "") == "0.5.1" ||
        throw(ArgumentError("Explicit EDG source must be SpaceAGORAEDG 0.5.1 with its declared UUID."))
    return path
end

function setup(environment; root=ROOT, with_edg=true,
               local_path=get(ENV, "SPACEAGORA_EDG_PATH", ""))
    root = realpath(root)
    destination = OptionalPackages.destination(environment; root)
    if !isempty(local_path)
        OptionalPackages.destination(destination; root=realpath(local_path))
    end
    with_edg || return CoreEnvironment.setup(destination; root)
    isempty(local_path) && OptionalPackages.checkout("edg"; root)
    provider = provider_path(root; local_path)
    run(OptionalPackages.companion_installation_command(provider, destination; root))
    return destination
end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 1 || error("Usage: julia scripts/setup_edg.jl <separate-environment>")
    EDGTestInstallation.setup(only(ARGS))
end

