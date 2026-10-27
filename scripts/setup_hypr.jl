"""Compatibility entrypoint for the HYPR companion's shared standalone installer."""
module HYPRInstallation
using Pkg, TOML, UUIDs
include("setup_core.jl")
include("optional_packages.jl")
using .CoreEnvironment: retire_loader_manifest!
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

function setup(environment; root=ROOT, with_hypr=true,
        local_path=get(ENV, "SPACEAGORA_HYPR_PATH", ""), migrate_from="")
    env = abspath(environment)
    ispath(env) && samefile(env, root) && error("Use a separate HYPR project, not the SpaceAGORA root project.")
    if !with_hypr
        isempty(migrate_from) || throw(ArgumentError("Use the HYPR migration option only with HYPR enabled."))
        return CoreEnvironment.setup(env; root)
    end
    # Destination checks precede source acquisition, exactly as for other companions.
    OptionalPackages.destination(env; root=realpath(root))
    source = source_spec(root; local_path)
    provider = isempty(local_path) ? OptionalPackages.checkout("hypr"; root) : realpath(source.path)
    command = OptionalPackages.companion_installation_command(provider, env; root, migrate_from)
    run(command)
    return env
end
end
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) == 1
        HYPRInstallation.setup(only(ARGS))
    elseif length(ARGS) == 3 && ARGS[1] == "--migrate-from"
        HYPRInstallation.setup(ARGS[3]; migrate_from=ARGS[2])
    else
        error("Usage: julia scripts/setup_hypr.jl [--migrate-from <old-project>] <separate-environment>")
    end
end
