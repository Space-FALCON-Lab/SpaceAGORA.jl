"""Prepare a core-only development or example project without acquiring companions."""
module CoreEnvironment
using Pkg, TOML
include("companion_support/legacy_environment.jl")
using .LegacyEnvironmentMigration: retire_loader_manifest!, migrated_records
const ROOT = normpath(joinpath(@__DIR__, ".."))

function setup(environment; root=ROOT)
    root = realpath(root)
    env = abspath(environment)
    ispath(env) && samefile(env, root) && error("Use a separate core project, not the SpaceAGORA root project.")
    project_path = joinpath(env, "Project.toml")
    manifest_path = joinpath(env, "Manifest.toml")
    source_manifest = isfile(manifest_path) ? manifest_path : joinpath(root, "Manifest.toml")
    records = migrated_records(project_path, source_manifest)
    project = records.project
    dependencies = get!(project, "deps", Dict{String,Any}())
    isfile(project_path) || merge!(dependencies, TOML.parsefile(joinpath(root, "Project.toml"))["deps"])
    dependencies["SpaceAGORA"] = "afbfb69f-5c0b-4832-b760-43725dff8540"
    get!(project, "sources", Dict{String,Any}())["SpaceAGORA"] = Dict("path" => root)
    mkpath(env)
    records.manifest === nothing || open(io -> TOML.print(io, records.manifest; sorted=true), manifest_path, "w")
    open(io -> TOML.print(io, project; sorted=true), project_path, "w")
    previous = Base.active_project()
    try
        Pkg.activate(env)
        Pkg.develop(PackageSpec(path=root); preserve=Pkg.PRESERVE_ALL)
        Pkg.instantiate()
    finally
        previous === nothing ? Pkg.activate() : Pkg.activate(dirname(previous))
    end
    return env
end
end
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: julia scripts/setup_core.jl <project-directory>")
    CoreEnvironment.setup(only(ARGS))
end
