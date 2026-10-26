using Test
using TOML
using Pkg

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const STDLIBS = Set([
    "Artifacts",
    "Base64",
    "Dates",
    "DelimitedFiles",
    "Distributed",
    "FileWatching",
    "InteractiveUtils",
    "LibCURL",
    "LibGit2",
    "Libdl",
    "LinearAlgebra",
    "Logging",
    "Markdown",
    "Mmap",
    "Pkg",
    "Printf",
    "Profile",
    "Random",
    "REPL",
    "SHA",
    "Serialization",
    "SharedArrays",
    "Sockets",
    "SparseArrays",
    "Statistics",
    "SuiteSparse",
    "TOML",
    "Tar",
    "Test",
    "UUIDs",
    "Unicode"
])

function missing_compat_entries(project_path::String)::Vector{String}
    project = TOML.parsefile(project_path)
    deps = Set(String.(keys(get(project, "deps", Dict()))))
    compat = Set(String.(keys(get(project, "compat", Dict()))))
    delete!(compat, "julia")
    return sort(collect(setdiff(setdiff(deps, STDLIBS), compat)))
end

for project_path in (
    joinpath(REPO_ROOT, "Project.toml"),
    joinpath(REPO_ROOT, "docs", "Project.toml"),
    joinpath(REPO_ROOT, "test", "Project.toml"),
    joinpath(REPO_ROOT, "examples", "Project.toml"),
)
    missing = missing_compat_entries(project_path)
    isempty(missing) || error("Missing compat entries in $(relpath(project_path, REPO_ROOT)): $(join(missing, ", "))")
end

# These environments run the checked-out core, so their bounds must admit its
# actual version. Starter bounds are placeholders rewritten by the generator;
# independently pinned companion and study environments have their own owners.
core_version = VersionNumber(TOML.parsefile(joinpath(REPO_ROOT, "Project.toml"))["version"])
for relative in ("docs", "test", "examples", "examples/odyssey_surrogate_env", "test/grid_atmosphere")
    project_path = joinpath(REPO_ROOT, relative, "Project.toml")
    project = TOML.parsefile(project_path)
    bound = get(get(project, "compat", Dict()), "SpaceAGORA", nothing)
    bound isa AbstractString || error("Missing SpaceAGORA compatibility bound in $relative/Project.toml")
    core_version in Pkg.Types.semver_spec(bound) ||
        error("$relative/Project.toml excludes checked-out SpaceAGORA $core_version (compat: $bound)")
end

println("project_compat_gate_ok")
