using Test
using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

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
    joinpath(REPO_ROOT, ".AGORA", "Project.toml")
)
    missing = missing_compat_entries(project_path)
    isempty(missing) || error("Missing compat entries in $(relpath(project_path, REPO_ROOT)): $(join(missing, ", "))")
end

println("project_compat_gate_ok")
