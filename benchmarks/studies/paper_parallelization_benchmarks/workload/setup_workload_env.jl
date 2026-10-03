# Create the paper harness's workload environment:
#
#   julia --startup-file=no setup_workload_env.jl <env_dir>
#
# The environment is the repository's exact Manifest.toml plus two path (develop)
# entries, SpaceAGORA (this checkout) and SpaceAGORAPaperWorkload (this
# directory), added with PRESERVE_ALL so no registered package moves off the
# repository's pinned version; the script fails if one did. The repository's own
# Project.toml and Manifest.toml are never written.
#
# The environment is used STACKED behind the repository project (JULIA_LOAD_PATH
# "@:<env_dir>:@v#.#:@stdlib" with --project=<repo>), never as the active
# project, so the only thing it contributes is where to find the workload
# package; everything else resolves from the repository exactly as it does
# without it. Instantiating does not precompile; build_workload.sh does that.
using Pkg

length(ARGS) == 1 || error("usage: setup_workload_env.jl <env_dir>")
env = abspath(ARGS[1])
workload_pkg = normpath(joinpath(@__DIR__, "SpaceAGORAPaperWorkload"))
repo = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
isfile(joinpath(repo, "Project.toml")) || error("repository root not found at $(repo)")

rm(env; recursive=true, force=true)
mkpath(env)
cp(joinpath(repo, "Manifest.toml"), joinpath(env, "Manifest.toml"); force=true)
write(joinpath(env, "Project.toml"), "")

ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"
# Offline: every package is already installed by the repository's own
# instantiate, and nothing here may resolve against a newer registry.
Pkg.offline(true)
Pkg.activate(env)
Pkg.develop([PackageSpec(path=repo), PackageSpec(path=workload_pkg)]; preserve=Pkg.PRESERVE_ALL)
Pkg.instantiate()

# Every registered package in the new manifest must carry the repository
# manifest's version. Standard libraries are left out: their versions are the
# ones the running Julia ships, so a host on another 1.12 patch release (TRX50
# runs 1.12.5, the manifest was written by 1.12.1) moves Pkg, Downloads and
# their jlls in any environment it resolves, the repository's own included.
is_stdlib_entry(entry) = haskey(entry, "uuid") &&
    Pkg.Types.is_stdlib(Base.UUID(entry["uuid"]))
parse_versions(path) = Dict(name => get(first(entries), "version", "stdlib/path")
                            for (name, entries) in Pkg.TOML.parsefile(path)["deps"]
                            if !is_stdlib_entry(first(entries)))
repo_v = parse_versions(joinpath(repo, "Manifest.toml"))
env_v = parse_versions(joinpath(env, "Manifest.toml"))
moved = [(n, repo_v[n], env_v[n]) for n in keys(repo_v) if haskey(env_v, n) && env_v[n] != repo_v[n]]
isempty(moved) || error("workload environment moved package versions off the repository manifest: $(moved)")
added = sort(collect(setdiff(keys(env_v), keys(repo_v))))
added == ["SpaceAGORAPaperWorkload"] ||
    error("workload environment added packages beyond the workload itself: $(added)")
@info "workload environment ready" env repo packages = length(env_v)
