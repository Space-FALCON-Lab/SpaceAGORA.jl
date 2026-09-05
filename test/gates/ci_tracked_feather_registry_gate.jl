# Every tracked feather file must be a registered truth input: either named in
# test/telemetry_benchmark_manifest.toml (graded by the telemetry regression)
# or under a relative_path declared in data/assets_manifest.toml. Large
# reference sets are fetched into gitignored paths, never committed; this gate
# keeps a force-added feather from slipping in silently (which is how ~480 MB
# of simulator output accumulated under data/telemetry before September 2026).
using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function _tracked_feathers()::Vector{String}
    output = read(`sh -lc "cd '$REPO_ROOT' && git ls-files '*.feather'"`, String)
    return [String(strip(l)) for l in split(chomp(output), '\n') if !isempty(strip(l))]
end

function _registered_paths()::Vector{String}
    paths = String[]
    manifest = TOML.parsefile(joinpath(REPO_ROOT, "test", "telemetry_benchmark_manifest.toml"))
    walk(x) = if x isa AbstractDict
        foreach(walk, values(x))
    elseif x isa AbstractVector
        foreach(walk, x)
    elseif x isa AbstractString && endswith(x, ".feather")
        push!(paths, normpath(x))
    end
    walk(manifest)
    assets = TOML.parsefile(joinpath(REPO_ROOT, "data", "assets_manifest.toml"))
    for entry in get(assets, "asset", Any[])
        rel = get(entry, "relative_path", "")
        isempty(rel) || rel == "." || push!(paths, normpath(rel))
    end
    return paths
end

registered = _registered_paths()
_covered(rel) = any(r -> r == rel || startswith(rel, r * "/"), registered)
unregistered = [rel for rel in _tracked_feathers() if !_covered(rel)]
isempty(unregistered) || error("Tracked feather files that no manifest registers (fetch them, do not commit them; see data/telemetry/PRIVATE_TELEMETRY.md):\n" * join(unregistered, '\n'))
println("tracked_feather_registry_gate_ok")
