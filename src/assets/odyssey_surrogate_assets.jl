module OdysseySurrogateAssets

import Artifacts
import Pkg
import SHA
import TOML

const _ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const _ARTIFACTS = joinpath(_ROOT, "Artifacts.toml")
const _MANIFEST = joinpath(_ROOT, "data", "odyssey_surrogate_assets.toml")
const _PACKAGE_UUID = Base.UUID("afbfb69f-5c0b-4832-b760-43725dff8540")

"""
    odyssey_surrogate_assets(; offline=false)

Retrieve and verify the public SPICE kernels and gravity coefficients for the
bounded Odyssey P20 exercise. Return explicit `kernels`, `gravity` paths and a
`provenance` record. Only this lazy artifact is fetched; native GRAM is unused.
With `offline=true`, an absent artifact fails before simulation starts. Every
payload is checked against its SHA256 even when Julia artifact overrides apply.
"""
function odyssey_surrogate_assets(; offline::Bool=false)
    manifest = TOML.parsefile(_MANIFEST)
    name = manifest["artifact_name"]
    meta = Artifacts.artifact_meta(name, _ARTIFACTS; pkg_uuid=_PACKAGE_UUID)
    meta === nothing && throw(ArgumentError("Odyssey scenario artifact binding is missing."))
    tree = meta["git-tree-sha1"]
    hash = Base.SHA1(tree)
    if !Artifacts.artifact_exists(hash)
        offline && throw(ArgumentError(
            "Odyssey scenario assets are not installed. Run odyssey_surrogate_assets() " *
            "once on a connected machine, then retry with offline=true."))
        try
            Pkg.Artifacts.ensure_artifact_installed(name, _ARTIFACTS; pkg_uuid=_PACKAGE_UUID)
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("Could not retrieve Odyssey scenario assets. Check connectivity and cache permissions, then retry. " * sprint(showerror, err)))
        end
    end
    root = Artifacts.artifact_path(hash)
    verified = Dict{String,String}()
    identities = Dict{String,Any}[]
    for entry in manifest["files"]
        relative = entry["path"]
        path = joinpath(root, relative)
        isfile(path) || throw(ArgumentError("Odyssey asset is missing: $path. Repair the managed artifact or its override."))
        filesize(path) == entry["bytes"] || throw(ArgumentError("Odyssey asset size mismatch: $path. Repair the managed artifact or its override."))
        digest = open(io -> bytes2hex(SHA.sha256(io)), path)
        digest == entry["sha256"] || throw(ArgumentError("Odyssey asset SHA256 mismatch: $path. Repair the managed artifact or its override; no fallback is attempted."))
        verified[entry["name"]] = path
        push!(identities, Dict("name" => entry["name"], "sha256" => digest, "source" => entry["source"]))
    end
    kernels = [verified[name] for name in manifest["kernel_order"]]
    provenance = Dict{String,Any}("id" => manifest["id"], "version" => manifest["version"],
        "git_tree_sha1" => tree, "files" => identities)
    return (; kernels, gravity=verified["Mars50c.csv"], provenance)
end

export odyssey_surrogate_assets
end
