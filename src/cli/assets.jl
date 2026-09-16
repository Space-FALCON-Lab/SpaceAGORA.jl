using TOML

"""
    AssetCheckItem

Single asset-root status entry produced by `check_assets()`.
"""
Base.@kwdef struct AssetCheckItem
    name::String
    scope::String
    path::String
    available::Bool
    required::Bool
    detail::String = ""
end

"""
    AssetCheckReport

Repository-local asset availability report returned by `check_assets()`.
"""
Base.@kwdef struct AssetCheckReport
    repo_root::String
    items::Vector{AssetCheckItem}
end

Base.@kwdef struct AssetManifestEntry
    name::String
    scope::String
    relative_path::String
    kind::String
    required::Bool
    licensing::String
    detail::String = ""
end

@inline function _asset_item(name::String, scope::String, path::String; available::Bool, required::Bool, detail::String="")
    return AssetCheckItem(
        name=name,
        scope=scope,
        path=path,
        available=available,
        required=required,
        detail=detail
    )
end

function load_asset_manifest(; repo_root::String=REPO_ROOT, manifest_path::String=joinpath(repo_root, "data", "assets_manifest.toml"))
    raw = TOML.parsefile(manifest_path)
    entries = AssetManifestEntry[]
    for entry in get(raw, "asset", Any[])
        push!(entries, AssetManifestEntry(
            name=String(entry["name"]),
            scope=String(entry["scope"]),
            relative_path=String(get(entry, "relative_path", ".")),
            kind=String(get(entry, "kind", "directory")),
            required=Bool(get(entry, "required", false)),
            licensing=String(get(entry, "licensing", "unspecified")),
            detail=String(get(entry, "detail", "")),
        ))
    end
    return entries
end

@inline function _manifest_entry_path(repo_root::String, entry::AssetManifestEntry)
    entry.relative_path in ("", ".") && return repo_root
    return joinpath(repo_root, splitpath(entry.relative_path)...)
end

@inline function _manifest_entry_available(repo_root::String, entry::AssetManifestEntry)
    path = _manifest_entry_path(repo_root, entry)
    entry.kind == "builtin" && return true
    entry.kind == "directory" && return isdir(path)
    entry.kind == "file" && return isfile(path)
    throw(ArgumentError("Unsupported asset manifest kind '$(entry.kind)' for $(entry.name)."))
end

"""
    check_assets(; repo_root=REPO_ROOT) -> AssetCheckReport

Inspect the standard SpaceAGORA asset roots and return a typed availability report.
"""
function check_assets(; repo_root::String=REPO_ROOT, manifest_path::String=joinpath(repo_root, "data", "assets_manifest.toml"))::AssetCheckReport
    entries = load_asset_manifest(; repo_root=repo_root, manifest_path=manifest_path)
    items = AssetCheckItem[]
    for entry in entries
        detail = isempty(entry.detail) ? "licensing=$(entry.licensing)" : "$(entry.detail) [licensing=$(entry.licensing)]"
        push!(items, _asset_item(
            entry.name,
            entry.scope,
            _manifest_entry_path(repo_root, entry);
            available=_manifest_entry_available(repo_root, entry),
            required=entry.required,
            detail=detail,
        ))
    end

    return AssetCheckReport(repo_root=repo_root, items=items)
end

"""
    render_asset_report(report; io=stdout)

Render a human-readable summary of an `AssetCheckReport`.
"""
function render_asset_report(report::AssetCheckReport; io::IO=stdout)
    println(io, "SpaceAGORA asset check")
    println(io, "repo_root=$(report.repo_root)")
    for item in report.items
        status = item.available ? "available" : (item.required ? "missing-required" : "missing-optional")
        println(io, "- $(item.name): $(status)")
        println(io, "  scope: $(item.scope)")
        println(io, "  path: $(item.path)")
        println(io, "  detail: $(item.detail)")
    end
    return nothing
end

function render_asset_manifest(entries::Vector{AssetManifestEntry}; io::IO=stdout)
    println(io, "SpaceAGORA asset manifest")
    println(io, "entries=$(length(entries))")
    for entry in entries
        println(io, "- $(entry.name)")
        println(io, "  scope: $(entry.scope)")
        println(io, "  kind: $(entry.kind)")
        println(io, "  relative_path: $(entry.relative_path)")
        println(io, "  required: $(entry.required)")
        println(io, "  licensing: $(entry.licensing)")
        !isempty(entry.detail) && println(io, "  detail: $(entry.detail)")
    end
    return nothing
end

function setup_open_assets(; repo_root::String=REPO_ROOT, io::IO=stdout)
    entries = load_asset_manifest(; repo_root=repo_root)
    report = check_assets(; repo_root=repo_root)
    manifest_path = joinpath(repo_root, "data", "assets_manifest.toml")
    println(io, "SpaceAGORA open-asset setup")
    println(io, "No downloads are required for baseline no-GRAM mode.")
    println(io, "Manifest path: $manifest_path")
    println(io)
    println(io, "Baseline/open entries:")
    for entry in entries
        entry.licensing == "licensed-external" && continue
        println(io, "- $(entry.name): $(entry.relative_path)")
    end
    println(io)
    println(io, "Licensed external entries remain user-provided:")
    for entry in entries
        entry.licensing == "licensed-external" || continue
        println(io, "- $(entry.name): $(entry.relative_path)")
    end
    println(io)
    render_asset_report(report; io=io)
    return report
end
