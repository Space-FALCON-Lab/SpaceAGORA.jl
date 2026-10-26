using TOML
import ..AssetCatalog

const ASSET_STATES = AssetCatalog.ASSET_STATES

"""
    AssetCheckItem

Single asset-root status entry produced by `check_assets()`.

`state` is one of `:usable`, `:missing`, `:empty`, `:pointer_only` or
`:invalid`. `available` is `true` exactly when `state === :usable`; it is kept
so existing readers keep working. `reason` explains any other state in one
sentence, including what to do about it.

The original six-argument positional constructor and keyword construction
without `state` derive the state from `available`. An explicit state must be
one of the supported values and agree with `available`.
"""
Base.@kwdef struct AssetCheckItem
    name::String
    scope::String
    path::String
    available::Bool
    required::Bool
    detail::String = ""
    state::Symbol = convert(Bool, available) ? :usable : :missing
    reason::String = ""

    function AssetCheckItem(name, scope, path, available, required, detail, state, reason)
        ready = convert(Bool, available)
        status = convert(Symbol, state)
        status in ASSET_STATES || throw(ArgumentError("Unknown asset readiness state $(repr(status)) for $(name)."))
        ready == (status === :usable) || throw(ArgumentError("Asset availability must agree with readiness state $(repr(status)) for $(name)."))
        return new(name, scope, path, ready, required, detail, status, reason)
    end
end

AssetCheckItem(name, scope, path, available, required, detail) =
    AssetCheckItem(name, scope, path, available, required, detail,
                   convert(Bool, available) ? :usable : :missing, "")

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

@inline function _asset_item(name::String, scope::String, path::String; state::Symbol, required::Bool, detail::String="", reason::String="")
    return AssetCheckItem(
        name=name,
        scope=scope,
        path=path,
        available=state === :usable,
        required=required,
        detail=detail,
        state=state,
        reason=reason,
    )
end

function load_asset_manifest(; repo_root::String=REPO_ROOT, manifest_path::String=joinpath(repo_root, "data", "assets_manifest.toml"))
    return AssetManifestEntry[AssetManifestEntry(; entry...) for entry in
        AssetCatalog.load_asset_manifest(; repo_root, manifest_path)]
end

@inline _manifest_entry_path(repo_root::String, entry::AssetManifestEntry) =
    AssetCatalog.manifest_entry_path(repo_root, entry)

@inline _manifest_entry_readiness(repo_root::String, entry::AssetManifestEntry) =
    AssetCatalog.manifest_entry_readiness(repo_root, entry)

@inline _manifest_entry_available(repo_root::String, entry::AssetManifestEntry) =
    AssetCatalog.manifest_entry_available(repo_root, entry)

"""
    check_assets(; repo_root=REPO_ROOT, manifest_path=joinpath(repo_root, "data", "assets_manifest.toml")) -> AssetCheckReport

Inspect the standard SpaceAGORA asset roots and return a typed readiness
report for baseline, optional, and high-fidelity assets. Each item carries a
`state` (`:usable`, `:missing`, `:empty`, `:pointer_only` or `:invalid`) and a
`reason` for any state other than `:usable`. The check reads directory
listings and file headers only; it never downloads. By default, inspect the
package repository root, not the current working directory.
"""
function check_assets(; repo_root::String=REPO_ROOT, manifest_path::String=joinpath(repo_root, "data", "assets_manifest.toml"))::AssetCheckReport
    report = AssetCatalog.check_assets(; repo_root, manifest_path)
    return AssetCheckReport(
        repo_root=report.repo_root,
        items=AssetCheckItem[AssetCheckItem(; item...) for item in report.items],
    )
end

_asset_status_label(item::AssetCheckItem) =
    item.state === :usable ? "usable" :
    replace(String(item.state), '_' => '-') * (item.required ? "-required" : "-optional")

"""
    asset_check_exit_status(report) -> Int

`0` when every required entry is usable, `1` otherwise. Optional entries never
change the exit status; their state is only reported.
"""
asset_check_exit_status(report::AssetCheckReport)::Int =
    any(item -> item.required && item.state !== :usable, report.items) ? 1 : 0

"""
    render_asset_report(report; io=stdout)

Render a human-readable asset readiness report. Each entry shows its status
(`usable`, or its state suffixed with `-required` or `-optional`), its scope,
path and licensing detail, and a `reason` line for any entry that is not
usable. The report ends with a one-line summary.
"""
function render_asset_report(report::AssetCheckReport; io::IO=stdout)
    println(io, "SpaceAGORA asset check")
    println(io, "repo_root=$(report.repo_root)")
    for item in report.items
        println(io, "- $(item.name): $(_asset_status_label(item))")
        println(io, "  scope: $(item.scope)")
        println(io, "  path: $(item.path)")
        println(io, "  detail: $(item.detail)")
        isempty(item.reason) || println(io, "  reason: $(item.reason)")
    end
    usable = count(item -> item.state === :usable, report.items)
    blocked = count(item -> item.required && item.state !== :usable, report.items)
    println(io, "summary: $(usable) of $(length(report.items)) entries usable; $(blocked) required entries not usable")
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

function _run_preset_assets(command::String, args::Vector{String}; io::IO=stdout)::Int
    if command == "list"
        isempty(args) || throw(ArgumentError("assets list takes no options"))
        for entry in available_surrogate_presets()
            status = entry["release_enabled"] ? "available" : "unpublished"
            println(io, entry["id"], "@", entry["version"], " (", entry["planet"], ", ", status, ")")
            println(io, "  ", entry["description"])
        end
        return 0
    end
    options = Dict{String,String}()
    offline = command == "check"
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--offline"
            offline = true
        elseif startswith(arg, "--")
            key_value = split(arg[3:end], '='; limit=2)
            key = first(key_value)
            key in ("preset", "version", "file") || throw(ArgumentError("Unknown preset option $arg"))
            haskey(options, key) && throw(ArgumentError("Duplicate --$key option"))
            if length(key_value) == 2
                options[key] = key_value[2]
            else
                i += 1
                i <= length(args) && !startswith(args[i], "--") || throw(ArgumentError("--$key needs a value"))
                options[key] = args[i]
            end
            isempty(options[key]) && throw(ArgumentError("--$key needs a nonempty value"))
        else
            throw(ArgumentError("Unexpected preset argument $arg"))
        end
        i += 1
    end
    all(key -> haskey(options, key), ("preset", "version")) || throw(ArgumentError(
        "assets $command requires --preset=<id> and --version=<version>; use assets list"))
    preset = options["preset"]
    version = options["version"]
    file = get(options, "file", "")
    # Fetch verifies the downloaded bytes. Check additionally loads and validates
    # the schema, axes and metadata; neither command initializes native GRAM.
    # A composite preset has no payload of its own: fetch verifies its components.
    components = command == "fetch" ? _composite_preset_components(preset, version) : nothing
    if components !== nothing
        isempty(file) || _composite_file_error(preset, version)
        println(io, "Composite ", preset, "@", version, ": verifying its components")
        for c in components
            resolved = resolve_surrogate_preset(c.id; version=c.version, offline)
            println(io, "Verified ", c.role, " component ", c.id, "@", c.version, ": ", resolved.file)
            println(io, "sha256=", resolved.expected_sha256)
        end
    elseif command == "fetch"
        resolved = resolve_surrogate_preset(preset; version, file, offline)
        println(io, "Verified ", preset, "@", version, ": ", resolved.file)
        println(io, "sha256=", resolved.expected_sha256)
    else
        model = surrogate_preset_model(preset; version, file, offline)
        provenance = atmosphere_provenance(model)
        println(io, "Validated ", preset, "@", version)
        TOML.print(io, provenance)
    end
    return 0
end
