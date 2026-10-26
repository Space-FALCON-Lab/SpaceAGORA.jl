# Named data selection only. GRAMSuite and the existing adapter own grid mathematics.
import Artifacts
import Pkg
import SHA
import TOML

const _SURROGATE_CATALOG = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "data", "surrogate_presets.toml"))
const _SURROGATE_ARTIFACTS = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "Artifacts.toml"))
const _SURROGATE_PACKAGE_UUID = Base.UUID("afbfb69f-5c0b-4832-b760-43725dff8540")
const _PRESET_METADATA_KEY = "spaceagora_named_surrogate"
_preset_error(message) = throw(ArgumentError(message))
# A preset's kind selects its model: a fixed grid (the default when absent), a near-surface payload, or a composite of
# two catalog entries by height, which has no payload of its own.
const _COMPOSITE_KIND = "composite_by_height"
const _PRESET_KINDS = ("gram_grid", "gram_near_surface_scalars", _COMPOSITE_KIND)
# A composite's only rule: the lower component below the handover height, the upper at or above it, nothing blended.
const _COMPOSITE_RULE = "lower_below_upper_at_or_above"
# The keys that identify a payload. A composite has none: its components' entries hold the payload identities.
const _PAYLOAD_KEYS = ("payload", "artifact_name", "distribution", "required_metadata", "axes", "generation_provenance")
# A near-surface payload's format says whether it stores winds: format 1 (published versions 1.0.0 and 1.1.0) stores
# none and returns a zero wind vector; format 2 adds a wind layer and returns its stored winds.
const _NEAR_SURFACE_FORMAT = "spaceagora_mars_near_surface_scalars_v1"
const _NEAR_SURFACE_WIND_FORMAT = "spaceagora_mars_near_surface_v2"
const _WIND_FIELDS = ("wind_ew_ms", "wind_ns_ms", "wind_up_ms")
_preset_kind(entry) = get(entry, "kind", "gram_grid")
_preset_digest(file) = open(io -> bytes2hex(SHA.sha256(io)), file)

function _preset_planet(planet::AbstractString)
    key = lowercase(strip(planet))
    key in ("earth", "mars", "venus", "jupiter", "uranus", "neptune", "titan") ||
        _preset_error("Unknown surrogate planet '$planet'. Choose an explicit catalog preset and supported planet.")
    return key
end
# A catalog axis is uniform (start, step, count) or an explicit, strictly increasing
# node list (nodes), in km or degrees. Longitude stays uniform so its periodic seam
# check keeps its meaning.
function _preset_axis_nodes(spec)
    if haskey(spec, "nodes")
        !any(key -> haskey(spec, key), ("start", "step", "count")) ||
            _preset_error("A preset axis gives either explicit nodes or start/step/count, not both.")
        nodes = spec["nodes"]
        nodes isa Vector && length(nodes) >= 2 && all(x -> x isa Real && !(x isa Bool) && isfinite(x), nodes) &&
            all(>(0), diff(Float64.(nodes))) || _preset_error("Invalid preset axis nodes: they must be at least two finite, strictly increasing values.")
        return Float64.(nodes)
    end
    spec["count"] isa Integer && !(spec["count"] isa Bool) && spec["count"] >= 2 || _preset_error("Invalid preset axis count.")
    all(x -> x isa Real && !(x isa Bool) && isfinite(x), (spec["start"], spec["step"])) && spec["step"] > 0 || _preset_error("Invalid preset axis spacing.")
    return collect(range(Float64(spec["start"]); step=Float64(spec["step"]), length=spec["count"]))
end
_preset_axis_bounds(spec) = haskey(spec, "nodes") ? [first(spec["nodes"]), last(spec["nodes"])] :
    [spec["start"], spec["start"] + (spec["count"]-1)*spec["step"]]
function _preset_catalog(path)
    isfile(path) || _preset_error("Surrogate preset catalog is missing: $path")
    bytes = read(path)
    catalog = TOML.parse(String(copy(bytes)))
    get(catalog, "schema_version", nothing) == 1 || _preset_error("Unsupported surrogate catalog schema.")
    revision = get(catalog, "catalog_revision", nothing)
    revision isa AbstractString && !isempty(revision) || _preset_error("The surrogate catalog needs a catalog_revision string.")
    entries = get(catalog, "presets", nothing)
    entries isa Vector || _preset_error("The surrogate catalog must contain preset entries.")
    all(entry -> entry isa AbstractDict, entries) || _preset_error("Every surrogate catalog preset entry must be a table.")
    seen = Set{Tuple{String,String}}()
    for (n, entry) in enumerate(entries)
        _preset_entry_checked(n, entry) do
            id = get(entry, "id", ""); version = get(entry, "version", "")
            occursin(r"^[a-z][a-z0-9_]*$", id) || _preset_error("Invalid surrogate preset ID '$id'.")
            occursin(r"^\d+\.\d+\.\d+$", version) || _preset_error("Preset '$id' needs an exact release version; moving versions are unsupported.")
            key = (id, version)
            key ∉ seen || _preset_error("Duplicate preset/version '$id@$version' in catalog.")
            push!(seen, key)
            kind = _preset_kind(entry)
            kind in _PRESET_KINDS || _preset_error("Preset '$id' has an unsupported kind '$kind'.")
            _preset_planet(entry["planet"])
            # A composite refers to other entries, which may come later in the file; it is checked once all are read.
            kind == _COMPOSITE_KIND || _preset_payload_entry(entry, kind)
        end
    end
    for (n, entry) in enumerate(entries)
        _preset_kind(entry) == _COMPOSITE_KIND && _preset_entry_checked(() -> _preset_composite_entry(entry, entries), n, entry)
    end
    return catalog, bytes2hex(SHA.sha256(bytes))
end

# Runs an entry's checks. A malformed entry, with a table or field missing or of the wrong type, fails with an
# ArgumentError naming the entry instead of the lower-level error the checks would raise.
function _preset_entry_checked(checks, n, entry)
    try
        checks()
    catch err
        err isa Union{ArgumentError, InterruptException} && rethrow()
        detail = first(split(sprint(showerror, err), '\n'))
        _preset_error("Surrogate catalog entry $n ('$(get(entry, "id", "?"))@$(get(entry, "version", "?"))') is malformed: $detail. " *
            "Each table and field the catalog requires must be present, with its expected type.")
    end
    return nothing
end

# A preset with a payload: a fixed grid or a near-surface payload.
function _preset_payload_entry(entry, kind)
    label = _preset_label(entry)
    payload = entry["payload"]
    occursin(r"^[0-9a-f]{64}$", payload["sha256"]) || _preset_error("Invalid preset payload SHA256.")
    payload["bytes"] isa Integer && payload["bytes"] > 0 || _preset_error("Invalid preset payload size.")
    file = payload["file"]
    file isa String && !isempty(file) && !isabspath(file) && !occursin('\\', file) &&
        all(x -> x ∉ ("", ".", ".."), split(file, '/')) || _preset_error("Unsafe preset payload filename.")
    get(entry, "release_enabled", nothing) isa Bool || _preset_error("Preset release_enabled must be explicit.")
    get(entry["atmosphere"], "query_elapsed_time_applied", nothing) === false ||
        _preset_error("Only explicitly frozen atmospheric presets are supported.")
    atmosphere = entry["atmosphere"]; domain = entry["domain"]
    required = entry["required_metadata"]; generation = required["generation_config"]
    _preset_planet(required["planet"]) == _preset_planet(entry["planet"]) || _preset_error("Preset metadata planet is inconsistent.")
    payload["format"] == required["format"] || _preset_error("Preset payload format is inconsistent.")
    kind == "gram_grid" ? _preset_grid_domain(entry) : _preset_near_surface_domain(entry)
    atmosphere["latitude"] == "geodetic" && generation["is_planetocentric"] === false && atmosphere["height"] == "ellipsoidal" || _preset_error("Preset coordinate datum is inconsistent or unsupported.")
    for (public_name, native_name) in (("equatorial_radius_m","equatorial_radius_km"),("polar_radius_m","polar_radius_km"))
        atmosphere[public_name] == 1000*generation[native_name] || _preset_error("Preset radii are inconsistent.")
    end
    epoch = match(r"^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)Z$", atmosphere["epoch_utc"])
    if kind == "gram_grid"
        initial = required["initial_time"]
        initial == generation["initial_time"] || _preset_error("Preset epoch metadata is inconsistent.")
        epoch !== nothing && all(parse(Float64,epoch.captures[i]) == initial[k] for (i,k) in enumerate(("year","month","day","hour","minute","second"))) || _preset_error("Preset frozen epoch differs from its generation metadata.")
    else
        epoch !== nothing && atmosphere["epoch_utc"] == required["epoch_utc"] || _preset_error("Preset frozen epoch differs from its payload metadata.")
    end
    # Resolution records the support statement in the provenance.
    get(entry, "support", nothing) isa AbstractDict || _preset_error("Preset $label needs its support table.")
    if entry["release_enabled"]
        name = get(entry, "artifact_name", nothing)
        name isa AbstractString && !isempty(name) || _preset_error("Released preset $label needs its artifact_name.")
        distribution = entry["distribution"]; urls = get(distribution, "urls", nothing)
        occursin(r"^[0-9a-f]{40}$", distribution["git_tree_sha1"]) && occursin(r"^[0-9a-f]{64}$", distribution["archive_sha256"]) &&
            urls isa AbstractVector && !isempty(urls) && all(url -> url isa AbstractString && !isempty(url), urls) ||
            _preset_error("Released presets require complete artifact publication identities.")
    else
        # An unpublished preset may have no distribution record yet.
        get(entry, "distribution", Dict{String,Any}()) isa AbstractDict || _preset_error("The distribution record of preset $label must be a table.")
    end
    return nothing
end

_preset_label(entry) = "'$(entry["id"])@$(entry["version"])'"
# An entry's declared output fields: a list of names.
function _preset_fields(entry)
    fields = get(entry["atmosphere"], "fields", Any[])
    fields isa AbstractVector && all(field -> field isa AbstractString, fields) ||
        _preset_error("The atmosphere fields of preset $(_preset_label(entry)) must be a list of field names.")
    return fields
end
# Whether an entry returns stored winds: a grid with stored local ENU winds, or a near-surface payload of format 2.
function _preset_stored_winds(entry)
    atmosphere = entry["atmosphere"]
    _preset_kind(entry) == "gram_near_surface_scalars" && return get(atmosphere, "winds_available", nothing) === true
    return get(atmosphere, "wind_convention", nothing) == "stored_nominal_local_ENU" &&
        get(atmosphere, "wind_boolean_suppresses_stored_components", nothing) === false && all(in(_preset_fields(entry)), _WIND_FIELDS)
end
# A composite preset names two catalog entries and the height that separates them; it has no payload of its own.
function _preset_composite_entry(entry, entries)
    name = "Composite preset $(_preset_label(entry))"; subject = "composite preset $(_preset_label(entry))"
    present = filter(key -> haskey(entry, key), _PAYLOAD_KEYS)
    isempty(present) || _preset_error("$name has no payload of its own, so it cannot declare $(join(present, ", ")); its components' entries hold the payload identities.")
    get(entry, "release_enabled", nothing) isa Bool || _preset_error("Preset release_enabled must be explicit.")
    get(entry, "components", nothing) isa AbstractVector || _preset_error("$name needs its components, a list of tables.")
    for key in ("handover", "atmosphere", "domain", "support")
        get(entry, key, nothing) isa AbstractDict || _preset_error("$name needs its $key table.")
    end
    components = entry["components"]
    length(components) == 2 && all(c -> c isa AbstractDict && get(c, "role", nothing) isa AbstractString, components) &&
        sort([c["role"] for c in components]) == ["lower", "upper"] ||
        _preset_error("$name needs exactly two components, one with role \"lower\" and one with role \"upper\".")
    side = Dict{String, Any}()
    for c in components
        id, version = get(c, "id", ""), get(c, "version", "")
        matches = filter(p -> p["id"] == id && p["version"] == version, entries)
        length(matches) == 1 || _preset_error("$name names the component '$id@$version', which is not in the catalog.")
        component = only(matches)
        _preset_kind(component) == _COMPOSITE_KIND && _preset_error("$name names the component '$id@$version', which is itself a composite; composites do not nest.")
        entry["release_enabled"] && !component["release_enabled"] && _preset_error("$name is released, but its component '$id@$version' is not.")
        _preset_planet(component["planet"]) == _preset_planet(entry["planet"]) ||
            _preset_error("$name is for $(entry["planet"]), but its component '$id@$version' is for $(component["planet"]).")
        side[c["role"]] = component
    end
    lower, upper = side["lower"], side["upper"]
    # The validated shape: terrain-following data below the handover and a fixed grid above it.
    _preset_kind(lower) == "gram_near_surface_scalars" && _preset_kind(upper) == "gram_grid" ||
        _preset_error("$name needs a near-surface preset as its lower component and a grid preset as its upper one, " *
            "not $(_preset_kind(lower)) $(_preset_label(lower)) and $(_preset_kind(upper)) $(_preset_label(upper)).")
    atmosphere = entry["atmosphere"]
    get(atmosphere, "query_elapsed_time_applied", nothing) === false || _preset_error("Only explicitly frozen atmospheric presets are supported.")
    get(atmosphere, "latitude", nothing) == "geodetic" && get(atmosphere, "height", nothing) == "ellipsoidal" ||
        _preset_error("Preset coordinate datum is inconsistent or unsupported.")
    for component in (lower, upper)
        theirs = component["atmosphere"]
        get(atmosphere, "epoch_utc", nothing) == theirs["epoch_utc"] ||
            _preset_error("$name declares the frozen instant $(get(atmosphere, "epoch_utc", nothing)), but its component $(_preset_label(component)) is frozen at $(theirs["epoch_utc"]).")
        all(key -> get(atmosphere, key, nothing) == theirs[key], ("equatorial_radius_m", "polar_radius_m")) ||
            _preset_error("$name declares reference-ellipsoid radii that differ from those of its component $(_preset_label(component)).")
        _preset_stored_winds(component) ||
            _preset_error("$name needs stored winds on both sides of its handover, but its component $(_preset_label(component)) stores none.")
    end
    get(atmosphere, "winds_available", nothing) === true && get(atmosphere, "wind_returned", nothing) == "stored_components" &&
        get(atmosphere, "wind_convention", nothing) == "stored_nominal_local_ENU" &&
        get(atmosphere, "wind_boolean_suppresses_stored_components", nothing) === false && all(in(_preset_fields(entry)), _WIND_FIELDS) ||
        _preset_error("$name returns its components' stored winds: it must declare winds_available = true, wind_returned = " *
            "\"stored_components\", wind_convention = \"stored_nominal_local_ENU\", wind_boolean_suppresses_stored_components = false, " *
            "and the wind fields $(join(_WIND_FIELDS, ", ")).")
    handover = entry["handover"]; h = get(handover, "height_m", nothing)
    h isa Real && !(h isa Bool) && isfinite(h) && h > 0 || _preset_error("$name needs a finite, positive handover height_m.")
    get(handover, "reference", nothing) == "ellipsoidal" ||
        _preset_error("The handover height of $subject must be measured above the reference ellipsoid (reference = \"ellipsoidal\").")
    get(handover, "rule", nothing) == _COMPOSITE_RULE || _preset_error("The handover rule of $subject must be \"$_COMPOSITE_RULE\": " *
        "the lower component below the handover height and the upper one at or above it, with nothing blended.")
    floor_m, ceiling_m = upper["domain"]["height_m"]
    floor_m <= h < ceiling_m || _preset_error("$name hands over at $h m, outside the range of its upper component $(_preset_label(upper)), " *
        "from $floor_m m (inclusive) to $ceiling_m m.")
    # Necessary only: the lower top is an areoid height, which the catalog cannot compare with an ellipsoidal one everywhere.
    h < lower["domain"]["top_areoid_height_m"] || _preset_error("$name hands over at $h m above the ellipsoid, not below the top of its lower " *
        "component $(_preset_label(lower)), $(lower["domain"]["top_areoid_height_m"]) m areoid height.")
    domain = entry["domain"]
    get(domain, "outside", nothing) == "error" && get(domain, "longitude_period_deg", nothing) == 360 ||
        _preset_error("Named presets require errors outside their supported domain.")
    get(domain, "height_top_m", nothing) == ceiling_m || _preset_error("The domain top of $subject differs from its upper component's ceiling, $ceiling_m m.")
    return nothing
end

function _preset_grid_domain(entry)
    domain = entry["domain"]; axes = entry["axes"]
    foreach(axis -> _preset_axis_nodes(axes[axis]), ("altitude", "latitude", "longitude"))
    haskey(axes["longitude"], "nodes") && _preset_error("Preset longitude must use start/step/count so the periodic seam is explicit.")
    bounds(axis, scale) = _preset_axis_bounds(axes[axis]) .* scale
    bounds("altitude",1000) == domain["height_m"] || _preset_error("Preset altitude domain differs from its axis.")
    bounds("latitude",1) == domain["latitude_deg"] || _preset_error("Preset latitude domain differs from its axis.")
    axes["longitude"]["start"] == 0 && axes["longitude"]["count"]*axes["longitude"]["step"] == 360 && domain["longitude_period_deg"] == 360 || _preset_error("Preset longitude must cover one periodic revolution without a duplicate seam.")
    domain["outside_height"] == domain["outside_latitude"] == "error" && domain["above_grid_vacuum"] === false || _preset_error("Named presets require errors outside their stored altitude and latitude domain.")
    return nothing
end
# A near-surface preset has no grid axes. Its domain is terrain-following: a planetocentric latitude limit, a
# minimum clearance above the local surface, an areoid-height top and a volcano-flank refusal, all errors outside.
function _preset_near_surface_domain(entry)
    domain = entry["domain"]; atmosphere = entry["atmosphere"]; required = entry["required_metadata"]
    format = required["format"]
    format in (_NEAR_SURFACE_FORMAT, _NEAR_SURFACE_WIND_FORMAT) ||
        _preset_error("Near-surface presets require payload format $_NEAR_SURFACE_FORMAT or $_NEAR_SURFACE_WIND_FORMAT.")
    haskey(entry, "axes") && _preset_error("A near-surface preset declares a terrain-following domain, not grid axes.")
    limits = ("planetocentric_latitude_max_deg", "minimum_clearance_m", "top_areoid_height_m", "surface_height_refused_at_or_above_m")
    all(key -> domain[key] isa Real && !(domain[key] isa Bool) && isfinite(domain[key]), limits) ||
        _preset_error("Near-surface preset domain limits must be finite numbers.")
    0 < domain["planetocentric_latitude_max_deg"] <= 90 && domain["minimum_clearance_m"] > 0 && domain["top_areoid_height_m"] > 0 ||
        _preset_error("Near-surface preset domain limits are out of range.")
    support = required["support"]
    support["phic_max_deg"] == domain["planetocentric_latitude_max_deg"] && 1000*support["min_clearance_km"] == domain["minimum_clearance_m"] &&
        1000*support["top_areoid_km"] == domain["top_areoid_height_m"] && 1000*support["zs_refuse_km"] == domain["surface_height_refused_at_or_above_m"] ||
        _preset_error("Near-surface preset domain differs from its payload support limits.")
    domain["longitude_period_deg"] == 360 || _preset_error("Near-surface preset longitude must be periodic over 360 degrees.")
    domain["outside"] == "error" || _preset_error("Named presets require errors outside their supported domain.")
    fields = _preset_fields(entry)
    if format == _NEAR_SURFACE_FORMAT
        get(atmosphere, "winds_available", nothing) === false && get(atmosphere, "wind_returned", nothing) == "zero_vector" && !any(in(fields), _WIND_FIELDS) ||
            _preset_error("A near-surface preset with payload format $format stores no winds: it must declare winds_available = false " *
                "and wind_returned = \"zero_vector\", and list no wind fields.")
    else
        get(atmosphere, "winds_available", nothing) === true && get(atmosphere, "wind_returned", nothing) == "stored_components" &&
            get(atmosphere, "wind_convention", nothing) == "stored_nominal_local_ENU" &&
            get(atmosphere, "wind_boolean_suppresses_stored_components", nothing) === false && all(in(fields), _WIND_FIELDS) ||
            _preset_error("A near-surface preset with payload format $format returns its stored winds: it must declare " *
                "winds_available = true, wind_returned = \"stored_components\", wind_convention = \"stored_nominal_local_ENU\", " *
                "wind_boolean_suppresses_stored_components = false, and the wind fields $(join(_WIND_FIELDS, ", ")).")
    end
    return nothing
end

"""
    available_surrogate_presets(; catalog_file=...)

List exact named surrogate versions and their publication status without downloading
or loading data. Planet-only and moving `latest` selection are intentionally absent.
"""
function available_surrogate_presets(; catalog_file::AbstractString=_SURROGATE_CATALOG)
    catalog, _ = _preset_catalog(catalog_file)
    return [Dict("id" => p["id"], "version" => p["version"], "planet" => p["planet"],
                 "release_enabled" => p["release_enabled"], "description" => get(p, "description", ""))
            for p in catalog["presets"]]
end

"""Exact-byte resolution of a named preset. Construct with `surrogate_preset_model`."""
struct SurrogatePresetResolution
    file::String
    expected_sha256::String
    planet::String
    provenance::Dict{String,Any}
    contract::Dict{String,Any}
end

function _verify_preset_file(file, payload)
    isfile(file) || _preset_error("Preset payload is missing or not a regular file: '$file'. An explicit file never falls back to a download.")
    filesize(file) > 0 || _preset_error("Preset payload is empty: '$file'. Retrieve the complete artifact.")
    header = open(io -> read(io, min(256, filesize(file))), file)
    startswith(String(header), "version https://git-lfs.github.com/spec/v1") &&
        _preset_error("Preset payload is an undownloaded Git LFS pointer: '$file'. Retrieve the actual payload bytes.")
    filesize(file) == payload["bytes"] || _preset_error("Preset payload size mismatch at '$file'. Expected $(payload["bytes"]) bytes; do not substitute a different grid.")
    digest = _preset_digest(file)
    digest == payload["sha256"] || _preset_error("Preset payload SHA256 mismatch at '$file'. Expected $(payload["sha256"]), got $digest. Remove a corrupt managed artifact or correct the explicit file/override; no fallback is attempted.")
    return file
end
function _preset_artifact(entry, artifacts_file; offline)
    name = entry["artifact_name"]
    isfile(artifacts_file) || _preset_error("Preset '$name' has no published artifact metadata.")
    meta = Artifacts.artifact_meta(name, String(artifacts_file); pkg_uuid=_SURROGATE_PACKAGE_UUID)
    meta === nothing && _preset_error("Preset '$name' has no published artifact binding.")
    get(meta, "lazy", false) === true || _preset_error("Preset '$name' must be a lazy artifact.")
    tree = get(meta, "git-tree-sha1", "")
    occursin(r"^[0-9a-f]{40}$", tree) || _preset_error("Invalid artifact tree identity for '$name'.")
    tree == entry["distribution"]["git_tree_sha1"] || _preset_error("Artifact tree identity differs from the catalog for '$name'.")
    downloads = get(meta, "download", Any[])
    expected = entry["distribution"]["archive_sha256"]
    !isempty(downloads) && all(d -> get(d, "sha256", "") == expected &&
        get(d, "url", "") in entry["distribution"]["urls"], downloads) ||
        _preset_error("Artifact download identities differ from the catalog for '$name'.")
    hash = Base.SHA1(tree)
    installed = false
    if !Artifacts.artifact_exists(hash)
        offline && _preset_error("Preset $(entry["id"])@$(entry["version"]) is not installed and offline=true. Fetch this exact preset on a connected machine with assets fetch --preset $(entry["id"]) --version $(entry["version"]), then retry offline.")
        label = "$(entry["id"])@$(entry["version"])"
        @info "Installing atmosphere preset $label ($(round(entry["payload"]["bytes"] / 1e6; digits=1)) MB payload) into the Julia artifact store from $(join((d["url"] for d in downloads), ", "))"
        # Pkg first asks its package server, which does not host this artifact, and
        # then uses the Artifacts.toml URL. Its "Downloading"/"Failure" status lines
        # go to this buffer instead of the terminal and are reported only on failure.
        pkg_output = IOBuffer()
        try
            Pkg.Artifacts.ensure_artifact_installed(name, String(artifacts_file); pkg_uuid=_SURROGATE_PACKAGE_UUID, io=pkg_output)
        catch err
            err isa InterruptException && rethrow()
            _preset_error("Could not install preset $label. No native fallback is available.\n$(_pkg_failure_detail(err, pkg_output))")
        end
        installed = true
    end
    return Artifacts.artifact_path(hash), tree, installed
end
function _pkg_failure_detail(err, pkg_output)
    detail = rstrip(sprint(showerror, err))
    output = rstrip(String(take!(pkg_output)))
    return isempty(output) ? detail : "$detail\nPkg output:\n$output"
end

# The keywords of resolve_surrogate_preset, which surrogate_preset_model forwards.
const _PRESET_RESOLUTION_KEYWORDS = (:version, :planet, :file, :offline, :allow_unreleased, :catalog_file, :artifacts_file)

"""
    resolve_surrogate_preset(id; version, planet="Mars", file="", offline=false,
                             allow_unreleased=false)

Resolve one exact preset and verify the payload size and SHA256 before any model is
constructed. Uses lazy Julia artifacts and honors artifact overrides, which must
still contain the exact pinned payload bytes. An explicit `file` never falls back.
An unreleased preset is available only with both an explicit file and
`allow_unreleased=true`; this marks local development, not public acceptance.
`catalog_file` and `artifacts_file` are advanced trusted-catalog overrides.
Resolution performs no native GRAM initialization and never selects native fallback.
A first installation logs its source and, once verified, its location.
A composite preset has no payload of its own and is refused here: build it with
`surrogate_preset_model`, or resolve each of its components.
"""
function resolve_surrogate_preset(id::AbstractString; version::AbstractString,
    planet::AbstractString="Mars", file::AbstractString="", offline::Bool=false,
    allow_unreleased::Bool=false, catalog_file::AbstractString=_SURROGATE_CATALOG,
    artifacts_file::AbstractString=_SURROGATE_ARTIFACTS)
    catalog, catalog_sha = _preset_catalog(catalog_file)
    candidates = filter(p -> p["id"] == id && p["version"] == version, catalog["presets"])
    length(candidates) == 1 || _preset_error("Unknown surrogate preset '$id@$version'. Inspect available_surrogate_presets(); select an exact listed version.")
    entry = only(candidates); key = _preset_planet(planet)
    key == _preset_planet(entry["planet"]) || _preset_error("Preset '$id@$version' is for $(entry["planet"]), not $planet.")
    _preset_kind(entry) == _COMPOSITE_KIND && _preset_error("Preset '$id@$version' is a composite of " *
        "$(_composite_component_text(entry)) and has no payload of its own. Build it with surrogate_preset_model, or resolve each component with resolve_surrogate_preset.")
    explicit = !isempty(strip(file))
    released = entry["release_enabled"]
    released || (allow_unreleased && explicit) || _preset_error("Preset '$id@$version' is not published. Public retrieval is disabled until its distribution identities are recorded; local development requires an explicit file and allow_unreleased=true.")
    source = "explicit_file"; tree = ""; installed = false
    resolved = if explicit
        abspath(expanduser(file))
    else
        directory, tree, installed = _preset_artifact(entry, artifacts_file; offline)
        source = "julia_artifact"
        joinpath(directory, entry["payload"]["file"])
    end
    _verify_preset_file(resolved, entry["payload"])
    installed && @info "Installed atmosphere preset $id@$version in $(dirname(resolved)); the archive and payload SHA256 checksums match the catalog."
    provenance = Dict{String,Any}(
        "backend" => _preset_kind(entry) == "gram_grid" ? "gram_grid_surrogate" : "gram_near_surface_surrogate",
        "preset_kind" => _preset_kind(entry), "preset_id" => String(id), "preset_version" => String(version),
        "planet" => key, "source_sha256" => entry["payload"]["sha256"], "payload_bytes" => entry["payload"]["bytes"],
        "catalog_sha256" => catalog_sha, "catalog_revision" => catalog["catalog_revision"],
        "resolution" => source, "artifact_git_tree_sha1" => tree, "resolved_file" => resolved,
        "public_release" => released, "local_development" => !released,
        "distribution" => deepcopy(get(entry, "distribution", Dict{String,Any}())),
        "generation_provenance" => deepcopy(get(entry, "generation_provenance", Dict{String,Any}())),
        "atmosphere" => deepcopy(entry["atmosphere"]), "domain" => deepcopy(entry["domain"]),
        "generation" => deepcopy(entry["required_metadata"]["generation_config"]),
        "support" => deepcopy(entry["support"]))
    return SurrogatePresetResolution(resolved, entry["payload"]["sha256"], key, provenance, deepcopy(entry))
end

function _preset_metadata_subset(actual, required, path="metadata")
    actual isa AbstractDict || _preset_error("Preset $path must be a metadata table.")
    for (key, expected) in required
        haskey(actual, key) || _preset_error("Preset $path.$key is missing.")
        if expected isa AbstractDict
            _preset_metadata_subset(actual[key], expected, "$path.$key")
        else
            isequal(actual[key], expected) || _preset_error("Preset $path.$key does not match the catalog contract.")
        end
    end
end
function _validate_preset_model(model, resolution)
    core = model.core; entry = resolution.contract
    core.source_sha256 == resolution.expected_sha256 || _preset_error("Constructed grid identity differs from the selected preset.")
    core.above_grid === :error || _preset_error("Named presets require errors outside altitude coverage.")
    _preset_metadata_subset(core.metadata, entry["required_metadata"])
    grid = core.surrogate; axes = entry["axes"]
    grid.planet_name == resolution.planet || _preset_error("Grid planet differs from the requested preset.")
    for (name, actual, scale) in (("altitude",grid.alt_nodes_m,1000.0), ("latitude",grid.lat_nodes_rad,pi/180), ("longitude",grid.lon_nodes_rad,pi/180))
        expected = _preset_axis_nodes(axes[name]) .* scale
        length(actual)==length(expected) && all(isapprox.(actual,expected;rtol=8eps(Float64),atol=0.0)) ||
            _preset_error("Preset $name grid axis differs from its declared domain and spacing.")
    end
    return model
end

function _validate_near_surface_preset_model(model, resolution)
    core = model.core; entry = resolution.contract
    core.source_sha256 == resolution.expected_sha256 || _preset_error("Constructed payload identity differs from the selected preset.")
    _preset_metadata_subset(core.metadata, entry["required_metadata"])
    evaluator = core.core; domain = entry["domain"]; atmosphere = entry["atmosphere"]
    evaluator.phic_max_deg == domain["planetocentric_latitude_max_deg"] && 1000*evaluator.min_clearance_km == domain["minimum_clearance_m"] &&
        1000*evaluator.top_km == domain["top_areoid_height_m"] && 1000*evaluator.zs_refuse_km == domain["surface_height_refused_at_or_above_m"] ||
        _preset_error("Near-surface payload support limits differ from the preset's declared domain.")
    1000*evaluator.a_km == atmosphere["equatorial_radius_m"] && 1000*evaluator.b_km == atmosphere["polar_radius_m"] ||
        _preset_error("Near-surface payload radii differ from the preset's declared datum.")
    # The payload format is checked through the required metadata above; the loaded model must also agree with the
    # declared winds (a wind layer exactly when the preset declares stored winds).
    _near_surface_model_has_winds(model) === atmosphere["winds_available"] ||
        _preset_error("The loaded near-surface payload $(atmosphere["winds_available"] ? "has no wind layer" : "has a wind layer"), " *
            "contrary to the preset's wind declaration.")
    return model
end

# The kind a request selects, read before any download so a missing model API fails first. An unknown
# id or version falls through to resolve_surrogate_preset, which reports it.
function _preset_requested_kind(id, kwargs)
    catalog, _ = _preset_catalog(get(kwargs, :catalog_file, _SURROGATE_CATALOG))
    version = get(kwargs, :version, nothing)
    for entry in catalog["presets"]
        entry["id"] == id && (version === nothing || entry["version"] == version) && return _preset_kind(entry)
    end
    return "gram_grid"
end

# The keyword constructors come from the GRAMSuite extension; without it a construction fails with a MethodError
# naming keywords the user never passed.
function _require_preset_constructor(near_surface::Bool)
    hasmethod(near_surface ? GRAMNearSurfaceAtmosphereModel : GRAMGridAtmosphereModel, Tuple{}) ||
        _preset_error("surrogate_preset_model needs the public GRAMSuite " *
        "Julia package, which provides the $(near_surface ? "near-surface" : "grid") atmosphere: run `import GRAMSuite` first (the Odyssey example " *
        "environment installs it). No native GRAM installation is used.")
    return nothing
end

# A composite's components, lower first, and their description in messages.
_composite_components(entry) = sort(entry["components"]; by = c -> c["role"])
_composite_component_text(entry) = join(["$(c["id"])@$(c["version"]) ($(c["role"]))" for c in _composite_components(entry)], " and ")
_composite_file_error(id, version) = _preset_error("Preset '$id@$version' is a composite and has no payload of its own, so an " *
    "explicit file cannot select it. To use local payload files, build each component from its own preset and file and " *
    "compose them with CombinedAtmosphereModel; that composition carries no named composite contract.")

# The components of a released composite entry as (role, id, version), lower first; an unreleased composite is refused.
function _composite_entry_components(entry)
    entry["release_enabled"] || _preset_error("Preset $(_preset_label(entry)) is not published. A composite preset has no " *
        "local-development path: build its components and compose them with CombinedAtmosphereModel.")
    return [(role = c["role"], id = c["id"], version = c["version"]) for c in _composite_components(entry)]
end
# The same for a preset named by id and version, or `nothing` for any other entry (an unknown id or version included,
# which resolution reports).
function _composite_preset_components(id, version; catalog_file::AbstractString=_SURROGATE_CATALOG)
    catalog, _ = _preset_catalog(catalog_file)
    for entry in catalog["presets"]
        entry["id"] == id && entry["version"] == version && _preset_kind(entry) == _COMPOSITE_KIND &&
            return _composite_entry_components(entry)
    end
    return nothing
end

# A composite preset: each component built and verified through its own entry, composed at the entry's handover.
function _composite_preset_model(id; version::AbstractString, planet::AbstractString="Mars", file::AbstractString="",
        offline::Bool=false, allow_unreleased::Bool=false, catalog_file::AbstractString=_SURROGATE_CATALOG,
        artifacts_file::AbstractString=_SURROGATE_ARTIFACTS)
    isempty(strip(file)) || _composite_file_error(id, version)
    catalog, catalog_sha = _preset_catalog(catalog_file)
    entry = only(filter(p -> p["id"] == id && p["version"] == version, catalog["presets"]))
    key = _preset_planet(planet)
    key == _preset_planet(entry["planet"]) || _preset_error("Preset '$id@$version' is for $(entry["planet"]), not $planet.")
    components = _composite_entry_components(entry)
    # The catalog fixes the kinds: a near-surface lower component and a grid upper one.
    _require_preset_constructor(true); _require_preset_constructor(false)
    lower, upper = (surrogate_preset_model(c.id; version=c.version, planet, offline, catalog_file, artifacts_file) for c in components)
    model = CombinedAtmosphereModel(lower, upper; handover_height_m=entry["handover"]["height_m"])
    model.metadata[_PRESET_METADATA_KEY] = Dict{String,Any}(
        "preset_kind" => _COMPOSITE_KIND, "preset_id" => String(id), "preset_version" => String(version), "planet" => key,
        "components" => [Dict{String,Any}("role" => c.role, "id" => c.id, "version" => c.version) for c in components],
        "catalog_sha256" => catalog_sha, "catalog_revision" => catalog["catalog_revision"],
        "public_release" => true, "local_development" => false,
        "atmosphere" => deepcopy(entry["atmosphere"]), "domain" => deepcopy(entry["domain"]), "support" => deepcopy(entry["support"]))
    return model
end

"""
    surrogate_preset_model(id; version, planet="Mars", file="", offline=false, ...)

Build the model a verified named preset declares: the existing `GRAMGridAtmosphereModel`
for a grid preset, or `GRAMNearSurfaceAtmosphereModel` for a near-surface preset. Load
`GRAMSuite` first. Validates the retained generation settings, datum and all grid
axes (or the near-surface domain limits) before use. Stored winds and the frozen-time
behavior are unchanged. A near-surface preset with payload format 1 (the published
versions 1.0.0 and 1.1.0) stores no winds and returns a zero wind vector; one with
format 2 returns its stored winds, and its catalog entry and payload must agree on
that. Queries outside the declared domain fail; longitude is periodic.
A composite preset (kind `composite_by_height`) builds each of its two components
through its own entry, as above, and composes them with `CombinedAtmosphereModel` at
the entry's handover height. It has no payload of its own, so it takes no `file`, and
an unreleased composite is refused.
Use `atmosphere_provenance(model)` to retain the selected contract in run outputs.
Without `GRAMSuite` loaded it raises an `ArgumentError` before resolving anything.
Accepts the keywords of `resolve_surrogate_preset` only. A named preset fixes its
domain policy, so grid options such as `above_grid` raise an `ArgumentError`;
construct the generic `GRAMGridAtmosphereModel` directly for another policy.
"""
function surrogate_preset_model(id::AbstractString; kwargs...)
    unsupported = [key for key in keys(kwargs) if key ∉ _PRESET_RESOLUTION_KEYWORDS]
    isempty(unsupported) || _preset_error("surrogate_preset_model does not accept the keyword(s) $(join(unsupported, ", ")). " *
        "A named preset fixes its domain policy: queries outside its stored altitude and latitude domain fail. " *
        "For another policy, such as above_grid=:vacuum, construct the generic GRAMGridAtmosphereModel directly, " *
        "for example with surrogate_file=resolve_surrogate_preset(id; version).file; that model carries no named preset contract. " *
        "Accepted keywords: $(join(_PRESET_RESOLUTION_KEYWORDS, ", ")).")
    kind = _preset_requested_kind(id, kwargs)
    kind == _COMPOSITE_KIND && return _composite_preset_model(id; kwargs...)
    near_surface = kind == "gram_near_surface_scalars"
    _require_preset_constructor(near_surface)
    resolved = resolve_surrogate_preset(id; kwargs...)
    model = if near_surface
        _validate_near_surface_preset_model(GRAMNearSurfaceAtmosphereModel(; planet=resolved.planet,
            surrogate_file=resolved.file, expected_sha256=resolved.expected_sha256), resolved)
    else
        _validate_preset_model(GRAMGridAtmosphereModel(; planet=resolved.planet, surrogate_file=resolved.file,
            expected_sha256=resolved.expected_sha256, above_grid=:error), resolved)
    end
    model.core.metadata[_PRESET_METADATA_KEY] = deepcopy(resolved.provenance)
    return model
end

"""
    atmosphere_provenance(model)

Return an independent metadata dictionary identifying an atmosphere backend. Named
surrogates include their immutable preset/version, payload and catalog identities,
frozen epoch, coordinates, support domain and limitations. A composition records its
handover rule and both components' records; a named composite preset adds its own
preset/version, components, catalog identity and support. Custom models are
identified by their actual type without inventing a preset or native-free claim.
"""
atmosphere_provenance(model::AbstractDensityModel) = Dict{String,Any}(
    "backend" => "configured_density_model", "model_type" => string(typeof(model)))
function atmosphere_provenance(model::GRAMNearSurfaceAtmosphereModel)
    metadata = model.core.metadata
    haskey(metadata, _PRESET_METADATA_KEY) && return deepcopy(metadata[_PRESET_METADATA_KEY])
    return Dict{String,Any}("backend" => "gram_near_surface_surrogate", "source_sha256" => model.core.source_sha256,
        "preset_status" => "user_supplied_payload_without_named_preset_contract")
end
function atmosphere_provenance(model::GRAMGridAtmosphereModel)
    metadata = model.core.metadata
    haskey(metadata, _PRESET_METADATA_KEY) && return deepcopy(metadata[_PRESET_METADATA_KEY])
    return Dict{String,Any}("backend" => "gram_grid_surrogate", "source_sha256" => model.core.source_sha256,
        "preset_status" => "user_supplied_grid_without_named_preset_contract", "above_grid" => string(model.core.above_grid))
end
function atmosphere_provenance(model::CombinedAtmosphereModel)
    record = Dict{String,Any}("backend" => "combined_native_free_snapshot",
        "handover_height_m" => model.handover_height_m, "handover_reference" => "height above the reference ellipsoid",
        "rule" => "lower component below the handover height, upper component at or above it; each keeps its own domain errors and winds",
        "lower" => atmosphere_provenance(model.lower), "upper" => atmosphere_provenance(model.upper))
    haskey(model.metadata, _PRESET_METADATA_KEY) && merge!(record, deepcopy(model.metadata[_PRESET_METADATA_KEY]))
    return record
end
