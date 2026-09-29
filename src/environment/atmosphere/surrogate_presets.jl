# Named data selection only. GRAMSuite and the existing adapter own grid mathematics.
import Artifacts
import Pkg
import SHA
import TOML

const _SURROGATE_CATALOG = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "surrogate_presets.toml"))
const _SURROGATE_ARTIFACTS = normpath(joinpath(@__DIR__, "..", "..", "..", "Artifacts.toml"))
const _SURROGATE_PACKAGE_UUID = Base.UUID("afbfb69f-5c0b-4832-b760-43725dff8540")
const _PRESET_METADATA_KEY = "spaceagora_named_surrogate"
_preset_error(message) = throw(ArgumentError(message))
_preset_digest(file) = open(io -> bytes2hex(SHA.sha256(io)), file)

function _preset_planet(planet::AbstractString)
    key = lowercase(strip(planet))
    key in ("earth", "mars", "venus", "jupiter", "uranus", "neptune", "titan") ||
        _preset_error("Unknown surrogate planet '$planet'. Choose an explicit catalog preset and supported planet.")
    return key
end
function _preset_catalog(path)
    isfile(path) || _preset_error("Surrogate preset catalog is missing: $path")
    bytes = read(path)
    catalog = TOML.parse(String(copy(bytes)))
    get(catalog, "schema_version", nothing) == 1 || _preset_error("Unsupported surrogate catalog schema.")
    entries = get(catalog, "presets", nothing)
    entries isa Vector || _preset_error("The surrogate catalog must contain preset entries.")
    seen = Set{Tuple{String,String}}()
    for entry in entries
        id = get(entry, "id", ""); version = get(entry, "version", "")
        occursin(r"^[a-z][a-z0-9_]*$", id) || _preset_error("Invalid surrogate preset ID '$id'.")
        occursin(r"^\d+\.\d+\.\d+$", version) || _preset_error("Preset '$id' needs an exact release version; moving versions are unsupported.")
        key = (id, version)
        key ∉ seen || _preset_error("Duplicate preset/version '$id@$version' in catalog.")
        push!(seen, key)
        _preset_planet(entry["planet"])
        payload = entry["payload"]
        occursin(r"^[0-9a-f]{64}$", payload["sha256"]) || _preset_error("Invalid preset payload SHA256.")
        payload["bytes"] isa Integer && payload["bytes"] > 0 || _preset_error("Invalid preset payload size.")
        file = payload["file"]
        file isa String && !isempty(file) && !isabspath(file) && !occursin('\\', file) &&
            all(x -> x ∉ ("", ".", ".."), split(file, '/')) || _preset_error("Unsafe preset payload filename.")
        get(entry, "release_enabled", nothing) isa Bool || _preset_error("Preset release_enabled must be explicit.")
        get(entry["atmosphere"], "query_elapsed_time_applied", nothing) === false ||
            _preset_error("Only explicitly frozen atmospheric presets are supported.")
        atmosphere = entry["atmosphere"]; domain = entry["domain"]; axes = entry["axes"]
        required = entry["required_metadata"]; generation = required["generation_config"]
        _preset_planet(required["planet"]) == _preset_planet(entry["planet"]) || _preset_error("Preset metadata planet is inconsistent.")
        payload["format"] == required["format"] || _preset_error("Preset payload format is inconsistent.")
        for axis in ("altitude", "latitude", "longitude")
            spec = axes[axis]
            spec["count"] isa Integer && !(spec["count"] isa Bool) && spec["count"] >= 2 || _preset_error("Invalid preset axis count.")
            all(x -> x isa Real && !(x isa Bool) && isfinite(x), (spec["start"], spec["step"])) && spec["step"] > 0 || _preset_error("Invalid preset axis spacing.")
        end
        bounds(axis, scale) = [axes[axis]["start"], axes[axis]["start"] + (axes[axis]["count"]-1)*axes[axis]["step"]] .* scale
        bounds("altitude",1000) == domain["height_m"] || _preset_error("Preset altitude domain differs from its axis.")
        bounds("latitude",1) == domain["latitude_deg"] || _preset_error("Preset latitude domain differs from its axis.")
        axes["longitude"]["start"] == 0 && axes["longitude"]["count"]*axes["longitude"]["step"] == 360 && domain["longitude_period_deg"] == 360 || _preset_error("Preset longitude must cover one periodic revolution without a duplicate seam.")
        domain["outside_height"] == domain["outside_latitude"] == "error" && domain["above_grid_vacuum"] === false || _preset_error("Named presets require errors outside their stored altitude and latitude domain.")
        atmosphere["latitude"] == "geodetic" && generation["is_planetocentric"] === false && atmosphere["height"] == "ellipsoidal" || _preset_error("Preset coordinate datum is inconsistent or unsupported.")
        for (public_name, native_name) in (("equatorial_radius_m","equatorial_radius_km"),("polar_radius_m","polar_radius_km"))
            atmosphere[public_name] == 1000*generation[native_name] || _preset_error("Preset radii are inconsistent.")
        end
        initial = required["initial_time"]
        initial == generation["initial_time"] || _preset_error("Preset epoch metadata is inconsistent.")
        epoch = match(r"^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)Z$", atmosphere["epoch_utc"])
        epoch !== nothing && all(parse(Float64,epoch.captures[i]) == initial[k] for (i,k) in enumerate(("year","month","day","hour","minute","second"))) || _preset_error("Preset frozen epoch differs from its generation metadata.")
        if entry["release_enabled"]
            distribution = entry["distribution"]
            occursin(r"^[0-9a-f]{40}$", distribution["git_tree_sha1"]) && occursin(r"^[0-9a-f]{64}$", distribution["archive_sha256"]) && !isempty(distribution["urls"]) || _preset_error("Released presets require complete artifact publication identities.")
        end
    end
    return catalog, bytes2hex(SHA.sha256(bytes))
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
        _preset_error("Preset payload is an undownloaded Git LFS pointer: '$file'. Retrieve the actual grid bytes.")
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
        @info "Installing atmosphere preset $label ($(round(entry["payload"]["bytes"] / 1e6; digits=1)) MB grid) into the Julia artifact store from $(join((d["url"] for d in downloads), ", "))"
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
    installed && @info "Installed atmosphere preset $id@$version in $(dirname(resolved)); the archive and grid SHA256 checksums match the catalog."
    provenance = Dict{String,Any}(
        "backend" => "gram_grid_surrogate", "preset_id" => String(id), "preset_version" => String(version),
        "planet" => key, "source_sha256" => entry["payload"]["sha256"], "payload_bytes" => entry["payload"]["bytes"],
        "catalog_sha256" => catalog_sha, "catalog_revision" => catalog["catalog_revision"],
        "resolution" => source, "artifact_git_tree_sha1" => tree, "resolved_file" => resolved,
        "public_release" => released, "local_development" => !released,
        "distribution" => deepcopy(entry["distribution"]),
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
        spec = axes[name]
        expected = collect(range(Float64(spec["start"]); step=Float64(spec["step"]), length=spec["count"])) .* scale
        length(actual)==length(expected) && all(isapprox.(actual,expected;rtol=8eps(Float64),atol=0.0)) ||
            _preset_error("Preset $name grid axis differs from its declared domain and spacing.")
    end
    return model
end

"""
    surrogate_preset_model(id; version, planet="Mars", file="", offline=false, ...)

Build the existing `GRAMGridAtmosphereModel` from a verified named preset. Load
`GRAMSuite` first. Validates the retained generation settings, datum and all grid
axes before use. Stored winds and the frozen-time behavior are unchanged. Queries
outside the declared altitude/latitude domain fail; the longitude axis is periodic.
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
    # The keyword constructor comes from the GRAMSuite extension; without it the
    # call below fails with a MethodError naming keywords the user never passed.
    hasmethod(GRAMGridAtmosphereModel, Tuple{}) || _preset_error("surrogate_preset_model needs the public GRAMSuite " *
        "Julia package, which provides the grid atmosphere: run `import GRAMSuite` first (the Odyssey example " *
        "environment installs it). No native GRAM installation is used.")
    resolved = resolve_surrogate_preset(id; kwargs...)
    model = GRAMGridAtmosphereModel(; planet=resolved.planet, surrogate_file=resolved.file,
        expected_sha256=resolved.expected_sha256, above_grid=:error)
    _validate_preset_model(model, resolved)
    model.core.metadata[_PRESET_METADATA_KEY] = deepcopy(resolved.provenance)
    return model
end

"""
    atmosphere_provenance(model)

Return an independent metadata dictionary identifying an atmosphere backend. Named
surrogates include their immutable preset/version, payload and catalog identities,
frozen epoch, coordinates, support domain and limitations. Custom models are
identified by their actual type without inventing a preset or native-free claim.
"""
atmosphere_provenance(model::AbstractDensityModel) = Dict{String,Any}(
    "backend" => "configured_density_model", "model_type" => string(typeof(model)))
function atmosphere_provenance(model::GRAMGridAtmosphereModel)
    metadata = model.core.metadata
    haskey(metadata, _PRESET_METADATA_KEY) && return deepcopy(metadata[_PRESET_METADATA_KEY])
    return Dict{String,Any}("backend" => "gram_grid_surrogate", "source_sha256" => model.core.source_sha256,
        "preset_status" => "user_supplied_grid_without_named_preset_contract", "above_grid" => string(model.core.above_grid))
end
