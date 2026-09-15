---
id: cli.assets_assetmanifestentry
label: AssetManifestEntry
kind: struct
source:
  file: src/cli/assets.jl
  symbol: AssetManifestEntry
  lines:
  - 27
  - 27
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Field `name`.
- id: scope
  type: String
  units: n/a
  required: true
  description: Field `scope`.
- id: relative_path
  type: String
  units: n/a
  required: true
  description: Field `relative_path`.
- id: kind
  type: String
  units: n/a
  required: true
  description: Field `kind`.
- id: required
  type: Bool
  units: n/a
  required: true
  description: Field `required`.
- id: licensing
  type: String
  units: n/a
  required: true
  description: Field `licensing`.
- id: detail
  type: String
  units: n/a
  required: false
  description: Field `detail` (default `""`).
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: AssetManifestEntry
  units: n/a
  description: Constructed `AssetManifestEntry` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- cli
charts:
- cli
origin: agent
---

# AssetManifestEntry

## Purpose
Typed representation of one `[[asset]]` table from `data/assets_manifest.toml`. It separates the static declaration of an asset (where it should live, what kind of object it is, whether it is required, and how it is licensed) from the runtime availability check captured by `AssetCheckItem`.

## Design & Implementation
Declared with `Base.@kwdef struct AssetManifestEntry`. Fields: `name::String`, `scope::String`, `relative_path::String` (relative to the repo root, with `""` or `"."` meaning the root itself), `kind::String` (one of `"directory"`, `"file"`, `"builtin"`), `required::Bool`, `licensing::String` (for example `"licensed-external"`, which `setup_open_assets` uses to split entries into open and user-provided groups) and `detail::String` defaulting to `""`. `load_asset_manifest` constructs it from TOML with defaults `relative_path="."`, `kind="directory"`, `required=false`, `licensing="unspecified"`. `render_asset_manifest` prints every field.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Field `name`. |
| in | `scope` | String | n/a | yes | Field `scope`. |
| in | `relative_path` | String | n/a | yes | Field `relative_path`. |
| in | `kind` | String | n/a | yes | Field `kind`. |
| in | `required` | Bool | n/a | yes | Field `required`. |
| in | `licensing` | String | n/a | yes | Field `licensing`. |
| in | `detail` | String | n/a | no | Field `detail` (default `""`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AssetManifestEntry | n/a | — | Constructed `AssetManifestEntry` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.assets_load_asset_manifest|load_asset_manifest]] · `callees` → `callers` · call · `src/cli/assets.jl:52-52`
- [[module.cli|SpaceAGORACLI]] · `api` → `module_api` · call · `src/cli/assets.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`kind` and `licensing` are free strings rather than enums, so a typo in the manifest is only caught when `_manifest_entry_available` throws `ArgumentError` for an unknown kind, and an unknown licensing value silently lands in the open group. `relative_path` is split with `splitpath`, which assumes the manifest uses the host path separator conventions. No validation rejects absolute paths or `..` segments that escape the repository.

## Provenance
Mapped from `src/cli/assets.jl` line 27.
