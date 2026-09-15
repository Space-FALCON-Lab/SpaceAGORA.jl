---
id: cli.assets_load_asset_manifest
label: load_asset_manifest
kind: function
source:
  file: src/cli/assets.jl
  symbol: load_asset_manifest
  lines:
  - 48
  - 48
inputs:
- id: repo_root
  type: String
  units: n/a
  required: false
  description: Keyword argument `repo_root` (default `REPO_ROOT`).
- id: manifest_path
  type: String
  units: n/a
  required: false
  description: Keyword argument `manifest_path` (default `joinpath(repo_root, "data",
    "assets_manifest.toml")`).
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
  type: Any
  units: n/a
  description: Return value of `load_asset_manifest`. Returns `entries`.
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

# load_asset_manifest

## Purpose
Parses `data/assets_manifest.toml` into a `Vector{AssetManifestEntry}`, applying defaults for optional keys, so that `check_assets`, `render_asset_manifest` and `setup_open_assets` all share one source of truth for which asset roots the repository expects.

## Design & Implementation
Keyword arguments are `repo_root::String=REPO_ROOT` and `manifest_path::String=joinpath(repo_root, "data", "assets_manifest.toml")`. It calls `TOML.parsefile(manifest_path)`, then iterates `get(raw, "asset", Any[])`, the array of `[[asset]]` tables. For each table it requires `name` and `scope` (indexing with `entry["name"]` throws `KeyError` if absent) and reads `relative_path` (default `"."`), `kind` (default `"directory"`), `required` (default `false`, converted with `Bool`), `licensing` (default `"unspecified"`) and `detail` (default `""`), converting each with `String` or `Bool` before pushing a new `AssetManifestEntry`. The vector is returned in manifest order.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `repo_root` | String | n/a | no | Keyword argument `repo_root` (default `REPO_ROOT`). |
| in | `manifest_path` | String | n/a | no | Keyword argument `manifest_path` (default `joinpath(repo_root, "data", "assets_manifest.toml")`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `load_asset_manifest`. Returns `entries`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.assets_setup_open_assets|setup_open_assets]] · `callees` → `callers` · call · `src/cli/assets.jl:135-135`
- [[cli.run_cli|run_cli]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:197-197`
- [[misc.assets_check_assets|check_assets]] · `callees` → `callers` · call · `src/cli/assets.jl:84-84`

**Downstream**

- `callees` → [[cli.assets_assetmanifestentry|AssetManifestEntry]] · `callers` · call · `src/cli/assets.jl:52-52`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/cli/assets.jl:52-52`
<!-- vulcan:connections:end -->

## Limitations
A missing manifest file raises the `SystemError` from `TOML.parsefile` uncaught, and a table lacking `name` or `scope` raises `KeyError` with no indication of which entry failed. Values of the wrong TOML type (for example `required = "yes"`) throw `MethodError` from the `Bool` conversion. Unknown keys are ignored silently, so misspelt optional keys revert to defaults without warning. The `repo_root` keyword only feeds the default manifest path; it does not validate that the manifest actually belongs to that root.

## Provenance
Mapped from `src/cli/assets.jl` line 48.
