---
id: cli.assets_assetcheckitem
label: AssetCheckItem
kind: struct
source:
  file: src/cli/assets.jl
  symbol: AssetCheckItem
  lines:
  - 8
  - 8
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
- id: path
  type: String
  units: n/a
  required: true
  description: Field `path`.
- id: available
  type: Bool
  units: n/a
  required: true
  description: Field `available`.
- id: required
  type: Bool
  units: n/a
  required: true
  description: Field `required`.
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
  type: AssetCheckItem
  units: n/a
  description: Constructed `AssetCheckItem` (keyword constructor via @kwdef).
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

# AssetCheckItem

## Purpose
Immutable record describing the availability of one asset root (a data directory, file or built-in resource) as determined by `check_assets()`. A vector of these forms the body of an `AssetCheckReport`, which the CLI renders for users deciding whether GRAM or other licensed data must be installed.

## Design & Implementation
Defined with `Base.@kwdef struct AssetCheckItem` so instances are built by keyword. Fields are `name::String` (manifest entry name), `scope::String` (which subsystem needs it), `path::String` (absolute resolved path under the repository root), `available::Bool` (result of the filesystem probe), `required::Bool` (whether absence blocks baseline use) and `detail::String` defaulting to `""`, which `check_assets` fills with the licensing tag and any manifest detail text. `_asset_item` is a thin keyword-forwarding constructor for it, and `render_asset_report` reads every field to print a status line of `available`, `missing-required` or `missing-optional`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Field `name`. |
| in | `scope` | String | n/a | yes | Field `scope`. |
| in | `path` | String | n/a | yes | Field `path`. |
| in | `available` | Bool | n/a | yes | Field `available`. |
| in | `required` | Bool | n/a | yes | Field `required`. |
| in | `detail` | String | n/a | no | Field `detail` (default `""`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AssetCheckItem | n/a | — | Constructed `AssetCheckItem` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.assets__asset_item|_asset_item]] · `callees` → `callers` · call · `src/cli/assets.jl:38-38`
- [[module.cli|SpaceAGORACLI]] · `api` → `module_api` · call · `src/cli/assets.jl`

**Downstream**

- `callees` → [[misc.assets_check_assets|check_assets]] · `callers` · call · `src/cli/assets.jl:20-20`
- `callees` → [[spaceagora.spaceagora_check_assets|check_assets]] · `callers` · call · `src/cli/assets.jl:20-20`
<!-- vulcan:connections:end -->

## Limitations
The struct carries no timestamp, so a report can go stale if assets are installed after it was produced. `available` is a plain Bool: a partially populated directory (present but missing files inside) counts as available. `detail` is free text mixing licensing and description, so consumers cannot filter on licensing without string parsing; the separate `AssetManifestEntry` retains the structured `licensing` field.

## Provenance
Mapped from `src/cli/assets.jl` line 8.
