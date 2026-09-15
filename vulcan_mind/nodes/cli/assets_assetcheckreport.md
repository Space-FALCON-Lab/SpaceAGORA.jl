---
id: cli.assets_assetcheckreport
label: AssetCheckReport
kind: struct
source:
  file: src/cli/assets.jl
  symbol: AssetCheckReport
  lines:
  - 22
  - 22
inputs:
- id: repo_root
  type: String
  units: n/a
  required: true
  description: Field `repo_root`.
- id: items
  type: Vector{AssetCheckItem}
  units: n/a
  required: true
  description: Field `items`.
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
  type: AssetCheckReport
  units: n/a
  description: Constructed `AssetCheckReport` (keyword constructor via @kwdef).
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

# AssetCheckReport

## Purpose
Container returned by `check_assets()` that pairs the repository root against which asset paths were resolved with the list of per-asset availability results. It is the typed object that `render_asset_report` prints and that `setup_open_assets` returns to its caller.

## Design & Implementation
A `Base.@kwdef struct AssetCheckReport` with two fields: `repo_root::String`, the absolute directory used to resolve every manifest `relative_path`, and `items::Vector{AssetCheckItem}` in the same order as the `[[asset]]` tables in `data/assets_manifest.toml`. `check_assets` builds it after iterating `load_asset_manifest` output, calling `_manifest_entry_path` and `_manifest_entry_available` for each entry. The struct is immutable but its `items` vector is a mutable container, so callers can append or filter in place.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `repo_root` | String | n/a | yes | Field `repo_root`. |
| in | `items` | Vector{AssetCheckItem} | n/a | yes | Field `items`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AssetCheckReport | n/a | — | Constructed `AssetCheckReport` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[misc.assets_check_assets|check_assets]] · `callees` → `callers` · call · `src/cli/assets.jl:98-98`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
There is no summary field (count of missing required assets), so callers must scan `items` themselves to decide whether the repository is usable; `render_asset_report` also performs no aggregate. Ordering is defined only by manifest order and nothing deduplicates repeated `name` values. Equality and hashing fall back to the default struct behaviour, which compares vectors element-wise but is not documented as part of the API.

## Provenance
Mapped from `src/cli/assets.jl` line 22.
