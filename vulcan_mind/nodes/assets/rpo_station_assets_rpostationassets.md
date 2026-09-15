---
id: assets.rpo_station_assets_rpostationassets
label: RPOStationAssets
kind: module
source:
  file: src/assets/rpo_station_assets.jl
  symbol: RPOStationAssets
  lines:
  - 1
  - 1
inputs:
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
  description: Value produced by this symbol.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- assets
charts:
- assets
origin: agent
---

# RPOStationAssets

## Purpose
The module that locates and loads target-station geometry for RPO scenarios, either as a demo point cloud or as sampled points from a Gateway CAD mesh.

## Design & Implementation
Anchors `_STATION_GEOMETRY_ROOT` at `data/rpo/station_geometry` relative to the source file with `@__DIR__`, so asset resolution does not depend on the working directory. It exports the two path resolvers, the demo point-cloud loader, and the CAD triangle and point-cloud loaders, and keeps STL parsing, centring and area-weighted surface sampling as private functions. Only `LinearAlgebra` and `Random` are imported; STL is parsed by hand rather than through a mesh package.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.assets|RPOStationAssets]] · `api` → `module_api` · call · `src/assets/rpo_station_assets.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only two asset kinds are wired in, `:demo` and `:gateway`, and every resolver ends in an `ArgumentError` telling the caller to add an artifact-backed loader; the module has no registry to extend, so a new station means editing the branch tables.

## Provenance
Mapped from `src/assets/rpo_station_assets.jl` line 1.
