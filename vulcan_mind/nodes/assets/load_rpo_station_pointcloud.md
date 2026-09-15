---
id: assets.load_rpo_station_pointcloud
label: load_rpo_station_pointcloud
kind: function
source:
  file: src/assets/rpo_station_assets.jl
  symbol: load_rpo_station_pointcloud
  lines:
  - 27
  - 41
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: RPOStationAssets namespace supplying station path resolution and geometry
    conventions.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: pointcloud
  type: Matrix{Float64}
  units: m
  description: Three-by-N station point matrix loaded from the selected repository
    geometry asset.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- assets
- geometry
charts:
- assets
origin: agent
---

# load_rpo_station_pointcloud

## Purpose
`load_rpo_station_pointcloud` loads the text point-cloud representation used by rendezvous-station examples and geometry utilities. It resolves a named station asset, parses rows into numeric coordinates, validates their dimensionality, and returns a matrix convenient for plotting and surface-target selection.

## Theory & Math
The loader performs no geometric fitting. It maps each file row `(x,y,z)` into a column of a `3 × N` matrix. Downstream targeting can compute distances or nearest points from this representation, but the loader preserves the asset coordinates without smoothing or resampling.

## Model & Assumptions
Each nonempty data row is expected to contain exactly three numeric values in a common length unit. The selected asset is assumed to exist under the repository data layout returned by `station_geometry_path`. Coordinate orientation and scale are inherited from the file and are not inferred from metadata.

## Design & Implementation
The function obtains the path from the station-kind helper, reads the file line by line, converts tokens to floating-point values, and rejects rows whose length is not three. It accumulates values into a matrix and returns the transposed shape used by the robotics and plotting callers. No random number generator is involved in this path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | RPOStationAssets namespace supplying station path resolution and geometry conventions. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `pointcloud` | Matrix{Float64} | m | — | Three-by-N station point matrix loaded from the selected repository geometry asset. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[assets.rpo_station_assets_station_geometry_path|station_geometry_path]] · `callers` · call · `src/assets/rpo_station_assets.jl:29-29`
<!-- vulcan:connections:end -->

## Limitations
Missing files and nonnumeric tokens fail during path resolution or parsing. The function does not verify units, orientation, duplicate points, surface closure, or physical scale. Very large assets are loaded in memory, and malformed rows are rejected rather than repaired, so callers need a data-preparation step for third-party geometry.

## Provenance
Mapped from `src/assets/rpo_station_assets.jl:27-41`.
