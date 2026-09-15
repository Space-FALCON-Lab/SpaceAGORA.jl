---
id: input.station_cad
label: Station geometry (STL / point cloud)
kind: external
inputs: []
outputs:
- id: geometry_file
  type: STL / CSV
  units: n/a
  description: Target station mesh or sampled point cloud for RPO clearance checks.
tags:
- master-flow
charts:
- master
origin: agent
---

# Station geometry (STL / point cloud)

## Purpose
The target-station geometry an RPO scenario plans around: the Gateway core STL mesh, or a demo point cloud, under `data/rpo/station_geometry`.

## Design & Implementation
Resolved and parsed by `RPOStationAssets` (`src/assets/rpo_station_assets.jl`), which reads binary or ASCII STL, centres the mesh on its bounding box, and samples an area-weighted point cloud with a fixed seed so clearance statistics are reproducible. The HYPR planner's clearance cost consumes the cloud.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `geometry_file` | STL / CSV | n/a | — | Target station mesh or sampled point cloud for RPO clearance checks. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `geometry_file` → [[flow.gnc|Guidance, navigation & control]] · `geometry_file` · dataflow · `src/assets/rpo_station_assets.jl`
<!-- vulcan:connections:end -->

## Limitations
Only two asset kinds are wired in and the loaders re-parse the STL on every call; the keep-out model is a sphere plus the chaser's bounding sphere, so the mesh detail matters only through the point cloud's density.
