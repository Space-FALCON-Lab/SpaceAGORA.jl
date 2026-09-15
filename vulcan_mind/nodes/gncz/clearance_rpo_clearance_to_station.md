---
id: gncz.clearance_rpo_clearance_to_station
label: rpo_clearance_to_station
kind: function
source:
  file: src/gnc/navigation/rpo_nav/distances/clearance.jl
  symbol: rpo_clearance_to_station
  lines:
  - 2
  - 6
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: NavigationHooks namespace providing the reference geometry types and
    the nearest station point query.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: clearance_report
  type: NamedTuple
  units: m
  description: Signed clearance, raw nearest-point distance, the nearest station point,
    and its index in the station point cloud.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncz
origin: agent
---

# rpo_clearance_to_station

## Purpose
`rpo_clearance_to_station` answers the central safety question of proximity operations, namely how much room remains between the chaser and the target station surface at a given body-frame point. It returns not only the margin but also the witness point that produced it, which lets a planner explain and visualise a violation instead of only detecting one.

## Theory & Math
Clearance is the nearest-surface distance reduced by two inflation terms, $c = d - r_{keepout} - \max(h_x, h_y, h_z)$, where $d$ is the Euclidean distance to the closest station point-cloud sample, $r_{keepout}$ is the configured keepout radius around the station, and the half-extent maximum is the circumscribing radius of the chaser box. Subtracting the largest half extent is the conservative choice, treating the chaser as a sphere that encloses its own bounding box.

## Model & Assumptions
Both the query point and the station geometry live in the same body frame, and the station is represented as a point cloud rather than a closed surface, so clearance is measured to samples and never becomes negative through penetration of a face. The sign of the result is meaningful, with negative values indicating an incursion into the inflated keepout volume.

## Design & Implementation
The companion routine in the same file skips the witness point and returns only the scalar margin using the squared-distance query, and the path statistics routine loops that scalar form over a three-by-N path to report the minimum clearance, the count of samples below a margin, and the violated fraction.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | NavigationHooks namespace providing the reference geometry types and the nearest station point query. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `clearance_report` | NamedTuple | m | — | Signed clearance, raw nearest-point distance, the nearest station point, and its index in the station point cloud. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.path_retiming_rpo_retime_path|rpo_retime_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:166-166`

**Downstream**

- `callees` → [[gnc.mesh_distance_nearest_station_point|nearest_station_point]] · `callers` · call · `src/gnc/navigation/rpo_nav/distances/clearance.jl:3-3`
<!-- vulcan:connections:end -->

## Limitations
Sphere inflation of the chaser discards attitude, so a slender chaser is penalised by its longest dimension in every orientation. Point-cloud sampling means the reported distance is only as accurate as the cloud density, and a coarse cloud can miss a thin protruding structure entirely.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/distances/clearance.jl:1-32`.
