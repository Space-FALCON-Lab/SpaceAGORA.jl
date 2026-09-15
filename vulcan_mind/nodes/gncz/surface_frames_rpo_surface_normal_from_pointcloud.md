---
id: gncz.surface_frames_rpo_surface_normal_from_pointcloud
label: rpo_surface_normal_from_pointcloud
kind: function
source:
  file: src/gnc/navigation/rpo_nav/distances/surface_frames.jl
  symbol: rpo_surface_normal_from_pointcloud
  lines:
  - 2
  - 9
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: NavigationHooks namespace providing the combined reference geometry
    and the nearest station point query.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: normal
  type: SVector{3, Float64}
  units: dimensionless
  description: Unit outward surface normal estimated at the station point nearest
    the query, falling back to the body x axis when the query sits on the surface.
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

# rpo_surface_normal_from_pointcloud

## Purpose
`rpo_surface_normal_from_pointcloud` estimates the local outward normal of the target station surface without a mesh, using only the point cloud that the navigation geometry already carries. The normal defines the approach direction for a standoff goal and the orientation of any surface-relative frame a planner constructs.

## Theory & Math
The estimate is the normalised vector from the nearest station sample to the query point, $\hat{n} = (p - q)/\lVert p - q \rVert$. For a query outside a locally smooth surface this direction converges to the true outward normal as the sample spacing shrinks, because the nearest-point map is the projection onto the surface and the residual is orthogonal to it in the limit.

## Model & Assumptions
The query point is assumed to lie outside the structure; a point inside would produce an inward-pointing estimate with no warning. When the query coincides with the nearest sample to within machine epsilon, the residual carries no direction, so the routine returns the body frame x axis as a defined but arbitrary fallback rather than dividing by zero or propagating a not-a-number result.

## Design & Implementation
The nearest-point query returns the witness point, its distance, and its index, of which only the witness point is used here. Working in static three-vectors keeps the estimate allocation free so it can be evaluated at every sample of a candidate approach path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | NavigationHooks namespace providing the combined reference geometry and the nearest station point query. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `normal` | SVector{3, Float64} | dimensionless | — | Unit outward surface normal estimated at the station point nearest the query, falling back to the body x axis when the query sits on the surface. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[gnc.mesh_distance_nearest_station_point|nearest_station_point]] · `callers` · call · `src/gnc/navigation/rpo_nav/distances/surface_frames.jl:3-3`
<!-- vulcan:connections:end -->

## Limitations
The estimate is first order and degrades near edges, corners, and thin structures, where the nearest sample may sit on a different face than the one being approached. It uses a single sample rather than a local neighbourhood fit, so it is sensitive to point-cloud noise and to the cloud spacing, and the epsilon fallback silently hides a degenerate query.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/distances/surface_frames.jl:1-10`.
