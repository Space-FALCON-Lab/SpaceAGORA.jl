---
id: gnc.replanning_rpo_geometry_with_replanning_spheres
label: rpo_geometry_with_replanning_spheres
kind: function
source:
  file: src/gnc/guidance/rpo/hypr_planning/replanning.jl
  symbol: rpo_geometry_with_replanning_spheres
  lines:
  - 151
  - 151
inputs:
- id: geometry
  type: RPOReferenceGeometry
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: spheres
  type: AbstractVector{RPOReplanningSphere}
  units: n/a
  required: true
  description: Positional argument `spheres`.
- id: sphere_surface_samples
  type: Integer
  units: n/a
  required: false
  description: Keyword argument `sphere_surface_samples` (default `96`).
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
  type: RPOReferenceGeometry
  units: n/a
  description: Return value of `rpo_geometry_with_replanning_spheres`. Returns `RPOReferenceGeometry(station;
    chaser=geometry.chaser)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# rpo_geometry_with_replanning_spheres

## Purpose
Builds an augmented `RPOReferenceGeometry` whose station point cloud includes surface samples of every active replanning sphere, so that the existing HYPR clearance and planning code treats dynamic obstacles exactly like fixed station structure.

## Design & Implementation
Signature `rpo_geometry_with_replanning_spheres(geometry::RPOReferenceGeometry, spheres::AbstractVector{RPOReplanningSphere}; sphere_surface_samples::Integer=96)`. With no spheres it returns `geometry` unchanged. Otherwise it allocates a `3 x (n_station + n_samples * length(spheres))` matrix, copies `geometry.station.points_body` into the leading columns, then appends `rpo_sphere_surface_points(sphere; n_points=sphere_surface_samples)` for each sphere at a running column offset. A new `RPOStationGeometry` is constructed with the same `keepout_radius_m` and the name suffixed with `"+replanning"`, and wrapped in `RPOReferenceGeometry(station; chaser=geometry.chaser)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `geometry` | RPOReferenceGeometry | n/a | yes | Positional argument `geometry`. |
| in | `spheres` | AbstractVector{RPOReplanningSphere} | n/a | yes | Positional argument `spheres`. |
| in | `sphere_surface_samples` | Integer | n/a | no | Keyword argument `sphere_surface_samples` (default `96`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOReferenceGeometry | n/a | — | Return value of `rpo_geometry_with_replanning_spheres`. Returns `RPOReferenceGeometry(station; chaser=geometry.chaser)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.replanning_rpo_replanning_decision|rpo_replanning_decision]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:216-216`

**Downstream**

- `callees` → [[gnc.replanning_rpo_sphere_surface_points|rpo_sphere_surface_points]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:159-159`
- `callees` → [[gncz.rpo_reference_geometry_rporeferencegeometry|RPOReferenceGeometry]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:168-168`
- `callees` → [[gncz.station_geometry_rpostationgeometry|RPOStationGeometry]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:163-163`
<!-- vulcan:connections:end -->

## Limitations
Sphere points are placed in the station body frame even though sphere centres are defined in RTN; the two coincide only when the station body frame is aligned with RTN, an assumption not checked here. The station `keepout_radius_m` is applied uniformly to sphere samples too, inflating dynamic obstacles by that radius. Any other fields of `RPOStationGeometry` beyond points, keep-out radius and name are dropped by the reconstruction. The matrix is rebuilt on every decision call, allocating `3 x N` floats each time.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr_planning/replanning.jl` line 151.
