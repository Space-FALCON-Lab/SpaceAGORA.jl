---
id: vehicle.robotics_closest_surface_target
label: closest_surface_target
kind: function
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: closest_surface_target
  lines:
  - 295
  - 295
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
- id: reference_position
  type: Any
  units: n/a
  required: true
  description: Positional argument `reference_position`.
- id: standoff_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `standoff_m` (default `0.0`).
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
  description: Return value of `closest_surface_target`. Returns `(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# closest_surface_target

## Purpose
Selects the point on a sampled surface nearest a reference position and derives a standoff target along the outward direction, used to aim the arm at cloth or structure.

## Design & Implementation
Validates the points as a three-by-N matrix with at least one column, then scans every column for the minimum squared distance to `reference_position`. The surface normal is approximated as the unit vector from the surface point toward the reference, falling back to the x axis if they coincide. Returns a named tuple with the standoff `target`, the `surface_point`, the `surface_normal` and the chosen `index`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `reference_position` | Any | n/a | yes | Positional argument `reference_position`. |
| in | `standoff_m` | Real | n/a | no | Keyword argument `standoff_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `closest_surface_target`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/robotics/robotics.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/vehicle/robotics/robotics.jl:315-315`
<!-- vulcan:connections:end -->

## Limitations
The normal is the direction to the reference, not the surface's geometric normal, so for a reference far off-axis the standoff is applied obliquely; the linear scan copies each column into an `SVector`, which is fine for thousands of points but not for dense scans.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 295.
