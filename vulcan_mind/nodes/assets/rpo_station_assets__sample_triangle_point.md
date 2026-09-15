---
id: assets.rpo_station_assets__sample_triangle_point
label: _sample_triangle_point
kind: function
source:
  file: src/assets/rpo_station_assets.jl
  symbol: _sample_triangle_point
  lines:
  - 106
  - 106
inputs:
- id: v1
  type: Any
  units: n/a
  required: true
  description: Positional argument `v1`.
- id: v2
  type: Any
  units: n/a
  required: true
  description: Positional argument `v2`.
- id: v3
  type: Any
  units: n/a
  required: true
  description: Positional argument `v3`.
- id: rng
  type: Any
  units: n/a
  required: true
  description: Positional argument `rng`.
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
  description: Return value of `_sample_triangle_point`. Returns `v1 .+ u .* (v2 .-
    v1) .+ v .* (v3 .- v1)`.
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

# _sample_triangle_point

## Purpose
Draws one point uniformly distributed over the surface of a triangle.

## Theory & Math
With $u, v \sim \mathcal{U}(0,1)$ and the reflection $(u, v) \mapsto (1-u, 1-v)$ applied when $u + v > 1$, the point

$$
p = v_1 + u\,(v_2 - v_1) + v\,(v_3 - v_1)
$$

is uniformly distributed over the triangle with vertices $v_1, v_2, v_3$.

## Design & Implementation
Samples two uniforms `u` and `v` from `rng`, and if their sum exceeds one reflects both through `1 - u`, `1 - v`, folding the unit square onto the lower-left triangle. The point is then `v1 + u (v2 - v1) + v (v3 - v1)`. The fold is what makes the distribution uniform over area rather than concentrated toward `v1`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `v1` | Any | n/a | yes | Positional argument `v1`. |
| in | `v2` | Any | n/a | yes | Positional argument `v2`. |
| in | `v3` | Any | n/a | yes | Positional argument `v3`. |
| in | `rng` | Any | n/a | yes | Positional argument `rng`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_sample_triangle_point`. Returns `v1 .+ u .* (v2 .- v1) .+ v .* (v3 .- v1)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[assets.rpo_station_assets__stl_pointcloud|_stl_pointcloud]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:141-141`
- [[module.assets|RPOStationAssets]] · `api` → `module_api` · call · `src/assets/rpo_station_assets.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Each call allocates the broadcast result as a fresh vector; in the point-cloud sampler this is called `n_points` times with column slices that are themselves copies, so the sampling loop allocates heavily relative to the arithmetic performed.

## Provenance
Mapped from `src/assets/rpo_station_assets.jl` line 106.
