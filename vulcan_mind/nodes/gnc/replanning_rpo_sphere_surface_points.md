---
id: gnc.replanning_rpo_sphere_surface_points
label: rpo_sphere_surface_points
kind: function
source:
  file: src/gnc/guidance/rpo/hypr_planning/replanning.jl
  symbol: rpo_sphere_surface_points
  lines:
  - 137
  - 137
inputs:
- id: sphere
  type: RPOReplanningSphere
  units: n/a
  required: true
  description: Positional argument `sphere`.
- id: n_points
  type: Integer
  units: n/a
  required: false
  description: Keyword argument `n_points` (default `96`).
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
  description: Return value of `rpo_sphere_surface_points`. Returns `pts`.
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

# rpo_sphere_surface_points

## Purpose
Generates an approximately uniform set of points on the surface of a replanning sphere so it can be appended to the station point cloud that the HYPR clearance checks treat as obstacle geometry.

## Theory & Math
Fibonacci sphere sampling: for $k = 0,\dots,n-1$, $z_k = 1 - \frac{2(k + 1/2)}{n}$, $\rho_k = \sqrt{1 - z_k^2}$, $\theta_k = k\,\pi(3 - \sqrt{5})$, and the surface point is $\mathbf{p}_k = \mathbf{c} + R\,(\rho_k\cos\theta_k,\ \rho_k\sin\theta_k,\ z_k)$ where $\mathbf{c}$ is `center_rtn` and $R$ is `radius_m`.

## Design & Implementation
Signature `rpo_sphere_surface_points(sphere::RPOReplanningSphere; n_points::Integer=96)` returning a `3 x n` `Matrix{Float64}`. It uses the Fibonacci (golden-angle) spiral: for `k = 0..n-1`, `z = 1 - 2(k + 0.5)/n`, `r = sqrt(max(0, 1 - z^2))`, `θ = k * π(3 - sqrt 5)`, and the point is `center_rtn + radius_m * (r cos θ, r sin θ, z)`. Each column is written with a broadcast assignment into the preallocated `zeros(3, n)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sphere` | RPOReplanningSphere | n/a | yes | Positional argument `sphere`. |
| in | `n_points` | Integer | n/a | no | Keyword argument `n_points` (default `96`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_sphere_surface_points`. Returns `pts`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.replanning_rpo_geometry_with_replanning_spheres|rpo_geometry_with_replanning_spheres]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:159-159`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Surface sampling leaves gaps of order `R sqrt(4π/n)` between points (about 0.36 R for n = 96), so a path can slip between samples and report a larger clearance than the true sphere allows; the `safe_distance_m` margin must absorb this. Only the shell is sampled, so a path passing through the interior of a large sphere far from any sample point is not penalised. `n_points < 1` produces an empty matrix without warning.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr_planning/replanning.jl` line 137.
