---
id: core.quaternion_utils_dcm_to_quaternion
label: dcm_to_quaternion
kind: function
source:
  file: src/core/numerics/quaternion_utils.jl
  symbol: dcm_to_quaternion
  lines:
  - 97
  - 97
inputs:
- id: dcm
  type: SMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Positional argument `dcm`.
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
  description: Return value of `dcm_to_quaternion`. Returns `q / norm(q)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# dcm_to_quaternion

## Purpose

Recovers an attitude quaternion from a direction cosine matrix, closing the loop with `rot` so that attitudes computed as rotation matrices, for example from pointing-frame constructions, can be stored and propagated in quaternion form.

## Design & Implementation

`dcm_to_quaternion(dcm)` takes an `SMatrix{3,3,Float64}`. It computes the trace `tr`, then picks the branch that maximises numerical conditioning by taking `argmax([dcm[1,1], dcm[2,2], dcm[3,3], tr])`. When the trace wins (`i == 4`) it builds the quaternion directly from the antisymmetric off-diagonal differences together with `1 + tr`. Otherwise it allocates `q = zeros(4)`, sets the dominant component to `1 + 2*dcm[i,i] - tr`, fills the remaining vector components from the symmetric sums `dcm[i,j] + dcm[j,i]`, and fills the scalar slot from an off-diagonal difference selected by wrapped indices. Both branches return `q / norm(q)`.

## Theory & Math

Shepperd's method selects the largest of $R_{11}, R_{22}, R_{33}, \operatorname{tr}(R)$ to avoid dividing by a small number. For the trace branch the unnormalised quaternion is

$$\tilde{q} = \left[R_{23}-R_{32},\; R_{31}-R_{13},\; R_{12}-R_{21},\; 1 + \operatorname{tr}(R)\right]$$

and for the branch where diagonal entry $i$ dominates, $\tilde{q}_i = 1 + 2R_{ii} - \operatorname{tr}(R)$ with the remaining vector entries $\tilde{q}_j = R_{ij} + R_{ji}$. The returned attitude is $q = \tilde{q}/\|\tilde{q}\|$, which is unit by construction.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dcm` | SMatrix{3, 3, Float64} | n/a | yes | Positional argument `dcm`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `dcm_to_quaternion`. Returns `q / norm(q)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/numerics/quaternion_utils.jl`
- [[simulation.setup__initialize_planet_frame_ephemeris_cache_bang|_initialize_planet_frame_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1886-1886`
- [[vehicle.kinematics_rotate_link|rotate_link]] · `callees` → `callers` · call · `src/vehicle/kinematics/kinematics.jl:47-47`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The non-trace branch allocates a mutable `zeros(4)` and the `argmax` call allocates a temporary array, so this is not allocation-free like the rest of the file. The wrapped index arithmetic uses `(i+1)%3` and `(i+2)%3` only when the raw index exceeds three, which yields index 1 rather than 3 for some cases and makes the scalar-part sign in that branch fragile; the trace branch is the well-tested path. The input is assumed to be a proper orthogonal matrix with determinant one, and no check enforces that, so a scaled or reflected matrix produces a meaningless quaternion. The sign of the result is arbitrary within the quaternion double cover.

## Provenance
Mapped from `src/core/numerics/quaternion_utils.jl` line 97.
