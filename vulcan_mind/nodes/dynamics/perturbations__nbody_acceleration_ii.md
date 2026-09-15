---
id: dynamics.perturbations__nbody_acceleration_ii
label: _nbody_acceleration_ii
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _nbody_acceleration_ii
  lines:
  - 976
  - 976
inputs:
- id: model
  type: NBodyGravityModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: x
  type: StateSample
  units: n/a
  required: true
  description: Positional argument `x`.
- id: third_bodies
  type: ThirdBodyEphemerisSample
  units: n/a
  required: true
  description: Positional argument `third_bodies`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_nbody_acceleration_ii`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _nbody_acceleration_ii

## Purpose
Sums the third-body perturbing accelerations on a spacecraft from precomputed body positions, in the wrench-based effector path.

## Theory & Math
$$
\vec{a} = \sum_k \mu_k \left( \frac{\vec{r}_k - \vec{r}}{|\vec{r}_k - \vec{r}|^3} - \frac{\vec{r}_k}{|\vec{r}_k|^3} \right)
$$

## Design & Implementation
For each body, forms the spacecraft-to-body vector and adds `μ_k` times the difference between the direct term and the indirect term on the primary. Accumulates in an `MVector` and returns an `SVector`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | NBodyGravityModel | n/a | yes | Positional argument `model`. |
| in | `x` | StateSample | n/a | yes | Positional argument `x`. |
| in | `third_bodies` | ThirdBodyEphemerisSample | n/a | yes | Positional argument `third_bodies`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_nbody_acceleration_ii`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_gravity_backbone_kick_acceleration_ii|gravity_backbone_kick_acceleration_ii]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1013-1013`
- [[dynamics.perturbations_wrench|wrench]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1001-1001`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The direct and indirect terms nearly cancel for distant bodies, losing precision; no Encke-style reformulation is used.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 976.
