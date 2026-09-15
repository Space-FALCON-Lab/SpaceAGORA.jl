---
id: dynamics.perturbations__nbody_body_force_ii
label: _nbody_body_force_ii
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _nbody_body_force_ii
  lines:
  - 117
  - 117
inputs:
- id: pos_primary_k
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_primary_k`.
- id: pos_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_ii`.
- id: mass
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass`.
- id: mu_k
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mu_k`.
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
  description: Return value of `_nbody_body_force_ii`.
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

# _nbody_body_force_ii

## Purpose
The per-body force contribution in the legacy N-body `calcForceTorque` path, identical in form to the acceleration kernel used by the wrench path but multiplied by spacecraft mass.

## Theory & Math
$$
\vec{F}_k = m\,\mu_k \left( \frac{\vec{r}_k - \vec{r}}{|\vec{r}_k - \vec{r}|^3} - \frac{\vec{r}_k}{|\vec{r}_k|^3} \right)
$$

## Design & Implementation
Forms the spacecraft-to-body vector and its norm, the primary-to-body norm, and returns `mass μ_k` times the difference between the direct term `(r_k - r) / |r_k - r|³` and the indirect term `r_k / |r_k|³`. Declared `@inline` so the threaded per-body loop has no call overhead.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos_primary_k` | SVector{3, Float64} | n/a | yes | Positional argument `pos_primary_k`. |
| in | `pos_ii` | SVector{3, Float64} | n/a | yes | Positional argument `pos_ii`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `mu_k` | Float64 | n/a | yes | Positional argument `mu_k`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_nbody_body_force_ii`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:951-951`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
For a distant perturber such as Jupiter the two terms are nearly equal and their difference loses several significant digits; no Encke or Battin reformulation is applied.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 117.
