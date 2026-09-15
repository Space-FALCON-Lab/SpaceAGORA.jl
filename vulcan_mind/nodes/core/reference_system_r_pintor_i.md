---
id: core.reference_system_r_pintor_i
label: r_pintor_i
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: r_pintor_i
  lines:
  - 36
  - 36
inputs:
- id: r_p
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `r_p`.
- id: v_p
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `v_p`.
- id: planet
  type: T
  units: n/a
  required: true
  description: Positional argument `planet`.
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
  type: Tuple{SVector{3,
  units: n/a
  description: 'Return value of `r_pintor_i`. Type parameters: `T`.'
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

# r_pintor_i

## Purpose
Rotates a planet-fixed position and velocity back into the J2000 inertial frame, the inverse of `r_intor_p!`.

## Theory & Math
$$
r_i = L_{PI}^\top r_p,\qquad v_i = L_{PI}^\top \left( v_p + \omega \times r_p \right)
$$

## Design & Implementation
Two methods. The legacy form applies the transpose of `L_PI` after adding the transport term `ω × r_p` to the velocity. The `et` form delegates to `_body_fixed_to_j2000_state` for a SPICE transform and splits the six-vector result into position and velocity.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_p` | SVector{3, Float64} | n/a | yes | Positional argument `r_p`. |
| in | `v_p` | SVector{3, Float64} | n/a | yes | Positional argument `v_p`. |
| in | `planet` | T | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{SVector{3, | n/a | — | Return value of `r_pintor_i`. Type parameters: `T`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system__planet_flattening|_planet_flattening]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:116-116`
- [[core.reference_system_latlongtooe|latlongtoOE]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:324-324`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
There is no `ephemerides_model` overload matching the inertial-to-fixed direction, so a caller using the analytic ephemerides has to reconstruct the inverse rotation itself.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 36.
