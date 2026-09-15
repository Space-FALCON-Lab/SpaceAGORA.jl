---
id: core.reference_system__body_fixed_to_j2000_state
label: _body_fixed_to_j2000_state
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: _body_fixed_to_j2000_state
  lines:
  - 70
  - 70
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
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
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
  type: SVector{6,
  units: n/a
  description: Return value of `_body_fixed_to_j2000_state`.
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

# _body_fixed_to_j2000_state

## Purpose
Transforms a body-fixed position-velocity pair into J2000 with SPICE, the inverse of `_j2000_to_body_fixed_state` with the same Earth fallback.

## Design & Implementation
Packs the inputs, tries `ITRF93 → J2000` for Earth with an `IAU_EARTH` fallback, and uses the planet's body-fixed frame for other bodies. Returns the transformed six-vector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_p` | SVector{3, Float64} | n/a | yes | Positional argument `r_p`. |
| in | `v_p` | SVector{3, Float64} | n/a | yes | Positional argument `v_p`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{6, | n/a | — | Return value of `_body_fixed_to_j2000_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system__planet_flattening|_planet_flattening]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:117-117`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- `callees` → [[core.reference_system__body_fixed_state_xform|_body_fixed_state_xform]] · `callers` · call · `src/core/interfaces/reference_system.jl:79-79`
- `callees` → [[core.reference_system__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callers` · call · `src/core/interfaces/reference_system.jl:84-84`
- `callees` → [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callers` · call · `src/core/interfaces/reference_system.jl:84-84`
- `callees` → [[environment.simple_ephemerides__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callers` · call · `src/core/interfaces/reference_system.jl:84-84`
<!-- vulcan:connections:end -->

## Limitations
Identical fallback masking as the forward transform; a mixed run where the forward transform succeeded with `ITRF93` and the inverse fell back to `IAU_EARTH` would introduce a small frame inconsistency with no diagnostic.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 70.
