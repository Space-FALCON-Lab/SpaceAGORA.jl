---
id: core.reference_system__j2000_to_body_fixed_state
label: _j2000_to_body_fixed_state
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: _j2000_to_body_fixed_state
  lines:
  - 53
  - 53
inputs:
- id: r_i
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `r_i`.
- id: v_i
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `v_i`.
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
  description: Return value of `_j2000_to_body_fixed_state`.
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

# _j2000_to_body_fixed_state

## Purpose
Transforms a J2000 position-velocity pair into the body-fixed frame with SPICE, using the high-precision Earth frame when its kernel is available.

## Design & Implementation
Packs the inputs into a six-vector. For Earth it tries the `ITRF93` transform and, on any exception, retries with `IAU_EARTH`; for other bodies it uses `_spice_body_fixed_frame`. Returns the transformed six-vector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_i` | SVector{3, Float64} | n/a | yes | Positional argument `r_i`. |
| in | `v_i` | SVector{3, Float64} | n/a | yes | Positional argument `v_i`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{6, | n/a | — | Return value of `_j2000_to_body_fixed_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system__planet_flattening|_planet_flattening]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:92-92`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- `callees` → [[core.reference_system__body_fixed_state_xform|_body_fixed_state_xform]] · `callers` · call · `src/core/interfaces/reference_system.jl:62-62`
- `callees` → [[core.reference_system__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callers` · call · `src/core/interfaces/reference_system.jl:67-67`
- `callees` → [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callers` · call · `src/core/interfaces/reference_system.jl:67-67`
- `callees` → [[environment.simple_ephemerides__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callers` · call · `src/core/interfaces/reference_system.jl:67-67`
<!-- vulcan:connections:end -->

## Limitations
The Earth fallback catches every exception, so a genuine kernel-loading error unrelated to `ITRF93` is masked as a silent precision downgrade; the try is also paid on every call when the high-precision frame is missing.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 53.
