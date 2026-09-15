---
id: core.reference_system__body_fixed_state_xform
label: _body_fixed_state_xform
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: _body_fixed_state_xform
  lines:
  - 47
  - 47
inputs:
- id: from_frame
  type: String
  units: n/a
  required: true
  description: Positional argument `from_frame`.
- id: to_frame
  type: String
  units: n/a
  required: true
  description: Positional argument `to_frame`.
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
  type: SMatrix{6,
  units: n/a
  description: Return value of `_body_fixed_state_xform`.
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

# _body_fixed_state_xform

## Purpose
Fetches the six-by-six state transformation matrix between two SPICE frames at an ephemeris time, under the SPICE lock.

## Design & Implementation
Acquires `_spice_lock()` in a `do` block, calls `sxform(from_frame, to_frame, et)` and wraps the result as an `SMatrix{6,6,Float64}`. Holding the lock is mandatory because CSPICE is not thread-safe and the GRAM extension shares the same statically linked symbols.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `from_frame` | String | n/a | yes | Positional argument `from_frame`. |
| in | `to_frame` | String | n/a | yes | Positional argument `to_frame`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{6, | n/a | — | Return value of `_body_fixed_state_xform`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system__body_fixed_to_j2000_state|_body_fixed_to_j2000_state]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:79-79`
- [[core.reference_system__j2000_to_body_fixed_state|_j2000_to_body_fixed_state]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:62-62`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- `callees` → [[core.reference_system__spice_lock|_spice_lock]] · `callers` · call · `src/core/interfaces/reference_system.jl:48-48`
- `callees` → [[dynamics.perturbations__spice_lock|_spice_lock]] · `callers` · call · `src/core/interfaces/reference_system.jl:48-48`
<!-- vulcan:connections:end -->

## Limitations
Each call performs a kernel lookup and takes a global lock, so calling it per RHS evaluation serialises the whole simulation on that lock; callers should compute once per step.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 47.
