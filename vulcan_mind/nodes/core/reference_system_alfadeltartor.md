---
id: core.reference_system_alfadeltartor
label: alfadeltartor
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: alfadeltartor
  lines:
  - 254
  - 254
inputs:
- id: R_RA_DEC
  type: Any
  units: n/a
  required: true
  description: Positional argument `R_RA_DEC`.
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
  type: AbstractArray
  units: n/a
  description: Return value of `alfadeltartor`. Returns `[x, y, z]`.
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

# alfadeltartor

## Purpose
Converts range, right ascension and declination back to an inertial Cartesian position.

## Design & Implementation
Reads `R, RA, DEC` from the input and returns `[R cos DEC cos RA, R cos DEC sin RA, R sin DEC]` as a three-element `Vector`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `R_RA_DEC` | Any | n/a | yes | Positional argument `R_RA_DEC`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractArray | n/a | — | Return value of `alfadeltartor`. Returns `[x, y, z]`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Returns a heap-allocated vector; no validation that declination lies within plus or minus pi over two.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 254.
