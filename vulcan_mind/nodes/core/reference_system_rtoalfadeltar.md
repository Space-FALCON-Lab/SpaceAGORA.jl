---
id: core.reference_system_rtoalfadeltar
label: rtoalfadeltar
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: rtoalfadeltar
  lines:
  - 235
  - 235
inputs:
- id: r
  type: Any
  units: n/a
  required: true
  description: Positional argument `r`.
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
  description: Return value of `rtoalfadeltar`. Returns `[r, RA, dec]`.
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

# rtoalfadeltar

## Purpose
Converts an inertial Cartesian position to range, right ascension and declination.

## Design & Implementation
Normalises the vector into direction cosines `l, m, n`, takes declination as `asin(n)`, and right ascension as `acos(l / cos(dec))` with the result reflected to `2π - acos` when `m ≤ 0` to resolve the quadrant. Returns a three-element `Vector`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r` | Any | n/a | yes | Positional argument `r`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractArray | n/a | — | Return value of `rtoalfadeltar`. Returns `[r, RA, dec]`. |
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
At the poles `cos(dec)` is zero and the division produces `NaN` or infinity; `atan(m, l)` would be both stabler and simpler. Returns a heap vector rather than an `SVector`.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 235.
