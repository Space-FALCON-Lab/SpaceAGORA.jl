---
id: dynamics.perturbations__srp_sun_position_from_spice_direct_j2000_m
label: _srp_sun_position_from_spice_direct_j2000_m
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _srp_sun_position_from_spice_direct_j2000_m
  lines:
  - 1415
  - 1415
inputs:
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
- id: primary_body_name
  type: String
  units: n/a
  required: true
  description: Positional argument `primary_body_name`.
- id: counter
  type: Base.Threads.Atomic{Int64}
  units: n/a
  required: true
  description: Positional argument `counter`.
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
  description: Return value of `_srp_sun_position_from_spice_direct_j2000_m`.
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

# _srp_sun_position_from_spice_direct_j2000_m

## Purpose
Performs a live SPICE query for the Sun's position relative to the primary body and records it in the SRP runtime call counter.

## Design & Implementation
Increments `counter` atomically and returns `spice_position_j2000_m("sun", et, primary_body_name)`. Declared `@inline`. It is the fallback when the SRP ephemeris cache does not cover the requested time and the memo is disabled.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `primary_body_name` | String | n/a | yes | Positional argument `primary_body_name`. |
| in | `counter` | Base.Threads.Atomic{Int64} | n/a | yes | Positional argument `counter`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_srp_sun_position_from_spice_direct_j2000_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__srp_sun_position_from_spice_j2000_m|_srp_sun_position_from_spice_j2000_m]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1412-1412`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[environment.simple_ephemerides_spice_position_j2000_m|spice_position_j2000_m]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1421-1421`
<!-- vulcan:connections:end -->

## Limitations
Takes the process-wide SPICE lock inside the callee, serialising with every other SPICE user including GRAM.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1415.
