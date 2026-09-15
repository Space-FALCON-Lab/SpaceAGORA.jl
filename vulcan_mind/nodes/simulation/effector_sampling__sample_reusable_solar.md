---
id: simulation.effector_sampling__sample_reusable_solar
label: _sample_reusable_solar
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: _sample_reusable_solar
  lines:
  - 275
  - 275
inputs:
- id: req
  type: EffectorEnvironmentRequirements
  units: n/a
  required: true
  description: Positional argument `req`.
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  description: Return value of `_sample_reusable_solar`. Returns `(p.shared_buffers.rhs_solar_prefilled[]
    && p.shared_buffers.rhs_flat_solar_t[] =`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _sample_reusable_solar

## Purpose
Chooses between the prefilled solar position and a fresh ephemeris lookup for an effector needing the Sun.

## Design & Implementation
Returns `nothing` when solar data is not required. It reuses `rhs_flat_solar_pos_ii[]` only when both `rhs_solar_prefilled[]` is set and the recorded `rhs_flat_solar_t[]` equals `t` exactly; otherwise it calls `sample_solar_ephemeris`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `req` | EffectorEnvironmentRequirements | n/a | yes | Positional argument `req`. |
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_sample_reusable_solar`. Returns `(p.shared_buffers.rhs_solar_prefilled[] && p.shared_buffers.rhs_flat_solar_t[] =`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.effector_sampling_sample_environment_with_reusable_buffers|sample_environment_with_reusable_buffers]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:296-296`

**Downstream**

- `callees` → [[core.effector_sampling_solarephemerissample|SolarEphemerisSample]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:278-278`
- `callees` → [[simulation.effector_sampling_sample_solar_ephemeris|sample_solar_ephemeris]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:279-279`
<!-- vulcan:connections:end -->

## Limitations
Unlike the atmosphere path, there is no freeze-per-step relaxation here, so the exact-time test means intermediate RHS stages at times other than the prefill time always pay a fresh lookup.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 275.
