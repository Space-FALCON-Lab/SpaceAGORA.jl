---
id: simulation.runtime__buffered_stage_environment_state
label: _buffered_stage_environment_state
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/runtime.jl
  symbol: _buffered_stage_environment_state
  lines:
  - 72
  - 72
inputs:
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
  description: Return value of `_buffered_stage_environment_state`. Returns `merge(kin,
    (rho=atmosphere.rho_kg_m3, T=atmosphere.temperature_k, wind=atmospher`.
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

# _buffered_stage_environment_state

## Purpose
Produces the full environment tuple for one satellite — kinematics plus density, temperature and wind — reading the atmosphere from the shared buffers rather than re-evaluating it.

## Design & Implementation
Merges the kinematics tuple with the result of the engine's `sample_buffered_atmosphere`, renaming `rho_kg_m3`, `temperature_k` and `wind_pp` to the shorter `rho`, `T` and `wind` keys that guidance and control code expects. `@inline` so the merge of two small named tuples is resolved at compile time.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_buffered_stage_environment_state`. Returns `merge(kin, (rho=atmosphere.rho_kg_m3, T=atmosphere.temperature_k, wind=atmospher`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/runtime.jl`

**Downstream**

- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:74-74`
- `callees` → [[simulation.runtime__stage_environment_kinematics|_stage_environment_kinematics]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:73-73`
<!-- vulcan:connections:end -->

## Limitations
It reads whatever the buffers currently hold, so if the density callback has not yet run at this `t` the returned atmosphere is from the previous sample without any staleness marker.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/runtime.jl` line 72.
