---
id: simulation.runtime__stage_environment_state
label: _stage_environment_state
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/runtime.jl
  symbol: _stage_environment_state
  lines:
  - 204
  - 204
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
- id: write_buffers
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `write_buffers` (default `true`).
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
  description: Return value of `_stage_environment_state`. Returns `merge(kin, (rho=atmosphere.rho_kg_m3,
    T=atmosphere.temperature_k, wind=atmospher`.
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

# _stage_environment_state

## Purpose
The general form of the environment tuple, letting the caller choose between reading buffered atmosphere and sampling fresh without disturbing the buffers.

## Design & Implementation
Computes kinematics as the buffered variant does, then branches on `write_buffers`: true calls `sample_buffered_atmosphere`, false calls `sample_atmosphere` with `write_buffers=false` so a diagnostic or predictive query does not overwrite the state the RHS is integrating against. The result is merged into the same `rho`, `T`, `wind` keys either way.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `write_buffers` | Bool | n/a | no | Keyword argument `write_buffers` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_stage_environment_state`. Returns `merge(kin, (rho=atmosphere.rho_kg_m3, T=atmosphere.temperature_k, wind=atmospher`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:706-706`
- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/runtime.jl`

**Downstream**

- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:207-207`
- `callees` → [[simulation.runtime__stage_environment_kinematics|_stage_environment_kinematics]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:205-205`
<!-- vulcan:connections:end -->

## Limitations
The non-writing path evaluates the density model directly, so with a native GRAM model it takes the GRAM lock and pays a full sample per call; a caller looping over many satellites with `write_buffers=false` can be far slower than the callback.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/runtime.jl` line 204.
