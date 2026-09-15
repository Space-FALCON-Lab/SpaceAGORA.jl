---
id: simulation.setup__rhs_env_config_from_buffers
label: _rhs_env_config_from_buffers
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_env_config_from_buffers
  lines:
  - 889
  - 889
inputs:
- id: shared_buffers
  type: Any
  units: n/a
  required: true
  description: Positional argument `shared_buffers`.
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
  type: SimulationModel.RhsPlanEnvConfig
  units: n/a
  description: Return value of `_rhs_env_config_from_buffers`.
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

# _rhs_env_config_from_buffers

## Purpose
Returns the run-scoped RHS plan snapshot from shared buffers, or builds a fresh one from the environment when none was captured.

## Design & Implementation
If buffers exist and their `rhs_env_config[]` is non-`nothing`, returns it; otherwise calls `_snapshot_rhs_plan_env_config()`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.RhsPlanEnvConfig | n/a | — | Return value of `_rhs_env_config_from_buffers`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.dynamics_rhs__prepare_rhs_flat_work_packets_bang|_prepare_rhs_flat_work_packets!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:518-518`
- [[simulation.dynamics_rhs__rhs_flat_use_packet_scheduler|_rhs_flat_use_packet_scheduler]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:571-571`
- [[simulation.dynamics_rhs__update_rhs_flat_packet_overhead_model_bang|_update_rhs_flat_packet_overhead_model!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:599-599`
- [[simulation.setup__rhs_effector_observed_cost_ns|_rhs_effector_observed_cost_ns]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:594-594`
- [[simulation.setup__update_effector_cost_model_bang|_update_effector_cost_model!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:546-546`
- [[simulation.setup__update_rhs_effector_cost_model_bang|_update_rhs_effector_cost_model!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:614-614`

**Downstream**

- `callees` → [[simulation.setup__effector_shared_buffers|_effector_shared_buffers]] · `callers` · call · `src/simulation/engine/setup.jl:898-898`
- `callees` → [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callers` · call · `src/simulation/engine/setup.jl:894-894`
<!-- vulcan:connections:end -->

## Limitations
The fallback parses over thirty environment variables, so a hand-built parameter set that never captured a snapshot pays that on every RHS call.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 889.
