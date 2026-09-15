---
id: simulation.setup__profile_forces_serial_rhs
label: _profile_forces_serial_rhs
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _profile_forces_serial_rhs
  lines:
  - 375
  - 375
inputs:
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
  type: Bool
  units: n/a
  description: Return value of `_profile_forces_serial_rhs`.
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

# _profile_forces_serial_rhs

## Purpose
Detects the benchmarking profiles that demand a strictly serial RHS so parallel-versus-serial comparisons in performance studies are not contaminated by automatic threading.

## Design & Implementation
Reads `SPACEAGORA_PARALLEL_PROFILE` through `_engine_env_get` with default `""`, lowercases and strips it, and returns `true` when it is one of `"r0"`, `"serial"`, `"r0_true_serial"`, or `"true_serial"`. The result becomes `RhsPlanEnvConfig.profile_forces_serial`, which short-circuits `_rhs_batch_parallel_enabled` to `false` regardless of mode.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_profile_forces_serial_rhs`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:852-852`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/setup.jl:376-376`
<!-- vulcan:connections:end -->

## Limitations
The list of profile names is hard-coded and must be kept in sync with the benchmarking scripts that set the variable. A profile value with a suffix (`"r0_v2"`) does not match. The same variable also selects the persistent-hint file name in `ParallelPolicy`, coupling two unrelated concerns.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 375.
