---
id: simulation.dynamics_rhs__robot_arm_present
label: _robot_arm_present
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _robot_arm_present
  lines:
  - 1668
  - 1668
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  description: Return value of `_robot_arm_present`.
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

# _robot_arm_present

## Purpose
Returns whether any robot arm exists in the run, caching the answer on shared buffers so the coupling scan is skipped on every RHS call in the common no-arm case.

## Design & Implementation
Returns true conservatively if the parameters lack shared buffers or the cache field, otherwise returns the cached flag, or computes it through `_any_robot_arm_effector` and stores it on first use. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_robot_arm_present`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__apply_coupled_robot_arm_rhs_bang|_apply_coupled_robot_arm_rhs!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1697-1697`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__any_robot_arm_effector|_any_robot_arm_effector]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1673-1673`
<!-- vulcan:connections:end -->

## Limitations
The conservative true for hand-built parameters forces the scan in unit tests; effector tuples are fixed per run so no invalidation is needed.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1668.
