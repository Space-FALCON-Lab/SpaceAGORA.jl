---
id: simulation.solver_policy__auto_stiff_smooth_gravity_effector
label: _auto_stiff_smooth_gravity_effector
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _auto_stiff_smooth_gravity_effector
  lines:
  - 113
  - 113
inputs:
- id: effector
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector`.
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
  description: Return value of `_auto_stiff_smooth_gravity_effector`.
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

# _auto_stiff_smooth_gravity_effector

## Purpose
Identifies the gravity effector types whose right-hand side is smooth enough to run on plain `Tsit5` without stiffness auto-switching.

## Design & Implementation
Returns `true` when `effector isa` one of `InverseSquaredGravityModel`, `InverseSquaredJ2GravityModel`, `GravitationalHarmonicsModel`, or `NBodyGravityModel` from `SimulationModel`. Used inside `_auto_stiff_smooth_gravity_reject_reason` for every dynamic effector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_auto_stiff_smooth_gravity_effector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`
- [[simulation.solver_policy__auto_stiff_smooth_gravity_reject_reason|_auto_stiff_smooth_gravity_reject_reason]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:130-130`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The allowlist is closed; new conservative effectors must be added here manually or they force the `AutoTsit5` path. Being smooth does not guarantee non-stiff behaviour near very low periapsis with high-degree harmonics; that regime is still routed to `Tsit5`.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 113.
