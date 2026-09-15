---
id: simulation.dynamics_rhs__effector_value
label: _effector_value
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _effector_value
  lines:
  - 1632
  - 1632
inputs:
- id: effector
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector`.
- id: name
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: Any
  units: n/a
  required: true
  description: Positional argument `default`.
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
  description: Return value of `_effector_value`.
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

# _effector_value

## Purpose
Reads an optional property from a robot-arm effector with a default, letting several effector types share the coupling extractor without a common interface.

## Design & Implementation
Returns `getproperty(effector, name)` if `hasproperty` holds, else `default`. Declared `@inline`. Used for the four stiffness and damping parameters in `_robot_arm_coupling_from_effector`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `name` | Symbol | n/a | yes | Positional argument `name`. |
| in | `default` | Any | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_effector_value`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__robot_arm_coupling_from_effector|_robot_arm_coupling_from_effector]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1640-1640`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No type check on the returned value, so a property present with the wrong type surfaces later inside the multibody dynamics.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1632.
