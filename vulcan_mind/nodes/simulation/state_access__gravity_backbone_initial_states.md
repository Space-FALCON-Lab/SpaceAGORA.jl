---
id: simulation.state_access__gravity_backbone_initial_states
label: _gravity_backbone_initial_states
kind: function
source:
  file: src/simulation/engine/state_access.jl
  symbol: _gravity_backbone_initial_states
  lines:
  - 94
  - 94
inputs:
- id: u0
  type: Any
  units: n/a
  required: true
  description: Positional argument `u0`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  description: Return value of `_gravity_backbone_initial_states`. Returns `deepcopy(_gravity_backbone_position_state(u0)),
    deepcopy(_gravity_backbone_veloc` or `q0, dq0`.
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

# _gravity_backbone_initial_states

## Purpose
Builds the `(q0, dq0)` pair of position and velocity ComponentVectors required to start the gravity-backbone second-order integrator from an arbitrary initial state `u0`.

## Design & Implementation
If `u0` is already a backbone `ArrayPartition`, it returns `deepcopy` of the position partition and of the velocity partition, isolating the new run from the source state. Otherwise it constructs fresh shapes by mapping over `args.dynamics_model.spacecraft` to make one `(pos = zeros(3),)` entry per spacecraft for `q0` and one `(vel = zeros(3),)` entry for `dq0`, wraps them as `ComponentVector(sc = ...)`, then copies `u0.sc[i].pos` and `u0.sc[i].vel` elementwise under `@inbounds`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u0` | Any | n/a | yes | Positional argument `u0`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gravity_backbone_initial_states`. Returns `deepcopy(_gravity_backbone_position_state(u0)), deepcopy(_gravity_backbone_veloc` or `q0, dq0`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/state_access.jl`
- [[simulation.execution__build_typed_solver_problem|_build_typed_solver_problem]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:27-27`

**Downstream**

- `callees` → [[simulation.state_access__gravity_backbone_position_state|_gravity_backbone_position_state]] · `callers` · call · `src/simulation/engine/state_access.jl:96-96`
- `callees` → [[simulation.state_access__gravity_backbone_velocity_state|_gravity_backbone_velocity_state]] · `callers` · call · `src/simulation/engine/state_access.jl:96-96`
- `callees` → [[simulation.state_access__is_gravity_backbone_state|_is_gravity_backbone_state]] · `callers` · call · `src/simulation/engine/state_access.jl:95-95`
<!-- vulcan:connections:end -->

## Limitations
Only three-element position and velocity are carried across, so mass, heat loads and attitude present in `u0` are dropped without warning, which is why backbone runs cannot resolve those channels afterwards. The `@inbounds` copy trusts that `args.dynamics_model.spacecraft` and `u0.sc` have identical length and ordering; a mismatch reads out of bounds rather than raising. Units are inherited from `u0` and never checked.

## Provenance
Mapped from `src/simulation/engine/state_access.jl` line 94.
