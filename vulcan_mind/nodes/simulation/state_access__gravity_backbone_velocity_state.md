---
id: simulation.state_access__gravity_backbone_velocity_state
label: _gravity_backbone_velocity_state
kind: function
source:
  file: src/simulation/engine/state_access.jl
  symbol: _gravity_backbone_velocity_state
  lines:
  - 12
  - 12
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
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
  description: Return value of `_gravity_backbone_velocity_state`. Returns `u.x[1]`.
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

# _gravity_backbone_velocity_state

## Purpose
Extracts the velocity half of a partitioned gravity-backbone integrator state, giving callers the ComponentVector that holds `sc[i].vel` for every spacecraft.

## Design & Implementation
Guards with `_is_gravity_backbone_state(u)` and throws `ArgumentError("Expected gravity-backbone second-order state.")` when the state is not the two-partition layout. On success it returns `u.x[1]`, the first member of the `ArrayPartition`, by convention the first derivative block of the second-order problem. No copy is made: the returned object aliases the live integrator state, so writes through it mutate the solver's working vector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gravity_backbone_velocity_state`. Returns `u.x[1]`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/state_access.jl`
- [[simulation.dynamics_rhs__gravity_backbone_half_kick_bang|_gravity_backbone_half_kick!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1527-1527`
- [[simulation.dynamics_rhs__gravity_backbone_state_sample|_gravity_backbone_state_sample]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1484-1484`
- [[simulation.state_access__gravity_backbone_initial_states|_gravity_backbone_initial_states]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:96-96`
- [[simulation.state_access__state_velocity_ii|_state_velocity_ii]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:48-48`

**Downstream**

- `callees` → [[simulation.state_access__is_gravity_backbone_state|_is_gravity_backbone_state]] · `callers` · call · `src/simulation/engine/state_access.jl:13-13`
<!-- vulcan:connections:end -->

## Limitations
The ordering convention (velocity first, position second) is an unchecked invariant shared with `_gravity_backbone_position_state`; swapping the partitions at construction time would silently produce positions here. Because the result aliases solver memory, retaining it across a step yields values that change underneath the holder.

## Provenance
Mapped from `src/simulation/engine/state_access.jl` line 12.
