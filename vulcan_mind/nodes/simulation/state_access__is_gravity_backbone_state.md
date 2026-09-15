---
id: simulation.state_access__is_gravity_backbone_state
label: _is_gravity_backbone_state
kind: function
source:
  file: src/simulation/engine/state_access.jl
  symbol: _is_gravity_backbone_state
  lines:
  - 5
  - 5
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
  type: Bool
  units: n/a
  description: Return value of `_is_gravity_backbone_state`.
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

# _is_gravity_backbone_state

## Purpose
Type predicate that decides whether an integrator state `u` is the partitioned second-order layout used by the gravity-backbone solver, as opposed to the ordinary flat ComponentVector state.

## Design & Implementation
Returns true only when `u` is a `RecursiveArrayTools.ArrayPartition` whose `u.x` tuple has exactly two members and both members expose a `:sc` property. The first partition holds velocities, the second positions. Marked `@inline` because every state accessor in this file calls it on the hot right-hand-side path; the check is purely structural, using `hasproperty`, so it costs nothing at runtime for concretely typed ComponentVectors.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_is_gravity_backbone_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.state_access__gravity_backbone_initial_states|_gravity_backbone_initial_states]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:95-95`
- [[simulation.state_access__gravity_backbone_position_state|_gravity_backbone_position_state]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:18-18`
- [[simulation.state_access__gravity_backbone_spacecraft_state|_gravity_backbone_spacecraft_state]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:23-23`
- [[simulation.state_access__gravity_backbone_velocity_state|_gravity_backbone_velocity_state]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:13-13`
- [[simulation.state_access__state_has_heat_loads|_state_has_heat_loads]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:80-80`
- [[simulation.state_access__state_heat_loads|_state_heat_loads]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:73-73`
- [[simulation.state_access__state_mass_kg|_state_mass_kg]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:61-61`
- [[simulation.state_access__state_quaternion|_state_quaternion]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:84-84`
- [[simulation.state_access__state_velocity_ii|_state_velocity_ii]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:47-47`
- [[simx.engine_state_access_state_position_ii|_state_position_ii]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:33-33`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Structural duck-typing means any two-part ArrayPartition carrying `.sc` fields is accepted, even if it is not a backbone state. It cannot distinguish the velocity partition from the position partition, and a three-partition state (for example if an attitude block were added) silently reports false and falls through to the flat-state branches.

## Provenance
Mapped from `src/simulation/engine/state_access.jl` line 5.
