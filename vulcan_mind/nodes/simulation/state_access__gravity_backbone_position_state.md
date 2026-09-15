---
id: simulation.state_access__gravity_backbone_position_state
label: _gravity_backbone_position_state
kind: function
source:
  file: src/simulation/engine/state_access.jl
  symbol: _gravity_backbone_position_state
  lines:
  - 17
  - 17
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
  description: Return value of `_gravity_backbone_position_state`. Returns `u.x[2]`.
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

# _gravity_backbone_position_state

## Purpose
Extracts the position half of a partitioned gravity-backbone integrator state, the ComponentVector carrying `sc[i].pos` for each modelled spacecraft.

## Design & Implementation
Mirrors the velocity accessor: it asserts `_is_gravity_backbone_state(u)` and otherwise throws `ArgumentError("Expected gravity-backbone second-order state.")`, then returns `u.x[2]`, the second member of the `ArrayPartition`. The function is `@inline` and allocation-free, returning a view onto the integrator's own storage rather than a copy, which keeps it usable inside saving callbacks and derivative evaluations.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gravity_backbone_position_state`. Returns `u.x[2]`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/state_access.jl`
- [[simulation.dynamics_rhs__gravity_backbone_state_sample|_gravity_backbone_state_sample]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1483-1483`
- [[simulation.state_access__gravity_backbone_initial_states|_gravity_backbone_initial_states]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:96-96`
- [[simulation.state_access__gravity_backbone_spacecraft_state|_gravity_backbone_spacecraft_state]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:24-24`

**Downstream**

- `callees` → [[simulation.state_access__is_gravity_backbone_state|_is_gravity_backbone_state]] · `callers` · call · `src/simulation/engine/state_access.jl:18-18`
<!-- vulcan:connections:end -->

## Limitations
The partition index is hard-coded, so the layout contract with the second-order problem construction is implicit rather than verified. Calling it on a flat ComponentVector state raises rather than degrading gracefully, so callers that accept both layouts must branch on `_is_gravity_backbone_state` first. The alias is mutable and unprotected.

## Provenance
Mapped from `src/simulation/engine/state_access.jl` line 17.
