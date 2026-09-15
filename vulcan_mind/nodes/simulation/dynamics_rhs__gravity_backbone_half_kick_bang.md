---
id: simulation.dynamics_rhs__gravity_backbone_half_kick_bang
label: _gravity_backbone_half_kick!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _gravity_backbone_half_kick!
  lines:
  - 1524
  - 1524
inputs:
- id: u_state
  type: Any
  units: n/a
  required: true
  description: Positional argument `u_state`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: half_dt
  type: Float64
  units: n/a
  required: true
  description: Positional argument `half_dt`.
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
  type: Nothing
  units: n/a
  description: Return value of `_gravity_backbone_half_kick!`; mutates `u_state` in
    place. Returns `nothing`.
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

# _gravity_backbone_half_kick!

## Purpose
Applies the explicit velocity half-kick from kick-structured effectors — third-body gravity — between backbone steps of the split solver.

## Design & Implementation
Loops active satellites, batched with Polyester when RHS batching is enabled or serially otherwise, building a state sample, computing `_gravity_backbone_kick_acceleration` and adding `half_dt` times it to the satellite's velocity in `vel_state`. Inactive satellites are skipped.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u_state` | Any | n/a | yes | Positional argument `u_state`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `half_dt` | Float64 | n/a | yes | Positional argument `half_dt`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_gravity_backbone_half_kick!`; mutates `u_state` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.solver_policy__solve_with_gravity_backbone_solver|_solve_with_gravity_backbone_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:581-581`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__gravity_backbone_kick_acceleration|_gravity_backbone_kick_acceleration]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1538-1538`
- `callees` → [[simulation.dynamics_rhs__gravity_backbone_state_sample|_gravity_backbone_state_sample]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1537-1537`
- `callees` → [[simulation.setup__rhs_batch_parallel_enabled|_rhs_batch_parallel_enabled]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1529-1529`
- `callees` → [[simulation.state_access__gravity_backbone_velocity_state|_gravity_backbone_velocity_state]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1527-1527`
<!-- vulcan:connections:end -->

## Limitations
Mutates the velocity state in place, so it must be called exactly once per half step; calling it twice doubles the kick with no guard.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1524.
