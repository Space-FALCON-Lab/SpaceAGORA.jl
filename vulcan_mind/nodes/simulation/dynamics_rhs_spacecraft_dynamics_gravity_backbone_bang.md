---
id: simulation.dynamics_rhs_spacecraft_dynamics_gravity_backbone_bang
label: spacecraft_dynamics_gravity_backbone!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: spacecraft_dynamics_gravity_backbone!
  lines:
  - 1554
  - 1554
inputs:
- id: ddu
  type: Any
  units: n/a
  required: true
  description: Positional argument `ddu`.
- id: dq
  type: Any
  units: n/a
  required: true
  description: Positional argument `dq`.
- id: q
  type: Any
  units: n/a
  required: true
  description: Positional argument `q`.
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
  description: Return value of `spacecraft_dynamics_gravity_backbone!`; mutates `ddu`
    in place. Returns `nothing`.
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

# spacecraft_dynamics_gravity_backbone!

## Purpose
The second-order backbone RHS for the split solver: given position and velocity states, writes the acceleration from position-only static gravity into the second-derivative output.

## Design & Implementation
Sets the current time, then loops active satellites — batched with Polyester using a `minbatch` derived from core count when RHS batching is enabled, serially otherwise — building a state sample from the separate `q` and `dq` arrays and writing `_gravity_backbone_core_acceleration` into `ddu.sc[i].vel`. Inactive satellites receive zero acceleration.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ddu` | Any | n/a | yes | Positional argument `ddu`. |
| in | `dq` | Any | n/a | yes | Positional argument `dq`. |
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `spacecraft_dynamics_gravity_backbone!`; mutates `ddu` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__gravity_backbone_core_acceleration|_gravity_backbone_core_acceleration]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1569-1569`
- `callees` → [[simulation.dynamics_rhs__gravity_backbone_state_sample|_gravity_backbone_state_sample]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1568-1568`
- `callees` → [[simulation.setup__rhs_batch_parallel_enabled|_rhs_batch_parallel_enabled]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1560-1560`
<!-- vulcan:connections:end -->

## Limitations
Only static gravity contributes; every other force reaches the state through the explicit remainder and kick steps, so the backbone alone is not a complete dynamics model.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1554.
