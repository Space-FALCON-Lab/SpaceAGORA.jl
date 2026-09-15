---
id: simulation.dynamics_rhs__gravity_backbone_state_sample
label: _gravity_backbone_state_sample
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _gravity_backbone_state_sample
  lines:
  - 1473
  - 1473
inputs:
- id: q_state
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_state`.
- id: dq_state
  type: Any
  units: n/a
  required: true
  description: Positional argument `dq_state`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  type: StateSample
  units: n/a
  description: Return value of `_gravity_backbone_state_sample`.
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

# _gravity_backbone_state_sample

## Purpose
Builds a `StateSample` from the split solver's separate position and velocity arrays, which do not carry the packed per-satellite layout the ordinary RHS sees.

## Design & Implementation
Extracts position from `q_state` and velocity from `dq_state` through `_gravity_backbone_xyz_chunk`, takes mass as `dry_mass + prop_mass` from the spacecraft model, and attaches the spacecraft. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q_state` | Any | n/a | yes | Positional argument `q_state`. |
| in | `dq_state` | Any | n/a | yes | Positional argument `dq_state`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | StateSample | n/a | — | Return value of `_gravity_backbone_state_sample`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__gravity_backbone_half_kick_bang|_gravity_backbone_half_kick!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1537-1537`
- [[simulation.dynamics_rhs_spacecraft_dynamics_gravity_backbone_bang|spacecraft_dynamics_gravity_backbone!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1568-1568`

**Downstream**

- `callees` → [[core.effector_sampling_statesample|StateSample]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1478-1478`
- `callees` → [[simulation.dynamics_rhs__gravity_backbone_xyz_chunk|_gravity_backbone_xyz_chunk]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1475-1475`
- `callees` → [[simulation.state_access__gravity_backbone_position_state|_gravity_backbone_position_state]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1483-1483`
- `callees` → [[simulation.state_access__gravity_backbone_velocity_state|_gravity_backbone_velocity_state]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1484-1484`
<!-- vulcan:connections:end -->

## Limitations
Mass is the configured total rather than the integrated mass state, so during a burn the backbone sees a slightly wrong mass — irrelevant for gravity, which is mass-independent, but wrong for any mass-dependent kick.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1473.
