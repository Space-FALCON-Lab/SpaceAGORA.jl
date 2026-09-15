---
id: simulation.dynamics_rhs__prefill_rhs_flat_state_samples_bang
label: _prefill_rhs_flat_state_samples!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _prefill_rhs_flat_state_samples!
  lines:
  - 292
  - 292
inputs:
- id: shared_buffers
  type: Any
  units: n/a
  required: true
  description: Positional argument `shared_buffers`.
- id: sc_state
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_state`.
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
  type: Nothing
  units: n/a
  description: Return value of `_prefill_rhs_flat_state_samples!`; mutates `shared_buffers`
    in place.
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

# _prefill_rhs_flat_state_samples!

## Purpose
Copies each active satellite's position, velocity, mass and — when attitude is simulated — quaternion and body rate into the flat state buffers, so worker threads read plain vectors instead of component-tree views.

## Design & Implementation
Reads the five buffer vectors from shared buffers, then loops satellites under `@inbounds`, skipping inactive ones, extracting position and velocity through `_extract_sample_pos_vel` and mass through `_extract_sample_mass_kg`, and copying attitude only when `orientation_sim`. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `sc_state` | Any | n/a | yes | Positional argument `sc_state`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_prefill_rhs_flat_state_samples!`; mutates `shared_buffers` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1029-1029`

**Downstream**

- `callees` → [[simulation.effector_sampling__extract_sample_mass_kg|_extract_sample_mass_kg]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:305-305`
- `callees` → [[simulation.effector_sampling__extract_sample_pos_vel|_extract_sample_pos_vel]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:302-302`
<!-- vulcan:connections:end -->

## Limitations
Buffers must already be sized by `_ensure_rhs_flat_effector_scratch!`; inactive satellites' entries are left stale rather than zeroed.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 292.
