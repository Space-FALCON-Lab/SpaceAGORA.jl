---
id: simulation.effector_sampling__extract_sample_pos_vel
label: _extract_sample_pos_vel
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: _extract_sample_pos_vel
  lines:
  - 4
  - 4
inputs:
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
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
  type: SVector
  units: n/a
  description: Return value of `_extract_sample_pos_vel`. Returns `x.pos_ii, x.vel_ii`
    or `SVector{3, Float64}(x[1], x[2], x[3]), SVector{3, Float64}(x[4], x[5], x[6])`.
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

# _extract_sample_pos_vel

## Purpose
Pulls inertial position and velocity out of whichever state representation the effector layer was handed — a component-tree view, a labelled vector, or a raw vector.

## Design & Implementation
Three cases in priority order: a view with `pos_ii` and `vel_ii` properties returns them directly; one with `pos` and `vel` properties, and the bare fallback, both read elements one to six into two `SVector{3,Float64}`. The middle and last branches are identical in effect, differing only in the property test that documents the intent. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `_extract_sample_pos_vel`. Returns `x.pos_ii, x.vel_ii` or `SVector{3, Float64}(x[1], x[2], x[3]), SVector{3, Float64}(x[4], x[5], x[6])`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.dynamics_rhs__prefill_rhs_flat_state_samples_bang|_prefill_rhs_flat_state_samples!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:302-302`
- [[simulation.effector_sampling__sample_atmosphere_from_planet_frame|_sample_atmosphere_from_planet_frame]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:91-91`
- [[simulation.effector_sampling_build_state_sample|build_state_sample]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:25-25`
- [[simulation.effector_sampling_sample_planet_frame|sample_planet_frame]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:43-43`
- [[simulation.effector_sampling_sample_planet_frame_with_lpi|sample_planet_frame_with_lpi]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:54-54`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The labelled-property branch does not actually use the properties it detects, so a state type whose `pos` is not stored in slots one to three would be misread without error.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 4.
