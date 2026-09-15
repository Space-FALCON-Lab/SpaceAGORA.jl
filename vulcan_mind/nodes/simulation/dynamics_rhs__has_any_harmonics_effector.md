---
id: simulation.dynamics_rhs__has_any_harmonics_effector
label: _has_any_harmonics_effector
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _has_any_harmonics_effector
  lines:
  - 687
  - 687
inputs:
- id: effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `effectors`.
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
  description: Return value of `_has_any_harmonics_effector`.
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

# _has_any_harmonics_effector

## Purpose
Tests whether any configured effector is a spherical-harmonics prepass effector, gating whether the flat path runs the harmonics batch kernel and whether the planet-frame rotation must be shared.

## Design & Implementation
Loops the effector tuple under `@inbounds` and returns true on the first effector for which `_harmonics_prepass_effector` holds. Declared `@inline` with a `::Bool` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effectors` | Tuple | n/a | yes | Positional argument `effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_has_any_harmonics_effector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1056-1056`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
None beyond the predicate it depends on; a harmonics model wrapped in another type would not be recognised.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 687.
