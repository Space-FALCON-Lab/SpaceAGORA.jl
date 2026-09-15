---
id: simulation.dynamics_rhs__has_any_batchable_effector
label: _has_any_batchable_effector
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _has_any_batchable_effector
  lines:
  - 667
  - 667
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
  description: Return value of `_has_any_batchable_effector`.
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

# _has_any_batchable_effector

## Purpose
Tests whether any configured effector has a vectorised all-satellites batch kernel, so the flat path can skip its batch phase entirely when none does.

## Design & Implementation
Loops the effector tuple under `@inbounds` and returns true on the first effector for which `_batchable_effector` holds — currently N-body, SRP, inverse-square and J2 gravity. Declared `@inline` with a `::Bool` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effectors` | Tuple | n/a | yes | Positional argument `effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_has_any_batchable_effector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1024-1024`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The batchable set is defined by the predicate it calls; if a new batch kernel is added to `_accumulate_batchable_effector_flat!` without updating that predicate, this test says false and the kernel is never reached.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 667.
