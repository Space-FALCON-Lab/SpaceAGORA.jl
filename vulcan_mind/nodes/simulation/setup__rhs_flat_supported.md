---
id: simulation.setup__rhs_flat_supported
label: _rhs_flat_supported
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_supported
  lines:
  - 919
  - 919
inputs:
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
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
  description: Return value of `_rhs_flat_supported`.
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

# _rhs_flat_supported

## Purpose
Decides whether the effector set can run through the flat constellation queue at all.

## Design & Implementation
For a single effector, true if it is a thread-safe harmonics model with batching enabled or an inverse-square model meeting `_rhs_single_invsq_flat_supported`. For several effectors, true if `_dynamic_effectors_parallel_supported` holds. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_rhs_flat_supported`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.rhs_calibration__rhs_plan_candidates|_rhs_plan_candidates]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:229-229`
- [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1104-1104`

**Downstream**

- `callees` → [[simulation.setup__dynamic_effector_threadsafe|_dynamic_effector_threadsafe]] · `callers` · call · `src/simulation/engine/setup.jl:922-922`
- `callees` → [[simulation.setup__dynamic_effectors_parallel_supported|_dynamic_effectors_parallel_supported]] · `callers` · call · `src/simulation/engine/setup.jl:926-926`
- `callees` → [[simulation.setup__rhs_harmonics_batch_enabled|_rhs_harmonics_batch_enabled]] · `callers` · call · `src/simulation/engine/setup.jl:923-923`
- `callees` → [[simulation.setup__rhs_single_invsq_flat_supported|_rhs_single_invsq_flat_supported]] · `callers` · call · `src/simulation/engine/setup.jl:924-924`
<!-- vulcan:connections:end -->

## Limitations
This is the older environment-reading form; the planner calls a two-argument variant taking the env snapshot, and the two must stay in agreement.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 919.
