---
id: simulation.setup__rhs_flat_has_batch_privileged_effector
label: _rhs_flat_has_batch_privileged_effector
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_has_batch_privileged_effector
  lines:
  - 739
  - 739
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
  description: Return value of `_rhs_flat_has_batch_privileged_effector`.
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

# _rhs_flat_has_batch_privileged_effector

## Purpose
Folds `_rhs_flat_batch_privileged_effector` over an effector tuple to report whether the flat queue should be eligible regardless of the `flat_min_effectors` count.

## Design & Implementation
`_rhs_flat_has_batch_privileged_effector(dynamic_effectors::Tuple)::Bool` loops with `@inbounds` and returns `true` on the first privileged effector, `false` otherwise. Used in `_rhs_execution_plan_uncached` as an alternative to `n_effectors >= env.flat_min_effectors` when deciding flat-queue eligibility.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_rhs_flat_has_batch_privileged_effector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1242-1242`

**Downstream**

- `callees` → [[simulation.setup__rhs_flat_batch_privileged_effector|_rhs_flat_batch_privileged_effector]] · `callers` · call · `src/simulation/engine/setup.jl:741-741`
<!-- vulcan:connections:end -->

## Limitations
Presence of a privileged effector does not by itself guarantee the flat queue is faster; the satellite-count and thread-budget gates still apply, and an inverse-square model on a tiny constellation would still be routed to serial by them. The loop is over a tuple so heterogeneous element types compile to a static unrolled check.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 739.
