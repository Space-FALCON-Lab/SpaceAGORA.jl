---
id: simulation.solver_policy__is_partitioned_second_order_problem
label: _is_partitioned_second_order_problem
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _is_partitioned_second_order_problem
  lines:
  - 397
  - 397
inputs:
- id: prob
  type: Any
  units: n/a
  required: true
  description: Positional argument `prob`.
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
  description: Return value of `_is_partitioned_second_order_problem`.
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

# _is_partitioned_second_order_problem

## Purpose
Tests whether `prob` is a `SecondOrderODEProblem` (position/velocity partitioned), which the symplectic and gravity-backbone modes require.

## Design & Implementation
Returns `false` if `prob` lacks a `problem_type` property; otherwise returns `getproperty(prob, :problem_type) isa SecondOrderODEProblem`. Used as a guard in `_solve_with_solver_policy` before selecting `KahanLi8`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `prob` | Any | n/a | yes | Positional argument `prob`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_is_partitioned_second_order_problem`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_solver_policy_solve_with_solver_policy|_solve_with_solver_policy]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:639-639`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Relies on the SciMLBase convention that `problem_type` carries the marker; a `DynamicalODEProblem` built differently may be rejected. Does not verify that the partition corresponds to translational states specifically.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 397.
