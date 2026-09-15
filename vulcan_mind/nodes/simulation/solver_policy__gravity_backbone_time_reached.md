---
id: simulation.solver_policy__gravity_backbone_time_reached
label: _gravity_backbone_time_reached
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _gravity_backbone_time_reached
  lines:
  - 178
  - 178
inputs:
- id: t_now
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t_now`.
- id: t_target
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t_target`.
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
  description: Return value of `_gravity_backbone_time_reached`.
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

# _gravity_backbone_time_reached

## Purpose
Tolerant equality test for whether the backbone integrator's cursor time reached the intended segment end, guarding against floating-point drift in accumulated steps.

## Design & Implementation
Computes `tol = max(1e-9, 64 * eps(Float64) * max(1.0, abs(t_target)))` and returns `abs(t_now - t_target) <= tol`. The absolute floor of 1 ns covers small times; the relative term scales with `t_target` for long simulations.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t_now` | Float64 | n/a | yes | Positional argument `t_now`. |
| in | `t_target` | Float64 | n/a | yes | Positional argument `t_target`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_gravity_backbone_time_reached`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.solver_policy__solve_with_gravity_backbone_solver|_solve_with_gravity_backbone_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:618-618`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:388-388`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The 64-ulp factor and 1e-9 s floor are hard-coded. For `t_target` around 1e7 s the tolerance is roughly 1.4e-7 s, which is looser than the fixed step in some fine configurations and could accept an early stop as complete.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 178.
