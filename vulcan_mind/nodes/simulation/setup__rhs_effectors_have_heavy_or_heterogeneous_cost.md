---
id: simulation.setup__rhs_effectors_have_heavy_or_heterogeneous_cost
label: _rhs_effectors_have_heavy_or_heterogeneous_cost
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_effectors_have_heavy_or_heterogeneous_cost
  lines:
  - 966
  - 966
inputs:
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
- id: heterogeneity_threshold
  type: Float64
  units: n/a
  required: true
  description: Positional argument `heterogeneity_threshold`.
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
  description: Return value of `_rhs_effectors_have_heavy_or_heterogeneous_cost`.
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

# _rhs_effectors_have_heavy_or_heterogeneous_cost

## Purpose
Decides whether the effector set is expensive or uneven enough to justify the flat constellation effector queue.

## Design & Implementation
Scans the effectors for their static cost estimate, tracking the minimum and maximum and flagging N-body, aerodynamic and harmonics effectors. Returns true if any of those three is present or if the max-to-min cost ratio, with the minimum floored at one nanosecond, reaches `heterogeneity_threshold`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `heterogeneity_threshold` | Float64 | n/a | yes | Positional argument `heterogeneity_threshold`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_rhs_effectors_have_heavy_or_heterogeneous_cost`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1237-1237`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__rhs_effector_static_cost_ns|_rhs_effector_static_cost_ns]] · `callers` · call · `src/simulation/engine/setup.jl:976-976`
<!-- vulcan:connections:end -->

## Limitations
The static cost table is a hard-coded estimate per type rather than a measurement, so an unusually cheap harmonics model still triggers the heavy path.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 966.
