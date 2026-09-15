---
id: simulation.dynamics_rhs__rhs_effector_static_cost_ns
label: _rhs_effector_static_cost_ns
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _rhs_effector_static_cost_ns
  lines:
  - 350
  - 350
inputs:
- id: effector
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector`.
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
  type: Float64
  units: n/a
  description: Return value of `_rhs_effector_static_cost_ns`.
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

# _rhs_effector_static_cost_ns

## Purpose
Converts an effector's cost rank into an estimated nanoseconds per item for the cost model's cold start.

## Design & Implementation
Returns the default per-item cost from the environment snapshot times the rank squared, so rank 4 is sixteen times rank 1. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_rhs_effector_static_cost_ns`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__rhs_effector_estimated_cost_ns|_rhs_effector_estimated_cost_ns]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:356-356`
- [[simulation.setup__rhs_effectors_have_heavy_or_heterogeneous_cost|_rhs_effectors_have_heavy_or_heterogeneous_cost]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:976-976`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:352-352`
- `callees` → [[simulation.dynamics_rhs__rhs_effector_cost_rank|_rhs_effector_cost_rank]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:351-351`
- `callees` → [[simulation.setup__effector_cost_ns_per_item_default|_effector_cost_ns_per_item_default]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:352-352`
<!-- vulcan:connections:end -->

## Limitations
The quadratic scaling is a heuristic chosen to spread the ranks, not a measurement.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 350.
