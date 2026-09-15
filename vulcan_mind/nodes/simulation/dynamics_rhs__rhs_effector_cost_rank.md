---
id: simulation.dynamics_rhs__rhs_effector_cost_rank
label: _rhs_effector_cost_rank
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _rhs_effector_cost_rank
  lines:
  - 333
  - 333
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
  type: Int
  units: n/a
  description: Return value of `_rhs_effector_cost_rank`.
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

# _rhs_effector_cost_rank

## Purpose
Assigns a static relative-cost rank to each effector type, the fallback the cost model uses before any measured timings exist.

## Design & Implementation
Returns 4 for harmonics and N-body, 3 for SRP, 2 for aerodynamics and J2, and 1 for point-mass gravity or anything unrecognised. Declared `@inline` with an `::Int` return. `_rhs_effector_static_cost_ns` squares this rank.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_rhs_effector_cost_rank`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__rhs_effector_static_cost_ns|_rhs_effector_static_cost_ns]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:351-351`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A hand-tuned table that does not distinguish a degree-8 harmonics field from a degree-165 one, nor a two-body N-body model from a ten-body one.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 333.
