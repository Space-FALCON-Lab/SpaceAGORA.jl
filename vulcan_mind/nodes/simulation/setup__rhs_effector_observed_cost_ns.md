---
id: simulation.setup__rhs_effector_observed_cost_ns
label: _rhs_effector_observed_cost_ns
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_effector_observed_cost_ns
  lines:
  - 585
  - 585
inputs:
- id: shared_buffers
  type: Any
  units: n/a
  required: true
  description: Positional argument `shared_buffers`.
- id: eff_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `eff_idx`.
- id: fallback_ns
  type: Float64
  units: n/a
  required: true
  description: Positional argument `fallback_ns`.
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
  description: Return value of `_rhs_effector_observed_cost_ns`.
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

# _rhs_effector_observed_cost_ns

## Purpose
Returns the measured per-call cost of a specific effector index once enough samples exist, otherwise the caller's fallback, giving the flat-queue planner per-effector heterogeneity information.

## Design & Implementation
Signature `(shared_buffers, eff_idx::Int, fallback_ns::Float64)::Float64`. Returns `fallback_ns` when buffers are `nothing` or lack `:rhs_effector_cost_ns`/`:rhs_effector_cost_samples`. Dereferences both refs and, if `eff_idx` is within both lengths, reads `estimate = Float64(costs[eff_idx])` and returns it when `samples[eff_idx] >= rhs_effector_cost_min_samples` (from `_rhs_env_config_from_buffers`) and the estimate is finite and positive. Otherwise returns `fallback_ns`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `eff_idx` | Int | n/a | yes | Positional argument `eff_idx`. |
| in | `fallback_ns` | Float64 | n/a | yes | Positional argument `fallback_ns`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_rhs_effector_observed_cost_ns`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.dynamics_rhs__rhs_effector_estimated_cost_ns|_rhs_effector_estimated_cost_ns]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:357-357`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/setup.jl:593-593`
- `callees` → [[simulation.setup__rhs_env_config_from_buffers|_rhs_env_config_from_buffers]] · `callers` · call · `src/simulation/engine/setup.jl:594-594`
<!-- vulcan:connections:end -->

## Limitations
Calls `_rhs_env_config_from_buffers` on every invocation; when the snapshot is unset that path re-parses dozens of environment variables per effector per RHS call. Out-of-range `eff_idx` silently returns the fallback instead of signalling a mismatch between the plan and the buffers.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 585.
