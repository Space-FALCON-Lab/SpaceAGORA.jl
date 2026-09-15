---
id: simulation.setup__effector_observed_cost_ns_per_item
label: _effector_observed_cost_ns_per_item
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_observed_cost_ns_per_item
  lines:
  - 513
  - 513
inputs:
- id: env
  type: SimulationModel.RhsPlanEnvConfig
  units: n/a
  required: true
  description: Positional argument `env`.
- id: shared_buffers
  type: Any
  units: n/a
  required: true
  description: Positional argument `shared_buffers`.
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
  description: Return value of `_effector_observed_cost_ns_per_item`.
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

# _effector_observed_cost_ns_per_item

## Purpose
Returns the best available estimate of nanoseconds per effector evaluation, preferring the runtime EMA in `shared_buffers` once it has enough samples and otherwise the configured prior.

## Design & Implementation
Takes `env::RhsPlanEnvConfig` and `shared_buffers`. Starts with `default_cost = env.effector_cost_ns_per_item_default`. Returns it when `shared_buffers === nothing` or when either `:effector_cost_ns_per_item` or `:effector_cost_samples` is missing. Otherwise reads `samples = Int(...[])` and `estimate = Float64(...[])` from the two `Ref`s and returns `estimate` only if `samples >= env.effector_cost_min_samples && isfinite(estimate) && estimate > 0.0`. Pure; does not mutate the refs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `env` | SimulationModel.RhsPlanEnvConfig | n/a | yes | Positional argument `env`. |
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_effector_observed_cost_ns_per_item`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__dynamic_effector_thread_decision|_dynamic_effector_thread_decision]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:679-679`
- [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1232-1232`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/setup.jl:520-520`
<!-- vulcan:connections:end -->

## Limitations
The estimate is a single scalar averaged across all effectors and satellites, so it cannot express heterogeneity; the per-effector estimator `_rhs_effector_observed_cost_ns` exists for that. The `Int(...)` conversion on the samples ref would throw if a non-integer value were stored, which no writer does.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 513.
