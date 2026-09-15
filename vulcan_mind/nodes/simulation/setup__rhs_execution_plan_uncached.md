---
id: simulation.setup__rhs_execution_plan_uncached
label: _rhs_execution_plan_uncached
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_execution_plan_uncached
  lines:
  - 1039
  - 1039
inputs:
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
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
  type: SimulationModel.RhsExecutionPlan
  units: n/a
  description: Return value of `_rhs_execution_plan_uncached`.
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

# _rhs_execution_plan_uncached

## Purpose
The routing chain that decides, for one RHS evaluation, whether satellites and effectors are evaluated serially, in satellite batches, per-satellite with inner effector threading, or through the flat constellation effector queue.

## Design & Implementation
An override from calibration returns immediately. Otherwise it reads the env and policy snapshots, counts active satellites, resolves the thread budget and outer-parallel state, and computes the effector-level decision. A forced `execution_mode` short-circuits to `:serial`, `:satellite_batch`, `:per_satellite_effector_reduce` or the flat queue with viability checks. Heuristic routing then tries single-harmonics flat, single-inverse-square flat, falls back to satellite batch when satellites or budget are one or flat is unsupported or the budget is below the flat minimum, and enters the flat queue only when satellite count, effector count and estimated work — from observed per-item cost — all exceed thresholds and outer parallelism permits. If satellite batching saturates the pool it is chosen; otherwise per-satellite effector reduction with the inner decision. Returns an isbits `RhsExecutionPlan`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.RhsExecutionPlan | n/a | — | Return value of `_rhs_execution_plan_uncached`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__rhs_execution_plan|_rhs_execution_plan]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1032-1032`

**Downstream**

- `callees` → [[parallel.env_config_effective_inner_thread_budget|effective_inner_thread_budget]] · `callers` · call · `src/simulation/engine/setup.jl:1060-1060`
- `callees` → [[simulation.config__policy_env_config|_policy_env_config]] · `callers` · call · `src/simulation/engine/setup.jl:1051-1051`
- `callees` → [[simulation.setup__dynamic_effector_thread_decision|_dynamic_effector_thread_decision]] · `callers` · call · `src/simulation/engine/setup.jl:1065-1065`
- `callees` → [[simulation.setup__effector_observed_cost_ns_per_item|_effector_observed_cost_ns_per_item]] · `callers` · call · `src/simulation/engine/setup.jl:1232-1232`
- `callees` → [[simulation.setup__effector_shared_buffers|_effector_shared_buffers]] · `callers` · call · `src/simulation/engine/setup.jl:1231-1231`
- `callees` → [[simulation.setup__policy_env_config|_policy_env_config]] · `callers` · call · `src/simulation/engine/setup.jl:1051-1051`
- `callees` → [[simulation.setup__rhs_effectors_have_heavy_or_heterogeneous_cost|_rhs_effectors_have_heavy_or_heterogeneous_cost]] · `callers` · call · `src/simulation/engine/setup.jl:1237-1237`
- `callees` → [[simulation.setup__rhs_flat_has_batch_privileged_effector|_rhs_flat_has_batch_privileged_effector]] · `callers` · call · `src/simulation/engine/setup.jl:1242-1242`
- `callees` → [[simulation.setup__rhs_flat_supported|_rhs_flat_supported]] · `callers` · call · `src/simulation/engine/setup.jl:1104-1104`
- `callees` → [[simulation.setup__rhs_invsq_flat_min_sats|_rhs_invsq_flat_min_sats]] · `callers` · call · `src/simulation/engine/setup.jl:1177-1177`
- `callees` → [[simulation.setup__rhs_single_harmonics_flat_supported|_rhs_single_harmonics_flat_supported]] · `callers` · call · `src/simulation/engine/setup.jl:1135-1135`
- `callees` → [[simulation.setup__rhs_single_invsq_flat_supported|_rhs_single_invsq_flat_supported]] · `callers` · call · `src/simulation/engine/setup.jl:1176-1176`
- `callees` → [[simulation.setup__satellite_batch_saturates_pool|_satellite_batch_saturates_pool]] · `callers` · call · `src/simulation/engine/setup.jl:1268-1268`
- `callees` → [[simulation.setup__with_serial_effector_decision|_with_serial_effector_decision]] · `callers` · call · `src/simulation/engine/setup.jl:1077-1077`
<!-- vulcan:connections:end -->

## Limitations
Roughly 260 lines with eleven return sites, most of which duplicate the satellite-batch tuple; the ordering of heuristics is the policy, and there is no trace of which branch fired beyond the returned `dominant_axis`.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1039.
