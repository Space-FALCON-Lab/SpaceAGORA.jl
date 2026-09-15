---
id: simulation.setup__effector_cost_min_samples
label: _effector_cost_min_samples
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_cost_min_samples
  lines:
  - 421
  - 421
inputs:
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
  description: Return value of `_effector_cost_min_samples`.
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

# _effector_cost_min_samples

## Purpose
Number of timing samples the effector-cost EMA must accumulate before the runtime estimate replaces the static prior in threading decisions.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_EFFECTOR_COST_MIN_SAMPLES", 4)`, clamped to at least 1. Compared against `shared_buffers.effector_cost_samples[]` in `_effector_observed_cost_ns_per_item`; the sample counter is incremented by `_update_effector_cost_model!` once per RHS evaluation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_effector_cost_min_samples`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:861-861`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:422-422`
<!-- vulcan:connections:end -->

## Limitations
The first sample after a cold start includes JIT compilation time and can skew the EMA for many steps; four samples with α = 0.2 do not fully wash that out. There is no upper bound, and a very high value effectively pins the decision to the prior.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 421.
