---
id: parallel.env_config_effective_inner_thread_budget
label: effective_inner_thread_budget
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: effective_inner_thread_budget
  lines:
  - 122
  - 122
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
  description: Return value of `effective_inner_thread_budget`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# effective_inner_thread_budget

## Purpose
Determines how many threads inner parallel loops may use, capping an operator-supplied budget by the actual pool size and treating non-positive budgets as 'use everything'.

## Design & Implementation
Reads `SPACEAGORA_INNER_THREAD_BUDGET` with default `"0"`, parsing with `parse(Int, ...)` inside `try`/`catch` that rethrows `ArgumentError("SPACEAGORA_INNER_THREAD_BUDGET must be an integer, got '<raw>'")`. With `available = _default_thread_pool_size()`, a `budget <= 0` yields `max(1, available)`; otherwise `max(1, min(available, budget))`. Returns `Int`. The result is the first field of `PolicyDecisionEnvConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `effective_inner_thread_budget`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.thread_execution__thread_worker_count|_thread_worker_count]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:187-187`
- [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callees` → `callers` · call · `src/parallel/policy/adaptive_decision.jl:12-12`
- [[parcore.env_config_snapshot_policy_decision_env|snapshot_policy_decision_env]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:187-187`
- [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callees` → `callers` · call · `src/parallel/policy/observation_tracking.jl:15-15`
- [[simulation.rhs_calibration__rhs_calib_signature|_rhs_calib_signature]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:77-77`
- [[simulation.rhs_calibration__rhs_plan_candidates|_rhs_plan_candidates]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:221-221`
- [[simulation.setup__dynamic_effector_thread_decision|_dynamic_effector_thread_decision]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:670-670`
- [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1060-1060`
- [[simx.engine_rhs_calibration_calibrate_rhs_plan_if_needed__calibrate_rhs_plan_if_needed_bang|_calibrate_rhs_plan_if_needed!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:317-317`

**Downstream**

- `callees` → [[parallel.env_config__default_thread_pool_size|_default_thread_pool_size]] · `callers` · call · `src/parallel/policy/env_config.jl:129-129`
<!-- vulcan:connections:end -->

## Limitations
A budget larger than the pool is silently clipped with no log message. When an outer parallel layer is active this function does not subtract threads already in use; that coordination lives in the policy decision. Empty string input throws rather than defaulting to 0.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 122.
