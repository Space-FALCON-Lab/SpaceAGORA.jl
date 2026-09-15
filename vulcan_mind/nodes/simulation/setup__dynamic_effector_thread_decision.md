---
id: simulation.setup__dynamic_effector_thread_decision
label: _dynamic_effector_thread_decision
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _dynamic_effector_thread_decision
  lines:
  - 625
  - 625
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
- id: active_sats
  type: Int
  units: n/a
  required: false
  description: Keyword argument `active_sats` (default `num_sats`).
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
  type: Any
  units: n/a
  description: Return value of `_dynamic_effector_thread_decision`. Returns `_dynamic_effector_thread_decision(`.
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

# _dynamic_effector_thread_decision

## Purpose
Central policy that decides, for one RHS evaluation, whether the per-effector force loop should run on multiple threads and how many, combining effector thread-safety, active satellite count, thread budget sharing, cost estimates, and the global parallel policy.

## Design & Implementation
Three methods. The core `(env::RhsPlanEnvConfig, penv::Union{Nothing,PolicyDecisionEnvConfig}, args::SimulationConfiguration, p, dynamic_effectors::Tuple, num_sats::Int; active_sats::Int=num_sats)` returns a named tuple `(use_threads, allotment, mode, policy_applied)`. It answers serial with `policy_applied=false` when `n_effectors <= 1`, when `_dynamic_effectors_parallel_supported` is false, or when `active_sats <= 1 && mode != :on` (source comments cite 1.5–2× slowdowns for single-satellite threading). Otherwise it computes `outer_active` and `budget` from `penv` (falling back to live `ParallelPolicy` reads), `share_budget = _effector_satellite_share_budget(num_sats, budget)` halved if `outer_active && allow_with_outer`, an `inner_floor` of `min(2, budget)` when no outer layer and more than one satellite, and `max_allotment = min(env.effector_max_threads, budget, max(share_budget, inner_floor))`. Work is estimated as `per_effector_cost_ns * n_effectors / max_allotment` and compared with `effector_work_ns_per_worker_threshold × (outer_active ? effector_outer_work_scale : 1)` to set `heavy_work`. `ParallelPolicy.thread_policy_decision` is then called with the mode, threshold, heaviness flags, outer flags, `source=:dynamic_effectors`, and `env=penv`; the final allotment is `min(policy.allotment, max_allotment)` and `use_threads` requires `allotment > 1`. The two convenience methods derive `env`/`penv` from `p` via `_rhs_env_config` and `_policy_env_config`, or pass `p = nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `active_sats` | Int | n/a | no | Keyword argument `active_sats` (default `num_sats`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_dynamic_effector_thread_decision`. Returns `_dynamic_effector_thread_decision(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1065-1065`

**Downstream**

- `callees` → [[parallel.env_config_effective_inner_thread_budget|effective_inner_thread_budget]] · `callers` · call · `src/simulation/engine/setup.jl:670-670`
- `callees` → [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callers` · call · `src/simulation/engine/setup.jl:685-685`
- `callees` → [[simulation.config__policy_env_config|_policy_env_config]] · `callers` · call · `src/simulation/engine/setup.jl:633-633`
- `callees` → [[simulation.setup__dynamic_effectors_parallel_supported|_dynamic_effectors_parallel_supported]] · `callers` · call · `src/simulation/engine/setup.jl:652-652`
- `callees` → [[simulation.setup__effector_observed_cost_ns_per_item|_effector_observed_cost_ns_per_item]] · `callers` · call · `src/simulation/engine/setup.jl:679-679`
- `callees` → [[simulation.setup__effector_outer_parallel_hint|_effector_outer_parallel_hint]] · `callers` · call · `src/simulation/engine/setup.jl:667-667`
- `callees` → [[simulation.setup__effector_satellite_share_budget|_effector_satellite_share_budget]] · `callers` · call · `src/simulation/engine/setup.jl:671-671`
- `callees` → [[simulation.setup__effector_shared_buffers|_effector_shared_buffers]] · `callers` · call · `src/simulation/engine/setup.jl:678-678`
- `callees` → [[simulation.setup__policy_env_config|_policy_env_config]] · `callers` · call · `src/simulation/engine/setup.jl:633-633`
<!-- vulcan:connections:end -->

## Limitations
When `penv === nothing` it reads `ENV` twice per call, which on the hot path is expensive; run_simulation always supplies a snapshot but tests may not. The `inner_floor` heuristic can grant two threads even when the cost model says the work is light, relying on the policy's `heavy_only` gate to refuse. The `args` argument is accepted but unused in the core method, retained for signature compatibility.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 625.
