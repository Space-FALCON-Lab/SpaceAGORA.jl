---
id: simulation.adaptive_routing__run_campaign_with_route_env
label: _run_campaign_with_route_env
kind: function
source:
  file: src/simulation/campaigns/adaptive_routing.jl
  symbol: _run_campaign_with_route_env
  lines:
  - 172
  - 172
inputs:
- id: f
  type: Any
  units: n/a
  required: true
  description: Positional argument `f`.
- id: spec
  type: MonteCarloSpec
  units: n/a
  required: true
  description: Positional argument `spec`.
- id: plan
  type: Any
  units: n/a
  required: true
  description: Positional argument `plan`.
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
  type: MonteCarloResult
  units: n/a
  description: Return value of `_run_campaign_with_route_env`. Returns `MonteCarloResult(samples,
    elapsed_s, length(active_workers))` or `withenv(env_pairs...) do`.
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

# _run_campaign_with_route_env

## Purpose
Runs a `MonteCarloSpec` according to a route plan, dispatching to the process pool for `:process`, or wrapping `run_monte_carlo` in the environment variables that mark the outer split as active and cap inner thread usage for `:threads`.

## Design & Implementation
`worker_count = min(spec.threads, length(spec.seeds))`; if it is 1 or less the spec is run directly via `run_monte_carlo(f, spec)`. For `plan.route === :process`, `campaign_process_pool()` is obtained, `ensure_process_workers!(pool, worker_count; warmup_fn=() -> f(first(spec.seeds)))` grows the pool and pays JIT cost up-front, and `_run_monte_carlo_process` executes the seeds on the first `worker_count` worker ids while `time_ns()` brackets the elapsed time. `spec.fail_fast` triggers `_throw_first_monte_carlo_failure`. Otherwise `withenv` sets `SPACEAGORA_OUTER_PARALLEL_ACTIVE="1"` and, only if `SPACEAGORA_INNER_THREAD_BUDGET` is unset or blank, `plan.inner_thread_budget`, before calling `run_monte_carlo`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | Any | n/a | yes | Positional argument `f`. |
| in | `spec` | MonteCarloSpec | n/a | yes | Positional argument `spec`. |
| in | `plan` | Any | n/a | yes | Positional argument `plan`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | MonteCarloResult | n/a | — | Return value of `_run_campaign_with_route_env`. Returns `MonteCarloResult(samples, elapsed_s, length(active_workers))` or `withenv(env_pairs...) do`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/campaigns/adaptive_routing.jl`
- [[simulation.adaptive_routing__run_campaign_adaptive|_run_campaign_adaptive]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:245-245`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:184-184`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:197-197`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:184-184`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:184-184`
- `callees` → [[parallel.ensure_process_workers_ensure_process_workers_bang|ensure_process_workers!]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:184-184`
- `callees` → [[simulation.monte_carlo__run_monte_carlo_process|_run_monte_carlo_process]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:187-187`
- `callees` → [[simulation.monte_carlo__throw_first_monte_carlo_failure|_throw_first_monte_carlo_failure]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:189-189`
- `callees` → [[simulation.monte_carlo_result|MonteCarloResult]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:190-190`
- `callees` → [[simulation.run_monte_carlo|run_monte_carlo]] · `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:174-174`
<!-- vulcan:connections:end -->

## Limitations
The warm-up executes `f` on `first(spec.seeds)` once per newly added worker, which doubles the cost of that seed and assumes `f` is side-effect free. `withenv` mutates process-global `ENV`, which is not safe if other tasks read those variables concurrently. Timing for the process route includes pool warm-up only indirectly (excluded), while thread-route timing comes from `run_monte_carlo`, making the two routes' feedback not perfectly comparable. Fewer workers than requested may be returned by the pool without warning.

## Provenance
Mapped from `src/simulation/campaigns/adaptive_routing.jl` line 172.
