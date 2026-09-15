---
id: parallel.env_config_auto_thread_min_budget
label: auto_thread_min_budget
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: auto_thread_min_budget
  lines:
  - 136
  - 136
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
  description: Return value of `auto_thread_min_budget`.
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

# auto_thread_min_budget

## Purpose
Defines the minimum thread budget below which `:auto` mode refuses to parallelise, with higher floors for lock-guarded density and thermal callbacks where oversubscription is costly.

## Design & Implementation
The zero-argument method returns `parse_thread_threshold_env("SPACEAGORA_AUTO_THREAD_MIN_BUDGET", 4)`. The `(source::Symbol)` method first computes that default, then: `:density_callback` reads `SPACEAGORA_DENSITY_CALLBACK_AUTO_THREAD_MIN_BUDGET` with default `max(default_budget, 16)`; `:density_callback_lockfree` reads `SPACEAGORA_DENSITY_CALLBACK_LOCKFREE_AUTO_THREAD_MIN_BUDGET` with default `default_budget` (the 16-thread floor is for locked GRAM only, as the inline comment explains); `:thermal_callback` reads `SPACEAGORA_THERMAL_CALLBACK_AUTO_THREAD_MIN_BUDGET` with default `max(default_budget, 16)`. Every other source returns `default_budget`. All values are clamped to at least 1 by the parser.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `auto_thread_min_budget`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callees` → `callers` · call · `src/parallel/policy/adaptive_decision.jl:13-13`
- [[parcore.env_config_snapshot_policy_decision_env|snapshot_policy_decision_env]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:189-189`
- [[simulation.setup__rhs_flat_min_thread_budget|_rhs_flat_min_thread_budget]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:799-799`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/parallel/policy/env_config.jl:137-137`
<!-- vulcan:connections:end -->

## Limitations
The 16-thread floor is hard-coded and assumes a locked native GRAM model; machines with fewer cores can never auto-parallelise density callbacks unless the override is set. Unknown `source` symbols silently receive the general default rather than erroring. Each call performs up to two `ENV` lookups, which is why `snapshot_policy_decision_env` caches the four variants.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 136.
