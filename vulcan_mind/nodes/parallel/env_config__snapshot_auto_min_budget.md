---
id: parallel.env_config__snapshot_auto_min_budget
label: _snapshot_auto_min_budget
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: _snapshot_auto_min_budget
  lines:
  - 197
  - 197
inputs:
- id: env
  type: PolicyDecisionEnvConfig
  units: n/a
  required: true
  description: Positional argument `env`.
- id: source
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `source`.
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
  description: Return value of `_snapshot_auto_min_budget`.
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

# _snapshot_auto_min_budget

## Purpose
Looks up the pre-computed auto-mode minimum thread budget for a given policy source from a `PolicyDecisionEnvConfig` snapshot, avoiding `ENV` access on the hot decision path.

## Design & Implementation
Takes `env::PolicyDecisionEnvConfig` and `source::Symbol`. Returns `env.auto_min_budget_density` for `:density_callback`, `env.auto_min_budget_density_lockfree` for `:density_callback_lockfree`, `env.auto_min_budget_thermal` for `:thermal_callback`, and `env.auto_min_budget_default` otherwise. The branch order mirrors `auto_thread_min_budget(source)` exactly so the snapshot and live paths agree. `@inline`, allocation-free, never throws.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `env` | PolicyDecisionEnvConfig | n/a | yes | Positional argument `env`. |
| in | `source` | Symbol | n/a | yes | Positional argument `source`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_snapshot_auto_min_budget`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callees` → `callers` · call · `src/parallel/policy/adaptive_decision.jl:13-13`
- [[parcore.env_config_snapshot_policy_decision_env|snapshot_policy_decision_env]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:197-197`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/parallel/policy/env_config.jl:208-208`
- `callees` → [[parallel.env_config_parse_nonnegative_int_env|parse_nonnegative_int_env]] · `callers` · call · `src/parallel/policy/env_config.jl:210-210`
- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/parallel/policy/env_config.jl:209-209`
<!-- vulcan:connections:end -->

## Limitations
Correctness depends on `snapshot_policy_decision_env` having populated the four fields from the same `auto_thread_min_budget` calls; any drift between the two switch statements would go unnoticed at compile time. Sources added later fall through to the default without warning.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 197.
