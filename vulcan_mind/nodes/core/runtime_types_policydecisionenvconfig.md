---
id: core.runtime_types_policydecisionenvconfig
label: PolicyDecisionEnvConfig
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: PolicyDecisionEnvConfig
  lines:
  - 578
  - 578
inputs:
- id: inner_thread_budget
  type: Int
  units: n/a
  required: true
  description: Field `inner_thread_budget`.
- id: outer_parallel_active
  type: Bool
  units: n/a
  required: true
  description: Field `outer_parallel_active`.
- id: auto_min_budget_default
  type: Int
  units: n/a
  required: true
  description: Field `auto_min_budget_default`.
- id: auto_min_budget_density
  type: Int
  units: n/a
  required: true
  description: Field `auto_min_budget_density`.
- id: auto_min_budget_density_lockfree
  type: Int
  units: n/a
  required: true
  description: Field `auto_min_budget_density_lockfree`.
- id: auto_min_budget_thermal
  type: Int
  units: n/a
  required: true
  description: Field `auto_min_budget_thermal`.
- id: adaptive_enabled
  type: Bool
  units: n/a
  required: true
  description: Field `adaptive_enabled`.
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
  type: PolicyDecisionEnvConfig
  units: n/a
  description: Constructed `PolicyDecisionEnvConfig`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# PolicyDecisionEnvConfig

## Purpose
The run-scoped snapshot of the environment variables that feed the inner threading policy decision, captured once so the hot path never parses `ENV`.

## Design & Implementation
An immutable struct with the inner thread budget, whether outer parallelism is active, four auto-mode minimum-budget thresholds (default, density, lock-free density, thermal) and the adaptive-enabled flag. Built at `run_simulation` setup and stored in `SharedBuffers.policy_env_config`; hot-path accessors fall back to live parsing only when the slot is `nothing`, which keeps hand-built test parameters working.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `inner_thread_budget` | Int | n/a | yes | Field `inner_thread_budget`. |
| in | `outer_parallel_active` | Bool | n/a | yes | Field `outer_parallel_active`. |
| in | `auto_min_budget_default` | Int | n/a | yes | Field `auto_min_budget_default`. |
| in | `auto_min_budget_density` | Int | n/a | yes | Field `auto_min_budget_density`. |
| in | `auto_min_budget_density_lockfree` | Int | n/a | yes | Field `auto_min_budget_density_lockfree`. |
| in | `auto_min_budget_thermal` | Int | n/a | yes | Field `auto_min_budget_thermal`. |
| in | `adaptive_enabled` | Bool | n/a | yes | Field `adaptive_enabled`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PolicyDecisionEnvConfig | n/a | — | Constructed `PolicyDecisionEnvConfig`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[grp.src_parallel_policy|parallel/policy/]] · `members_out` → `callers` · call · `src/parallel/policy/env_config.jl:186-186`
- [[parcore.env_config_snapshot_policy_decision_env|snapshot_policy_decision_env]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:186-186`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Values are frozen at run start, so changing an environment variable mid-run has no effect on that run — deliberate, but a source of confusion when tuning interactively.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 578.
