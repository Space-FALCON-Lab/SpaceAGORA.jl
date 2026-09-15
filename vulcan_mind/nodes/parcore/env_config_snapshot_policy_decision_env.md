---
id: parcore.env_config_snapshot_policy_decision_env
label: snapshot_policy_decision_env
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: snapshot_policy_decision_env
  lines:
  - 185
  - 233
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ParallelPolicy namespace supplying PolicyDecisionEnvConfig and the
    per-variable parsers defined in this file.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: env_config
  type: PolicyDecisionEnvConfig
  units: n/a
  description: Immutable snapshot of every parsed SPACEAGORA_* threading variable,
    taken once and reused across many policy decisions.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parcore
origin: agent
---

# snapshot_policy_decision_env

## Purpose
`snapshot_policy_decision_env` reads the threading-related environment variables once and returns them as a `PolicyDecisionEnvConfig`. Call sites that make policy decisions in a tight loop hold the snapshot and pass it to `thread_policy_decision`, which removes the repeated dictionary lookup and string parsing from the hot path.

## Model & Assumptions
The snapshot is a point-in-time value. Changing an environment variable after a snapshot has been taken has no effect on holders of that snapshot, which is intentional: a decision sequence within one propagation stays internally consistent even if something else mutates `ENV`. Consumers that want live behaviour simply pass no snapshot at all.

## Design & Implementation
The file owns every parser for the threading configuration surface. `parse_bool_env` accepts 1/0, true/false, yes/no and on/off and throws an `ArgumentError` naming the variable and the offending value for anything else. `parse_parallel_mode_env` maps off/none/serial/0/false/no to `:off`, on/thread/threads/1/true/yes to `:on`, and the literal auto to `:auto`. `parse_thread_threshold_env` and `parse_nonnegative_int_env` parse integers and clamp to one and zero respectively, `parse_float_env` parses floats, and `inner_scheduler_mode` accepts static, strided or dynamic. The snapshot function calls these in sequence and packs the results into the record, so validation happens once at snapshot time rather than at each decision.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ParallelPolicy namespace supplying PolicyDecisionEnvConfig and the per-variable parsers defined in this file. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `env_config` | PolicyDecisionEnvConfig | n/a | — | Immutable snapshot of every parsed SPACEAGORA_* threading variable, taken once and reused across many policy decisions. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.env_config__telemetry_bucket|_telemetry_bucket]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:178-178`
- [[simulation.setup__initialize_runtime_env_config_bang|_initialize_runtime_env_config!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:913-913`

**Downstream**

- `callees` → [[core.runtime_types_policydecisionenvconfig|PolicyDecisionEnvConfig]] · `callers` · call · `src/parallel/policy/env_config.jl:186-186`
- `callees` → [[parallel.env_config__adaptive_desire_cap|_adaptive_desire_cap]] · `callers` · call · `src/parallel/policy/env_config.jl:231-231`
- `callees` → [[parallel.env_config__snapshot_auto_min_budget|_snapshot_auto_min_budget]] · `callers` · call · `src/parallel/policy/env_config.jl:197-197`
- `callees` → [[parallel.env_config_adaptive_delta|adaptive_delta]] · `callers` · call · `src/parallel/policy/env_config.jl:215-215`
- `callees` → [[parallel.env_config_adaptive_rho|adaptive_rho]] · `callers` · call · `src/parallel/policy/env_config.jl:223-223`
- `callees` → [[parallel.env_config_auto_thread_min_budget|auto_thread_min_budget]] · `callers` · call · `src/parallel/policy/env_config.jl:189-189`
- `callees` → [[parallel.env_config_effective_inner_thread_budget|effective_inner_thread_budget]] · `callers` · call · `src/parallel/policy/env_config.jl:187-187`
- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/parallel/policy/env_config.jl:208-208`
- `callees` → [[parallel.env_config_parse_float_env|parse_float_env]] · `callers` · call · `src/parallel/policy/env_config.jl:216-216`
- `callees` → [[parallel.env_config_parse_nonnegative_int_env|parse_nonnegative_int_env]] · `callers` · call · `src/parallel/policy/env_config.jl:210-210`
- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/parallel/policy/env_config.jl:209-209`
<!-- vulcan:connections:end -->

## Limitations
Strict parsing means a typo in a variable name is silent — the default is used — while a typo in a value throws, which is an asymmetry worth knowing when debugging a configuration. The snapshot cannot represent a variable that legitimately changes during a run, such as one toggled by an outer route as it switches strategies; those call sites must re-snapshot. Every field is eagerly computed, so taking a snapshot to read one value still parses all of them.

## Provenance
Mapped from `src/parallel/policy/env_config.jl:185-233`.
