---
id: parallel.env_config_parse_bool_env
label: parse_bool_env
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: parse_bool_env
  lines:
  - 1
  - 1
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: Bool
  units: n/a
  required: true
  description: Positional argument `default`.
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
  type: Bool
  units: n/a
  description: Return value of `parse_bool_env`.
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

# parse_bool_env

## Purpose
Reads a boolean feature flag from the process environment with a typed default, accepting the common spellings operators use in shell scripts so that `SPACEAGORA_*` switches behave consistently across the parallel policy layer.

## Design & Implementation
Takes `name::String` and `default::Bool`. `get(ENV, name, default ? "1" : "0")` supplies the default as a string, then the raw value is `strip`ped and lowercased. Membership in `("1","true","yes","on")` returns `true`; `("0","false","no","off")` returns `false`. Any other text throws `ArgumentError` listing the accepted spellings. Marked `@inline`; performs a `Dict` lookup on `ENV` every call, so callers on hot paths snapshot the result.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | Bool | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `parse_bool_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.env_config__persistent_hint_path|_persistent_hint_path]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:112-112`
- [[parallel.env_config__safe_token|_safe_token]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:73-73`
- [[parallel.env_config__snapshot_auto_min_budget|_snapshot_auto_min_budget]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:208-208`
- [[parallel.env_config_persistent_hints_persist_enabled|persistent_hints_persist_enabled]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:76-76`
- [[parallel.env_config_persistent_hints_state_reset_requested|persistent_hints_state_reset_requested]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:80-80`
- [[parallel.thread_execution_threaded_foreach_worker_persistent|threaded_foreach_worker_persistent]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:156-156`
- [[parcore.env_config_snapshot_policy_decision_env|snapshot_policy_decision_env]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:208-208`
- [[simulation.assembly__requires_staged_density_callback|_requires_staged_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:36-36`
- [[simulation.from_env__solver_config_from_env|_solver_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:125-125`
- [[simulation.setup__effector_allow_with_outer|_effector_allow_with_outer]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:410-410`
- [[simulation.setup__effector_heavy_only|_effector_heavy_only]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:414-414`
- [[simulation.setup__ephemeris_reuse_enabled|_ephemeris_reuse_enabled]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:233-233`
- [[simulation.setup__harmonics_batch_allow_with_outer|_harmonics_batch_allow_with_outer]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:837-837`
- [[simulation.setup__nbody_ephemeris_cache_enabled|_nbody_ephemeris_cache_enabled]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:205-205`
- [[simulation.setup__planet_frame_cache_enabled|_planet_frame_cache_enabled]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:217-217`
- [[simulation.setup__rhs_harmonics_batch_enabled|_rhs_harmonics_batch_enabled]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:808-808`
- [[simulation.setup__rhs_harmonics_flat_experimental_enabled|_rhs_harmonics_flat_experimental_enabled]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:804-804`
- [[simulation.setup__rhs_plan_step_cache_enabled|_rhs_plan_step_cache_enabled]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:997-997`
- [[simulation.setup__spice_rhs_memo_enabled|_spice_rhs_memo_enabled]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:229-229`
- [[simulation.setup__srp_ephemeris_cache_enabled|_srp_ephemeris_cache_enabled]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:193-193`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
An empty string set in the environment (`export FLAG=`) is not treated as the default and throws, since `""` matches neither list. Reads `ENV` directly, so it is not thread-safe against concurrent `ENV` mutation and reflects overrides applied by the engine at call time rather than at process start.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 1.
