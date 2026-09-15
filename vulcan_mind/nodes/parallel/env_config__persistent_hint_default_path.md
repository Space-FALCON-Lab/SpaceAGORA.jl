---
id: parallel.env_config__persistent_hint_default_path
label: _persistent_hint_default_path
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: _persistent_hint_default_path
  lines:
  - 92
  - 92
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
  type: String
  units: n/a
  description: Return value of `_persistent_hint_default_path`.
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

# _persistent_hint_default_path

## Purpose
Computes the default location of the persisted inner-policy state TOML, keyed by profile, machine label and thread count so hints learned under one configuration are not applied to another.

## Design & Implementation
Reads `SPACEAGORA_PARALLEL_PROFILE` and `SPACEAGORA_PERF_MACHINE_LABEL` (both defaulting to `"default"`), passes each through `_safe_token`, and takes `Threads.nthreads()` as `threads`. Returns `joinpath(pwd(), "output", "parallel_policy_state", "inner_policy_state_<profile>_<machine>_t<threads>.toml")`. The directory is not created here.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_persistent_hint_default_path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/env_config.jl`
- [[parallel.env_config__persistent_hint_path|_persistent_hint_path]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:107-107`

**Downstream**

- `callees` → [[parallel.env_config__safe_token|_safe_token]] · `callers` · call · `src/parallel/policy/env_config.jl:93-93`
<!-- vulcan:connections:end -->

## Limitations
Anchoring on `pwd()` means the same run launched from a different working directory reads a different state file. The thread count in the name only captures the default pool size, not `SPACEAGORA_INNER_THREAD_BUDGET`, so runs with different budgets on the same machine share hints.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 92.
