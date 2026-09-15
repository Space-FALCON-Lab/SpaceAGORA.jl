---
id: parallel.env_config_persistent_hints_state_reset_requested
label: persistent_hints_state_reset_requested
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: persistent_hints_state_reset_requested
  lines:
  - 79
  - 79
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
  type: Bool
  units: n/a
  description: Return value of `persistent_hints_state_reset_requested`.
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

# persistent_hints_state_reset_requested

## Purpose
Signals that any previously persisted thread-policy state file should be discarded before the run starts, letting an operator clear stale hints after a hardware or workload change.

## Design & Implementation
Returns `parse_bool_env("SPACEAGORA_PARALLEL_POLICY_STATE_RESET", false)`. The default is `false` so state is preserved unless explicitly requested. The function only reports the request; deletion or ignoring of the file at `_persistent_hint_path()` is performed by the policy state loader that calls it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `persistent_hints_state_reset_requested`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/env_config.jl`
- [[parallel.persistent_hints__load_persistent_hint_state_locked_bang|_load_persistent_hint_state_locked!]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:65-65`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:170-170`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/parallel/policy/env_config.jl:80-80`
<!-- vulcan:connections:end -->

## Limitations
The flag is read from `ENV` on each call, so a caller checking it more than once during a run could observe different answers if the environment is mutated. It does not verify that a state file exists or that the path is writable.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 79.
