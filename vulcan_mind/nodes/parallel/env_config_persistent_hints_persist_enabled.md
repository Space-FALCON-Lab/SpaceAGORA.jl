---
id: parallel.env_config_persistent_hints_persist_enabled
label: persistent_hints_persist_enabled
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: persistent_hints_persist_enabled
  lines:
  - 75
  - 75
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
  description: Return value of `persistent_hints_persist_enabled`.
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

# persistent_hints_persist_enabled

## Purpose
Decides whether learned thread-policy hints should be written back to disk at the end of a run, defaulting to the same value as the persistent-hints master switch.

## Design & Implementation
Calls `parse_bool_env("SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST", persistent_hints_enabled())`. Because the default argument is evaluated eagerly, `persistent_hints_enabled()` (itself a `parse_bool_env` on `SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS`) is consulted every call even when the persist variable is explicitly set. Returns `Bool`; throws `ArgumentError` if either variable contains an unrecognised spelling.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `persistent_hints_persist_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/env_config.jl`
- [[parallel.persistent_hints__save_persistent_hint_state_locked_bang|_save_persistent_hint_state_locked!]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:113-113`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/parallel/policy/env_config.jl:76-76`
<!-- vulcan:connections:end -->

## Limitations
Setting persist to true while hints are disabled is allowed but ineffective, since nothing accumulates to persist. Two environment reads per call; the result is not cached in `PolicyDecisionEnvConfig`.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 75.
