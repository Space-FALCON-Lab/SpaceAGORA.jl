---
id: parallel.env_config__persistent_hint_path
label: _persistent_hint_path
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: _persistent_hint_path
  lines:
  - 104
  - 104
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
  description: Return value of `_persistent_hint_path`.
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

# _persistent_hint_path

## Purpose
Resolves the final absolute path of the persisted policy state file, honouring an operator override before falling back to the derived default.

## Design & Implementation
Reads `SPACEAGORA_PARALLEL_POLICY_STATE_PATH` and strips it. When empty, returns `normpath(_persistent_hint_default_path())`. Otherwise, if `isabspath(override)` the override is normalised as-is; a relative override is joined onto `pwd()` first. Always returns a `String` that has passed through `normpath`, so `..` segments are collapsed.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_persistent_hint_path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/env_config.jl`
- [[parallel.persistent_hints__load_persistent_hint_state_locked_bang|_load_persistent_hint_state_locked!]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:64-64`

**Downstream**

- `callees` → [[parallel.env_config__persistent_hint_default_path|_persistent_hint_default_path]] · `callers` · call · `src/parallel/policy/env_config.jl:107-107`
- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/parallel/policy/env_config.jl:112-112`
<!-- vulcan:connections:end -->

## Limitations
No check that the parent directory exists or is writable; failures surface later at load or save time. Relative overrides are interpreted against the current working directory at call time, which can differ between the loading and saving calls if the process changes directory.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 104.
