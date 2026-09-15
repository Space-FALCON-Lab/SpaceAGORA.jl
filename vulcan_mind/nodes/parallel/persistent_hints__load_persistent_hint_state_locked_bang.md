---
id: parallel.persistent_hints__load_persistent_hint_state_locked_bang
label: _load_persistent_hint_state_locked!
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _load_persistent_hint_state_locked!
  lines:
  - 58
  - 58
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
  type: Nothing
  units: n/a
  description: Return value of `_load_persistent_hint_state_locked!`.
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

# _load_persistent_hint_state_locked!

## Purpose
Reads the on-disk hint history into the process-wide state exactly once, honouring the cold-start and disabled modes.

## Design & Implementation
Returns immediately if `state.loaded` is already set, then marks it loaded and records the resolved path. If a state reset was requested it empties the history and returns, deliberately ignoring the file; if hints are disabled or the file is absent it returns with an empty history. Otherwise it parses the TOML, and for each row with a non-empty `signature`, positive `allotment` and parseable `stats`, it merges the counts into the matching bucket with `get!`, summing rather than replacing so a file with duplicate rows still loads coherently. Parse failures are swallowed.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_load_persistent_hint_state_locked!`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- [[parallel.persistent_hints__ensure_persistent_hint_state_loaded_bang|_ensure_persistent_hint_state_loaded!]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:146-146`

**Downstream**

- `callees` → [[parallel.env_config__persistent_hint_path|_persistent_hint_path]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:64-64`
- `callees` → [[parallel.env_config_persistent_hints_state_reset_requested|persistent_hints_state_reset_requested]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:65-65`
- `callees` → [[parallel.persistent_hints__hint_payload_stats|_hint_payload_stats]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:94-94`
- `callees` → [[parallel.types_adaptivechoicestats|AdaptiveChoiceStats]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:100-100`
<!-- vulcan:connections:end -->

## Limitations
It must be called with `_persistent_hint_lock` held, which the name records but nothing enforces; `state.loaded` is set before parsing, so a file that fails to parse is never retried in this process.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 58.
