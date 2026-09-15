---
id: parallel.env_config__safe_token
label: _safe_token
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: _safe_token
  lines:
  - 67
  - 67
inputs:
- id: raw
  type: Any
  units: n/a
  required: true
  description: Positional argument `raw`.
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
  description: Return value of `_safe_token`.
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

# _safe_token

## Purpose
Sanitises an arbitrary label (profile name, machine label) into a filesystem-safe token used in persistent policy state filenames.

## Design & Implementation
Accepts any `raw` convertible via `String(raw)`, lowercases and strips it, then applies `replace(token, r"[^a-z0-9._-]+" => "_")` so every run of characters outside `[a-z0-9._-]` collapses to a single underscore. An empty result (including all-whitespace input) becomes the literal `"default"`. Returns a `String`; used twice in `_persistent_hint_default_path`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | Any | n/a | yes | Positional argument `raw`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_safe_token`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/env_config.jl`
- [[parallel.env_config__persistent_hint_default_path|_persistent_hint_default_path]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:93-93`
- [[parallel.persistent_hints__hint_workload_signature|_hint_workload_signature]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:195-195`
- [[parallel.persistent_hints_hint_layer_stats_snapshot|hint_layer_stats_snapshot]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:389-389`
- [[simulation.rhs_calibration__calib_machine_label|_calib_machine_label]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:51-51`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/parallel/policy/env_config.jl:73-73`
<!-- vulcan:connections:end -->

## Limitations
Distinct inputs can collide after sanitisation (`"my host"` and `"my_host"` both become `my_host`), so two machines could share a state file. Leading dots survive, allowing tokens like `..` fragments; the regex does not guard against path traversal beyond removing separators. Non-ASCII letters are replaced rather than transliterated.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 67.
