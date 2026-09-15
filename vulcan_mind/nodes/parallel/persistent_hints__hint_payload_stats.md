---
id: parallel.persistent_hints__hint_payload_stats
label: _hint_payload_stats
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _hint_payload_stats
  lines:
  - 19
  - 19
inputs:
- id: payload
  type: Any
  units: n/a
  required: true
  description: Positional argument `payload`.
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_hint_payload_stats`.
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

# _hint_payload_stats

## Purpose
Parses one stats dictionary from the hint file back into `AdaptiveChoiceStats`, tolerating missing, mistyped or inconsistent fields rather than failing the whole load.

## Design & Implementation
Returns `nothing` for any non-dictionary payload. Each of the five fields is read with `get` and converted inside its own `try` block, falling back to zero on any conversion error, and clamped at zero. A record with no samples is discarded. Successes are capped at `samples` and failures at `samples - successes`, so the three counts are always mutually consistent whatever the file said.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `payload` | Any | n/a | yes | Positional argument `payload`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_hint_payload_stats`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- [[parallel.persistent_hints__load_persistent_hint_state_locked_bang|_load_persistent_hint_state_locked!]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:94-94`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:37-37`
- `callees` → [[parallel.types_adaptivechoicestats|AdaptiveChoiceStats]] · `callers` · call · `src/parallel/policy/persistent_hints.jl:49-49`
<!-- vulcan:connections:end -->

## Limitations
A partially corrupt record is silently repaired to zeros rather than reported, so a hint file damaged in the field loses history without any diagnostic to the user.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 19.
