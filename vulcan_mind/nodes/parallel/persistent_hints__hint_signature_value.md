---
id: parallel.persistent_hints__hint_signature_value
label: _hint_signature_value
kind: function
source:
  file: src/parallel/policy/persistent_hints.jl
  symbol: _hint_signature_value
  lines:
  - 371
  - 371
inputs:
- id: signature
  type: String
  units: n/a
  required: true
  description: Positional argument `signature`.
- id: key
  type: String
  units: n/a
  required: true
  description: Positional argument `key`.
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
  description: Return value of `_hint_signature_value`.
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

# _hint_signature_value

## Purpose
Extracts one field's value from a serialised signature string, used when aggregating history by profile, machine or source.

## Design & Implementation
Splits the signature on `|`, finds the first token starting with `key=`, and returns the remainder after the equals sign — an empty string if there is nothing after it — or `nothing` if no token matches. String indexing uses `lastindex` so multi-byte characters in the key do not break the slice.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `signature` | String | n/a | yes | Positional argument `signature`. |
| in | `key` | String | n/a | yes | Positional argument `key`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_hint_signature_value`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- [[parallel.persistent_hints_hint_layer_stats_snapshot|hint_layer_stats_snapshot]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:404-404`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It matches by prefix only, so a key that is a prefix of another key, such as `heavy` against `heavy_only`, returns the wrong token if the shorter one appears later in the string; the current signature layout avoids that ordering but the function does not guard against it.

## Provenance
Mapped from `src/parallel/policy/persistent_hints.jl` line 371.
