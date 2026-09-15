---
id: simulation.config__gram_track_cache_mode
label: _gram_track_cache_mode
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/config.jl
  symbol: _gram_track_cache_mode
  lines:
  - 50
  - 50
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
  type: Symbol
  units: n/a
  description: Return value of `_gram_track_cache_mode`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _gram_track_cache_mode

## Purpose
`_gram_track_cache_mode()::Symbol` resolves the master on/off/auto switch for the GRAM ground-track cache from the environment. The cache is off by default for a measured reason recorded in the source: on a representative entry case, point-to-point GRAM queries cost about 0.45 s while the cache plus Allen-Eggers path cost about 42.8 s, so enabling it is an explicit opt-in rather than a default.

## Design & Implementation
It prefers `SPACEAGORA_GRAM_TRACK_CACHE` and falls back to the older `SPACEAGORA_GRAM_SEGMENT_CACHE`, defaulting to the string `"off"` when neither is set; the two-name lookup preserves compatibility with scripts written before the rename. The raw value is normalised with `lowercase(strip(raw))` and then matched against three sets of accepted spellings: `("off", "none", "0", "false", "no")` yields `:off`, `("on", "1", "true", "yes")` yields `:on`, and the literal `"auto"` yields `:auto`. Anything else throws `ArgumentError("Unsupported GRAM track-cache mode '$mode'. Use one of: off, auto, on.")`, so a typo disables nothing silently.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_gram_track_cache_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.gram_cache_config_gram_track_cache_config|_gram_track_cache_config]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/config.jl:70-70`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Precedence is by variable name only: if both `SPACEAGORA_GRAM_TRACK_CACHE` and `SPACEAGORA_GRAM_SEGMENT_CACHE` are exported, the older one is ignored with no warning, which is surprising for a user who set only the legacy name in a wrapper script and the new name in a CI job. `"auto"` has no accepted synonyms while `:on` and `:off` have five each, an asymmetry that makes `"automatic"` an error. The error message quotes the already-lowercased, stripped `mode` rather than the user's original text, hiding stray characters that whitespace trimming removed.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/config.jl` line 50.
