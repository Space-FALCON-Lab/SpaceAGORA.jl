---
id: simulation.config__parse_float_env
label: _parse_float_env
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/config.jl
  symbol: _parse_float_env
  lines:
  - 15
  - 15
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: Float64
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
  type: Float64
  units: n/a
  description: Return value of `_parse_float_env`.
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

# _parse_float_env

## Purpose
`_parse_float_env(name::String, default::Float64)::Float64` reads a required-with-default floating-point tuning knob out of the process environment for the GRAM track-cache configuration. It exists so that cache tolerances such as `SPACEAGORA_GRAM_SEGMENT_CACHE_TRANSITION_BAND_M` can be overridden per run without editing code, while a malformed value is reported as a clear configuration error rather than silently falling back.

## Design & Implementation
The body is three statements. `raw = strip(get(ENV, name, string(default)))` fetches the variable, substituting the stringified default when it is absent, and trims surrounding whitespace so a trailing newline from a shell export does not break parsing. `parse(Float64, raw)` runs inside a `try`; any parse failure is caught and rethrown as `ArgumentError("$name must be a floating-point value, got '$raw'")`, which names both the offending variable and the exact trimmed text. The function is `@inline` and return-type-annotated `::Float64`, keeping it free of dynamic dispatch at its call sites.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | Float64 | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_parse_float_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__gram_entry_target_cd|_gram_entry_target_cd]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:34-34`
- [[simulation.config__gram_entry_target_dt|_gram_entry_target_dt]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:35-35`
- [[simulation.vacuum_predicted_gram__vacuum_gram_cache_npoints|_vacuum_gram_cache_npoints]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:33-33`
- [[simulation_a.gram_cache_config_gram_track_cache_config|_gram_track_cache_config]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/config.jl:143-143`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Round-tripping the default through `string()` and back through `parse` is wasteful and, for values not exactly representable in the shortest decimal form, is only correct because Julia's `string(::Float64)` round-trips; a caller passing a default computed at run time still pays a formatting and parsing cost. The `catch` swallows every exception type, so an interrupt raised while parsing is misreported as a bad value. `ENV` is process-global and unsynchronised, so a concurrent write from another task can be observed mid-read. An empty-string override parses as an error rather than being treated as unset.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/config.jl` line 15.
