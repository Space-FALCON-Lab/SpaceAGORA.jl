---
id: simulation.config__parse_float_env_optional
label: _parse_float_env_optional
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/config.jl
  symbol: _parse_float_env_optional
  lines:
  - 25
  - 25
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
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
  description: Return value of `_parse_float_env_optional`.
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

# _parse_float_env_optional

## Purpose
`_parse_float_env_optional(name::String)::Union{Nothing, Float64}` distinguishes an absent environment variable from one that is present but set to zero. The GRAM track-cache configuration builder relies on that distinction to layer regime-specific overrides on top of backward-compatible global knobs via chained `something(...)` calls, where `nothing` means 'fall through to the next source'.

## Design & Implementation
It short-circuits with `haskey(ENV, name) || return nothing`, so the absent case allocates nothing and never parses. When present, `strip(ENV[name])` trims whitespace and `parse(Float64, raw)` runs under a `try`; on failure it throws `ArgumentError("$name must be a floating-point value, got '$raw'")`. The declared return type is the small `Union{Nothing, Float64}`, which Julia stores without boxing, so the layered `something(something(regime, compat), literal)` idiom in `_gram_track_cache_config` stays allocation-free.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_parse_float_env_optional`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.gram_cache_config_gram_track_cache_config|_gram_track_cache_config]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/config.jl:73-73`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Presence is judged by `haskey` alone, so an exported-but-empty variable (`FOO=`) is treated as present and then fails to parse, surfacing as a configuration error where the user almost certainly meant 'unset'. Failure is thrown rather than returned, so a caller cannot distinguish malformed from absent without a `try` of its own. As with every `ENV` read here, there is a race window if another thread mutates the environment between the `haskey` check and the `ENV[name]` lookup, which would raise a `KeyError` instead of returning `nothing`.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/config.jl` line 25.
