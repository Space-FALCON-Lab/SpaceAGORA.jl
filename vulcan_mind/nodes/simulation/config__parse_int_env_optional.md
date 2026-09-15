---
id: simulation.config__parse_int_env_optional
label: _parse_int_env_optional
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/config.jl
  symbol: _parse_int_env_optional
  lines:
  - 36
  - 36
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
  description: Return value of `_parse_int_env_optional`.
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

# _parse_int_env_optional

## Purpose
`_parse_int_env_optional(name::String)::Union{Nothing, Int}` is the integer counterpart used for the GRAM track-cache sample-count knobs — `SPACEAGORA_GRAM_TRACK_CACHE_NPOS`, its entry and orbit variants, and the legacy `SPACEAGORA_GRAM_SEGMENT_CACHE_POINTS`. Returning `nothing` for an unset variable lets the configuration builder cascade from the most specific override to the least before landing on the literal defaults of 16 entry points and 48 orbit points.

## Design & Implementation
Structurally identical to the float version: `haskey(ENV, name) || return nothing`, then `strip(ENV[name])` and `parse(Int, raw)` inside a `try`, with failures converted to `ArgumentError("$name must be an integer value, got '$raw'")`. `@inline` plus the concrete `Union{Nothing, Int}` return annotation keeps the result unboxed. `parse(Int, ...)` accepts an optional sign and rejects any decimal point, so `"48.0"` is a hard error rather than a truncation — deliberate, because a fractional sample count is always a mistake.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_parse_int_env_optional`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.gram_cache_config_gram_track_cache_config|_gram_track_cache_config]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/config.jl:76-76`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No range check is applied here; negative or absurdly large counts pass through and are only clamped downstream by `max(2, ...)` in `_gram_track_cache_config`, which means a user who sets a count of one or of a billion gets either a silent correction or an out-of-memory refresh rather than a diagnostic. Values exceeding `typemax(Int)` raise a parse error reported as a formatting problem. The blanket `catch` again converts unrelated exceptions into a misleading message.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/config.jl` line 36.
