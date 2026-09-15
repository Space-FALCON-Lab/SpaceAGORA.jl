---
id: parallel.profile_definitions__normalize_profile_token
label: _normalize_profile_token
kind: function
source:
  file: src/parallel/routing/profile_definitions.jl
  symbol: _normalize_profile_token
  lines:
  - 67
  - 67
inputs:
- id: raw
  type: AbstractString
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
  description: Return value of `_normalize_profile_token`.
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

# _normalize_profile_token

## Purpose
Reduces a user-written profile string to the one spelling the parser's alias tables are keyed on, so casing, padding and separator style do not matter at the command line.

## Design & Implementation
Converts to `String`, strips leading and trailing whitespace, lowercases, rewrites every hyphen to an underscore, and finally removes interior spaces outright. The order matters: stripping before the space removal means the interior-space pass only has to deal with separators a user typed inside a name, such as `outer only`. Declared `@inline` with a `::String` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | AbstractString | n/a | yes | Positional argument `raw`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_normalize_profile_token`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/profile_definitions.jl`
- [[parallel.profile_definitions_parse_parallel_profile|parse_parallel_profile]] · `callees` → `callers` · call · `src/parallel/routing/profile_definitions.jl:80-80`

**Downstream**

- `callees` → [[parallel.profile_definitions_parse_parallel_profile|parse_parallel_profile]] · `callers` · feedback · `src/parallel/routing/profile_definitions.jl:74-74`
<!-- vulcan:connections:end -->

## Limitations
Only hyphens and spaces are folded, so a token written with a dot or a slash separator survives normalisation unchanged and then fails the alias match with an unhelpfully literal error message.

## Provenance
Mapped from `src/parallel/routing/profile_definitions.jl` line 67.
