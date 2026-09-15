---
id: analysis.types__parse_scenario_list
label: _parse_scenario_list
kind: function
source:
  file: src/analysis/verification/telemetry_verification/types.jl
  symbol: _parse_scenario_list
  lines:
  - 245
  - 245
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
  type: Vector{String}
  units: n/a
  description: Return value of `_parse_scenario_list`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# _parse_scenario_list

## Purpose
Turns a comma-separated scenario selection string such as `"odyssey,vex"` into a normalised `Vector{String}` of scenario names for the telemetry verification CLI, where an empty result means every scenario in the manifest.

## Design & Implementation
An `@inline` one-liner: `split(String(raw), ',')` tokenises on commas, each token is `strip`ped, tokens that are empty after stripping are discarded, and the survivors are `lowercase`d into a typed `String[]` comprehension. The return annotation `::Vector{String}` fixes the element type even for an empty result. Accepting `AbstractString` lets `SubString` inputs from environment parsing pass without conversion.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | AbstractString | n/a | yes | Positional argument `raw`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{String} | n/a | — | Return value of `_parse_scenario_list`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[envana.ana_manifest_parsing_parse_cli|parse_cli]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:719-719`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the comma is a separator; semicolons or whitespace-separated lists are treated as a single token. Lower-casing uses Julia's Unicode-aware `lowercase`, so names must be compared with the same normalisation elsewhere. Duplicates are preserved, not de-duplicated. No validation is done against the manifest here, so unknown names surface only later when scenarios are filtered.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/types.jl` line 245.
