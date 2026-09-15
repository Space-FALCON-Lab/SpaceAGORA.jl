---
id: parallel.outer_route_selection__route_density_bucket
label: _route_density_bucket
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _route_density_bucket
  lines:
  - 60
  - 60
inputs:
- id: family
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `family`.
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
  description: Return value of `_route_density_bucket`.
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

# _route_density_bucket

## Purpose
Normalises the atmosphere density model family name into a short canonical token for the `dens=` signature field and for the native-GRAM detection in `_is_native_gram_point_density`.

## Design & Implementation
`@inline _route_density_bucket(family::AbstractString)::String` lowercases and strips the input, then maps aliases: ("none", "vacuum", "noatmosphere") -> "none"; ("gram_point", "gram") -> "gram_pt"; ("gram_surrogate", "gram_offline_surrogate") -> "gram_srg"; ("exponential", "exp") -> "exp"; ("polynomialfit", "polyfit", "poly") -> "poly"; ("nrlmsise00", "nrl") -> "nrl"; empty -> "unknown". Any other token is returned with `|` replaced by `_` so it cannot corrupt the pipe-delimited signature.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `family` | AbstractString | n/a | yes | Positional argument `family`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_route_density_bucket`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__is_native_gram_point_density|_is_native_gram_point_density]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:269-269`
- [[parallel.outer_route_selection__outer_route_signature_hierarchy|_outer_route_signature_hierarchy]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:194-194`
- [[parallel.outer_route_selection_outer_route_signature|outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:166-166`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Unknown family names pass through verbatim, so a typo creates a fresh signature bucket with no history rather than raising. Only `|` is sanitised; `=` or whitespace inside an unusual token would still survive into the signature string. Allocates a new string on every call.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 60.
