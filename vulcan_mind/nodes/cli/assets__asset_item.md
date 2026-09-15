---
id: cli.assets__asset_item
label: _asset_item
kind: function
source:
  file: src/cli/assets.jl
  symbol: _asset_item
  lines:
  - 37
  - 37
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: scope
  type: String
  units: n/a
  required: true
  description: Positional argument `scope`.
- id: path
  type: String
  units: n/a
  required: true
  description: Positional argument `path`.
- id: available
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `available`.
- id: required
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `required`.
- id: detail
  type: String
  units: n/a
  required: false
  description: Keyword argument `detail` (default `""`).
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
  type: AssetCheckItem
  units: n/a
  description: Return value of `_asset_item`. Returns `AssetCheckItem(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- cli
charts:
- cli
origin: agent
---

# _asset_item

## Purpose
Positional-plus-keyword convenience constructor for `AssetCheckItem`, used by `check_assets` so that the three identifying strings can be passed positionally while the status flags remain explicit keywords. It exists purely to keep the call site in `check_assets` readable.

## Design & Implementation
Declared `@inline` with signature `_asset_item(name::String, scope::String, path::String; available::Bool, required::Bool, detail::String="")`. The body forwards every argument to the `Base.@kwdef` keyword constructor `AssetCheckItem(name=..., scope=..., path=..., available=..., required=..., detail=...)` and returns the new struct. `available` and `required` have no defaults, so a caller must state both; `detail` defaults to the empty string exactly as in the struct definition. No validation or transformation of the inputs occurs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `scope` | String | n/a | yes | Positional argument `scope`. |
| in | `path` | String | n/a | yes | Positional argument `path`. |
| in | `available` | Bool | n/a | yes | Keyword argument `available`. |
| in | `required` | Bool | n/a | yes | Keyword argument `required`. |
| in | `detail` | String | n/a | no | Keyword argument `detail` (default `""`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AssetCheckItem | n/a | — | Return value of `_asset_item`. Returns `AssetCheckItem(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[misc.assets_check_assets|check_assets]] · `callees` → `callers` · feedback · `src/cli/assets.jl:88-88`

**Downstream**

- `callees` → [[cli.assets_assetcheckitem|AssetCheckItem]] · `callers` · call · `src/cli/assets.jl:38-38`
<!-- vulcan:connections:end -->

## Limitations
Because it adds no behaviour over the keyword constructor, it duplicates the field list and must be updated in lockstep if `AssetCheckItem` gains a field. Strict `String` typing rejects `SubString` or `AbstractString` inputs, which can surprise callers that slice paths. The function does not normalise or `abspath` the `path` argument, so relative paths pass through unchanged.

## Provenance
Mapped from `src/cli/assets.jl` line 37.
