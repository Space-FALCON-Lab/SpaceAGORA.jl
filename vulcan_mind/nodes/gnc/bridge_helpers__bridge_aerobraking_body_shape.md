---
id: gnc.bridge_helpers__bridge_aerobraking_body_shape
label: _bridge_aerobraking_body_shape
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _bridge_aerobraking_body_shape
  lines:
  - 73
  - 73
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: mission
  type: Any
  units: n/a
  required: false
  description: Positional argument `mission` (default `nothing`).
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
  description: Return value of `_bridge_aerobraking_body_shape`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _bridge_aerobraking_body_shape

## Purpose
Resolves the vehicle body-shape label used by the aerobraking aerodynamic and heating models, defaulting to `"Spacecraft"`.

## Design & Implementation
Calls `_bridge_optional_field(args, :body_shape, "Spacecraft")` and forces the result through `String(value)`, so `Symbol` or `SubString` inputs are normalised to a `String`. The `mission` argument is accepted but unused. Return type is annotated `::String`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `mission` | Any | n/a | no | Positional argument `mission` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_bridge_aerobraking_body_shape`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.bridge_helpers__make_aerobraking_runtime_settings|_make_aerobraking_runtime_settings]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:150-150`

**Downstream**

- `callees` → [[gnc.bridge_helpers__bridge_optional_field|_bridge_optional_field]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:74-74`
<!-- vulcan:connections:end -->

## Limitations
`String(value)` throws `MethodError` when the stored field is a number or `nothing`. No validation is performed against the set of shapes the downstream aerodynamic coefficient tables actually support, so a typo such as `"Spacecarft"` is only detected later when a lookup fails. The `mission` object is never consulted even if it carries a body definition.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 73.
