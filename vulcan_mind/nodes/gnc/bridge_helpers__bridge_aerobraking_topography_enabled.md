---
id: gnc.bridge_helpers__bridge_aerobraking_topography_enabled
label: _bridge_aerobraking_topography_enabled
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _bridge_aerobraking_topography_enabled
  lines:
  - 48
  - 48
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  type: Bool
  units: n/a
  description: Return value of `_bridge_aerobraking_topography_enabled`.
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

# _bridge_aerobraking_topography_enabled

## Purpose
Determines whether the aerobraking runtime should use a topographic (spherical-harmonic) planet surface instead of a spherical one, honouring the typed `environment_model.topography` flag first and a legacy string field second.

## Design & Implementation
If `args.environment_model` exists and has a `:topography` property, its value is converted with `Bool(...)` and returned. Otherwise `_bridge_optional_field(args, :topography_model, "")` is read and compared for exact equality with the string `"Spherical Harmonics"`; any other string, including the empty default, yields `false`. The return type is annotated `::Bool`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_bridge_aerobraking_topography_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.bridge_helpers__make_aerobraking_runtime_settings|_make_aerobraking_runtime_settings]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:147-147`

**Downstream**

- `callees` → [[gnc.bridge_helpers__bridge_optional_field|_bridge_optional_field]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:52-52`
<!-- vulcan:connections:end -->

## Limitations
The legacy comparison is case- and whitespace-sensitive, so `"spherical harmonics"` silently disables topography. `Bool(...)` throws `InexactError` for numeric flags other than 0 or 1. Presence of `environment_model` with a `topography` field takes precedence even if a legacy `topography_model` string is also present, with no warning about the conflict.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 48.
