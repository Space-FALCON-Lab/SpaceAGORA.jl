---
id: analysis.scenario_builders__gramofflinesurrogatefallbackbase
label: _GRAMOfflineSurrogateFallbackBase
kind: struct
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _GRAMOfflineSurrogateFallbackBase
  lines:
  - 273
  - 273
inputs:
- id: planet_name
  type: String
  units: n/a
  required: true
  description: Field `planet_name`.
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
  type: _GRAMOfflineSurrogateFallbackBase
  units: n/a
  description: Constructed `_GRAMOfflineSurrogateFallbackBase` (keyword constructor
    via @kwdef).
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

# _GRAMOfflineSurrogateFallbackBase

## Purpose
A stand-in base model that lets the GRAM offline surrogate run when the native GRAM shared library is absent.

## Design & Implementation
A `Base.@kwdef struct` holding only `planet_name`. It implements `_gram_point_density` by returning zero density, 200 K and zero wind, so any surrogate query that would fall back to the native point model instead sees vacuum. It is constructed by `_try_libraryless_gram_surrogate` and wrapped in `GRAMAtmosphereModelSurrogate` with the default surrogate file.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Field `planet_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | _GRAMOfflineSurrogateFallbackBase | n/a | — | Constructed `_GRAMOfflineSurrogateFallbackBase` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__try_libraryless_gram_surrogate|_try_libraryless_gram_surrogate]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:312-312`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[environment.density_models__gram_point_density|_gram_point_density]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:277-277`
<!-- vulcan:connections:end -->

## Limitations
Below the surrogate's point-fallback altitude the density silently becomes zero rather than failing, so a low-periapsis scenario under this fallback flies a physically wrong atmosphere while reporting success.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 273.
