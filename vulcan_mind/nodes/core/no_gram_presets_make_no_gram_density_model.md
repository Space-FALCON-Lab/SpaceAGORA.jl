---
id: core.no_gram_presets_make_no_gram_density_model
label: make_no_gram_density_model
kind: function
source:
  file: src/core/state/no_gram_presets.jl
  symbol: make_no_gram_density_model
  lines:
  - 46
  - 46
inputs:
- id: planet
  type: AbstractPlanet
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: density_model
  type: AbstractDensityModel
  units: n/a
  required: true
  description: Positional argument `density_model`.
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
  type: Any
  units: n/a
  description: Return value of `make_no_gram_density_model`. Returns `density_model`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# make_no_gram_density_model

## Purpose
Resolves an atmosphere designation into a concrete density model for the no-GRAM baseline, so a quickstart run needs no GRAM data files.

## Design & Implementation
Mirrors the planet resolver in shape. An `AbstractDensityModel` argument passes through as the identity so a hand-built model such as a piecewise exponential can be supplied directly. A `Symbol` maps `:none` to `NoAtmosphereModel()` and `:exponential` to `ExponentialAtmosphereModel(planet)`, the latter taking the resolved planet so the scale height matches. Anything else raises `ArgumentError` naming both the bad value and the two supported options. The string method strips, lowercases and delegates.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | AbstractPlanet | n/a | yes | Positional argument `planet`. |
| in | `density_model` | AbstractDensityModel | n/a | yes | Positional argument `density_model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `make_no_gram_density_model`. Returns `density_model`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.no_gram_presets_make_no_gram_planet|make_no_gram_planet]] · `callees` → `callers` · feedback · `src/core/state/no_gram_presets.jl:37-37`
- [[parcore.no_gram_presets_make_no_gram_environment|make_no_gram_environment]] · `callees` → `callers` · feedback · `src/core/state/no_gram_presets.jl:83-83`

**Downstream**

- `callees` → [[envana.env_simple_ephemerides_simpleephemeridesmodel|SimpleEphemeridesModel]] · `callers` · call · `src/core/state/no_gram_presets.jl:68-68`
- `callees` → [[environment.density_models_exponentialatmospheremodel|ExponentialAtmosphereModel]] · `callers` · call · `src/core/state/no_gram_presets.jl:55-55`
- `callees` → [[environment.density_models_noatmospheremodel|NoAtmosphereModel]] · `callers` · call · `src/core/state/no_gram_presets.jl:53-53`
- `callees` → [[parcore.no_gram_presets_make_no_gram_environment|make_no_gram_environment]] · `callers` · call · `src/core/state/no_gram_presets.jl:64-64`
<!-- vulcan:connections:end -->

## Limitations
The exponential model is constructed from the planet alone, so it carries that planet's single nominal scale height and cannot represent the diurnal or dust-driven density variation the GRAM path provides.

## Provenance
Mapped from `src/core/state/no_gram_presets.jl` line 46.
