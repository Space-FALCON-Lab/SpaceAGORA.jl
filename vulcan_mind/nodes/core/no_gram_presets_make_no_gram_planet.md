---
id: core.no_gram_presets_make_no_gram_planet
label: make_no_gram_planet
kind: function
source:
  file: src/core/state/no_gram_presets.jl
  symbol: make_no_gram_planet
  lines:
  - 18
  - 18
inputs:
- id: planet
  type: AbstractPlanet
  units: n/a
  required: true
  description: Positional argument `planet`.
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
  description: Return value of `make_no_gram_planet`. Returns `planet`.
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

# make_no_gram_planet

## Purpose
Resolves a caller's loose planet designation into a concrete planet object for the no-GRAM onboarding mode, without touching SPICE kernels or GRAM assets.

## Design & Implementation
Three methods of one generic. Given an `AbstractPlanet` it is the identity, marked `@inline`, so an already-constructed planet passes through untouched. Given a `Symbol` it lowercases the name and returns `Earth()`, `Mars()` or `Venus()`, raising `ArgumentError` listing the supported keys for anything else. Given an `AbstractString` it strips whitespace, lowercases, and delegates to the symbol method. The identity method is what lets configuration code accept either a preset key or a fully specified planet at the same argument position.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | AbstractPlanet | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `make_no_gram_planet`. Returns `planet`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.no_gram_presets_nogrampresets|NoGramPresets]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:13-13`
- [[parcore.no_gram_presets_make_no_gram_environment|make_no_gram_environment]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:82-82`

**Downstream**

- `callees` → [[core.no_gram_presets_make_no_gram_density_model|make_no_gram_density_model]] · `callers` · feedback · `src/core/state/no_gram_presets.jl:37-37`
- `callees` → [[environment.density_models_exponentialatmospheremodel|ExponentialAtmosphereModel]] · `callers` · call · `src/core/state/no_gram_presets.jl:42-42`
- `callees` → [[environment.density_models_noatmospheremodel|NoAtmosphereModel]] · `callers` · call · `src/core/state/no_gram_presets.jl:41-41`
- `callees` → [[environment.density_models_piecewiseexponentialatmospheremodel|PiecewiseExponentialAtmosphereModel]] · `callers` · call · `src/core/state/no_gram_presets.jl:44-44`
- `callees` → [[environment.planets_earth|Earth]] · `callers` · call · `src/core/state/no_gram_presets.jl:25-25`
- `callees` → [[environment.planets_mars|Mars]] · `callers` · call · `src/core/state/no_gram_presets.jl:27-27`
- `callees` → [[environment.planets_venus|Venus]] · `callers` · call · `src/core/state/no_gram_presets.jl:29-29`
<!-- vulcan:connections:end -->

## Limitations
Only three planets are recognised, and the supported set is a literal branch rather than a registry, so extending it means editing this function; the string method strips and lowercases but does not accept common aliases such as `terra`.

## Provenance
Mapped from `src/core/state/no_gram_presets.jl` line 18.
