---
id: core.no_gram_presets_nogrampresets
label: NoGramPresets
kind: module
source:
  file: src/core/state/no_gram_presets.jl
  symbol: NoGramPresets
  lines:
  - 1
  - 1
inputs:
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
  description: Value produced by this symbol.
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

# NoGramPresets

## Purpose
The module that packages the no-GRAM onboarding presets, giving a new user a runnable environment configuration without GRAM data files or SPICE kernels on disk.

## Design & Implementation
Imports the abstract planet and density types, the three concrete planets, the no-atmosphere and exponential atmosphere models, the simple analytic ephemerides, the Maxwellian heat model and the `EnvironmentModel` constructor, and exports exactly three functions: `make_no_gram_planet`, `make_no_gram_density_model` and `make_no_gram_environment`. Keeping the preset surface this narrow means the onboarding path is a documented three-function contract rather than a set of struct literals a newcomer must reconstruct.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/state/no_gram_presets.jl`

**Downstream**

- `callees` → [[core.no_gram_presets_make_no_gram_planet|make_no_gram_planet]] · `callers` · call · `src/core/state/no_gram_presets.jl:13-13`
<!-- vulcan:connections:end -->

## Limitations
The presets pin one ephemerides implementation and one default thermal model, so a run that needs kernel-accurate third-body positions has to leave this module entirely rather than override a single field.

## Provenance
Mapped from `src/core/state/no_gram_presets.jl` line 1.
