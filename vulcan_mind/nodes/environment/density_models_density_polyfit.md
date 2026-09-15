---
id: environment.density_models_density_polyfit
label: density_polyfit
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: density_polyfit
  lines:
  - 1081
  - 1081
inputs:
- id: h
  type: Float64
  units: n/a
  required: true
  description: Positional argument `h`.
- id: p
  type: params
  units: n/a
  required: true
  description: Positional argument `p`.
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
  type: Tuple{Float64,
  units: n/a
  description: 'Return value of `density_polyfit`. Type parameters: `params`.'
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# density_polyfit

## Purpose
Convenience evaluator building the planet's polynomial model on the fly, used by the GRAM `getDensity` methods above the entry interface on non-Keplerian runs.

## Design & Implementation
Constructs `PolynomialFitAtmosphereModel(planet)` from the run configuration and calls its seven-argument `getDensity` with zero latitude, longitude and time.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `h` | Float64 | n/a | yes | Positional argument `h`. |
| in | `p` | params | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `density_polyfit`. Type parameters: `params`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[ext.spaceagoragramsuiteext__gram_spice_ephemeris_state|_gram_spice_ephemeris_state]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:288-288`
- [[gnc.closed_form_solution_closed_form_calculation|closed_form_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:191-191`
- [[gnc.trajectory_predictor_closed_form_targeting|closed_form_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:424-424`
- [[gncx.energy_profile_solver_security_mode|security_mode]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/energy_profile_solver.jl:11-11`
- [[gncx.heat_rate_models_lambdas|lambdas]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl:23-23`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:191-191`
- [[simulation.model_selection__gram_isolated_pool_density_state|_gram_isolated_pool_density_state]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:101-101`

**Downstream**

- `callees` → [[environment.density_models_polynomialfitatmospheremodel|PolynomialFitAtmosphereModel]] · `callers` · call · `src/environment/atmosphere/density_models.jl:1082-1082`
- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `src/environment/atmosphere/density_models.jl:1083-1083`
<!-- vulcan:connections:end -->

## Limitations
It rebuilds the model — including `collect`ing the coefficient vector — on every call, which is allocation on the density hot path; the coefficients could be cached once per run.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 1081.
