---
id: simulation.setup__any_effector_consumes_atmosphere
label: _any_effector_consumes_atmosphere
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _any_effector_consumes_atmosphere
  lines:
  - 72
  - 72
inputs:
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
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
  description: Return value of `_any_effector_consumes_atmosphere`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _any_effector_consumes_atmosphere

## Purpose
Determines whether at least one dynamic effector in the model will read atmospheric density to produce force, covering both built-in aerodynamic models and user-defined effectors that declare the requirement through the public extension hook.

## Design & Implementation
Takes `dynamic_effectors::Tuple`. First delegates to `SimulationModel.SimulationCallbacks._uses_atmospheric_dynamic_effector(dynamic_effectors)`, which recognises the built-in `AerodynamicCoefficient*` types; a `true` short-circuits. Otherwise an `@inbounds` loop calls `SimulationModel.environment_requirements(effector).atmosphere` on each element and returns `true` on the first hit. Returns `false` when neither path matches. Pure and allocation-free for tuple inputs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_any_effector_consumes_atmosphere`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__warn_density_without_atmospheric_effector|_warn_density_without_atmospheric_effector]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:93-93`

**Downstream**

- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/setup.jl:75-75`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/setup.jl:75-75`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/setup.jl:75-75`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/setup.jl:75-75`
- `callees` → [[simulation.assembly__uses_atmospheric_dynamic_effector|_uses_atmospheric_dynamic_effector]] · `callers` · call · `src/simulation/engine/setup.jl:73-73`
<!-- vulcan:connections:end -->

## Limitations
Relies on every custom effector having an `environment_requirements` method returning a struct with an `atmosphere` field; a missing method raises `MethodError` here rather than a clear configuration error. Effectors with `A = 0` or otherwise disabled still count as consumers if their type is recognised.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 72.
