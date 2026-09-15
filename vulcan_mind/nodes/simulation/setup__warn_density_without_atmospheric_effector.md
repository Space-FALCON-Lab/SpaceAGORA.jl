---
id: simulation.setup__warn_density_without_atmospheric_effector
label: _warn_density_without_atmospheric_effector
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _warn_density_without_atmospheric_effector
  lines:
  - 87
  - 87
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
  type: Nothing
  units: n/a
  description: Return value of `_warn_density_without_atmospheric_effector`. Returns
    `nothing`.
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

# _warn_density_without_atmospheric_effector

## Purpose
Emits a single loud warning when a density model other than `NoAtmosphereModel` is configured but no effector will ever apply aerodynamic force, a configuration that has historically produced silently drag-free studies.

## Design & Implementation
Returns early if `_density_without_aero_warning_enabled()` is false, if `args.environment_model.density_model isa SimulationModel.EnvironmentModels.NoAtmosphereModel`, or if `_any_effector_consumes_atmosphere(args.dynamics_model.dynamic_effectors)` is true. Otherwise it calls `@warn` with a multi-line message naming the density model type and explaining both the built-in and `environment_requirements(model).atmosphere = true` routes, using `maxlog = 1` so repeated setups do not spam. Always returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_warn_density_without_atmospheric_effector`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:167-167`

**Downstream**

- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/setup.jl:100-100`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/setup.jl:100-100`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/setup.jl:100-100`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/setup.jl:100-100`
- `callees` → [[simulation.setup__any_effector_consumes_atmosphere|_any_effector_consumes_atmosphere]] · `callers` · call · `src/simulation/engine/setup.jl:93-93`
- `callees` → [[simulation.setup__density_without_aero_warning_enabled|_density_without_aero_warning_enabled]] · `callers` · call · `src/simulation/engine/setup.jl:88-88`
<!-- vulcan:connections:end -->

## Limitations
This is a warning, not an error, so an unintended configuration proceeds. `maxlog = 1` is per call-site for the whole Julia session, so a second scenario in the same process with the same mistake goes unreported. The check is type-based and cannot detect an aerodynamic effector whose reference area is zero.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 87.
