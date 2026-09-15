---
id: simulation.constellation_ensemble__ensemble_member_configuration
label: _ensemble_member_configuration
kind: function
source:
  file: src/simulation/campaigns/constellation_ensemble.jl
  symbol: _ensemble_member_configuration
  lines:
  - 32
  - 32
inputs:
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: spacecraft
  type: SpacecraftModel
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: member_tag
  type: String
  units: n/a
  required: true
  description: Positional argument `member_tag`.
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
  type: SimulationConfiguration
  units: n/a
  description: Return value of `_ensemble_member_configuration`.
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

# _ensemble_member_configuration

## Purpose
Builds the single-satellite `SimulationConfiguration` that one ensemble worker will hand to `run_simulation`. It is the point where a multi-spacecraft configuration is reduced to a configuration containing exactly one `SpacecraftModel` while keeping every other model (environment, guidance, navigation, control) shared by reference.

## Design & Implementation
Calls the `SimulationConfiguration` keyword constructor, forwarding `file_paths`, `mission_configuration`, `environment_model`, `guidance_model`, `navigation_model`, `control_model`, `initial_time`, `integration_tolerances` and `solver_config` unchanged from `args`. Two fields differ: `simulation_settings` is replaced by `_ensemble_member_settings(args.simulation_settings, member_tag)` to give the member private output directories, and `dynamics_model` becomes `DynamicsModel([spacecraft], args.dynamics_model.dynamic_effectors)`, a one-element spacecraft vector paired with the original dynamic effector list. No copying occurs here; isolation between concurrent workers is achieved later by `deepcopy` inside the worker task.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `spacecraft` | SpacecraftModel | n/a | yes | Positional argument `spacecraft`. |
| in | `member_tag` | String | n/a | yes | Positional argument `member_tag`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationConfiguration | n/a | — | Return value of `_ensemble_member_configuration`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.campaigns_constellation_ensemble_run_constellation_ensemble|run_constellation_ensemble]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:145-145`

**Downstream**

- `callees` → [[parcore.simulation_configuration_simulationconfiguration|SimulationConfiguration]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:37-37`
- `callees` → [[simulation.constellation_ensemble__ensemble_member_settings|_ensemble_member_settings]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:39-39`
- `callees` → [[vehicle.model_dynamicsmodel|DynamicsModel]] · `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:42-42`
<!-- vulcan:connections:end -->

## Limitations
Because the shared models are passed by reference, the function is only safe under concurrency when the caller deep-copies the result; it does not do so itself. `dynamic_effectors` are reused as-is, so an effector that internally indexes spacecraft by position in the original constellation will see index 1 for every member. Any `SimulationConfiguration` field added in future must be forwarded explicitly or it reverts to the constructor default.

## Provenance
Mapped from `src/simulation/campaigns/constellation_ensemble.jl` line 32.
