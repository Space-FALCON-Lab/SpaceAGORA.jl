---
id: spaceagora.spaceagora_spaceagora
label: SpaceAGORA
kind: module
source:
  file: src/SpaceAGORA.jl
  symbol: SpaceAGORA
  lines:
  - 3
  - 3
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
- spaceagora
charts:
- spaceagora
origin: agent
---

# SpaceAGORA

## Purpose
`SpaceAGORA` is the top-level Julia package module. It stitches together the parallel-routing, process-pool, runtime-services, core model, simulation engine, campaign, telemetry-verification, RPO asset, RPO visualization and CLI submodules into a single namespace, re-exports their public surface, and defines a small set of thin forwarding entrypoints (`run_simulation`, `prewarm_nbody_ephemeris_cache`, `load_nbody_ephemeris_cache!`, `check_assets`, `render_asset_report`, `run_cli`).

## Design & Implementation
The module is declared `__precompile__(true)` and imports `@compile_workload`/`@setup_workload` from PrecompileTools. Ten `include` calls (in dependency order, from `parallel/routing/parallel_profiles.jl` through `cli/spaceagora_cli.jl`) load the nested modules; ordering matters because later files `using` earlier ones. A long block of `using .Submodule: names` statements pulls specific symbols into scope, followed by matching `export` lists so downstream code can write `using SpaceAGORA`. The forwarding methods are one-line `f(args...; kwargs...) = Submodule.f(args...; kwargs...)` definitions; `run_simulation` has an extra method keyed on a leading `SimulationEngineConfig`. The final `include("precompile_workload.jl")` runs a representative workload at precompile time to reduce first-call latency. No state is held at module level beyond what the submodules own.

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

- [[module.spaceagora|SpaceAGORA]] · `api` → `module_api` · call · `src/SpaceAGORA.jl`

**Downstream**

- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/SpaceAGORA.jl:286-286`
- `callees` → [[core.effector_sampling_gravity_backbone_kick_structure|gravity_backbone_kick_structure]] · `callers` · call · `src/SpaceAGORA.jl:330-330`
- `callees` → [[core.effector_sampling_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/SpaceAGORA.jl:307-307`
- `callees` → [[core.effector_sampling_solver_partition|solver_partition]] · `callers` · call · `src/SpaceAGORA.jl:295-295`
- `callees` → [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/SpaceAGORA.jl:265-265`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/SpaceAGORA.jl:286-286`
- `callees` → [[dynamics.aerodynamic_wrench_models_solver_partition|solver_partition]] · `callers` · call · `src/SpaceAGORA.jl:295-295`
- `callees` → [[dynamics.aerodynamic_wrench_models_wrench|wrench]] · `callers` · call · `src/SpaceAGORA.jl:274-274`
- `callees` → [[dynamics.calc_force_torque|calcForceTorque]] · `callers` · call · `src/SpaceAGORA.jl:265-265`
- `callees` → [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callers` · call · `src/SpaceAGORA.jl:265-265`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/SpaceAGORA.jl:286-286`
- `callees` → [[dynamics.perturbations_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callers` · call · `src/SpaceAGORA.jl:319-319`
- `callees` → [[dynamics.perturbations_gravity_backbone_kick_acceleration_ii|gravity_backbone_kick_acceleration_ii]] · `callers` · call · `src/SpaceAGORA.jl:342-342`
- `callees` → [[dynamics.perturbations_gravity_backbone_kick_structure|gravity_backbone_kick_structure]] · `callers` · call · `src/SpaceAGORA.jl:330-330`
- `callees` → [[dynamics.perturbations_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/SpaceAGORA.jl:307-307`
- `callees` → [[dynamics.perturbations_wrench|wrench]] · `callers` · call · `src/SpaceAGORA.jl:274-274`
- `callees` → [[dynamics.robot_arm_reaction_effector_calcforcetorque|calcForceTorque]] · `callers` · call · `src/SpaceAGORA.jl:265-265`
- `callees` → [[envana.env_simple_ephemerides_simpleephemeridesmodel|SimpleEphemeridesModel]] · `callers` · call · `src/SpaceAGORA.jl:243-243`
- `callees` → [[environment.density_models_exponentialatmospheremodel|ExponentialAtmosphereModel]] · `callers` · call · `src/SpaceAGORA.jl:195-195`
- `callees` → [[environment.density_models_getdensitybatch_bang|getDensityBatch!]] · `callers` · call · `src/SpaceAGORA.jl:368-368`
- `callees` → [[environment.density_models_init_nrlmsise_space_indices_bang|init_nrlmsise_space_indices!]] · `callers` · call · `src/SpaceAGORA.jl:232-232`
- `callees` → [[environment.density_models_noatmospheremodel|NoAtmosphereModel]] · `callers` · call · `src/SpaceAGORA.jl:187-187`
- `callees` → [[environment.density_models_nrlmsise00atmospheremodel|NRLMSISE00AtmosphereModel]] · `callers` · call · `src/SpaceAGORA.jl:218-218`
- `callees` → [[environment.density_models_piecewiseexponentialatmospheremodel|PiecewiseExponentialAtmosphereModel]] · `callers` · call · `src/SpaceAGORA.jl:208-208`
- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `src/SpaceAGORA.jl:354-354`
- `callees` → [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/SpaceAGORA.jl:265-265`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/SpaceAGORA.jl:286-286`
- `callees` → [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callers` · call · `src/SpaceAGORA.jl:319-319`
- `callees` → [[environment.gravity_models_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/SpaceAGORA.jl:307-307`
- `callees` → [[environment.gravity_models_wrench|wrench]] · `callers` · call · `src/SpaceAGORA.jl:274-274`
- `callees` → [[gnc.momentum_manager_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/SpaceAGORA.jl:377-377`
- `callees` → [[gnc.momentum_manager_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/SpaceAGORA.jl:386-386`
- `callees` → [[gnc.propulsive_maneuvers_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/SpaceAGORA.jl:386-386`
- `callees` → [[gnc.propulsive_maneuvers_calccontrolmassflowrate|calcControlMassFlowRate]] · `callers` · call · `src/SpaceAGORA.jl:395-395`
- `callees` → [[gnc.robot_arm_control_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/SpaceAGORA.jl:377-377`
- `callees` → [[gnc.robot_arm_control_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/SpaceAGORA.jl:386-386`
- `callees` → [[gnc.robot_arm_control_calccontrolmassflowrate|calcControlMassFlowRate]] · `callers` · call · `src/SpaceAGORA.jl:395-395`
- `callees` → [[gnc.rpo_mpc_control_model_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/SpaceAGORA.jl:386-386`
- `callees` → [[gnc.rpo_mpc_control_model_calccontrolmassflowrate|calcControlMassFlowRate]] · `callers` · call · `src/SpaceAGORA.jl:395-395`
- `callees` → [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/SpaceAGORA.jl:377-377`
- `callees` → [[gnc.targeting_control_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/SpaceAGORA.jl:386-386`
- `callees` → [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/SpaceAGORA.jl:377-377`
- `callees` → [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/SpaceAGORA.jl:377-377`
- `callees` → [[simulation.run_simulation|run_simulation]] · `callers` · call · `src/SpaceAGORA.jl:492-492`
- `callees` → [[simx.engine_execution_run_simulation|run_simulation]] · `callers` · call · `src/SpaceAGORA.jl:492-492`
- `callees` → [[spaceagora.run_simulation|run_simulation]] · `callers` · call · `src/SpaceAGORA.jl:492-492`
<!-- vulcan:connections:end -->

## Limitations
Everything exported is decided statically here, so adding a public symbol in a submodule requires editing both the `using` and `export` blocks or it stays invisible to `using SpaceAGORA`. The `include` order is an implicit, undocumented dependency graph; reordering silently breaks loading. Because forwarding methods use `args...; kwargs...` they give no type-level signature or argument validation at this layer, and method errors surface from the inner module. The precompile workload is executed unconditionally, lengthening package precompilation.

## Provenance
Mapped from `src/SpaceAGORA.jl` line 3.
