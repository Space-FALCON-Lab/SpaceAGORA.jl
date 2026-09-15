---
id: simx.engine_adapters_from_env_simulation_engine_config_from_env
label: simulation_engine_config_from_env
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: simulation_engine_config_from_env
  lines:
  - 161
  - 197
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: engine_api
  type: Module
  units: n/a
  required: true
  description: RuntimeServices namespace and the SimulationEngine/SimulationCampaigns
    scope it aggregates, supplying the SPICE lock and the SimulationModel types used
    here.
- id: env
  type: AbstractDict
  units: n/a
  required: true
  description: Environment mapping to read, defaulting to the process ENV; tests pass
    a plain Dict to build a configuration without touching global state.
- id: solver_strict
  type: Bool
  units: n/a
  required: true
  description: When true, malformed solver environment values raise instead of falling
    back to defaults.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: engine_config
  type: SimulationEngineConfig
  units: n/a
  description: Typed configuration aggregating ParallelConfig, SolverConfig, RuntimePolicyConfig
    and ArtifactConfig.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simx
origin: agent
---
# simulation_engine_config_from_env

## Purpose
`simulation_engine_config_from_env` is the single adapter that converts the `SPACEAGORA_*` environment variables into the typed engine configuration tree. Every hot path in the engine reads struct fields instead of calling `get(ENV, ...)`, and this function is where that translation happens exactly once per run.

## Model & Assumptions
Defaults encode the shipping behaviour: parallel routing modes default to the string `"auto"`, the ephemeris caches (`srp_ephemeris_cache`, `nbody_ephemeris_cache`, `planet_frame_cache`) and `spice_rhs_memo` default to true, `warn_normalize` and `save_bundle` default to true, while `allow_typed_normalize`, `gram_per_sat_instances`, `outer_parallel_active` and `parallel_policy_adaptive` default to false. Boolean parsing goes through `_parse_bool`, which accepts a missing value and returns the supplied default.

## Design & Implementation
The four sub-configurations are constructed in order. `ParallelConfig` reads the profile name plus five per-subsystem routing modes covering effectors, RHS batching and the density, control and thermal callbacks. The solver branch wraps the incoming dictionary in a closure `env_get = (name, default) -> String(get(env, name, default))` and hands it to `_solver_config_from_env`, which is where solver-mode, multirate and split-IMEX symbols are parsed and where `strict` decides between raising and defaulting. `RuntimePolicyConfig` and `ArtifactConfig` are plain boolean reads. Because the parameter is any `AbstractDict`, the whole engine configuration is testable without mutating process environment state.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `env` | AbstractDict | n/a | yes | Environment mapping to read, defaulting to the process ENV; tests pass a plain Dict to build a configuration without touching global state. |
| in | `solver_strict` | Bool | n/a | yes | When true, malformed solver environment values raise instead of falling back to defaults. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `engine_config` | SimulationEngineConfig | n/a | — | Typed configuration aggregating ParallelConfig, SolverConfig, RuntimePolicyConfig and ArtifactConfig. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[misc.precompile_workload_run_spaceagora_precompile_workload|_run_spaceagora_precompile_workload]] · `callees` → `callers` · call · `src/precompile_workload.jl:43-43`
- [[simulation.from_env__solver_config_from_env|_solver_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:144-144`
- [[simulation.solver_policy__auto_stiff_switched|_auto_stiff_switched]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:82-82`

**Downstream**

- `callees` → [[simulation.from_env__parse_bool|_parse_bool]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:164-164`
- `callees` → [[simulation.from_env__solver_config_from_env|_solver_config_from_env]] · `callers` · feedback · `src/simulation/engine/adapters/from_env.jl:174-174`
- `callees` → [[simulation.simulation_engine_config|SimulationEngineConfig]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:191-191`
- `callees` → [[simx.engine_config_artifact_config_artifactconfig|ArtifactConfig]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:186-186`
- `callees` → [[simx.engine_config_parallel_config_parallelconfig|ParallelConfig]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:162-162`
- `callees` → [[simx.engine_config_runtime_policy_config_runtimepolicyconfig|RuntimePolicyConfig]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:176-176`
<!-- vulcan:connections:end -->

## Limitations
Only the variables named here are honoured; a typo in a `SPACEAGORA_*` name is indistinguishable from the variable being unset and silently yields the default. The function reads the environment once, so changing a variable after the configuration is built has no effect until the next call, and `_engine_active_overrides_ref` scoping applies only through `_with_engine_env_overrides`.

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl:161-197`, with the solver parser at line 78 and the override scope helper at line 273 of the same file.
