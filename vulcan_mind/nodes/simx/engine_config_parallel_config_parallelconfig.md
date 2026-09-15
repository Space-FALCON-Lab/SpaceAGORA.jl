---
id: simx.engine_config_parallel_config_parallelconfig
label: ParallelConfig
kind: struct
source:
  file: src/simulation/engine/config/parallel_config.jl
  symbol: ParallelConfig
  lines:
  - 8
  - 17
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
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: parallel_policy
  type: ParallelConfig
  units: n/a
  description: Eight-field record naming the parallel profile and the routing mode
    for effectors, RHS batching and the density, control and thermal callbacks.
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
# ParallelConfig

## Purpose
`ParallelConfig` is the typed description of how a run is allowed to use threads. It names the parallel profile and then carries one routing mode per parallelisable subsystem, so the scheduler decisions taken deep inside the right-hand side are traceable to a single immutable record built at run start.

## Model & Assumptions
The profile string defaults to empty, meaning no named profile is selected and the policy layer falls back to its own heuristics. `outer_parallel_active` records whether an outer campaign pool is already consuming cores, which is how a Monte Carlo or ensemble run tells the inner engine to stop competing with itself. `parallel_policy_adaptive` enables the adaptive policy. The five mode strings each default to `"auto"`, meaning the policy layer chooses between serial and threaded execution from workload size.

## Design & Implementation
The five modes are `effector_parallel_mode`, `rhs_batch_parallel_mode`, `density_callback_parallel_mode`, `control_callback_parallel_mode` and `thermal_callback_parallel_mode`. They are strings rather than symbols because they arrive directly from environment variables and are compared against a small fixed vocabulary. The struct is built by `simulation_engine_config_from_env` from `SPACEAGORA_PARALLEL_PROFILE`, `SPACEAGORA_OUTER_PARALLEL_ACTIVE`, `SPACEAGORA_PARALLEL_POLICY_ADAPTIVE`, `SPACEAGORA_EFFECTOR_PARALLEL`, `SPACEAGORA_RHS_BATCH_PARALLEL`, `SPACEAGORA_DENSITY_CALLBACK_PARALLEL`, `SPACEAGORA_CONTROL_CALLBACK_PARALLEL` and `SPACEAGORA_THERMAL_CALLBACK_PARALLEL`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `parallel_policy` | ParallelConfig | n/a | — | Eight-field record naming the parallel profile and the routing mode for effectors, RHS batching and the density, control and thermal callbacks. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.simulation_engine_config|SimulationEngineConfig]] · `callees` → `callers` · call · `src/simulation/engine/config/simulation_engine_config.jl:9-9`
- [[simx.engine_adapters_from_env_simulation_engine_config_from_env|simulation_engine_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:162-162`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The record carries no thread counts; the actual budget comes from `ParallelPolicy.effective_inner_thread_budget`. An unrecognised mode string is not rejected at construction, so a misspelling is only caught where the mode is interpreted. Nested outer and inner parallelism is coordinated through the single `outer_parallel_active` flag, which is coarser than a real core reservation.

## Provenance
Mapped from `src/simulation/engine/config/parallel_config.jl:8-17`; populated at `src/simulation/engine/adapters/from_env.jl:162-171`.
