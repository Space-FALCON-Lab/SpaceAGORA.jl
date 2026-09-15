---
id: simx.engine_config_runtime_policy_config_runtimepolicyconfig
label: RuntimePolicyConfig
kind: struct
source:
  file: src/simulation/engine/config/runtime_policy_config.jl
  symbol: RuntimePolicyConfig
  lines:
  - 7
  - 15
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
- id: runtime_policy
  type: RuntimePolicyConfig
  units: n/a
  description: Seven-field record controlling normalization policy, GRAM instancing
    and the three ephemeris caches plus the SPICE right-hand-side memo.
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
# RuntimePolicyConfig

## Purpose
`RuntimePolicyConfig` collects the execution-time switches that change how the engine caches and how strictly it polices legacy inputs. It is read once per run and then consulted by the setup phase, so a hot right-hand-side call never re-parses these decisions.

## Model & Assumptions
`warn_normalize` defaults to true and `allow_typed_normalize` to false, which together mean a configuration requesting normalized propagation is rejected outright by the typed pipeline unless the transition escape hatch is set. `gram_per_sat_instances` defaults to false so a single GRAM atmosphere instance is shared and serialised behind the shared lock. The three cache flags `srp_ephemeris_cache`, `nbody_ephemeris_cache` and `planet_frame_cache` default to true, as does `spice_rhs_memo`.

## Design & Implementation
All seven fields are `Bool` with `Base.@kwdef` defaults, and the struct is immutable so that a cache flag cannot change between the setup phase that allocates the cache buffer and the integration that reads it. The record is populated from `SPACEAGORA_WARN_NORMALIZE`, `SPACEAGORA_ALLOW_TYPED_NORMALIZE`, `SPACEAGORA_GRAM_PER_SAT_INSTANCES`, `SPACEAGORA_SRP_EPHEMERIS_CACHE`, `SPACEAGORA_NBODY_EPHEMERIS_CACHE`, `SPACEAGORA_PLANET_FRAME_CACHE` and `SPACEAGORA_SPICE_RHS_MEMO`. The cache flags gate the corresponding `_initialize_*_cache_buffer!` calls in the setup file.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `runtime_policy` | RuntimePolicyConfig | n/a | — | Seven-field record controlling normalization policy, GRAM instancing and the three ephemeris caches plus the SPICE right-hand-side memo. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.simulation_engine_config|SimulationEngineConfig]] · `callees` → `callers` · call · `src/simulation/engine/config/simulation_engine_config.jl:11-11`
- [[simx.engine_adapters_from_env_simulation_engine_config_from_env|simulation_engine_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:176-176`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Turning a cache off does not free memory already allocated by an earlier run in the same process, because the reuse caches are module level dictionaries keyed by mission window. The flags are all-or-nothing per cache type; there is no per-body or per-spacecraft granularity, and no field bounds the cache memory footprint.

## Provenance
Mapped from `src/simulation/engine/config/runtime_policy_config.jl:7-15`; populated at `src/simulation/engine/adapters/from_env.jl:176-184`.
