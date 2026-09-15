---
id: simulation.from_env__engine_env_overrides
label: _engine_env_overrides
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _engine_env_overrides
  lines:
  - 236
  - 236
inputs:
- id: config
  type: SimulationEngineConfig
  units: n/a
  required: true
  description: Positional argument `config`.
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
  type: Dict{String,
  units: n/a
  description: Return value of `_engine_env_overrides`.
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

# _engine_env_overrides

## Purpose
Serialises a typed `SimulationEngineConfig` back into the `SPACEAGORA_*` string environment contract, producing the override dictionary that `_with_engine_env_overrides` installs so legacy `ENV`-reading code observes the config.

## Design & Implementation
`_engine_env_overrides(config::SimulationEngineConfig)::Dict{String,String}` builds a dict with 23 fixed keys: boolean runtime-policy and artifact flags rendered with `_env_bool` as `"1"`/`"0"` (`WARN_NORMALIZE`, `ALLOW_TYPED_NORMALIZE`, `GRAM_PER_SAT_INSTANCES`, `SRP_EPHEMERIS_CACHE`, `NBODY_EPHEMERIS_CACHE`, `PLANET_FRAME_CACHE`, `SPICE_RHS_MEMO`, `SAVE_BUNDLE`, `WARN_DEPRECATED_CONFIG`, `OUTER_PARALLEL_ACTIVE`, `PARALLEL_POLICY_ADAPTIVE`, `AUTO_STIFF_GRAVITY_TSIT5`), the five parallel mode strings, and solver symbols/ints via `string(...)`. Optional fields are added only when set: `PARALLEL_PROFILE` if non-empty, and `SOLVER_MAXITERS`, `SYMPLECTIC_DT_S`, `GRAVITY_BACKBONE_DT_S`, `MULTIRATE_SLOW_DT_S` if not `nothing`. Finally `merge!(overrides, config.env_overrides)` lets free-form user overrides win.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | SimulationEngineConfig | n/a | yes | Positional argument `config`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Dict{String, | n/a | — | Return value of `_engine_env_overrides`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- [[simulation.from_env__with_engine_env_overrides|_with_engine_env_overrides]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:274-274`

**Downstream**

- `callees` → [[simulation.from_env__env_bool|_env_bool]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:238-238`
<!-- vulcan:connections:end -->

## Limitations
The key list must be kept in sync by hand with `simulation_engine_config_from_env`; a new config field that is not added here is silently dropped when round-tripping through the environment. `config.env_overrides` can overwrite any typed value without validation. `Float64` fields are rendered with `string`, so values like `1.0e-5` re-parse correctly but lose no precision only because both sides use shortest round-trip printing.

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 236.
