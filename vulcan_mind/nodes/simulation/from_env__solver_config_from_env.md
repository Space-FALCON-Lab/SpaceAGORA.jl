---
id: simulation.from_env__solver_config_from_env
label: _solver_config_from_env
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _solver_config_from_env
  lines:
  - 78
  - 78
inputs:
- id: env_get
  type: Any
  units: n/a
  required: false
  description: Positional argument `env_get` (default `_engine_env_get`).
- id: strict
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `strict` (default `true`).
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
  type: SolverConfig
  units: n/a
  description: Return value of `_solver_config_from_env`.
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

# _solver_config_from_env

## Purpose
Assembles a typed `SolverConfig` from the `SPACEAGORA_SOLVER_*`, `SPACEAGORA_SYMPLECTIC_*`, `SPACEAGORA_GRAVITY_BACKBONE_*`, `SPACEAGORA_SPLIT_IMEX_*`, `SPACEAGORA_MULTIRATE_*` and `SPACEAGORA_AUTO_STIFF_*` knobs, using an injectable lookup so it can read either the process environment or an explicit dictionary.

## Design & Implementation
`_solver_config_from_env(env_get=_engine_env_get; strict::Bool=true)::SolverConfig`. Each knob is read through `env_get(name, default)` and parsed inside `_parse_or_default(strict, default) do ... end`: `solver_mode` via `_parse_solver_mode_sym` (default `:tsit5`); `maxiters` as an optional positive `Int` (empty means `nothing`); `symplectic_dt_s`, `gravity_backbone_dt_s`, `multirate_slow_dt_s` via `_parse_float_opt`; `split_imex_solver` (default `:kencarp4`); `multirate_fast_substeps` as a positive `Int` (default 8); `multirate_slow_solver` (default `:tsit5`) and `multirate_fast_solver` (default `:auto_stiff`) via `_parse_multirate_solver_sym`. Two knobs bypass `env_get`: `auto_stiff_gravity_tsit5` uses `ParallelPolicy.parse_bool_env(..., true)` and `auto_stiff_switch_max` uses `parse_thread_threshold_env(..., 50)`, both reading `ENV` directly. The keyword constructor of `SolverConfig` receives all eleven fields.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `env_get` | Any | n/a | no | Positional argument `env_get` (default `_engine_env_get`). |
| in | `strict` | Bool | n/a | no | Keyword argument `strict` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SolverConfig | n/a | — | Return value of `_solver_config_from_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.from_env__parse_split_imex_solver_sym|_parse_split_imex_solver_sym]] · `callees` → `callers` · feedback · `src/simulation/engine/adapters/from_env.jl:51-51`
- [[simx.engine_adapters_from_env_simulation_engine_config_from_env|simulation_engine_config_from_env]] · `callees` → `callers` · feedback · `src/simulation/engine/adapters/from_env.jl:174-174`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:158-158`

**Downstream**

- `callees` → [[core.simulation_configuration_solverconfig|SolverConfig]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:128-128`
- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:125-125`
- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:126-126`
- `callees` → [[simulation.from_env__parse_float_opt|_parse_float_opt]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:96-96`
- `callees` → [[simulation.from_env__parse_multirate_solver_sym|_parse_multirate_solver_sym]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:119-119`
- `callees` → [[simulation.from_env__parse_or_default|_parse_or_default]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:79-79`
- `callees` → [[simulation.from_env__parse_solver_mode_sym|_parse_solver_mode_sym]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:80-80`
- `callees` → [[simulation.from_env__parse_split_imex_solver_sym|_parse_split_imex_solver_sym]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:103-103`
- `callees` → [[simx.engine_adapters_from_env_simulation_engine_config_from_env|simulation_engine_config_from_env]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:144-144`
- `callees` → [[simx.engine_config_solver_config_solverconfig|SolverConfig]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:128-128`
<!-- vulcan:connections:end -->

## Limitations
The two `AUTO_STIFF` knobs ignore the injected `env_get`, so passing a custom dictionary to `simulation_engine_config_from_env` does not control them and they leak from the real process environment. Positive-integer checks throw `ArgumentError` with hard-coded messages. No cross-field validation is performed (for example a symplectic dt with a non-symplectic mode is accepted).

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 78.
