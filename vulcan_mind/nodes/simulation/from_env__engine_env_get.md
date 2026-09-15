---
id: simulation.from_env__engine_env_get
label: _engine_env_get
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _engine_env_get
  lines:
  - 199
  - 199
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: String
  units: n/a
  required: false
  description: Positional argument `default` (default `""`).
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
  type: String
  units: n/a
  description: Return value of `_engine_env_get`.
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

# _engine_env_get

## Purpose
Primary environment lookup for the simulation engine: returns a `SPACEAGORA_*` knob from the active `SimulationEngineConfig` override scope when one is installed, otherwise from the process `ENV`, so engine code can read configuration without knowing which source is in effect.

## Design & Implementation
`_engine_env_get(name::String, default::String="")::String` is `@inline`. It reads the global `_engine_active_overrides_ref[]` (a `Ref{Union{Nothing, Dict{String,String}}}`). If non-`nothing`, the override dict is consulted exclusively with `get(active_overrides, name, default)`; only when no scope is active does it fall back to `get(ENV, name, default)`. The result is always converted to `String`. It is the default `env_get` callable for `_solver_config_from_env`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | String | n/a | no | Positional argument `default` (default `""`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_engine_env_get`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.rhs_calibration__calib_machine_label|_calib_machine_label]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:49-49`
- [[simulation.rhs_calibration__rhs_calib_path|_rhs_calib_path]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:95-95`
- [[simulation.rhs_calibration__rhs_calibrate_n_timed|_rhs_calibrate_n_timed]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:41-41`
- [[simulation.rhs_calibration__rhs_calibrate_n_warmup|_rhs_calibrate_n_warmup]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:36-36`
- [[simulation.rhs_calibration__rhs_calibration_mode|_rhs_calibration_mode]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:26-26`
- [[simulation.setup__density_without_aero_warning_enabled|_density_without_aero_warning_enabled]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:66-66`
- [[simulation.setup__gram_per_sat_instances_enabled|_gram_per_sat_instances_enabled]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:148-148`
- [[simulation.setup__parse_nonnegative_int_env|_parse_nonnegative_int_env]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:182-182`
- [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:160-160`
- [[simulation.setup__parse_unit_float_env|_parse_unit_float_env]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:171-171`
- [[simulation.setup__profile_forces_serial_rhs|_profile_forces_serial_rhs]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:376-376`
- [[simulation.setup__rhs_execution_mode_env|_rhs_execution_mode_env]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:759-759`
- [[simulation.setup__typed_allow_transition_normalize|_typed_allow_transition_normalize]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:31-31`
- [[simulation.setup__typed_normalize_warning_enabled|_typed_normalize_warning_enabled]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:30-30`
- [[simulation.setup__typed_save_bundle_enabled|_typed_save_bundle_enabled]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:32-32`
- [[simulation.solver_policy__solver_bool_env|_solver_bool_env]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:308-308`
- [[simulation.solver_policy__solver_save_everystep|_solver_save_everystep]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:303-303`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:204-204`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Inside an active override scope, plain `ENV` values for names absent from the override dict are ignored entirely, so a knob that `_engine_env_overrides` does not enumerate becomes invisible; the `_with_env_fallback` variants exist for that case. The override `Ref` is process-global, not task-local, so concurrent scopes on different tasks interfere.

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 199.
