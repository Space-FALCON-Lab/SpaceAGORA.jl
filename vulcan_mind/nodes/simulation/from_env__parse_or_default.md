---
id: simulation.from_env__parse_or_default
label: _parse_or_default
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _parse_or_default
  lines:
  - 72
  - 72
inputs:
- id: f
  type: Function
  units: n/a
  required: true
  description: Positional argument `f`.
- id: strict
  type: Bool
  units: n/a
  required: true
  description: Positional argument `strict`.
- id: default
  type: Any
  units: n/a
  required: true
  description: Positional argument `default`.
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
  description: 'Return value of `_parse_or_default`. Returns `strict ? f() : (try     f()
    catch e     e isa ArgumentError ? default : rethrow(`.'
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

# _parse_or_default

## Purpose
Strictness switch for individual solver knobs: in strict mode a parser's `ArgumentError` propagates, in lenient mode a malformed knob falls back to its default so one typo cannot prevent construction of the whole `SimulationEngineConfig`.

## Design & Implementation
`_parse_or_default(f::Function, strict::Bool, default)`. With `strict == true` it simply returns `f()`. Otherwise it wraps `f()` in `try`/`catch` and returns `default` only when the caught exception `isa ArgumentError`; every other exception type is re-raised with `rethrow()`. The `do`-block call form in `_solver_config_from_env` makes each knob read as `_parse_or_default(strict, default) do ... end`. `_active_solver_config` uses `strict=true`; `simulation_engine_config_from_env` defaults to `solver_strict=false`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | Function | n/a | yes | Positional argument `f`. |
| in | `strict` | Bool | n/a | yes | Positional argument `strict`. |
| in | `default` | Any | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_parse_or_default`. Returns `strict ? f() : (try     f() catch e     e isa ArgumentError ? default : rethrow(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.from_env__solver_config_from_env|_solver_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:79-79`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:72-72`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:72-72`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/simulation/engine/adapters/from_env.jl:72-72`
<!-- vulcan:connections:end -->

## Limitations
Lenient mode swallows the error silently with no logging, so a misconfigured solver knob is only discoverable by inspecting the resulting config. Only `ArgumentError` is treated as recoverable; a `tryparse` returning `nothing` must be converted to `ArgumentError` by the parser for the fallback to engage.

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 72.
