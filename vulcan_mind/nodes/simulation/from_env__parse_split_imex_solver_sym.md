---
id: simulation.from_env__parse_split_imex_solver_sym
label: _parse_split_imex_solver_sym
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _parse_split_imex_solver_sym
  lines:
  - 40
  - 40
inputs:
- id: raw
  type: String
  units: n/a
  required: true
  description: Positional argument `raw`.
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
  type: Symbol
  units: n/a
  description: Return value of `_parse_split_imex_solver_sym`.
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

# _parse_split_imex_solver_sym

## Purpose
Parses `SPACEAGORA_SPLIT_IMEX_SOLVER` into the KenCarp variant symbol used when the solver mode is `:split_imex`.

## Design & Implementation
`_parse_split_imex_solver_sym(raw::String)::Symbol` normalises with `lowercase(strip(raw))` and returns `:kencarp4` for `kencarp4|ken4|default`, `:kencarp47` for `kencarp47|ken47`, and `:kencarp58` for `kencarp58|ken58`. Anything else throws `ArgumentError("Unsupported SPACEAGORA_SPLIT_IMEX_SOLVER='<raw>'. Use one of: kencarp4, kencarp47, kencarp58.")`. `_solver_config_from_env` supplies `"kencarp4"` as the default and wraps the call in `_parse_or_default`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | String | n/a | yes | Positional argument `raw`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_parse_split_imex_solver_sym`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- [[simulation.from_env__solver_config_from_env|_solver_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:103-103`

**Downstream**

- `callees` → [[simulation.from_env__solver_config_from_env|_solver_config_from_env]] · `callers` · feedback · `src/simulation/engine/adapters/from_env.jl:51-51`
<!-- vulcan:connections:end -->

## Limitations
An empty string is not accepted as default, so an explicitly blank variable errors in strict mode. Only three IMEX schemes are exposed; other OrdinaryDiffEq IMEX methods cannot be selected without code changes.

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 40.
