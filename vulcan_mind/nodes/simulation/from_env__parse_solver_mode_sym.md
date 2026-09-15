---
id: simulation.from_env__parse_solver_mode_sym
label: _parse_solver_mode_sym
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _parse_solver_mode_sym
  lines:
  - 13
  - 13
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
  description: Return value of `_parse_solver_mode_sym`.
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

# _parse_solver_mode_sym

## Purpose
Translates the `SPACEAGORA_SOLVER_MODE` string into the canonical solver-mode `Symbol` consumed by `SolverConfig`, accepting several spellings per mode.

## Design & Implementation
`_parse_solver_mode_sym(raw::String)::Symbol` lowercases and strips the input and returns `:tsit5` for `tsit5|default|<empty>`, `:symplectic` for `symplectic|kahanli8|verlet`, `:gravity_backbone_split` for `gravity_backbone_split|gravity-backbone-split|gravity_backbone|gravity-backbone`, `:auto_stiff` for `auto_stiff|auto-stiff|autostiff|auto`, `:rodas5p` for `rodas5p|rodas|stiff`, `:split_imex` for `split_imex|split-imex|split|imex`, `:multirate` for `multirate|multirate_split|split_multirate|mr`, and `:dp8` for `dp8|dormandprince8|dop8`. Unknown values throw an `ArgumentError` listing the eight canonical names.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | String | n/a | yes | Positional argument `raw`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_parse_solver_mode_sym`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- [[simulation.from_env__solver_config_from_env|_solver_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:80-80`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The alias set is hard-coded; adding a solver requires editing this function, `_engine_env_overrides` and the engine dispatch together. `kahanli8` and `verlet` both collapse to `:symplectic`, so the specific symplectic scheme cannot be selected through this variable.

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 13.
