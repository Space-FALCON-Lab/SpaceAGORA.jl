---
id: simulation.from_env__parse_multirate_solver_sym
label: _parse_multirate_solver_sym
kind: function
source:
  file: src/simulation/engine/adapters/from_env.jl
  symbol: _parse_multirate_solver_sym
  lines:
  - 28
  - 28
inputs:
- id: raw
  type: String
  units: n/a
  required: true
  description: Positional argument `raw`.
- id: env_name
  type: String
  units: n/a
  required: true
  description: Positional argument `env_name`.
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
  description: Return value of `_parse_multirate_solver_sym`.
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

# _parse_multirate_solver_sym

## Purpose
Maps the string value of `SPACEAGORA_MULTIRATE_SLOW_SOLVER` or `SPACEAGORA_MULTIRATE_FAST_SOLVER` onto one of the solver symbols the multirate integrator understands, with the variable name threaded through for the error message.

## Design & Implementation
`_parse_multirate_solver_sym(raw::String, env_name::String)::Symbol` lowercases and strips `raw`, then matches alias groups: `tsit5|tsit|default -> :tsit5`, `auto_stiff|auto-stiff|autostiff|auto -> :auto_stiff`, `rodas5p|rodas|stiff -> :rodas5p`, `kencarp4|ken4 -> :kencarp4`, `dp8|dormandprince8|dop8 -> :dp8`. Any other value throws `ArgumentError("Unsupported <env_name>='<raw>'. Use one of: tsit5, dp8, auto_stiff, rodas5p, kencarp4.")`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | String | n/a | yes | Positional argument `raw`. |
| in | `env_name` | String | n/a | yes | Positional argument `env_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_parse_multirate_solver_sym`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- [[simulation.from_env__solver_config_from_env|_solver_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:119-119`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
An empty string is not aliased to the default here (unlike `_parse_solver_mode_sym`), so an explicitly blank variable raises. The alias table is separate from the top-level solver-mode table, so the two can drift. `:symplectic` and `:gravity_backbone_split` are intentionally not selectable for multirate sub-solvers.

## Provenance
Mapped from `src/simulation/engine/adapters/from_env.jl` line 28.
