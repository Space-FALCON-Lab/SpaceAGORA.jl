---
id: simulation.solver_policy__solver_bool_env
label: _solver_bool_env
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _solver_bool_env
  lines:
  - 307
  - 307
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: Bool
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
  type: Bool
  units: n/a
  description: Return value of `_solver_bool_env`.
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

# _solver_bool_env

## Purpose
Generic boolean parser for solver environment switches (`SPACEAGORA_SOLVER_SAVE_ON`, `_SAVE_START`, `_SAVE_END`) with a caller-supplied default.

## Design & Implementation
Calls `_engine_env_get(name, default ? "true" : "false")`, normalises with `lowercase(strip(...))`, and returns `raw in ("1", "true", "yes", "on")`. Used by `_solve_with_explicit_solver` and `_solve_with_fixed_step_solver`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | Bool | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_solver_bool_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`
- [[simulation.solver_policy__solve_with_explicit_solver|_solve_with_explicit_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:330-330`
- [[simulation.solver_policy__solve_with_fixed_step_solver|_solve_with_fixed_step_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:387-387`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/solver_policy.jl:308-308`
<!-- vulcan:connections:end -->

## Limitations
Values such as `"0"`, `"no"`, `"off"`, or typos all map to `false` indistinguishably, so a mis-spelled `"ture"` disables saving silently. The function is called on every solve, allocating strings each time.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 307.
