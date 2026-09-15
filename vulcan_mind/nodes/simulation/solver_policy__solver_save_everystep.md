---
id: simulation.solver_policy__solver_save_everystep
label: _solver_save_everystep
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _solver_save_everystep
  lines:
  - 302
  - 302
inputs:
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
  description: Return value of `_solver_save_everystep`.
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

# _solver_save_everystep

## Purpose
Reads the `SPACEAGORA_SOLVER_SAVE_EVERYSTEP` switch through the engine environment adapter and parses it as a boolean.

## Design & Implementation
Fetches `_engine_env_get("SPACEAGORA_SOLVER_SAVE_EVERYSTEP", "true")`, lowercases and strips it, and returns membership in `("1", "true", "yes", "on")`. Callers only consult it when `_engine_env_haskey_with_env_fallback` says the key is set; otherwise `needs_full_solution` decides.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_solver_save_everystep`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`
- [[simulation.solver_policy__solve_with_explicit_solver|_solve_with_explicit_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:328-328`
- [[simulation.solver_policy__solve_with_fixed_step_solver|_solve_with_fixed_step_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:385-385`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/solver_policy.jl:303-303`
<!-- vulcan:connections:end -->

## Limitations
Any unrecognised value (including `"false "` variants with unusual characters) is treated as `false` rather than raising. The default string `"true"` is only used when the key is present but empty. Allocates a lowercase copy on every call.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 302.
