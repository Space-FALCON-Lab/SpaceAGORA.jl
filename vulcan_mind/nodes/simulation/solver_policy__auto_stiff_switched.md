---
id: simulation.solver_policy__auto_stiff_switched
label: _auto_stiff_switched
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _auto_stiff_switched
  lines:
  - 66
  - 66
inputs:
- id: sol
  type: Any
  units: n/a
  required: true
  description: Positional argument `sol`.
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
  description: Return value of `_auto_stiff_switched`.
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

# _auto_stiff_switched

## Purpose
Detects whether an `AutoTsit5` solution ever switched between its explicit and implicit algorithms by inspecting `sol.alg_choice`.

## Design & Implementation
Returns `false` if `sol` lacks `alg_choice` or the vector is empty. Otherwise records `first(choices)` and scans with `@inbounds` for any element that differs, returning `true` at the first difference. This is used by `_solve_with_solver_policy` to populate `fallback_used` and by the multirate driver to count `auto_switch_events`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sol` | Any | n/a | yes | Positional argument `sol`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_auto_stiff_switched`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.solver_policy__rodas5p_alg|_rodas5p_alg]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:700-700`
- [[simulation.solver_policy__solve_with_multirate_solver|_solve_with_multirate_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:458-458`
- [[simx.engine_solver_policy_solve_with_solver_policy|_solve_with_solver_policy]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:700-700`

**Downstream**

- `callees` → [[simx.engine_adapters_from_env_simulation_engine_config_from_env|simulation_engine_config_from_env]] · `callers` · call · `src/simulation/engine/solver_policy.jl:82-82`
<!-- vulcan:connections:end -->

## Limitations
Requires per-step saving: when `save_everystep=false` the `alg_choice` vector holds only endpoints and a mid-run switch that switched back is invisible. A solution that started on the implicit solver and stayed there reports `false`. The scan is O(number of steps).

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 66.
