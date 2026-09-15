---
id: simulation.solver_policy__requires_componentwise_tolerances
label: _requires_componentwise_tolerances
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _requires_componentwise_tolerances
  lines:
  - 8
  - 8
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  description: Return value of `_requires_componentwise_tolerances`.
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

# _requires_componentwise_tolerances

## Purpose
Decides whether the solver must build per-component tolerance vectors rather than using the scalar orbit tolerances.

## Design & Implementation
Reads `args.integration_tolerances` and returns `true` when `args.mission_configuration.orientation_sim` is set or when any of `reltol_mass`, `abstol_mass`, `reltol_heat_load`, `abstol_heat_load` is non-zero. Quaternion and angular-rate tolerances are not checked directly because they only matter when `orientation_sim` is true.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_requires_componentwise_tolerances`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`
- [[simulation.solver_policy__build_solver_tolerances|_build_solver_tolerances]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:17-17`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Non-zero `reltol_angular_rate` or `abstol_angular_rate` with `orientation_sim=false` is silently ignored. The decision is boolean-only, so a mixed configuration still pays the full `ComponentVector` copy cost in `_build_solver_tolerances` even if only one component differs.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 8.
