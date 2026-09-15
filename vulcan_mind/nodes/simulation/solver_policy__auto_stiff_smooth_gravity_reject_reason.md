---
id: simulation.solver_policy__auto_stiff_smooth_gravity_reject_reason
label: _auto_stiff_smooth_gravity_reject_reason
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _auto_stiff_smooth_gravity_reject_reason
  lines:
  - 120
  - 120
inputs:
- id: cfg
  type: SolverConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_auto_stiff_smooth_gravity_reject_reason`.
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

# _auto_stiff_smooth_gravity_reject_reason

## Purpose
Explains, as a `String`, why an `:auto_stiff` run cannot be routed to plain `Tsit5`, or returns `nothing` when the run is smooth-gravity-only.

## Design & Implementation
Sequential guards each returning a human-readable reason: `cfg.auto_stiff_gravity_tsit5` false; `orientation_sim` true; non-empty control, guidance, or navigation effectors; empty dynamic effectors. Then, for each effector, it requires `_auto_stiff_smooth_gravity_effector` and calls `SimulationModel.environment_requirements(effector)`, rejecting if `req.atmosphere` or `req.solar` is set. `_auto_stiff_smooth_gravity_eligible` wraps this as `isnothing(...)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | SolverConfig | n/a | yes | Positional argument `cfg`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_auto_stiff_smooth_gravity_reject_reason`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`

**Downstream**

- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/solver_policy.jl:131-131`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/solver_policy.jl:131-131`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/solver_policy.jl:131-131`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/solver_policy.jl:131-131`
- `callees` → [[simulation.solver_policy__auto_stiff_smooth_gravity_effector|_auto_stiff_smooth_gravity_effector]] · `callers` · call · `src/simulation/engine/solver_policy.jl:130-130`
<!-- vulcan:connections:end -->

## Limitations
Reasons are English strings rather than an enum, so callers cannot dispatch programmatically. The reason text is built with string interpolation on every call even when only the boolean is needed. Third-body ephemeris dependence (`req.third_body_names`) is not checked here, unlike the gravity-backbone gate.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 120.
