---
id: simulation.solver_policy__multirate_solver_spec_from_sym
label: _multirate_solver_spec_from_sym
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _multirate_solver_spec_from_sym
  lines:
  - 237
  - 237
inputs:
- id: mode
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `mode`.
- id: field_name
  type: String
  units: n/a
  required: true
  description: Positional argument `field_name`.
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
  description: Return value of `_multirate_solver_spec_from_sym`.
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

# _multirate_solver_spec_from_sym

## Purpose
Translates a solver symbol into an algorithm instance, label, and `auto_switch_capable` flag for the slow or fast stage of the multirate integrator.

## Design & Implementation
Supported symbols: `:tsit5` -> `Tsit5()`, `:auto_stiff` -> `AutoTsit5(Rodas5P(autodiff=AutoFiniteDiff()))` with `auto_switch_capable=true`, `:rodas5p`, `:kencarp4`, and `:dp8`. Anything else throws `ArgumentError` naming `field_name` ("multirate_slow_solver" or "multirate_fast_solver") and the allowed set. `_multirate_slow_solver_spec` and `_multirate_fast_solver_spec` are thin wrappers over the two config fields.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `mode` | Symbol | n/a | yes | Positional argument `mode`. |
| in | `field_name` | String | n/a | yes | Positional argument `field_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_multirate_solver_spec_from_sym`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.solver_policy__multirate_fast_solver_spec|_multirate_fast_solver_spec]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:253-253`
- [[simulation.solver_policy__multirate_slow_solver_spec|_multirate_slow_solver_spec]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:252-252`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The `:auto_stiff` variant here does not pass `switch_max` from `cfg.auto_stiff_switch_max`, so multirate stages use OrdinaryDiffEq's default of 5 rather than the tuned value used by the main `:auto_stiff` mode. The dense `Rodas5P` never receives a sparse linear solver.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 237.
