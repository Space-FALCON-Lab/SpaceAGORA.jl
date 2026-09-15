---
id: simulation.solver_policy__multirate_slow_dt_s
label: _multirate_slow_dt_s
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _multirate_slow_dt_s
  lines:
  - 229
  - 229
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
  type: Float64
  units: n/a
  description: Return value of `_multirate_slow_dt_s`.
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

# _multirate_slow_dt_s

## Purpose
Computes the macro (slow) step in seconds for Strang-split multirate integration, defaulting to `min(dt_max_orbit, 2.0)` and never exceeding `dt_max_orbit`.

## Design & Implementation
`default_dt = min(args.integration_tolerances.dt_max_orbit, 2.0)`; picks `cfg.multirate_slow_dt_s` when set; throws `ArgumentError` unless `dt > 0.0`; returns `min(dt, dt_max_orbit)`. A one-argument overload uses `_active_solver_config()`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | SolverConfig | n/a | yes | Positional argument `cfg`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_multirate_slow_dt_s`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`
- [[simulation.solver_policy__solve_with_multirate_solver|_solve_with_multirate_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:429-429`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The 2.0 s default cap is a hard-coded literal chosen for aerobraking passes and may be far too small for high orbits, inflating macro-step counts. A configured value above `dt_max_orbit` is silently clamped without a log message.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 229.
