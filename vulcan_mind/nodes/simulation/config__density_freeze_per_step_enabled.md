---
id: simulation.config__density_freeze_per_step_enabled
label: _density_freeze_per_step_enabled
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _density_freeze_per_step_enabled
  lines:
  - 20
  - 20
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
  description: Return value of `_density_freeze_per_step_enabled`.
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

# _density_freeze_per_step_enabled

## Purpose
Reports whether atmospheric density should be sampled once per integrator step and held constant across all right-hand-side evaluations within that step.

## Design & Implementation
Returns `_parse_bool_env("SPACEAGORA_DENSITY_FREEZE_PER_STEP", false)`. The source comment defers the rationale to the `CallbackEnvConfig.density_freeze_per_step` docstring. Defaulting to `false` keeps the physically exact behaviour of re-evaluating density at every stage of the integrator.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_density_freeze_per_step_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:183-183`

**Downstream**

- `callees` → [[simulation.config__parse_bool_env|_parse_bool_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:21-21`
<!-- vulcan:connections:end -->

## Limitations
When enabled, the frozen value makes the right-hand side inconsistent between Runge-Kutta stages, which can degrade the effective order of the integrator and confuse adaptive error control during rapid density changes near periapsis. The knob is global: it cannot be enabled per satellite or per density model.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 20.
