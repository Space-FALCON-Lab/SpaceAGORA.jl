---
id: dynamics.aerodynamic_wrench_models__simulation_model_module_for_aero
label: _simulation_model_module_for_aero
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _simulation_model_module_for_aero
  lines:
  - 168
  - 168
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
  type: Any
  units: n/a
  description: Return value of `_simulation_model_module_for_aero`. Returns `mod`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _simulation_model_module_for_aero

## Purpose
Walks up the module ancestry from the including module to find the `SimulationModel` module that defines `planet_frame_lpi`.

## Design & Implementation
Starts at `@__MODULE__` and loops: if `isdefined(mod, :planet_frame_lpi)` returns `mod`; otherwise moves to `parentmodule(mod)` until the parent equals itself (top level), then calls `error(...)` with a descriptive message.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_simulation_model_module_for_aero`. Returns `mod`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Not called anywhere in this file; it is a lookup utility for other code. The search stops at the first ancestor defining the symbol, which could be a shadowing definition. `error` produces a generic `ErrorException` rather than a typed exception.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 168.
