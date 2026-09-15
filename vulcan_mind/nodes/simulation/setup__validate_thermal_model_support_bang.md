---
id: simulation.setup__validate_thermal_model_support_bang
label: _validate_thermal_model_support!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _validate_thermal_model_support!
  lines:
  - 55
  - 55
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
  type: Nothing
  units: n/a
  description: Return value of `_validate_thermal_model_support!`; mutates `args`
    in place. Returns `nothing`.
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

# _validate_thermal_model_support!

## Purpose
Ensures the configured thermal model implements the heat-rate interface the engine will call, so a missing method surfaces at setup rather than deep inside the RHS.

## Design & Implementation
Fetches `args.environment_model.thermal_model` and checks `hasmethod(SimulationModel.getHeatRate, Tuple{typeof(thermal_model), Float64, Float64, Float64, Float64, Float64})`. The five `Float64` slots correspond to the documented signature `getHeatRate(model, S, T, ρ, v, α)`, meaning heat-flux area, temperature, density, velocity magnitude, and angle of attack. On failure it throws `ArgumentError` naming the model type and the required signature. Returns `nothing`; no state is mutated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_validate_thermal_model_support!`; mutates `args` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:165-165`

**Downstream**

- `callees` → [[vehx.thermal_models_getheatrate|getHeatRate]] · `callers` · call · `src/simulation/engine/setup.jl:60-60`
<!-- vulcan:connections:end -->

## Limitations
`hasmethod` with concrete `Float64` argument types accepts methods declared on abstract `Real`, but would reject a model whose method is specialised on `Float32` even if it would work through promotion. The return type annotation `::Float64` mentioned in the error is not verified. Models implementing the method via a fallback on `Any` also pass, so the check cannot detect a placeholder.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 55.
