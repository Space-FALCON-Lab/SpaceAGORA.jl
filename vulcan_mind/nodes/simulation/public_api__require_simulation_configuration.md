---
id: simulation.public_api__require_simulation_configuration
label: _require_simulation_configuration
kind: function
source:
  file: src/simulation/engine/public_api.jl
  symbol: _require_simulation_configuration
  lines:
  - 5
  - 5
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
  type: Any
  units: n/a
  description: Return value of `_require_simulation_configuration`. Returns `args`.
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

# _require_simulation_configuration

## Purpose
Type gate at the `run_simulation` entry boundary. It asserts that the caller supplied a `SimulationConfiguration` and returns that same object unchanged so downstream dispatch sees a concrete, fully typed configuration.

## Design & Implementation
Declared `@inline` and written as a single short-circuit: `args isa SimulationConfiguration || throw(ArgumentError(...))`, with the message interpolating `typeof(args)` so the rejected type is named in the error text. The returned value is the argument itself, which lets callers write `typed_args = _require_simulation_configuration(args)` and then invoke the typed `run_simulation` method without a second conversion step. It performs no coercion and no field validation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_require_simulation_configuration`. Returns `args`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/public_api.jl:25-25`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The check is nominal only: any value of the right struct type passes even if its inner fields are empty, inconsistent, or dimensionally wrong, so configuration errors surface much later inside the solver. Subtypes or look-alike configuration objects from other packages are rejected outright with no conversion path offered.

## Provenance
Mapped from `src/simulation/engine/public_api.jl` line 5.
