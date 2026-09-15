---
id: simulation.save_fields__save_mass
label: _save_mass
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_mass
  lines:
  - 152
  - 152
inputs:
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
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
  description: Return value of `_save_mass`. Returns `masses`.
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

# _save_mass

## Purpose
Save-time getter for current spacecraft mass in kilograms, tracking propellant depletion and any staging events over the run.

## Design & Implementation
Marked `@inline`. It allocates a `Vector{Float64}` of length `num_sats` and fills it in an `@inbounds` loop with `_simulation_engine_module()._state_mass_kg(u, integrator.p.args, i)`. Passing `integrator.p.args` alongside the state lets the accessor resolve whether mass is an integrated state component, as it is when a thruster model burns propellant, or a static configuration value, and return a single consistent scalar either way. No unit conversion is applied, since the state already carries kilograms.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_save_mass`. Returns `masses`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:178-178`

**Downstream**

- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:155-155`
<!-- vulcan:connections:end -->

## Limitations
The single scalar hides the composition of the mass: dry mass, propellant, and any jettisoned elements are not separable in the output. Where mass is integrated, the saved value inherits any solver drift and can fall below the configured dry mass if a thruster model does not clamp its own depletion, and nothing in this getter detects or flags a non-positive result.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 152.
