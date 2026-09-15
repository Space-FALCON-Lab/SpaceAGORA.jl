---
id: simulation.runtime__extract_mass_kg
label: _extract_mass_kg
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/runtime.jl
  symbol: _extract_mass_kg
  lines:
  - 5
  - 5
inputs:
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
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
  description: Return value of `_extract_mass_kg`.
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

# _extract_mass_kg

## Purpose
Recovers the current spacecraft mass from a state vector that may or may not carry one, so the track cache can use it for drag prediction without every configuration being mass-augmented.

## Design & Implementation
If the state has a `mass` property it reads slot seven as `Float64`. Otherwise it returns slot seven when the vector has at least seven elements and `NaN` when it does not. The `NaN` sentinel lets the downstream cache refresh detect the absent case with a finiteness test instead of a separate flag. `@inline` with a `::Float64` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_extract_mass_kg`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime_update_density_sat_bang|update_density_sat!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:232-232`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:232-232`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:7-7`
<!-- vulcan:connections:end -->

## Limitations
Both branches read the same slot seven, so the property test only changes whether a short vector yields `NaN` or an error; a seven-element state whose seventh slot is not mass is misinterpreted without diagnostic.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/runtime.jl` line 5.
