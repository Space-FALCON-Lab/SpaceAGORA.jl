---
id: simulation.save_fields__save_latitude_deg
label: _save_latitude_deg
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_latitude_deg
  lines:
  - 99
  - 99
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
  description: Return value of `_save_latitude_deg`. Returns `latitudes_deg`.
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

# _save_latitude_deg

## Purpose
Save-time getter for planetodetic latitude in degrees for each spacecraft, giving the ground track its north-south coordinate in the saved output.

## Design & Implementation
Marked `@inline`. It computes the ephemeris time once as `et_start[] + Float64(t)`, then per spacecraft converts the inertial state into the planet-fixed frame through `r_intor_p!` and takes the second element of `rtolatlong(rp, planet, ephemerides_model)`, which is latitude in radians, converting it with `rad2deg`. Degrees are stored rather than radians because this field is an output-boundary quantity meant for plotting and comparison against mission products, while the internal callbacks work in radians throughout.

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
| out | `result` | Any | n/a | — | Return value of `_save_latitude_deg`. Returns `latitudes_deg`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:176-176`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:101-101`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:107-107`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:108-108`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:105-105`
<!-- vulcan:connections:end -->

## Limitations
The conversion is repeated per save field, so enabling altitude, latitude and longitude together triples the frame transformation work; nothing caches the intermediate `rp`. Whether the reported value is geodetic or geocentric is decided inside `rtolatlong` by the planet model and is not restated here. Near the poles, longitude becomes ill-conditioned while latitude saturates, and no special handling is applied.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 99.
