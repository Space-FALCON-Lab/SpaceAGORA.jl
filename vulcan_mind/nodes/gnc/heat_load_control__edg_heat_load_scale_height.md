---
id: gnc.heat_load_control__edg_heat_load_scale_height
label: _edg_heat_load_scale_height
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_heat_load_scale_height
  lines:
  - 3
  - 3
inputs:
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
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
  description: 'Return value of `_edg_heat_load_scale_height`. Returns `max(Float64(density_model.H),
    1.0)` or `hasproperty(planet, :H) ? max(Float64(planet.H), 1.0) : 7_000.0`.'
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _edg_heat_load_scale_height

## Purpose
Returns the atmospheric density scale height in metres used by the heat-load costate equations, taken from the density model if it exposes one, else the planet, else a fixed fallback.

## Design & Implementation
`@inline` accessor over `p::ODEParams`. It reads `p.args.environment_model.density_model` and `.planet`; if `hasproperty(density_model, :H)` it returns `max(Float64(density_model.H), 1.0)`, otherwise if the planet has `H` it returns `max(Float64(planet.H), 1.0)`, and otherwise the literal `7_000.0` m. The 1.0 m floor prevents division by zero in `_edg_heat_load_lambdas`, which divides by `scale_height`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_heat_load_scale_height`. Returns `max(Float64(density_model.H), 1.0)` or `hasproperty(planet, :H) ? max(Float64(planet.H), 1.0) : 7_000.0`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_heat_load_profile_for_k|_edg_heat_load_profile_for_k]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:559-559`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/heat_load_control.jl:7-7`
<!-- vulcan:connections:end -->

## Limitations
A single constant scale height is assumed across the whole drag pass even though real density profiles (and GRAM-based models) have altitude-dependent scale heights. The 7 km fallback is Earth-like and is wrong for Mars (~11 km) or Titan (~40 km) if neither the density model nor the planet carries `H`.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 3.
