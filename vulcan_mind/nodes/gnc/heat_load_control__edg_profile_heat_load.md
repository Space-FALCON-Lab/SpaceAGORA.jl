---
id: gnc.heat_load_control__edg_profile_heat_load
label: _edg_profile_heat_load
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_profile_heat_load
  lines:
  - 422
  - 422
inputs:
- id: config
  type: Any
  units: n/a
  required: true
  description: Positional argument `config`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: track
  type: Any
  units: n/a
  required: true
  description: Positional argument `track`.
- id: alpha_profile
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `alpha_profile`.
- id: heat_rate_control
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `heat_rate_control`.
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
  description: Return value of `_edg_profile_heat_load`. Returns `_edg_integrate_series(track.time,
    qdot)`.
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

# _edg_profile_heat_load

## Purpose
Convenience composition that computes the total heat load over a track for one alpha profile by integrating the heat-rate history.

## Design & Implementation
Signature `(config, p, track, alpha_profile; heat_rate_control::Bool)`. Calls `_edg_profile_heat_rates` with the same keyword and passes the result together with `track.time` to `_edg_integrate_series`. Returns a `Float64` in J/cm^2. It is the objective evaluated inside the `residual` closure of `_edg_solve_heat_load_switches`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `track` | Any | n/a | yes | Positional argument `track`. |
| in | `alpha_profile` | Vector{Float64} | n/a | yes | Positional argument `alpha_profile`. |
| in | `heat_rate_control` | Bool | n/a | yes | Keyword argument `heat_rate_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_profile_heat_load`. Returns `_edg_integrate_series(track.time, qdot)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control_residual|residual]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:677-677`
- [[gnc.targeting_control__edg_targeting_outcome_with_heat_load|_edg_targeting_outcome_with_heat_load]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:873-873`
- [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:677-677`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_integrate_series|_edg_integrate_series]] · `callers` · call · `src/gnc/control/heat_load_control.jl:424-424`
- `callees` → [[gnc.heat_load_control__edg_profile_heat_rates|_edg_profile_heat_rates]] · `callers` · call · `src/gnc/control/heat_load_control.jl:423-423`
<!-- vulcan:connections:end -->

## Limitations
Recomputes the full heat-rate vector on every call, allocating `length(track.time)` doubles; the root solver calls it dozens of times per switch solve. Inherits the single-node, single-scalar limitations of `_edg_profile_heat_rates`.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 422.
