---
id: gnc.heat_load_control__edg_integrated_heat_load_trajectory
label: _edg_integrated_heat_load_trajectory
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_integrated_heat_load_trajectory
  lines:
  - 173
  - 173
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
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: pos0
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos0`.
- id: vel0
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel0`.
- id: mass
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: times
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `times`.
- id: alpha_profile
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `alpha_profile`.
- id: altitude_offset_m
  type: Float64
  units: n/a
  required: false
  description: Positional argument `altitude_offset_m` (default `0.0`).
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
  description: Return value of `_edg_integrated_heat_load_trajectory`. Returns `gravity`
    or `gravity + drag` or `(`.
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

# _edg_integrated_heat_load_trajectory

## Purpose
Numerically propagates the spacecraft through the drag pass under point-mass gravity and angle-of-attack-dependent drag, producing a higher-fidelity track than the closed-form approximation for the TPBVP switch solver.

## Design & Implementation
Signature `(config, p, spacecraft, pos0, vel0, mass, t, times, alpha_profile, altitude_offset_m = 0.0)`. It allocates `positions` and `velocities` as `Vector{SVector{3,Float64}}` and steps a classical fixed-step RK4 across `times`, using the nested `acceleration(r, v, tau, alpha)` closure with the alpha for interval `j` taken from `alpha_profile[min(j, end)]`. After propagation it recomputes altitude (`norm(r) - Rp_e + altitude_offset_m`), flight-path angle, speed, density, temperature and speed ratio at every node. Returns a NamedTuple with the seven track vectors plus `positions` and `velocities`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `pos0` | SVector{3, Float64} | n/a | yes | Positional argument `pos0`. |
| in | `vel0` | SVector{3, Float64} | n/a | yes | Positional argument `vel0`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `times` | Vector{Float64} | n/a | yes | Positional argument `times`. |
| in | `alpha_profile` | Vector{Float64} | n/a | yes | Positional argument `alpha_profile`. |
| in | `altitude_offset_m` | Float64 | n/a | no | Positional argument `altitude_offset_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_integrated_heat_load_trajectory`. Returns `gravity` or `gravity + drag` or `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_heat_load_profile_for_k|_edg_heat_load_profile_for_k]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:575-575`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_total_ref_area|_edg_total_ref_area]] · `callers` · call · `src/gnc/control/heat_load_control.jl:186-186`
<!-- vulcan:connections:end -->

## Limitations
The integrator is fixed-step RK4 on the caller's grid (about 1 s), with no error control, and evaluates the full panel aerodynamic model four times per step. Gravity is spherical point-mass with no J2, and lift is ignored entirely. `altitude_offset_m` is a constant correction that assumes the geodetic-minus-spherical altitude difference is unchanged along the pass.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 173.
