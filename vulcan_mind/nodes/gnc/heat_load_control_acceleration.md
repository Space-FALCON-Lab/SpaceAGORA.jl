---
id: gnc.heat_load_control_acceleration
label: acceleration
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: acceleration
  lines:
  - 193
  - 193
inputs:
- id: r
  type: Any
  units: n/a
  required: true
  description: Positional argument `r`.
- id: v
  type: Any
  units: n/a
  required: true
  description: Positional argument `v`.
- id: tau
  type: Any
  units: n/a
  required: true
  description: Positional argument `tau`.
- id: alpha
  type: Any
  units: n/a
  required: true
  description: Positional argument `alpha`.
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
  description: Return value of `acceleration`. Returns `gravity` or `gravity + drag`.
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

# acceleration

## Purpose
Nested closure inside `_edg_integrated_heat_load_trajectory` that evaluates the translational acceleration (gravity plus drag) at a trial position, velocity, elapsed time and angle of attack for each RK4 stage.

## Theory & Math
$$\mathbf{a} = -\frac{\mu}{|\mathbf{r}|^3}\mathbf{r} - \frac{\rho\,|\mathbf{v}|\,C_D(\alpha, S, T)\,A}{2m}\,\mathbf{v}$$ where $\rho$ is density (kg/m^3), $A$ the total reference area (m^2), $m$ the mass (kg), and $C_D$ the area-weighted drag coefficient.

## Design & Implementation
Captures `planet`, `p`, `spacecraft`, `area`, `mass`, `t`, and `altitude_offset_m`. Computes `gravity = -μ r / |r|^3`, samples density and temperature at `altitude = |r| - Rp_e + altitude_offset_m` via `_edg_sample_prediction_atmosphere(p, altitude, t + tau)`, and returns gravity alone when density or speed is non-positive. Otherwise it forms the molecular speed ratio `S = sqrt(γ/2) |v| / sqrt(γ R T)`, calls `_edg_weighted_aero_coefficients(spacecraft, T, S, alpha)` for CD, and adds `drag = -0.5 ρ |v|^2 CD A / m * (v / |v|)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r` | Any | n/a | yes | Positional argument `r`. |
| in | `v` | Any | n/a | yes | Positional argument `v`. |
| in | `tau` | Any | n/a | yes | Positional argument `tau`. |
| in | `alpha` | Any | n/a | yes | Positional argument `alpha`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `acceleration`. Returns `gravity` or `gravity + drag`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_integrated_max_energy_depletion_trajectory|_edg_integrated_max_energy_depletion_trajectory]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:585-585`
- [[gnc.targeting_control_max_energy_alpha|max_energy_alpha]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:624-624`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_sample_prediction_atmosphere|_edg_sample_prediction_atmosphere]] · `callers` · call · `src/gnc/control/heat_load_control.jl:197-197`
- `callees` → [[gnc.heat_load_control__edg_weighted_aero_coefficients|_edg_weighted_aero_coefficients]] · `callers` · call · `src/gnc/control/heat_load_control.jl:204-204`
<!-- vulcan:connections:end -->

## Limitations
Drag is applied along the inertial velocity vector, ignoring atmospheric co-rotation and winds. Lift from `_edg_weighted_aero_coefficients` is discarded. Because the closure captures `mass` as a constant, propellant use during the pass is not modelled.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 193.
