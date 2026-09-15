---
id: environment.gravity_models_aerobraking_gravity_force_ii
label: aerobraking_gravity_force_ii
kind: function
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: aerobraking_gravity_force_ii
  lines:
  - 143
  - 143
inputs:
- id: gm_code
  type: Integer
  units: n/a
  required: true
  description: Positional argument `gm_code`.
- id: mass
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass`.
- id: pos_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_ii`.
- id: vel_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel_ii`.
- id: pos_pp
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_pp`.
- id: lat
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lat`.
- id: lon
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lon`.
- id: alt
  type: Float64
  units: n/a
  required: true
  description: Positional argument `alt`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: gram_atmosphere
  type: Any
  units: n/a
  required: true
  description: Positional argument `gram_atmosphere`.
- id: gram
  type: Any
  units: n/a
  required: true
  description: Positional argument `gram`.
- id: n_bodies_list
  type: Any
  units: n/a
  required: true
  description: Positional argument `n_bodies_list`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `aerobraking_gravity_force_ii`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# aerobraking_gravity_force_ii

## Purpose
Gravity force assembly for the aerobraking guidance predictors (`f_ctrl!` in `eom_predictor.jl` and related control code), which use a legacy integer gravity-model code and loosely typed `args` rather than the effector model types. It combines central or J2 gravity with optional N-body perturbations and spherical-harmonics acceleration.

## Design & Implementation
Selects `_inverse_squared_j2_gravity_accel(pos_ii, planet)` when `gm_code == 2`, otherwise `_inverse_squared_gravity_accel` (the comment notes legacy codes 0, 1 and 3 all reduce to inverse-square). Multiplies by `mass` (kg) to get force. If `_gravity_runtime_field(args, :n_bodies, ())` is non-empty it adds `mass * PerturbationEffectors.gravity_n_bodies(et, pos_ii, planet, n_bodies_list[k])` for every `k`. If `_gravity_runtime_field(args, :gravity_harmonics, 0) == 1` it reads degree `L` and order `M`, evaluates `acc_gravity_pines!(pos_pp, planet.Clm, planet.Slm, L, M, μ, Rp_e, planet)` in the planet-fixed frame and rotates with `planet.L_PI'`. Arguments `vel_ii`, `lat`, `lon`, `alt`, `gram_atmosphere` and `gram` are accepted but unused. Returns an `SVector{3,Float64}` in newtons.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `gm_code` | Integer | n/a | yes | Positional argument `gm_code`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `pos_ii` | SVector{3, Float64} | n/a | yes | Positional argument `pos_ii`. |
| in | `vel_ii` | SVector{3, Float64} | n/a | yes | Positional argument `vel_ii`. |
| in | `pos_pp` | SVector{3, Float64} | n/a | yes | Positional argument `pos_pp`. |
| in | `lat` | Float64 | n/a | yes | Positional argument `lat`. |
| in | `lon` | Float64 | n/a | yes | Positional argument `lon`. |
| in | `alt` | Float64 | n/a | yes | Positional argument `alt`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `gram_atmosphere` | Any | n/a | yes | Positional argument `gram_atmosphere`. |
| in | `gram` | Any | n/a | yes | Positional argument `gram`. |
| in | `n_bodies_list` | Any | n/a | yes | Positional argument `n_bodies_list`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `aerobraking_gravity_force_ii`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.constraint_tracking_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:203-203`
- [[gnc.control_commands_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:223-223`
- [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:236-236`
- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:286-286`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:203-203`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:223-223`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:236-236`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:286-286`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:203-203`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:236-236`

**Downstream**

- `callees` → [[environment.gravity_models__gravity_runtime_field|_gravity_runtime_field]] · `callers` · call · `src/environment/gravity/gravity_models.jl:168-168`
- `callees` → [[environment.gravity_models__inverse_squared_gravity_accel|_inverse_squared_gravity_accel]] · `callers` · call · `src/environment/gravity/gravity_models.jl:163-163`
- `callees` → [[environment.gravity_models__inverse_squared_j2_gravity_accel|_inverse_squared_j2_gravity_accel]] · `callers` · call · `src/environment/gravity/gravity_models.jl:160-160`
<!-- vulcan:connections:end -->

## Limitations
The J2 branch is evaluated on the inertial position `pos_ii`, so the oblateness term is aligned with inertial z, not the spin axis, unlike the sampled `wrench` path; for a planet whose pole is not along J2000 z this is an error. `planet.L_PI` must have been updated to the current epoch by the caller (the predictor mutates it each step). The harmonics loop iterates `eachindex(n_bodies)` but indexes `n_bodies_list[k]`, assuming both collections align. Six of the fourteen arguments are dead.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 143.
