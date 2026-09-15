---
id: gnc.eom_predictor_f_ctrl_bang
label: f_ctrl!
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl
  symbol: f_ctrl!
  lines:
  - 67
  - 67
inputs:
- id: y_dot
  type: Any
  units: n/a
  required: true
  description: Positional argument `y_dot`.
- id: in_cond
  type: Any
  units: n/a
  required: true
  description: Positional argument `in_cond`.
- id: param
  type: Any
  units: n/a
  required: true
  description: Positional argument `param`.
- id: t0
  type: Any
  units: n/a
  required: true
  description: Positional argument `t0`.
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
  description: Return value of `f_ctrl!`; mutates `y_dot` in place.
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

# f_ctrl!

## Purpose
In-place right-hand side of the 10-state augmented aerobraking system integrated by `asim_ctrl_targeting_plot`. The state `in_cond` holds inertial position (1:3, m), inertial velocity (4:6, m/s), three costates `lambdav`, `lambdagamma`, `lambdah` (7:9) and accumulated heat load (10). It writes translational dynamics, costate dynamics and heat rate into `y_dot` so the shooting solver can propagate an optimal-control bang-bang angle-of-attack profile through one drag pass.

## Theory & Math
State $y = [r_{ii}, v_{ii}, \lambda_v, \lambda_\gamma, \lambda_h, Q]$. With density $\rho$, reference area $A$, mass $m$, gravity $g = |F_g|/m$, inertial speed $v$, radius $r$, flight-path angle $\gamma$, scale height $H$ and angle of attack $\alpha$:
$$\dot\lambda_v = -\frac{3\rho v^2 \alpha}{\pi} + \lambda_v\frac{\rho A C_D v}{m} - \lambda_\gamma\left(\frac{\rho A C_L}{2m} + \frac{g}{v^2} + \frac{1}{r}\right) - \lambda_h \gamma$$
$$\dot\lambda_\gamma = \lambda_v g - \lambda_h v$$
$$\dot\lambda_h = \frac{\rho v^3 \alpha}{\pi H} - \lambda_v\left(\frac{\rho A C_D v^2}{2mH} + \frac{2 g \gamma}{r}\right) + \lambda_\gamma\left(\frac{\rho A C_L v}{2mH} - \frac{2g}{r v} + \frac{v}{r^2}\right)$$
Switching: $\alpha = 0$ if $\lambda_v < \lambda_{sw} = \dfrac{2 m v}{A\, C_{D,slope}\, \pi}$, else $\alpha = \alpha_{max}$, where $C_{D,slope} = (C_D(\pi/2) - C_D(0))/(\pi/2)$.

## Design & Implementation
`param` is the aerobraking runtime context (mission, `ip`, `args`, GRAM handles, `control_gain`, `settings`). Each call rebuilds the epoch from `date_initial + t0*seconds`, calls SPICE `utc2et` and `pxform("J2000", "IAU_<planet>")`, and mutates `m.planet.L_PI` with that rotation. Inertial state is converted to planet-relative via `r_intor_p!`, then to altitude/lat/lon via `rtolatlong`. Density model selection follows `ip.dm` (0 constant, 1 exponential, 2 none, 3 GRAM via `density_gram` with `pyconvert`). The angle of attack is a bang-bang switch: `aoa = 0` when `lambdav_ii < lambda_switch` where `lambda_switch = 2 m v_ii / (A CD_slope pi)`, otherwise `m.aerodynamics.α`. When `heat_rate_control` is true and the Maxwellian convective heat rate exceeds `settings.max_heat_rate`, `_control_solarpanels_heatrate` (and optionally `control_struct_load`) reduce `aoa`. Forces are drag and lift from `aerodynamic_coefficient_fM` at bank angle 0, gravity from `aerobraking_gravity_force_ii`, and SRP (4.56e-6 N/m^2 at 1 AU) if `settings.srp_enabled`. Costate derivatives use the closed-form expressions for an exponential atmosphere with scale height `m.planet.H`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `y_dot` | Any | n/a | yes | Positional argument `y_dot`. |
| in | `in_cond` | Any | n/a | yes | Positional argument `in_cond`. |
| in | `param` | Any | n/a | yes | Positional argument `param`. |
| in | `t0` | Any | n/a | yes | Positional argument `t0`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `f_ctrl!`; mutates `y_dot` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:51-51`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:55-55`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:67-67`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:238-238`
- `callees` → [[core.reference_system_config_clock|clock]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:91-91`
- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:155-155`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:116-116`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:149-149`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:264-264`
- `callees` → [[dynx.coupled_perturbations_srp|srp]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:256-256`
- `callees` → [[environment.gravity_models_aerobraking_gravity_force_ii|aerobraking_gravity_force_ii]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:236-236`
- `callees` → [[gnc.guidance_hooks__control_solarpanels_heatrate|_control_solarpanels_heatrate]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:211-211`
- `callees` → [[gnc.tracking_executor_control_struct_load|control_struct_load]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:214-214`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:121-121`
<!-- vulcan:connections:end -->

## Limitations
This closure captures `CD_slope`, `MonteCarlo`, `wind_m`, `heat_rate_control` and `cnf_state` from the enclosing scope, so it is not reusable outside `asim_ctrl_targeting_plot`. `m.planet.L_PI` is mutated on every derivative evaluation, which is a shared-state hazard if the mission object is used concurrently. The costate equations assume an exponential atmosphere with constant scale height `m.planet.H`, which is inconsistent with the GRAM (`ip.dm == 3`) density branch. The relative wind norm is never guarded against zero, so a stationary planet-relative state would produce NaNs. Bank angle is hard-coded to 0 and thrust is omitted entirely.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl` line 67.
