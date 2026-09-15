---
id: gnc.constraint_tracking_f_ctrl_bang
label: f_ctrl!
kind: function
source:
  file: src/gnc/control/aerobraking/constraint_tracking.jl
  symbol: f_ctrl!
  lines:
  - 51
  - 51
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
  description: Return value of `f_ctrl!`; mutates `y_dot` in place. Returns `y_dot`.
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
In-place right-hand side for the optimal-control aerobraking pass: it propagates the inertial translational state together with three costates and the accumulated heat load, applying the bang-bang panel-angle switching law that the costates imply.

## Theory & Math
$$\lambda_{\mathrm{switch}} = \frac{2 k_{cf} m \lVert v_{ii} \rVert}{A_{\mathrm{tot}}\, C_{D}' \pi}, \qquad C_D' = \frac{C_D(\pi/2) - C_D(0)}{\pi/2}$$ with $k_{cf}$ the control gain, $m$ the vehicle mass in $\mathrm{kg}$, $A_{\mathrm{tot}}$ the total reference area in $\mathrm{m^2}$. The costate dynamics integrated are $$\dot\lambda_v = -\frac{3 k_{cf} \rho v^2 \alpha}{\pi} + \lambda_v \frac{\rho A C_D v}{m} - \lambda_\gamma\left(\frac{\rho A C_L}{2m} + \frac{g}{v^2} + \frac{1}{r}\right) - \lambda_h \gamma,$$ $$\dot\lambda_\gamma = \lambda_v g - \lambda_h v, \qquad \dot\lambda_h = \frac{k_{cf}\rho v^3 \alpha}{\pi H} - \lambda_v\left(\frac{\rho A C_D v^2}{2 m H} + \frac{2 g \gamma}{r}\right) + \lambda_\gamma\left(\frac{\rho A C_L v}{2 m H} - \frac{2 g}{r v} + \frac{v}{r^2}\right),$$ where $\rho$ is density in $\mathrm{kg/m^3}$, $H$ the atmospheric scale height in $\mathrm{m}$, $g$ the local gravitational acceleration magnitude in $\mathrm{m/s^2}$, $r$ the inertial radius in $\mathrm{m}$ and $\gamma$ the inertial flight path angle in radians.

## Design & Implementation
Unpacks the runtime context from `param`, builds the clock from `date_initial + t0*seconds`, and splits `in_cond` into position 1:3, velocity 4:6, costates `lambdav_ii`, `lambdagamma_ii`, `lambdah_ii` at indices 7 to 9. It transforms to the planet-relative frame with `r_intor_p`, computes latitude, longitude and altitude, and selects the atmosphere model by `ip.dm`, dispatching to constant, exponential, none, or GRAM density. The Mach number gives the molecular speed ratio. The switching threshold `lambda_switch` is compared against `lambdav_ii`, and depending on `settings.heat_load_solution` being 0 or 1 the angle of attack is set to either 0.0001 radians or `m.aerodynamics.α`. If the resulting convective heat rate exceeds `settings.max_heat_rate`, the angle is replaced by the solar-panel heat-rate solve and the rate is pinned to the limit. Lift, drag, gravity and optional solar radiation pressure sum into `force_ii`, and the function writes `y_dot[1:3] = vel_ii`, `y_dot[4:6] = force_ii / mass`, the three costate derivatives into slots 7 to 9, and the instantaneous heat rate into slot 10, mutating `y_dot` in place.

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
| out | `result` | Any | n/a | — | Return value of `f_ctrl!`; mutates `y_dot` in place. Returns `y_dot`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:51-51`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:55-55`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:67-67`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:205-205`
- `callees` → [[core.reference_system_config_clock|clock]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:69-69`
- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:126-126`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:120-120`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:233-233`
- `callees` → [[dynx.coupled_perturbations_srp|srp]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:224-224`
- `callees` → [[environment.gravity_models_aerobraking_gravity_force_ii|aerobraking_gravity_force_ii]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:203-203`
- `callees` → [[gncx.tracking_executor_control_solarpanels_heatrate|control_solarpanels_heatrate]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:180-180`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:92-92`
<!-- vulcan:connections:end -->

## Limitations
The closure references `et` and `el_time` which are not bound in its own scope, so the SPICE rotation `pxform` call and the GRAM density branch depend on outer-scope definitions and will raise `UndefVarError` where they are absent. The local `heat_rate_control = true` assignment overrides the argument of the same name, making the limiter unconditional. `CD_slope` is captured from the enclosing function, evaluated once at the initial temperature and speed ratio rather than at the current state, so the switching threshold is only locally valid. The bang-bang angle of 0.0001 radians is a magic value standing in for zero to avoid a coefficient-model singularity, and the near-discontinuous switch makes the adaptive step controller work hard without an explicit switching callback.

## Provenance
Mapped from `src/gnc/control/aerobraking/constraint_tracking.jl` line 51.
