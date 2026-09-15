---
id: gnc.trajectory_predictor_f_ctrl_rf_bang
label: f_ctrl_rf!
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl
  symbol: f_ctrl_rf!
  lines:
  - 25
  - 25
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
  description: Return value of `f_ctrl_rf!`; mutates `y_dot` in place. Returns `y_dot`.
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

# f_ctrl_rf!

## Purpose
Non-dimensionalised right-hand side of the T-EDG targeting propagation. It integrates inertial position and velocity, vehicle mass and accumulated heat load through a drag passage while applying the heat-rate-limited solar-panel attitude law, and is the dynamics the switch-time targeting shooting method is built on.

## Theory & Math
Canonical scaling uses length unit $DU$, time unit $TU$ and mass unit $MU$, so $\dot{\bar r} = v\,TU/DU$, $\dot{\bar v} = (F/m)\,TU^2/DU$ and $\dot{\bar m} = -\|T\|/(g_e I_{sp})\cdot TU/MU$ with $g_e$ standard gravity and $I_{sp}$ specific impulse in seconds. The ideal velocity increment is $\Delta v = g_e I_{sp}\ln(m_0/m)$. Knudsen number is estimated as $Kn = 1.26\sqrt{\gamma}\,M/(Re + 10^{-5})$ to flag departure from free-molecular flow below $Kn = 0.1$.

## Design & Implementation
`f_ctrl_rf!(y_dot, in_cond, param, t0)` first rescales the independent variable as `t0 = t0 * cnf_state.TU`, so the solver works in canonical time units while the physics is evaluated in seconds. It increments three counters in `cnf_state` (`count_aerobraking`, `count_dori`, `count_phase`) on every derivative evaluation. Atmosphere selection follows `ip.dm`, aerodynamic coefficients follow `ip.am` across `aerodynamic_coefficient_constant`, `aerodynamic_coefficient_fM` and `aerodynamic_coefficient_no_ballistic_flight`, and the thrust profile follows `ip.tc` across `no_maneuver`, `abms` and `deceleration_drag_passage`. The attitude law is the switch test `t0 >= t_switch`: past the switch it sets `cnf_state.α = 0`, before it calls `_control_solarpanels_heatrate` with the elapsed time `t0 - cnf_state.t_switch_targeting`; when `settings.struct_control_enabled` the result is further reduced by `min(cnf_state.α, α_struct)` from `control_struct_load`. Thrust direction uses Rodrigues rotation of the drag unit vector about the drag-lift normal through `settings.thrust_phi`. Derivatives are written non-dimensionally: `y_dot[1:3] = vel_ii*(TU/DU)`, `y_dot[4:6] = force_ii/mass*(TU^2/DU)`, `y_dot[7] = -norm(thrust_ii)/(g_e*Isp)*TU/MU` for mass depletion, and `y_dot[8] = heat_rate*TU^3/MU`. Propellant exhaustion is detected when `mass - settings.dry_mass <= 0.5` kg, which zeroes `m.engines.T`.

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
| out | `result` | Any | n/a | — | Return value of `f_ctrl_rf!`; mutates `y_dot` in place. Returns `y_dot`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:25-25`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:288-288`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:204-204`
- `callees` → [[core.reference_system_config_clock|clock]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:56-56`
- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:160-160`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:81-81`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:129-129`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_constant|aerodynamic_coefficient_constant]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:316-316`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:318-318`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_no_ballistic_flight|aerodynamic_coefficient_no_ballistic_flight]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:320-320`
- `callees` → [[dynx.coupled_perturbations_srp|srp]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:306-306`
- `callees` → [[environment.gravity_models_aerobraking_gravity_force_ii|aerobraking_gravity_force_ii]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:286-286`
- `callees` → [[gnc.bridge_helpers__bridge_verbose_enabled|_bridge_verbose_enabled]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:203-203`
- `callees` → [[gnc.guidance_hooks__control_solarpanels_heatrate|_control_solarpanels_heatrate]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:264-264`
- `callees` → [[gnc.tracking_executor_control_struct_load|control_struct_load]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:254-254`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:90-90`
- `callees` → [[vehicle.thermal_models_heatrate_convective_radiative|heatrate_convective_radiative]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:228-228`
<!-- vulcan:connections:end -->

## Limitations
The function writes `cnf_state.α` and three call counters on every derivative evaluation, so the attitude actually integrated depends on the solver's internal stage ordering and rejected steps, and the counters count evaluations rather than steps. That same shared state makes concurrent propagation of two vehicles unsafe. Setting `m.engines.T = 0` on propellant exhaustion permanently mutates the mission object, so a subsequent pass starts with a dead engine even after a state reset. The `energy` quantity is computed and discarded. Zeroing the angle of attack to exactly `0` after the switch differs from the `0.0001` rad floor used elsewhere in the aerobraking control path.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl` line 25.
