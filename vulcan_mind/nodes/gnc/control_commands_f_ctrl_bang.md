---
id: gnc.control_commands_f_ctrl_bang
label: f_ctrl!
kind: function
source:
  file: src/gnc/control/aerobraking/control_commands.jl
  symbol: f_ctrl!
  lines:
  - 55
  - 55
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
In-place right-hand side of the ten-state aerobraking optimal-control problem integrated by `asim_ctrl`. It advances inertial position and velocity, the three costates of the energy-depletion problem, and accumulated heat load, applying the bang-bang angle-of-attack law at every evaluation.

## Theory & Math
The costate dynamics integrated are $\dot\lambda_v = -3k\rho v^2\alpha/\pi + \lambda_v \rho A C_D v/m - \lambda_\gamma\left(\frac{\rho A C_L}{2m} + \frac{g}{v^2} + \frac1r\right) - \lambda_h\gamma$, $\dot\lambda_\gamma = \lambda_v g - \lambda_h v$, and $\dot\lambda_h = k\rho v^3\alpha/(\pi H) - \lambda_v\left(\frac{\rho A C_D v^2}{2mH} + \frac{2g\gamma}{r}\right) + \lambda_\gamma\left(\frac{\rho A C_L v}{2mH} - \frac{2g}{rv} + \frac{v}{r^2}\right)$, with $k$ the control gain, $\rho$ density, $A$ reference area, $m$ mass, $H$ the atmospheric scale height, $r$ inertial radius, $v$ inertial speed and $g$ the local gravity magnitude. The bang-bang switching surface is $\lambda_v = 2 k m v/(A\,\pi\,dC_D/d\alpha)$.

## Design & Implementation
`f_ctrl!(y_dot, in_cond, param, t0)` writes into `y_dot` and reads the closure-captured `cnf_state`, `CD_slope`, `time_switch_eval`, `time_switch_2` and `MonteCarlo` flags plus the `param` context (mission, `ip`, `args`, `date_initial`, `control_gain`, `settings`). It rebuilds the epoch as `date_initial + t0*seconds`, calls `pxform("J2000", "IAU_"*uppercase(planet.name), et)` and mutates `m.planet.L_PI` in place with the resulting 3x3 matrix. State slices are `in_cond[1:3]` position (m), `in_cond[4:6]` velocity (m/s), and `in_cond[7:9]` the costates `lambdav`, `lambdagamma`, `lambdah`. Density, temperature and wind come from one of `density_constant`, `density_exp`, `density_no` or `density_gram` selected by `ip.dm`, with the GRAM branch converting Python values through `pyconvert`. Molecular speed ratio is `S = sqrt(gamma/2)*Mach`. Angle of attack is chosen by bang-bang: in switch-evaluation mode against the threshold `lambda_switch = k_cf*2*mass*vel_ii_mag/(area_tot*CD_slope*pi)`, otherwise by whether `t0` lies inside `[cnf_state.time_switch_1, time_switch_2]`; the sense of the comparison flips with `settings.heat_load_solution`, and the low-drag attitude is the literal `0.0001` rad rather than zero. When `heat_rate_control` is set and the Maxwellian heat rate exceeds `settings.max_heat_rate`, `control_solarpanels_heatrate` overrides the attitude and the heat rate is pinned to the limit. Forces are drag and lift from `aerodynamic_coefficient_fM`, gravity from `aerobraking_gravity_force_ii`, and optional SRP at `p_srp_unscaled = 4.56e-6` N/m^2. The final writes are `y_dot[1:3] = vel_ii`, `y_dot[4:6] = force_ii/mass`, the three costate rates, and `y_dot[10] = heat_rate`.

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

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:225-225`
- `callees` → [[core.reference_system_config_clock|clock]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:74-74`
- `callees` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:138-138`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:99-99`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:132-132`
- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:253-253`
- `callees` → [[dynx.coupled_perturbations_srp|srp]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:244-244`
- `callees` → [[environment.gravity_models_aerobraking_gravity_force_ii|aerobraking_gravity_force_ii]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:223-223`
- `callees` → [[gncx.tracking_executor_control_solarpanels_heatrate|control_solarpanels_heatrate]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:200-200`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:104-104`
<!-- vulcan:connections:end -->

## Limitations
The function mutates `m.planet.L_PI` in place on every derivative evaluation, so the mission object is shared mutable state and two satellites integrated concurrently on the same mission corrupt each other's rotation matrix; the SPICE `pxform` call is additionally repeated a few lines later into a local `L_PI`, doubling the kernel lookups per step. The low-drag attitude is `0.0001` rad rather than `0.0`, a magic value that leaks a small residual drag into the supposedly minimum-drag arc. `settings.heat_load_solution` values `2` and `3` are unhandled in the switch-evaluation branch, leaving `aoa` undefined and raising an `UndefVarError`. Pinning `heat_rate = settings.max_heat_rate` after the attitude override records the limit rather than the achieved rate, biasing the integrated heat load. The GRAM path calls into Python from inside the ODE right-hand side, which is neither thread-safe nor cheap.

## Provenance
Mapped from `src/gnc/control/aerobraking/control_commands.jl` line 55.
