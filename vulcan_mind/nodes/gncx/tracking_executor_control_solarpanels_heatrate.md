---
id: gncx.tracking_executor_control_solarpanels_heatrate
label: control_solarpanels_heatrate
kind: function
source:
  file: src/gnc/control/aerobraking/tracking_executor.jl
  symbol: control_solarpanels_heatrate
  lines:
  - 84
  - 180
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace supplying the heat-rate tracking executor together
    with the mission, argument dictionary, drag-passage indicator, and sampled thermal
    state.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: alpha_cmd
  type: Float64
  units: rad
  description: Commanded solar-panel angle of attack that holds the free-molecular
    heat rate at the configured limit.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncx
origin: agent
---
# control_solarpanels_heatrate

## Purpose
`control_solarpanels_heatrate` is the executor that converts a heat-rate limit into a commanded solar-panel angle of attack. Given the local temperature, density, and molecular speed ratio sampled from the environment, it returns the angle at which the free-molecular convective heat rate equals the mission heat-rate limit, so the panels shed as much energy as possible without exceeding the thermal constraint.

## Theory & Math
The commanded angle solves the free-molecular heat-rate constraint

$$\dot q(\alpha) = L\left[\left(S^2 + \frac{\gamma}{\gamma-1} - \frac{\gamma+1}{2(\gamma-1)}\frac{T_w}{T_p}\right)\left(e^{-(S\sin\alpha)^2} + \sqrt{\pi}\,S\sin\alpha\,\big(1+\mathrm{erf}(S\sin\alpha)\big)\right) - \tfrac12 e^{-(S\sin\alpha)^2}\right] = \dot q_{\max}$$

with $L = 10^{-4}\,\tau\,\rho R T_p \sqrt{R T_p / 2\pi}$ and $S$ the molecular speed ratio. The analytic derivative with respect to $\alpha$ is supplied to Newton's method rather than differenced.

## Model & Assumptions
The heat-rate model is the Maxwellian free-molecular expression built from the thermal accommodation factor, the gas constant, the specific-heat ratio, the wall-to-freestream temperature ratio, and the speed ratio projected onto the panel through `sin(alpha)`. Wall temperature is set equal to the freestream temperature. The routine first brackets the problem by evaluating the heat rate at the maximum and minimum admissible angles: if even the maximum angle stays under the limit the maximum is commanded, and if the minimum angle already exceeds the limit the minimum is commanded.

## Design & Implementation
Between those bounds the constraint is solved as a scalar root problem. Both the residual `f` and its analytic derivative `df` are formed in closed form and handed to `Roots.Newton()`, warm-started from the previous commanded angle held in the configuration state. Because the residual can admit multiple roots, a failed or out-of-range solve is retried from a seed chosen by comparing which bracket endpoint lies closer to the limit, and a final failure falls through `_control_exception_fallback` to the minimum angle. The result is clamped so an angle outside the physical range collapses to zero, and outside a drag passage the previously held angle is returned unchanged.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace supplying the heat-rate tracking executor together with the mission, argument dictionary, drag-passage indicator, and sampled thermal state. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `alpha_cmd` | Float64 | rad | — | Commanded solar-panel angle of attack that holds the free-molecular heat rate at the configured limit. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.constraint_tracking_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:180-180`
- [[gnc.control_commands_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:200-200`
- [[gnc.tracking_executor__control_solarpanels_openloop_impl|_control_solarpanels_openloop_impl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:278-278`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:180-180`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:200-200`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:155-155`
- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:85-85`
- `callees` → [[gnc.bridge_helpers__bridge_verbose_enabled|_bridge_verbose_enabled]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:134-134`
- `callees` → [[gnc.heat_rate_control_df|df]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:122-122`
- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:112-112`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:112-112`
- `callees` → [[gnc.tracking_executor__control_exception_fallback|_control_exception_fallback]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:147-147`
- `callees` → [[gnc.tracking_executor_df|df]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:122-122`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:112-112`
- `callees` → [[gnc.tracking_executor_heat_rate_calc|heat_rate_calc]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:104-104`
<!-- vulcan:connections:end -->

## Limitations
Newton iteration on a multi-rooted residual is sensitive to the warm start; the retry heuristic reduces but does not eliminate convergence to the wrong branch. Wall temperature is not modelled as a separate thermal state, so radiative equilibrium and soak-back are ignored. The commanded angle is returned but the panel rotation itself is applied elsewhere, so a caller that discards the return value silently loses the command.

## Provenance
Mapped from `src/gnc/control/aerobraking/tracking_executor.jl:84-180`.
