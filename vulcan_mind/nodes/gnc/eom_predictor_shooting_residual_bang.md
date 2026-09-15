---
id: gnc.eom_predictor_shooting_residual_bang
label: shooting_residual!
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl
  symbol: shooting_residual!
  lines:
  - 735
  - 735
inputs:
- id: resid
  type: Any
  units: n/a
  required: true
  description: Positional argument `resid`.
- id: z
  type: Any
  units: n/a
  required: true
  description: Positional argument `z`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: param
  type: Any
  units: n/a
  required: true
  description: Positional argument `param`.
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
  type: Nothing
  units: n/a
  description: Return value of `shooting_residual!`; mutates `resid` in place. Returns
    `nothing`.
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

# shooting_residual!

## Purpose
Residual function for the indirect single-shooting problem solved by `nlsolve` in `asim_ctrl_targeting_plot`. Given a guess `z` for the three initial costates it propagates the augmented state through one drag pass with `f_ctrl!` and returns, in `resid`, the mismatch between the terminal costates and the transversality conditions for the terminal velocity and altitude weights.

## Theory & Math
Terminal transversality conditions with weight $v_E$ on final energy:
$$\lambda_v(t_f) = v_E\, v_f, \qquad \lambda_\gamma(t_f) = 0, \qquad \lambda_h(t_f) = v_E \frac{\mu}{r_f^2}$$
where $v_f = |v(t_f)|$, $r_f = |r(t_f)|$ and $\mu$ is `m.planet.μ`.

## Design & Implementation
`p` is the tuple `(r0, v0, hf, vf, γf, v_E)`; `param` is the runtime context. The initial condition is `[r0; v0; z[1:3]; 0.0]` (position m, velocity m/s, costates, heat load). It builds an `ODEProblem(f_ctrl!, in_cond, (time_0, 1500), param)` and solves with `Tsit5()`, `abstol=reltol=1e-9`, and the `out_drag_pass` callback. From the terminal state it computes `r_fin`, `v_fin`, `gamma_fin = asin(r·v/(|r||v|))` and the terminal costates. Residuals written in place are `resid[1] = λv_f - v_E v_f`, `resid[2] = λγ_f`, `resid[3] = λh_f - v_E μ / r_f^2`. `hf`, `vf`, `γf` and `gamma_fin` are unpacked but not used in the residuals.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `resid` | Any | n/a | yes | Positional argument `resid`. |
| in | `z` | Any | n/a | yes | Positional argument `z`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `param` | Any | n/a | yes | Positional argument `param`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `shooting_residual!`; mutates `resid` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:735-735`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1041-1041`
- `callees` → [[gnc.bridge_helpers__bridge_verbose_enabled|_bridge_verbose_enabled]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1037-1037`
- `callees` → [[gnc.bridge_helpers__make_aerobraking_runtime_context|_make_aerobraking_runtime_context]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1014-1014`
- `callees` → [[gnc.bridge_helpers__with_control_gain|_with_control_gain]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1013-1013`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1072-1072`
<!-- vulcan:connections:end -->

## Limitations
The final time is hard-coded to `1500` (not `time_0 + 1500` as in the final propagation later in the file), so if `time_0 > 1500` the time span is reversed. Each residual evaluation runs a full ODE solve, and `nlsolve` uses finite-difference Jacobians, so each Newton step costs roughly four propagations. Convergence depends on the hard-coded initial guess `z0 = [-500, 80000, 9]`. The `hf`, `vf`, `γf` targets are ignored, so the solver only enforces the energy-weighted transversality conditions. No check is made on `sol.retcode`.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl` line 735.
