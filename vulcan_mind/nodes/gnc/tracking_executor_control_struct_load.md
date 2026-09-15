---
id: gnc.tracking_executor_control_struct_load
label: control_struct_load
kind: function
source:
  file: src/gnc/control/aerobraking/tracking_executor.jl
  symbol: control_struct_load
  lines:
  - 31
  - 31
inputs:
- id: ip
  type: Any
  units: n/a
  required: true
  description: Positional argument `ip`.
- id: m
  type: Any
  units: n/a
  required: true
  description: Positional argument `m`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: S
  type: Any
  units: n/a
  required: true
  description: Positional argument `S`.
- id: T_p
  type: Any
  units: n/a
  required: true
  description: Positional argument `T_p`.
- id: q
  type: Any
  units: n/a
  required: true
  description: Positional argument `q`.
- id: MonteCarlo
  type: Any
  units: n/a
  required: false
  description: Positional argument `MonteCarlo` (default `false`).
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
  description: Return value of `control_struct_load`. Returns `α`.
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

# control_struct_load

## Purpose
Structural-load limiter for the solar panels: chooses the panel angle of attack so that aerodynamic drag does not exceed the force implied by the configured maximum dynamic pressure, using bisection on the drag-versus-angle curve.

## Design & Implementation
Arguments are `ip`, the mission model `m`, `args`, speed ratio `S`, plasma/atmosphere temperature `T_p`, dynamic pressure `q` (Pa) and a `MonteCarlo` flag forwarded to `aerodynamic_coefficient_fM`. It sets `max_α = m.aerodynamics.α`, `min_α = 0.0001` rad, total area from `config.get_spacecraft_reference_area(m.body)`, and the 90-degree drag coefficient `CD90`. `drag_limit = max_dyn_press * CD90 * area_tot` with `max_dyn_press` from `_bridge_required_field(args, :max_dyn_press)`. If `drag_max < drag_limit` the full angle is returned; if `drag_min > drag_limit` the minimum; otherwise `find_zero(f, (0, pi/2), Roots.Bisection())` solves `q*CD(x)*area - drag_limit = 0`, with failures routed through `_control_exception_fallback` to `min_α`. A final guard sets `α = 0` if outside `[0, max_α]`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ip` | Any | n/a | yes | Positional argument `ip`. |
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `S` | Any | n/a | yes | Positional argument `S`. |
| in | `T_p` | Any | n/a | yes | Positional argument `T_p`. |
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `MonteCarlo` | Any | n/a | no | Positional argument `MonteCarlo` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `control_struct_load`. Returns `α`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:214-214`
- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:254-254`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:214-214`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:254-254`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:38-38`
- `callees` → [[gnc.bridge_helpers__bridge_required_field|_bridge_required_field]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:42-42`
<!-- vulcan:connections:end -->

## Limitations
The bisection bracket is hard-coded to `(0, pi/2)` regardless of `max_α`, so a root above `max_α` is clamped to 0, not to `max_α`. The `else` branch after the three conditions is unreachable for finite values but leaves `α` unassigned for NaN inputs, causing an `UndefVarError`. `min_α` is a hard-coded 1e-4 rad. The panel-rotation block is commented out, so this function only returns the angle.

## Provenance
Mapped from `src/gnc/control/aerobraking/tracking_executor.jl` line 31.
