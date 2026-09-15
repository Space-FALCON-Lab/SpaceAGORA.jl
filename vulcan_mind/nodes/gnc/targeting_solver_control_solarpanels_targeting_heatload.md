---
id: gnc.targeting_solver_control_solarpanels_targeting_heatload
label: control_solarpanels_targeting_heatload
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl
  symbol: control_solarpanels_targeting_heatload
  lines:
  - 139
  - 139
inputs:
- id: energy_f
  type: Any
  units: n/a
  required: true
  description: Positional argument `energy_f`.
- id: param
  type: Any
  units: n/a
  required: true
  description: Positional argument `param`.
- id: OE
  type: Any
  units: n/a
  required: true
  description: Positional argument `OE`.
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
  description: Return value of `control_solarpanels_targeting_heatload`. Returns `(energy_fin
    - energy_f) / 1e6` or `v_E_fin`.
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

# control_solarpanels_targeting_heatload

## Purpose
Solves for the heat-rate-constrained control parameter `v_E` such that a pass propagated by `_control_asim_ctrl_rf` finishes at the target specific energy `energy_f`. This variant of solar-panel targeting is used when the angle-of-attack profile is shaped by a heat-load limit rather than a single switch time.

## Design & Implementation
Defines the closure `func_targeting_heatload(v_E)` which unpacks `m = param.mission`, `ip = param.ip`, `time_0`, `args` and `gram_atmosphere` from `param`, then calls `_control_asim_ctrl_rf(ip, m, time_0, OE, args, v_E, 1.0, false, gram_atmosphere; cnf=_bridge_get_cnf(param))`, discarding the second return value. The final energy `|v|^2/2 - μ/|r|` is taken from the last column of `sol`, and the residual `(energy_fin - energy_f)/1e6` is returned. `find_zero` with `Roots.Brent()` on the bracket `[1, 1000]` and `rtol=1e-8` yields `v_E_fin`, which is returned. Verbose printing of each `(v_E, energy_fin)` pair is gated on `_bridge_verbose_enabled(param.args)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `energy_f` | Any | n/a | yes | Positional argument `energy_f`. |
| in | `param` | Any | n/a | yes | Positional argument `param`. |
| in | `OE` | Any | n/a | yes | Positional argument `OE`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `control_solarpanels_targeting_heatload`. Returns `(energy_fin - energy_f) / 1e6` or `v_E_fin`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl`

**Downstream**

- `callees` → [[gnc.bridge_helpers__bridge_verbose_enabled|_bridge_verbose_enabled]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:140-140`
<!-- vulcan:connections:end -->

## Limitations
The bracket `[1, 1000]` and the fixed positional arguments `1.0, false` to `_control_asim_ctrl_rf` are hard-coded with no documentation of their meaning or units. A non-bracketing residual causes `find_zero` to throw. The tolerance `rtol=1e-8` on `v_E` is far tighter than the `1e-5` used by the numerical-integration variant, so this solver runs more full trajectory propagations per call. `OE` is passed through untyped and must match what `_control_asim_ctrl_rf` expects.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl` line 139.
