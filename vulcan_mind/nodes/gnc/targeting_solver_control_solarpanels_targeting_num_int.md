---
id: gnc.targeting_solver_control_solarpanels_targeting_num_int
label: control_solarpanels_targeting_num_int
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl
  symbol: control_solarpanels_targeting_num_int
  lines:
  - 116
  - 116
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
- id: time_0
  type: Any
  units: n/a
  required: true
  description: Positional argument `time_0`.
- id: in_cond
  type: Any
  units: n/a
  required: true
  description: Positional argument `in_cond`.
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
  description: Return value of `control_solarpanels_targeting_num_int`. Returns `(energy_fin
    - energy_f) / 1e6` or `t_switch`.
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

# control_solarpanels_targeting_num_int

## Purpose
Finds the solar-panel switch time `t_switch` (seconds after `time_0`) at which the drag configuration must change so that a numerically integrated aerobraking pass ends at the requested specific energy `energy_f`. It is the numerical-integration variant of the targeting family, contrasted with the heat-load and closed-form variants.

## Design & Implementation
Reads `log_enabled = _bridge_verbose_enabled(param.args)` and defines the closure `func_targeting_num_int(t_switch)`, which propagates a pass with `asim_ctrl_targeting(t_switch, param, time_0, in_cond; cnf=_bridge_get_cnf(param))`, computes the final energy `|v|^2/2 - μ/|r|` from `sol[4:6,end]` and `sol[1:3,end]` using `param.mission.planet.μ`, and returns `(energy_fin - energy_f) / 1e6` so the residual is O(1) in MJ/kg. The root is bracketed on `[0, 600]` seconds and solved with `find_zero(..., Roots.Brent(), verbose=log_enabled, rtol=1e-5)`. Returns the scalar `t_switch`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `energy_f` | Any | n/a | yes | Positional argument `energy_f`. |
| in | `param` | Any | n/a | yes | Positional argument `param`. |
| in | `time_0` | Any | n/a | yes | Positional argument `time_0`. |
| in | `in_cond` | Any | n/a | yes | Positional argument `in_cond`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `control_solarpanels_targeting_num_int`. Returns `(energy_fin - energy_f) / 1e6` or `t_switch`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl`

**Downstream**

- `callees` → [[gnc.bridge_helpers__bridge_verbose_enabled|_bridge_verbose_enabled]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:117-117`
<!-- vulcan:connections:end -->

## Limitations
The bracket `[0, 600]` s is hard-coded; if the residual has the same sign at both ends `find_zero` throws an `ArgumentError` rather than returning a saturated value. Each Brent iteration runs a full trajectory integration, so cost is roughly a dozen ODE solves. The `1e6` scaling constant is an implicit unit assumption (J/kg). `in_cond` and `param` are untyped, and `param` must expose `.args` and `.mission` fields.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl` line 116.
