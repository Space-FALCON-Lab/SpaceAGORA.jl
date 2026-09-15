---
id: gnc.guidance_hooks__control_solarpanels_heatrate
label: _control_solarpanels_heatrate
kind: function
source:
  file: src/gnc/guidance/guidance_hooks.jl
  symbol: _control_solarpanels_heatrate
  lines:
  - 35
  - 35
inputs:
- id: args
  type: Vararg{Any}
  units: n/a
  required: false
  description: Positional argument `args` (variadic).
- id: kwargs
  type: Vararg{Any}
  units: n/a
  required: false
  description: Keyword argument `kwargs` (variadic).
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
  description: Return value of `_control_solarpanels_heatrate`. Returns `_control_module().control_solarpanels_heatrate(args...;
    kwargs...)`.
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

# _control_solarpanels_heatrate

## Purpose
Forwards to `ControlHooks.control_solarpanels_heatrate`, the routine that solves for the solar-panel angle of attack which holds convective heating at the configured maximum heat rate during a drag passage.

## Design & Implementation
Written as `_control_module().control_solarpanels_heatrate(args...; kwargs...)` with `@inline`. Guidance calls it from inside derivative evaluations, for example in the constraint-tracking right-hand side where an exceeded `settings.max_heat_rate` triggers a panel-angle solve using the local atmospheric state vector of temperature, density and molecular speed ratio.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vararg{Any} | n/a | no | Positional argument `args` (variadic). |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_control_solarpanels_heatrate`. Returns `_control_module().control_solarpanels_heatrate(args...; kwargs...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:211-211`
- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:264-264`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:211-211`
- [[gncy.guidance_hooks_guidancehooks|GuidanceHooks]] · `callees` → `callers` · call · `src/gnc/guidance/guidance_hooks.jl:35-35`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:264-264`

**Downstream**

- `callees` → [[gnc.guidance_hooks__control_module|_control_module]] · `callers` · call · `src/gnc/guidance/guidance_hooks.jl:36-36`
<!-- vulcan:connections:end -->

## Limitations
Being on the derivative path, the dynamic module lookup and variadic splat happen at every solver stage where the heat-rate limit is active, which is measurable inside a stiff pass. No result caching or memoisation exists at this level. Any error raised by the underlying panel solve propagates out through the integrator, aborting the pass rather than saturating the command.

## Provenance
Mapped from `src/gnc/guidance/guidance_hooks.jl` line 35.
