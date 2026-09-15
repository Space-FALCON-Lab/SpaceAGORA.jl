---
id: gnc.targeting_control_max_energy_alpha
label: max_energy_alpha
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: max_energy_alpha
  lines:
  - 592
  - 592
inputs:
- id: r
  type: Any
  units: n/a
  required: true
  description: Positional argument `r`.
- id: v
  type: Any
  units: n/a
  required: true
  description: Positional argument `v`.
- id: tau
  type: Any
  units: n/a
  required: true
  description: Positional argument `tau`.
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
  description: Return value of `max_energy_alpha`. Returns `alpha`.
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

# max_energy_alpha

## Purpose
The closure that selects the predicted angle at each step of the max-energy-depletion integration.

## Design & Implementation
Computes absolute time, returns `min_alpha_rad` if the heat-load sub-mode is on and the time lies within the finite `heat_load_switches` window, and otherwise samples the environment and constrains `max_alpha_rad` through `_edg_targeting_constrained_alpha` with the captured `alpha_past`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r` | Any | n/a | yes | Positional argument `r`. |
| in | `v` | Any | n/a | yes | Positional argument `v`. |
| in | `tau` | Any | n/a | yes | Positional argument `tau`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `max_energy_alpha`. Returns `alpha`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.heat_load_control_acceleration|acceleration]] · `callers` · call · `src/gnc/control/targeting_control.jl:624-624`
- `callees` → [[gnc.targeting_control__edg_targeting_constrained_alpha|_edg_targeting_constrained_alpha]] · `callers` · call · `src/gnc/control/targeting_control.jl:602-602`
- `callees` → [[gnc.targeting_control__edg_targeting_prediction_environment|_edg_targeting_prediction_environment]] · `callers` · call · `src/gnc/control/targeting_control.jl:601-601`
- `callees` → [[gnc.targeting_control_acceleration|acceleration]] · `callers` · call · `src/gnc/control/targeting_control.jl:624-624`
<!-- vulcan:connections:end -->

## Limitations
It reads `alpha_past` from the enclosing scope but the enclosing loop updates that variable after the call, so the warm start is always one step stale.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 592.
