---
id: gnc.targeting_control_acceleration
label: acceleration
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: acceleration
  lines:
  - 479
  - 479
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
- id: alpha
  type: Any
  units: n/a
  required: true
  description: Positional argument `alpha`.
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
  description: Return value of `acceleration`. Returns `gravity + aero`.
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

# acceleration

## Purpose
The closure that gives total acceleration — point-mass gravity plus aerodynamics — for the prediction integrators' RK4 stages.

## Design & Implementation
Defined identically inside both `_edg_integrated_targeting_trajectory` and `_edg_integrated_max_energy_depletion_trajectory`. It returns `-μ r / |r|³` plus `_edg_targeting_aero_acceleration` at `t + tau` with the given `alpha`, capturing `planet`, `config`, `p`, `spacecraft`, `mass` and `t`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r` | Any | n/a | yes | Positional argument `r`. |
| in | `v` | Any | n/a | yes | Positional argument `v`. |
| in | `tau` | Any | n/a | yes | Positional argument `tau`. |
| in | `alpha` | Any | n/a | yes | Positional argument `alpha`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `acceleration`. Returns `gravity + aero`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_integrated_max_energy_depletion_trajectory|_edg_integrated_max_energy_depletion_trajectory]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:585-585`
- [[gnc.targeting_control_max_energy_alpha|max_energy_alpha]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:624-624`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.targeting_control__edg_targeting_aero_acceleration|_edg_targeting_aero_acceleration]] · `callers` · call · `src/gnc/control/targeting_control.jl:482-482`
- `callees` → [[gnc.targeting_control__edg_targeting_constrained_alpha|_edg_targeting_constrained_alpha]] · `callers` · call · `src/gnc/control/targeting_control.jl:492-492`
- `callees` → [[gnc.targeting_control__edg_targeting_prediction_environment|_edg_targeting_prediction_environment]] · `callers` · call · `src/gnc/control/targeting_control.jl:491-491`
<!-- vulcan:connections:end -->

## Limitations
Gravity is the point-mass term only, so J2 and higher harmonics are absent from the prediction even when the truth propagator includes them; over a single pass this introduces a small but systematic apoapsis error.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 479.
