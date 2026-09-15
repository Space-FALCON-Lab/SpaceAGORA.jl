---
id: gnc.targeting_control__edg_control_pos_vel_mass
label: _edg_control_pos_vel_mass
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_control_pos_vel_mass
  lines:
  - 57
  - 57
inputs:
- id: sc
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc`.
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
  description: Return value of `_edg_control_pos_vel_mass`. Returns `pos, vel, mass`.
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

# _edg_control_pos_vel_mass

## Purpose
Extracts position, velocity and mass from a satellite state under either the labelled or the positional layout.

## Design & Implementation
Reads `sc.pos`, `sc.vel` and `sc.mass` when those properties exist, otherwise slots one to three, four to six and seven, returning `NaN` for mass if the vector is shorter than seven. Returns a tuple of two `SVector{3}` and a `Float64`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sc` | Any | n/a | yes | Positional argument `sc`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_control_pos_vel_mass`. Returns `pos, vel, mass`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_environment_state|_edg_environment_state]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:66-66`
- [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:262-262`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/targeting_control.jl:60-60`
<!-- vulcan:connections:end -->

## Limitations
A `NaN` mass propagates into the drag-passage prediction, which `_edg_predict_mass` must guard; this function itself gives no diagnostic.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 57.
