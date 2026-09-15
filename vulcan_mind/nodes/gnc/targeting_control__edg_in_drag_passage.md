---
id: gnc.targeting_control__edg_in_drag_passage
label: _edg_in_drag_passage
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_in_drag_passage
  lines:
  - 91
  - 91
inputs:
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: env
  type: Any
  units: n/a
  required: true
  description: Positional argument `env`.
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
  type: Bool
  units: n/a
  description: Return value of `_edg_in_drag_passage`.
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

# _edg_in_drag_passage

## Purpose
Tests whether the vehicle is currently below the entry interface altitude, which gates when switch times are solved.

## Design & Implementation
Converts `environment_model.EI` from kilometres to metres and returns true when both it and the sampled altitude are finite and altitude is at or below it. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `env` | Any | n/a | yes | Positional argument `env`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_edg_in_drag_passage`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_recompute_switches_bang|_edg_recompute_switches!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:112-112`
- [[gncx.struct_load_control__edg_structural_alpha|_edg_structural_alpha]] · `callees` → `callers` · call · `src/gnc/control/struct_load_control.jl:105-105`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/targeting_control.jl:92-92`
<!-- vulcan:connections:end -->

## Limitations
The entry interface is a single altitude with no hysteresis, so a vehicle skimming the interface toggles the passage flag on every crossing.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 91.
