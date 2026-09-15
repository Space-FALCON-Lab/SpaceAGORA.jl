---
id: gnc.targeting_control_apoapsis_residual
label: apoapsis_residual
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: apoapsis_residual
  lines:
  - 1032
  - 1032
inputs:
- id: t_switch
  type: Any
  units: n/a
  required: true
  description: Positional argument `t_switch`.
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
  description: Return value of `apoapsis_residual`. Returns `outcome.apoapsis_radius_m
    - target_apoapsis`.
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

# apoapsis_residual

## Purpose
The residual whose root is the switch time that achieves the target apoapsis radius, preferred over the energy residual when a target apoapsis is configured.

## Design & Implementation
Evaluates the candidate at `t_switch` and returns `apoapsis_radius_m - target_apoapsis`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t_switch` | Any | n/a | yes | Positional argument `t_switch`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `apoapsis_residual`. Returns `outcome.apoapsis_radius_m - target_apoapsis`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.targeting_control_evaluate_candidate|evaluate_candidate]] · `callers` · call · `src/gnc/control/targeting_control.jl:1033-1033`
<!-- vulcan:connections:end -->

## Limitations
Apoapsis is `Inf` for an unbound predicted orbit, so this residual is non-finite for candidates that fail to capture, which the caller's sign check catches by falling back to the energy solve.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 1032.
