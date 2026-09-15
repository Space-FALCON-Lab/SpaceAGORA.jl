---
id: gnc.targeting_control_energy_residual
label: energy_residual
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: energy_residual
  lines:
  - 1027
  - 1027
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
  description: Return value of `energy_residual`. Returns `outcome.energy_jkg - target_energy`.
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

# energy_residual

## Purpose
The residual whose root is the switch time that achieves the target energy, for the Brent solve.

## Design & Implementation
Evaluates the candidate at `t_switch` and returns `energy_jkg - target_energy`, capturing `target_energy` from the enclosing solve.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t_switch` | Any | n/a | yes | Positional argument `t_switch`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `energy_residual`. Returns `outcome.energy_jkg - target_energy`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.targeting_control_evaluate_candidate|evaluate_candidate]] · `callers` · call · `src/gnc/control/targeting_control.jl:1028-1028`
<!-- vulcan:connections:end -->

## Limitations
Each evaluation is a full prediction; Brent at `rtol = 1e-7` typically needs a dozen or more.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 1027.
