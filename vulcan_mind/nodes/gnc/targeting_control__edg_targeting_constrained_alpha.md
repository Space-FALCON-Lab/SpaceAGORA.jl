---
id: gnc.targeting_control__edg_targeting_constrained_alpha
label: _edg_targeting_constrained_alpha
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_targeting_constrained_alpha
  lines:
  - 365
  - 365
inputs:
- id: config
  type: AerobrakingEnergyDepletionConfig
  units: n/a
  required: true
  description: Positional argument `config`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: env
  type: Any
  units: n/a
  required: true
  description: Positional argument `env`.
- id: base_alpha
  type: Float64
  units: n/a
  required: true
  description: Positional argument `base_alpha`.
- id: alpha_past
  type: Float64
  units: n/a
  required: true
  description: Positional argument `alpha_past`.
- id: heat_rate_control
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `heat_rate_control`.
- id: structural_control
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `structural_control`.
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
  type: Float64
  units: n/a
  description: Return value of `_edg_targeting_constrained_alpha`.
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

# _edg_targeting_constrained_alpha

## Purpose
Applies the heat-rate and structural constraints to a predicted base angle inside the prediction integrator, mirroring what the live controller would do.

## Design & Implementation
Clamps the base angle and returns the minimum immediately if the base is already there. With `heat_rate_control` it runs the Maxwellian root solve with the thermal accommodation factor and warm start `alpha_past`; with `structural_control` it runs the structural-load solve. Returns the clamped minimum of the two.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | AerobrakingEnergyDepletionConfig | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `env` | Any | n/a | yes | Positional argument `env`. |
| in | `base_alpha` | Float64 | n/a | yes | Positional argument `base_alpha`. |
| in | `alpha_past` | Float64 | n/a | yes | Positional argument `alpha_past`. |
| in | `heat_rate_control` | Bool | n/a | yes | Keyword argument `heat_rate_control`. |
| in | `structural_control` | Bool | n/a | yes | Keyword argument `structural_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_edg_targeting_constrained_alpha`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control_acceleration|acceleration]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:492-492`
- [[gnc.targeting_control_max_energy_alpha|max_energy_alpha]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:602-602`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/targeting_control.jl:382-382`
- `callees` → [[gnc.heat_rate_control__energy_depletion_heatrate_root_alpha|_energy_depletion_heatrate_root_alpha]] · `callers` · call · `src/gnc/control/targeting_control.jl:384-384`
- `callees` → [[gnc.struct_load_control__energy_depletion_struct_load_root_alpha|_energy_depletion_struct_load_root_alpha]] · `callers` · call · `src/gnc/control/targeting_control.jl:398-398`
<!-- vulcan:connections:end -->

## Limitations
The early return at minimum angle skips both solves, which is correct but means the prediction never reports the constrained values for those samples.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 365.
