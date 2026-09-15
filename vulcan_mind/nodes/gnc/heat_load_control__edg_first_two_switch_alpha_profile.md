---
id: gnc.heat_load_control__edg_first_two_switch_alpha_profile
label: _edg_first_two_switch_alpha_profile
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_first_two_switch_alpha_profile
  lines:
  - 439
  - 439
inputs:
- id: config
  type: Any
  units: n/a
  required: true
  description: Positional argument `config`.
- id: alpha_profile
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `alpha_profile`.
- id: high_profile
  type: Vector{Float64}
  units: n/a
  required: false
  description: Positional argument `high_profile` (default `fill(config.max_alpha_rad,
    length(alpha_profile))`).
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
  description: Return value of `_edg_first_two_switch_alpha_profile`. Returns `applied_profile`.
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

# _edg_first_two_switch_alpha_profile

## Purpose
Projects an arbitrary bang-bang alpha profile onto a two-switch profile: high alpha everywhere except the first low-alpha interval, matching what the flight controller can execute.

## Design & Implementation
Signature `(config, alpha_profile, high_profile = fill(config.max_alpha_rad, length(alpha_profile)))`. Copies `high_profile`, locates the first low interval with `_edg_first_low_alpha_interval_indices`, and if one exists overwrites `applied_profile[first_low:last_low] .= config.min_alpha_rad`. Returns the new vector. Passing a constrained `high_profile` lets the high segments carry heat-rate and structural limits.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `alpha_profile` | Vector{Float64} | n/a | yes | Positional argument `alpha_profile`. |
| in | `high_profile` | Vector{Float64} | n/a | no | Positional argument `high_profile` (default `fill(config.max_alpha_rad, length(alpha_profile))`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_first_two_switch_alpha_profile`. Returns `applied_profile`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control_residual|residual]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:674-674`
- [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:674-674`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_first_low_alpha_interval_indices|_edg_first_low_alpha_interval_indices]] · `callers` · call · `src/gnc/control/heat_load_control.jl:441-441`
<!-- vulcan:connections:end -->

## Limitations
Any low-alpha nodes after the first interval are silently promoted to the high profile, so the predicted heat load can exceed the raw optimum. The function does not verify that `high_profile` and `alpha_profile` have the same length.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 439.
