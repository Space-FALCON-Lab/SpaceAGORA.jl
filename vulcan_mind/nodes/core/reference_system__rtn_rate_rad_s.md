---
id: core.reference_system__rtn_rate_rad_s
label: _rtn_rate_rad_s
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: _rtn_rate_rad_s
  lines:
  - 639
  - 639
inputs:
- id: r_target_ii
  type: Any
  units: n/a
  required: true
  description: Positional argument `r_target_ii`.
- id: v_target_ii
  type: Any
  units: n/a
  required: true
  description: Positional argument `v_target_ii`.
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
  description: Return value of `_rtn_rate_rad_s`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# _rtn_rate_rad_s

## Purpose
Computes the instantaneous angular rate of the RTN frame about its cross-track axis, needed to convert between inertial and rotating-frame relative velocities.

## Theory & Math
$$
n = \frac{|\vec{r} \times \vec{v}|}{|\vec{r}|^2}
$$

## Design & Implementation
Returns `|r × v| / |r|²`, which equals the true-anomaly rate for any conic. Raises `ArgumentError` on a zero-length position. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_target_ii` | Any | n/a | yes | Positional argument `r_target_ii`. |
| in | `v_target_ii` | Any | n/a | yes | Positional argument `v_target_ii`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_rtn_rate_rad_s`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.reference_system_inertial_to_rtn_relative_state|inertial_to_rtn_relative_state]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:660-660`
- [[core.reference_system_rtn_to_inertial_relative_state|rtn_to_inertial_relative_state]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:684-684`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/interfaces/reference_system.jl`

**Downstream**

- `callees` → [[core.reference_system_inertial_to_rtn_relative_state|inertial_to_rtn_relative_state]] · `callers` · feedback · `src/core/interfaces/reference_system.jl:648-648`
<!-- vulcan:connections:end -->

## Limitations
This is the exact instantaneous rate, which for an eccentric target differs from the mean motion the HCW matrices assume; callers mixing this with `rpo_hcw_continuous_mats` on eccentric references introduce a modelling inconsistency.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 639.
