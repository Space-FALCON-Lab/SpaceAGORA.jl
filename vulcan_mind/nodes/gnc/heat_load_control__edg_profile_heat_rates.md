---
id: gnc.heat_load_control__edg_profile_heat_rates
label: _edg_profile_heat_rates
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_profile_heat_rates
  lines:
  - 388
  - 388
inputs:
- id: config
  type: Any
  units: n/a
  required: true
  description: Positional argument `config`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: track
  type: Any
  units: n/a
  required: true
  description: Positional argument `track`.
- id: alpha_profile
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `alpha_profile`.
- id: heat_rate_control
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `heat_rate_control`.
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
  description: Return value of `_edg_profile_heat_rates`. Returns `qdot`.
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

# _edg_profile_heat_rates

## Purpose
Evaluates the stagnation heat-rate history along a track for a given alpha profile, optionally capping it at the configured heat-rate limit where active control would hold alpha down.

## Design & Implementation
Signature `(config, p, track, alpha_profile; heat_rate_control::Bool)`. It reads `taf` from the thermal model (default 1.0) and `heat_rate_limit = config.heat_rate_limit_w_cm2`. For each node it calls `_energy_depletion_heat_rate_calc(taf, max(0, rho), T, T, planet.R, planet.γ, max(0, speed_ratio), alpha)`, passing the atmospheric temperature as both gas and wall temperature. When `heat_rate_control` is true, the node's alpha exceeds `min_alpha + 1e-8`, and the limit is finite and positive, the rate is clipped to the limit. Non-finite or negative rates become 0. Returns `qdot::Vector{Float64}` in W/cm^2.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `track` | Any | n/a | yes | Positional argument `track`. |
| in | `alpha_profile` | Vector{Float64} | n/a | yes | Positional argument `alpha_profile`. |
| in | `heat_rate_control` | Bool | n/a | yes | Keyword argument `heat_rate_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_profile_heat_rates`. Returns `qdot`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_balanced_tpbvp_heat_load_window|_edg_balanced_tpbvp_heat_load_window]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:486-486`
- [[gnc.heat_load_control__edg_profile_heat_load|_edg_profile_heat_load]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:423-423`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/heat_load_control.jl:390-390`
- `callees` → [[gnc.heat_rate_control__energy_depletion_heat_rate_calc|_energy_depletion_heat_rate_calc]] · `callers` · call · `src/gnc/control/heat_load_control.jl:395-395`
<!-- vulcan:connections:end -->

## Limitations
Using the atmospheric temperature as the wall temperature ignores surface heating. The clip to `heat_rate_limit` assumes the controller can always achieve the limit exactly, which is optimistic when the root solver saturates at `min_alpha`. No accounting for view factors or per-panel differences; a single scalar rate per node.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 388.
