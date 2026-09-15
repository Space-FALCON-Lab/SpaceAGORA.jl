---
id: gnc.heat_load_control__edg_first_low_alpha_interval_indices
label: _edg_first_low_alpha_interval_indices
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_first_low_alpha_interval_indices
  lines:
  - 427
  - 427
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
  type: Tuple
  units: n/a
  description: Return value of `_edg_first_low_alpha_interval_indices`. Returns `(first_low,
    last_low)`.
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

# _edg_first_low_alpha_interval_indices

## Purpose
Finds the first contiguous run of minimum-alpha nodes in a bang-bang profile, which defines the single low-drag window the controller will actually fly.

## Design & Implementation
Takes `config` and `alpha_profile::Vector{Float64}`. Builds the boolean mask `low = alpha <= config.min_alpha_rad + 1e-8`, locates `first_low = findfirst(identity, low)` (returning `nothing` when there is no low node), then extends `last_low` forward while the next entry is also low. Returns the tuple `(first_low, last_low)` of 1-based indices.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `alpha_profile` | Vector{Float64} | n/a | yes | Positional argument `alpha_profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_edg_first_low_alpha_interval_indices`. Returns `(first_low, last_low)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_first_two_switch_alpha_profile|_edg_first_two_switch_alpha_profile]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:441-441`
- [[gnc.heat_load_control__edg_low_alpha_switch_window|_edg_low_alpha_switch_window]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:450-450`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the first run is returned; if the optimal profile has two or more separated low-alpha intervals the later ones are discarded by design. The `1e-8` rad tolerance is hard-coded. Allocates a `BitVector`-like mask through `map` each call.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 427.
