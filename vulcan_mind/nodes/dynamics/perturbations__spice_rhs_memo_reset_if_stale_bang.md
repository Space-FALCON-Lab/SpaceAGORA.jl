---
id: dynamics.perturbations__spice_rhs_memo_reset_if_stale_bang
label: _spice_rhs_memo_reset_if_stale!
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _spice_rhs_memo_reset_if_stale!
  lines:
  - 1016
  - 1016
inputs:
- id: memo
  type: SpiceRhsMemo
  units: n/a
  required: true
  description: Positional argument `memo`.
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
- id: primary_body_name
  type: String
  units: n/a
  required: true
  description: Positional argument `primary_body_name`.
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
  type: Nothing
  units: n/a
  description: Return value of `_spice_rhs_memo_reset_if_stale!`; mutates `memo` in
    place. Returns `nothing`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _spice_rhs_memo_reset_if_stale!

## Purpose
Clears the single-epoch SPICE memo when the requested epoch or primary body differs from the memoised one, so a stale position is never served across integrator stages.

## Design & Implementation
Compares `memo.et` with `et` and `memo.primary_body_name` with the request, and on any mismatch overwrites both and empties the position dictionary. It must be called with `memo.lock` held, which both memoised lookups do. Returns `nothing`. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `memo` | SpiceRhsMemo | n/a | yes | Positional argument `memo`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `primary_body_name` | String | n/a | yes | Positional argument `primary_body_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_spice_rhs_memo_reset_if_stale!`; mutates `memo` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__nbody_body_position_from_spice_memoized_j2000_m|_nbody_body_position_from_spice_memoized_j2000_m]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1060-1060`
- [[dynamics.perturbations__srp_sun_position_from_spice_memoized_j2000_m|_srp_sun_position_from_spice_memoized_j2000_m]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1432-1432`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Exact floating-point comparison on `et`, so evaluations at stage times differing by an ulp discard the memo; and clearing on primary change means a multi-planet configuration thrashes it.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1016.
