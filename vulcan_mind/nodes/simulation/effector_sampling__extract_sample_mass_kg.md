---
id: simulation.effector_sampling__extract_sample_mass_kg
label: _extract_sample_mass_kg
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: _extract_sample_mass_kg
  lines:
  - 14
  - 14
inputs:
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
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
  description: Return value of `_extract_sample_mass_kg`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _extract_sample_mass_kg

## Purpose
Recovers current mass from a state sample under any of three layouts, returning `NaN` when no mass is carried.

## Design & Implementation
Prefers a `mass_kg` property, then a `mass` property read from slot seven, then slot seven of any vector with at least seven elements, and otherwise `NaN`. The `NaN` sentinel lets downstream cache code use a finiteness check rather than a separate flag. `@inline` with a `::Float64` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_extract_sample_mass_kg`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- [[simulation.dynamics_rhs__prefill_rhs_flat_state_samples_bang|_prefill_rhs_flat_state_samples!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:305-305`
- [[simulation.effector_sampling__sample_atmosphere_from_planet_frame|_sample_atmosphere_from_planet_frame]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:92-92`
- [[simulation.effector_sampling_build_state_sample|build_state_sample]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:31-31`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:16-16`
<!-- vulcan:connections:end -->

## Limitations
Duplicates the logic of the density callback's own mass extractor with one extra branch, so the two can drift apart; both silently misread a seven-element state whose seventh slot is not mass.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 14.
