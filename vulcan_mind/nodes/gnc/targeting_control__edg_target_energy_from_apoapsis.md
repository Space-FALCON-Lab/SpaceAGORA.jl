---
id: gnc.targeting_control__edg_target_energy_from_apoapsis
label: _edg_target_energy_from_apoapsis
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_target_energy_from_apoapsis
  lines:
  - 307
  - 307
inputs:
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: target_apoapsis_radius_m
  type: Float64
  units: n/a
  required: true
  description: Positional argument `target_apoapsis_radius_m`.
- id: periapsis_radius_m
  type: Float64
  units: n/a
  required: true
  description: Positional argument `periapsis_radius_m`.
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
  description: Return value of `_edg_target_energy_from_apoapsis`. Returns `NaN` or
    `-planet.μ / (target_apoapsis_radius_m + periapsis_radius_m)`.
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

# _edg_target_energy_from_apoapsis

## Purpose
Converts a target apoapsis radius and the current periapsis radius into the specific orbital energy the targeting solve should aim for.

## Theory & Math
$$
\epsilon_{\text{target}} = -\frac{\mu}{r_{a,\text{target}} + r_p}
$$

## Design & Implementation
Returns `NaN` unless both radii are finite and positive, otherwise `-μ / (r_a + r_p)`, which is `-μ / 2a` with `2a` the sum of the apsides.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `target_apoapsis_radius_m` | Float64 | n/a | yes | Positional argument `target_apoapsis_radius_m`. |
| in | `periapsis_radius_m` | Float64 | n/a | yes | Positional argument `periapsis_radius_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_target_energy_from_apoapsis`. Returns `NaN` or `-planet.μ / (target_apoapsis_radius_m + periapsis_radius_m)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It assumes periapsis is unchanged by the pass, which is approximately true for a shallow drag pass but not for a deep one; the targeting solve corrects for this with an apoapsis-based refinement afterwards.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 307.
