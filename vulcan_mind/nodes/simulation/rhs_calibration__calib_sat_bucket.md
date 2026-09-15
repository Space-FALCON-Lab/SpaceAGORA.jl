---
id: simulation.rhs_calibration__calib_sat_bucket
label: _calib_sat_bucket
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _calib_sat_bucket
  lines:
  - 62
  - 62
inputs:
- id: n
  type: Int
  units: n/a
  required: true
  description: Positional argument `n`.
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
  type: String
  units: n/a
  description: Return value of `_calib_sat_bucket`.
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

# _calib_sat_bucket

## Purpose
Maps an active-satellite count to a coarse bucket label so that calibration signatures group constellations of similar size rather than requiring an exact count match.

## Design & Implementation
Pure `@inline` function taking `n::Int` and returning one of nine literal strings via a cascade of `<=` comparisons: `"1"` for n <= 1, then the power-of-two ranges `"2_4"`, `"5_8"`, `"9_16"`, `"17_32"`, `"33_64"`, `"65_128"`, `"129_256"`, and `"257p"` for anything larger. The string becomes the `sats=` field inside `_rhs_calib_signature`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n` | Int | n/a | yes | Positional argument `n`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_calib_sat_bucket`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- [[simulation.rhs_calibration__rhs_calib_signature|_rhs_calib_signature]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:86-86`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Negative or zero counts fall into the `"1"` bucket without error. Bucket boundaries are hard-coded; a constellation of 5 satellites and one of 8 share a plan even though the optimal flat allotment can differ between them. Because the bucket is embedded in the signature, editing these ranges silently invalidates all previously persisted calibrations.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 62.
