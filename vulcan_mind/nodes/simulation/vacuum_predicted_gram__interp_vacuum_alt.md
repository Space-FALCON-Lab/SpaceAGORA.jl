---
id: simulation.vacuum_predicted_gram__interp_vacuum_alt
label: _interp_vacuum_alt
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl
  symbol: _interp_vacuum_alt
  lines:
  - 158
  - 158
inputs:
- id: cache
  type: VacuumPredictedGRAMCache
  units: n/a
  required: true
  description: Positional argument `cache`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  description: Return value of `_interp_vacuum_alt`.
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

# _interp_vacuum_alt

## Purpose
Linearly interpolates the vacuum-predicted altitude stored at the cache knots to an arbitrary query time, giving callers a cheap estimate of where the drag-free reference trajectory would be at time `t` without re-propagating.

## Design & Implementation
Signature `_interp_vacuum_alt(cache::VacuumPredictedGRAMCache, t::Float64)::Float64`, `@inline`. With `n = length(cache.vac_alts)`, it computes the segment `idx = clamp(floor(Int, (t - cache.t0)/cache.h) + 1, 1, n - 1)` and the normalised offset `x = (t - (cache.t0 + (idx-1) cache.h))/cache.h`, then returns `(1 - x) vac_alts[idx] + x vac_alts[idx+1]` under `@inbounds`. Altitudes were computed at build time by `rtolatlong` on the planet-fixed position, so they are geodetic altitudes in metres.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | VacuumPredictedGRAMCache | n/a | yes | Positional argument `cache`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_interp_vacuum_alt`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Times outside the cached horizon are extrapolated linearly from the end segment because of the clamp, with `x` outside `[0, 1]`. If the cache is invalid or `vac_alts` is empty, `n - 1` is zero and the `@inbounds` read is undefined behaviour; the function does not check `cache.valid`. Linear interpolation between knots ~30 s apart under-resolves altitude near periapsis where the trajectory curves fastest.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl` line 158.
