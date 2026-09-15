---
id: simulation.setup__planet_transform_key
label: _planet_transform_key
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _planet_transform_key
  lines:
  - 244
  - 244
inputs:
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
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
  type: NTuple{9,
  units: n/a
  description: Return value of `_planet_transform_key`.
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

# _planet_transform_key

## Purpose
Reduces the physical constants that define a planet's rotation and shape to a 9-tuple of integers so that planet-frame caches are only reused between runs whose planet definitions are numerically identical.

## Design & Implementation
Takes any `planet` with fields `ω` (3-vector angular velocity, rad/s), `α` and `δ` (pole right ascension and declination, rad), `Rp_e`, `Rp_p`, `Rp_m` (equatorial, polar, mean radii, km), and `μ` (gravitational parameter, km³/s²). Returns `NTuple{9, Int64}` with `ω` components, `α`, `δ` scaled by `1e12`; radii scaled by `1e6`; and `μ` scaled by `1e-6`. Each is passed through `round(Int64, ...)`. The tuple is embedded in `PlanetFrameEphemerisReuseKey`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NTuple{9, | n/a | — | Return value of `_planet_transform_key`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__planet_frame_ephemeris_reuse_key|_planet_frame_ephemeris_reuse_key]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:293-293`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The scale factors are hard-coded: 1e12 on rad/s gives picoradian resolution, 1e6 on km gives millimetre resolution, but 1e-6 on μ gives only 1e6 km³/s² resolution, so bodies whose μ differ by less than that (about 0.25% of Earth's) collide. Overflow throws `InexactError` for planets with `μ` above roughly 9e24. Fields other than these nine (for example a J2 coefficient) are not part of the key.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 244.
