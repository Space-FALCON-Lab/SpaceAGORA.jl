---
id: analysis.scenario_builders__period_seconds
label: _period_seconds
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _period_seconds
  lines:
  - 1
  - 1
inputs:
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: ra
  type: Float64
  units: n/a
  required: true
  description: Positional argument `ra`.
- id: rp
  type: Float64
  units: n/a
  required: true
  description: Positional argument `rp`.
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
  description: Return value of `_period_seconds`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# _period_seconds

## Purpose
Computes the Keplerian orbital period from apoapsis and periapsis radii, used to convert a target orbit count into a mission duration.

## Theory & Math
$$
T = 2\pi \sqrt{\frac{a^3}{\mu}},\qquad a = \tfrac{1}{2}(r_a + r_p)
$$

## Design & Implementation
Forms the semi-major axis as half the sum of `ra` and `rp`, then returns `2π sqrt(a³/μ)` using the planet's `μ`. `@inline` with a `::Float64` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `ra` | Float64 | n/a | yes | Positional argument `ra`. |
| in | `rp` | Float64 | n/a | yes | Positional argument `rp`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_period_seconds`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_orbit_args|_make_orbit_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:569-569`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It is the two-body period; drag, oblateness and third bodies shorten or lengthen real orbits, so the mission time derived from it drifts from the actual orbit count over a long campaign.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 1.
