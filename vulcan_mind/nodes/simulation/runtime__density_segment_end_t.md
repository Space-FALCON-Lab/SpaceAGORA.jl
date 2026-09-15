---
id: simulation.runtime__density_segment_end_t
label: _density_segment_end_t
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/runtime.jl
  symbol: _density_segment_end_t
  lines:
  - 44
  - 44
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: cache_cfg
  type: Any
  units: n/a
  required: true
  description: Positional argument `cache_cfg`.
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
  description: Return value of `_density_segment_end_t`.
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

# _density_segment_end_t

## Purpose
Chooses how far ahead a GRAM track cache should predict, bounding the prediction by the end of the current solve segment or the mission rather than always using the full cache horizon.

## Design & Implementation
Returns `p.shared_buffers.solve_segment_end_time[]` if it is finite and later than `t`; else the mission time from `p.args.mission_configuration` under the same test; else `t + cache_cfg.orbit_horizon_s`. The cascade means a cache built near a segment boundary does not waste GRAM calls predicting past a point where the integrator will restart anyway.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `cache_cfg` | Any | n/a | yes | Positional argument `cache_cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_density_segment_end_t`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/runtime.jl`
- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:173-173`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:49-49`
<!-- vulcan:connections:end -->

## Limitations
A finite but past segment end is treated as absent and falls through, so a stale `solve_segment_end_time` produces a full-horizon prediction rather than an error.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/runtime.jl` line 44.
