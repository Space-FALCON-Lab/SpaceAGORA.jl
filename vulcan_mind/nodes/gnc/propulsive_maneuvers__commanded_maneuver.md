---
id: gnc.propulsive_maneuvers__commanded_maneuver
label: _commanded_maneuver
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _commanded_maneuver
  lines:
  - 100
  - 100
inputs:
- id: controlModel
  type: BaseThrusterModel
  units: n/a
  required: true
  description: Positional argument `controlModel`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
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
  description: Return value of `_commanded_maneuver`. Returns `(` or `nothing`.
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

# _commanded_maneuver

## Purpose
Resolves the delta-v magnitude and thrust direction to execute for spacecraft `i`, preferring a guidance command and falling back to the thruster model's own arrays, and normalises negative delta-v into a retrograde direction.

## Design & Implementation
First tries `_guidance_maneuver_command(p, i)`. When a command exists it converts `delta_v_mps` and `direction_rad` to `Float64`; if the delta-v is finite and negative it is replaced by its absolute value and the direction is forced to `π`. It then writes both back into `controlModel.Δv[i]` and `controlModel.direction[i]`, mutating the model, and returns a named tuple carrying `command.source_orbit`. With no guidance command it bounds-checks `i` against both `controlModel.Δv` and `controlModel.direction`, returning `nothing` if either is too short; otherwise it applies the same sign normalisation — positive delta-v forces direction `0.0`, negative forces magnitude and `π` — writes the values back, and returns the tuple with `source_orbit = -1`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `controlModel` | BaseThrusterModel | n/a | yes | Positional argument `controlModel`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_commanded_maneuver`. Returns `(` or `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:496-496`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:103-103`
- `callees` → [[gnc.propulsive_maneuvers__guidance_maneuver_command|_guidance_maneuver_command]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:101-101`
<!-- vulcan:connections:end -->

## Limitations
The function mutates the caller's thruster model as a side effect of what reads like a query, so a guidance command permanently overwrites the model's configured direction. The model-array path overrides any non-zero commanded direction with exactly `0.0` or `π`, discarding off-axis pointing. A non-finite delta-v skips normalisation entirely and is returned unchanged for downstream validation to reject.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 100.
