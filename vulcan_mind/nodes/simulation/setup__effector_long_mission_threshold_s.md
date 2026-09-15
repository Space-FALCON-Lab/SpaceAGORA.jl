---
id: simulation.setup__effector_long_mission_threshold_s
label: _effector_long_mission_threshold_s
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_long_mission_threshold_s
  lines:
  - 469
  - 469
inputs:
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
  description: Return value of `_effector_long_mission_threshold_s`.
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

# _effector_long_mission_threshold_s

## Purpose
Mission duration, in seconds, above which a time-based mission is considered long enough for effector threading setup costs to amortise.

## Design & Implementation
Returns `_parse_positive_float_env("SPACEAGORA_EFFECTOR_LONG_MISSION_THRESHOLD_S", 5400.0)`, a 90-minute default that corresponds to roughly one low-Earth orbit. `_mission_is_long_for_effector_threads` compares `mission_cfg.mission_time` against it when `mission_type != MissionOrbits`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_effector_long_mission_threshold_s`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__mission_is_long_for_effector_threads|_mission_is_long_for_effector_threads]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:503-503`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:470-470`
- `callees` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/simulation/engine/setup.jl:470-470`
<!-- vulcan:connections:end -->

## Limitations
Read via `_parse_positive_float_env` on each call rather than snapshotted, so it is one of the few knobs in this family that still touches the environment inside `_mission_is_long_for_effector_threads`. The threshold ignores step count, which is what actually determines how many times threading overhead is paid.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 469.
