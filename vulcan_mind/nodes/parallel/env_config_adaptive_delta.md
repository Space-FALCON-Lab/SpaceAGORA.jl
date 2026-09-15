---
id: parallel.env_config_adaptive_delta
label: adaptive_delta
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: adaptive_delta
  lines:
  - 215
  - 215
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
  description: Return value of `adaptive_delta`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# adaptive_delta

## Purpose
Supplies the discount factor δ used by the adaptive thread-policy controller to weight older reward observations less than recent ones, validated to lie in (0, 1].

## Design & Implementation
Parses `SPACEAGORA_PARALLEL_POLICY_DELTA` with `parse_float_env` (default `0.85`). If `!(0.0 < δ <= 1.0)` it throws `ArgumentError("SPACEAGORA_PARALLEL_POLICY_DELTA must satisfy 0 < δ <= 1, got '<δ>'")`. Returns `Float64`. A value of exactly 1.0 disables discounting; values near 0 make the controller react almost solely to the latest sample.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `adaptive_delta`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.env_config_snapshot_policy_decision_env|snapshot_policy_decision_env]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:215-215`
- [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callees` → `callers` · call · `src/parallel/policy/observation_tracking.jl:70-70`

**Downstream**

- `callees` → [[parallel.env_config_parse_float_env|parse_float_env]] · `callers` · call · `src/parallel/policy/env_config.jl:216-216`
<!-- vulcan:connections:end -->

## Limitations
`NaN` fails the compound comparison and is rejected, but there is no dedicated message for it. The variable is read on each call rather than captured in `PolicyDecisionEnvConfig`, so it is the caller's responsibility to read it once per controller instance.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 215.
