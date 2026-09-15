---
id: analysis.example_support__example_smoke_mission_time
label: _example_smoke_mission_time
kind: function
source:
  file: src/analysis/verification/telemetry_verification/example_support.jl
  symbol: _example_smoke_mission_time
  lines:
  - 10
  - 10
inputs:
- id: default_time
  type: Float64
  units: n/a
  required: true
  description: Positional argument `default_time`.
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
  description: Return value of `_example_smoke_mission_time`. Returns `min(default_time,
    120.0)` or `min(default_time, parsed)`.
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

# _example_smoke_mission_time

## Purpose
Caps an example's mission duration in smoke mode so the run finishes in seconds, while letting a harness raise or lower the cap.

## Design & Implementation
Reads `SPACEAGORA_EXAMPLE_SMOKE_MISSION_TIME` with a default of `"120.0"` and parses it with `tryparse`. If parsing fails or the value is not strictly positive it falls back to `min(default_time, 120.0)`; otherwise it returns `min(default_time, parsed)`. Taking the minimum against the example's own `default_time` means the cap only ever shortens a mission and can never extend one beyond what the example asked for.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `default_time` | Float64 | n/a | yes | Positional argument `default_time`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_example_smoke_mission_time`. Returns `min(default_time, 120.0)` or `min(default_time, parsed)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.example_support__example_smoke_args|_example_smoke_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:31-31`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A malformed override is silently replaced by the 120-second default rather than reported, so a typo in the variable produces a plausible-looking but unintended run length; the units are seconds by convention and not validated.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/example_support.jl` line 10.
