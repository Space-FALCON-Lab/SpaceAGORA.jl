---
id: analysis.example_support__link_q
label: _link_q
kind: function
source:
  file: src/analysis/verification/telemetry_verification/example_support.jl
  symbol: _link_q
  lines:
  - 98
  - 98
inputs:
- id: q
  type: Any
  units: n/a
  required: true
  description: Positional argument `q`.
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
  description: Return value of `_link_q`. Returns `q === nothing ?`.
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

# _link_q

## Purpose
Turns an optional attitude quaternion tuple into the mutable four-vector a `Link` stores, defaulting to identity when none was given.

## Design & Implementation
A local closure inside `make_three_body_spacecraft`. When `q` is `nothing` it returns `MVector{4,Float64}(0, 0, 0, 1)`, the scalar-last identity that the `Link` constructor also defaults to, so callers omitting attitude keys get bit-identical models to the historical behaviour. Otherwise it splats the tuple into an `MVector`. It is applied to the bus and both panel quaternions.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_link_q`. Returns `q === nothing ?`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`

**Downstream**

- `callees` → [[vehx.spacecraft_model_spacecraftmodel|SpacecraftModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:104-104`
<!-- vulcan:connections:end -->

## Limitations
The tuple is not normalised or checked for unit length, so a hand-typed quaternion with a small error is stored as-is and skews the link's rotation; the scalar-last convention is assumed without validation.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/example_support.jl` line 98.
