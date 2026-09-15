---
id: analysis.scenario_builders__make_spacecraft
label: _make_spacecraft
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _make_spacecraft
  lines:
  - 350
  - 350
inputs:
- id: cfg
  type: SpacecraftConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: ic
  type: AbstractInitialCondition
  units: n/a
  required: true
  description: Positional argument `ic`.
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
  description: Return value of `_make_spacecraft`. Returns `make_three_body_spacecraft(`.
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

# _make_spacecraft

## Purpose
Builds the scenario's three-body vehicle from its manifest spacecraft block and an initial condition.

## Design & Implementation
An `@inline` forward to `make_three_body_spacecraft`, mapping every `SpacecraftConfig` field — bus and panel dimensions, masses, panel offset, propellant mass, id, ram-face convention and the three optional attitude quaternions — onto the corresponding keyword.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | SpacecraftConfig | n/a | yes | Positional argument `cfg`. |
| in | `ic` | AbstractInitialCondition | n/a | yes | Positional argument `ic`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_make_spacecraft`. Returns `make_three_body_spacecraft(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_orbit_args|_make_orbit_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:559-559`
- [[analysis.scenario_builders__make_time_aligned_args|_make_time_aligned_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:602-602`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[analysis.example_support_make_three_body_spacecraft|make_three_body_spacecraft]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:351-351`
<!-- vulcan:connections:end -->

## Limitations
`reflection_coefficient` is not forwarded and so takes the builder's default of one for every telemetry scenario; SRP reflectivity is instead carried by the separate `srp_cr` effector parameter.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 350.
