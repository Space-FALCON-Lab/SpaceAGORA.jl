---
id: analysis.types_abstractscenarioconfig
label: AbstractScenarioConfig
kind: struct
source:
  file: src/analysis/verification/telemetry_verification/types.jl
  symbol: AbstractScenarioConfig
  lines:
  - 88
  - 88
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
  type: AbstractScenarioConfig
  units: n/a
  description: Abstract supertype `AbstractScenarioConfig`; no fields.
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

# AbstractScenarioConfig

## Purpose
Common abstract supertype for the two telemetry verification scenario descriptions, `OrbitEventsScenarioConfig` and `TimeAlignedScenarioConfig`, so that manifest parsing, scenario filtering, and the study runner can dispatch on one type while each subtype carries its own comparison-specific fields.

## Design & Implementation
Declared as `abstract type AbstractScenarioConfig end` with no fields or interface methods of its own. Both concrete subtypes are `Base.@kwdef struct`s that share a core of fields (`name`, `planet_name`, `units_x`, `units_y`, `tolerances_quick`, `tolerances_full`, `initial_time::InitialTime`, `spacecraft::SpacecraftConfig`, gravity, n-body, SRP, drag, and `atmosphere_truth`/`calibration` sub-configs, `EI_km`). Code that needs the shared fields accesses them by name, relying on convention rather than an enforced interface.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractScenarioConfig | n/a | — | Abstract supertype `AbstractScenarioConfig`; no fields. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No abstract accessor functions are defined, so a new subtype that omits any of the conventionally shared fields will fail at runtime with a `FieldError` in the runner rather than at definition time. The type carries no parameters, so containers of scenarios are `Vector{AbstractScenarioConfig}` with dynamic dispatch on each element. Field overlap between subtypes is duplicated by hand, so a new shared option must be added in both places.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/types.jl` line 88.
