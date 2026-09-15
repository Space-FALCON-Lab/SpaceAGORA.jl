---
id: gnc.bridge_helpers__bridge_aerobraking_max_heat_rate
label: _bridge_aerobraking_max_heat_rate
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _bridge_aerobraking_max_heat_rate
  lines:
  - 84
  - 84
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: mission
  type: Any
  units: n/a
  required: false
  description: Positional argument `mission` (default `nothing`).
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
  description: Return value of `_bridge_aerobraking_max_heat_rate`.
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

# _bridge_aerobraking_max_heat_rate

## Purpose
Returns the maximum allowed convective heat-rate limit (W/m^2 as used by the aerobraking guidance) with a mission-derived default.

## Design & Implementation
Builds `default_limit` from `mission.aerodynamics.heat_rate_limit` when `mission` is not `nothing` and both nested properties exist, else `Inf`. Then reads `_bridge_optional_field(args, :max_heat_rate, default_limit)` and converts with `Float64(...)`. An explicit `max_heat_rate` in `args` therefore overrides the mission value. Return type `::Float64`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `mission` | Any | n/a | no | Positional argument `mission` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_bridge_aerobraking_max_heat_rate`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.bridge_helpers__make_aerobraking_runtime_settings|_make_aerobraking_runtime_settings]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:152-152`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:86-86`
- `callees` → [[gnc.bridge_helpers__bridge_optional_field|_bridge_optional_field]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:90-90`
<!-- vulcan:connections:end -->

## Limitations
An `Inf` default effectively disables the constraint when neither source provides it, and no warning is logged. Units are not checked or converted, so a mission file in kW/m^2 and a runtime override in W/m^2 mix silently. Negative limits are accepted.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 84.
