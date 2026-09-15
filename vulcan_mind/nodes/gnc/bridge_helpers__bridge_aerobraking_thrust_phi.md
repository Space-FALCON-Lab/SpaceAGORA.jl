---
id: gnc.bridge_helpers__bridge_aerobraking_thrust_phi
label: _bridge_aerobraking_thrust_phi
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _bridge_aerobraking_thrust_phi
  lines:
  - 120
  - 120
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
  description: Return value of `_bridge_aerobraking_thrust_phi`.
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

# _bridge_aerobraking_thrust_phi

## Purpose
Returns the thrust pointing angle `phi` (radians, as consumed by the aerobraking maneuver planner) with a default sourced from the mission engine definition.

## Design & Implementation
Computes `default_phi` as `Float64(mission.engines.ϕ)` when `mission` exposes that nested field, else `0.0`. Reads `_bridge_optional_field(args, :phi, default_phi)` so a runtime `phi` overrides the mission value, and coerces with `Float64(...)`. The mission field uses the Unicode identifier `ϕ` while the runtime field uses ASCII `phi`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `mission` | Any | n/a | no | Positional argument `mission` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_bridge_aerobraking_thrust_phi`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.bridge_helpers__make_aerobraking_runtime_settings|_make_aerobraking_runtime_settings]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:157-157`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:122-122`
- `callees` → [[gnc.bridge_helpers__bridge_optional_field|_bridge_optional_field]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:126-126`
<!-- vulcan:connections:end -->

## Limitations
The two spellings (`ϕ` versus `phi`) are easy to confuse; a runtime container keyed by `:ϕ` is ignored and the default applies. No angle wrapping or unit validation is applied, so degrees passed by mistake propagate. A `0.0` default represents a specific (tangential) thrust direction rather than an absence of data.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 120.
