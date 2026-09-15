---
id: gnc.bridge_helpers__bridge_aerobraking_dry_mass
label: _bridge_aerobraking_dry_mass
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _bridge_aerobraking_dry_mass
  lines:
  - 107
  - 107
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
  description: Return value of `_bridge_aerobraking_dry_mass`.
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

# _bridge_aerobraking_dry_mass

## Purpose
Resolves the spacecraft dry mass in kilograms for the aerobraking bridge from the typed dynamics model, the mission body definition, or a raw `dry_mass` runtime field, in that priority order.

## Design & Implementation
If `args.dynamics_model.spacecraft` exists and is non-empty, returns `Float64(first(spacecraft).dry_mass)`; only the first spacecraft is consulted. Otherwise, if `mission.body.dry_mass` exists it is returned. Finally `_bridge_required_field(args, :dry_mass)` is used, which throws `ArgumentError` when no source supplies a value. Return type `::Float64`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `mission` | Any | n/a | no | Positional argument `mission` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_bridge_aerobraking_dry_mass`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.bridge_helpers__make_aerobraking_runtime_settings|_make_aerobraking_runtime_settings]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:156-156`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:110-110`
- `callees` → [[gnc.bridge_helpers__bridge_required_field|_bridge_required_field]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:115-115`
<!-- vulcan:connections:end -->

## Limitations
Multi-spacecraft configurations silently use spacecraft 1's dry mass for every vehicle. An empty `spacecraft` vector falls through rather than erroring, which hides misconfigured dynamics models. `first(spacecraft)` assumes the element has a `dry_mass` property and throws otherwise. No positivity check is performed.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 107.
