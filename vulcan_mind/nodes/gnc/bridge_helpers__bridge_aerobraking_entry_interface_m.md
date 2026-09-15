---
id: gnc.bridge_helpers__bridge_aerobraking_entry_interface_m
label: _bridge_aerobraking_entry_interface_m
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _bridge_aerobraking_entry_interface_m
  lines:
  - 57
  - 57
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
  description: Return value of `_bridge_aerobraking_entry_interface_m`.
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

# _bridge_aerobraking_entry_interface_m

## Purpose
Returns the atmospheric entry-interface altitude in metres for the aerobraking bridge, converting from the kilometre value stored in configuration.

## Design & Implementation
Prefers `args.environment_model.EI` when that nested property exists; otherwise requires a top-level `:EI` via `_bridge_required_field`, which throws `ArgumentError` if absent. The kilometre value is converted with `Float64(...)` and multiplied by `1e3`. The `mission` positional argument (default `nothing`) is accepted for signature symmetry with the other `_bridge_aerobraking_*` accessors but is not consulted. Return type is `::Float64`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `mission` | Any | n/a | no | Positional argument `mission` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_bridge_aerobraking_entry_interface_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.bridge_helpers__bridge_aerobraking_exit_interface_m|_bridge_aerobraking_exit_interface_m]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:67-67`
- [[gncz.bridge_helpers__make_aerobraking_runtime_settings|_make_aerobraking_runtime_settings]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:148-148`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:59-59`
- `callees` → [[gnc.bridge_helpers__bridge_required_field|_bridge_required_field]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:61-61`
<!-- vulcan:connections:end -->

## Limitations
There is no fallback to `mission`, so a runtime lacking `EI` in both locations aborts with an `ArgumentError` even when the mission definition could supply one. No range check is applied: a negative or zero altitude propagates silently into guidance interface tests. The unit convention (config in km, output in m) is implicit in the `1e3` literal.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 57.
