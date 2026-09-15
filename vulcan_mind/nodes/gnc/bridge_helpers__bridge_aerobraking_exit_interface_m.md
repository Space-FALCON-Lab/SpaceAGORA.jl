---
id: gnc.bridge_helpers__bridge_aerobraking_exit_interface_m
label: _bridge_aerobraking_exit_interface_m
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _bridge_aerobraking_exit_interface_m
  lines:
  - 66
  - 66
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
  description: Return value of `_bridge_aerobraking_exit_interface_m`.
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

# _bridge_aerobraking_exit_interface_m

## Purpose
Returns the atmospheric exit-interface altitude in metres, defaulting to the entry interface when no separate exit altitude `AE` is configured.

## Design & Implementation
Computes `default_m = _bridge_aerobraking_entry_interface_m(args, mission)`, then reads `_bridge_optional_field(args, :AE, default_m / 1e3)`; the default is divided back to kilometres so both branches share the trailing `* 1e3` conversion. The result is coerced with `Float64(...)`. Because the entry accessor is called unconditionally, a missing `EI` field throws `ArgumentError` even when `AE` is present.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `mission` | Any | n/a | no | Positional argument `mission` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_bridge_aerobraking_exit_interface_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.bridge_helpers__make_aerobraking_runtime_settings|_make_aerobraking_runtime_settings]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:149-149`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:68-68`
- `callees` → [[gnc.bridge_helpers__bridge_aerobraking_entry_interface_m|_bridge_aerobraking_entry_interface_m]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:67-67`
- `callees` → [[gnc.bridge_helpers__bridge_optional_field|_bridge_optional_field]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:68-68`
<!-- vulcan:connections:end -->

## Limitations
The default is always computed, so the entry-interface lookup cost and its failure mode apply regardless of whether `AE` exists. Nothing enforces `AE <= EI` or positivity. Rounding through `/ 1e3` then `* 1e3` can introduce a last-bit floating-point difference between the default exit altitude and the entry altitude, which matters for exact equality comparisons downstream.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 66.
