---
id: gnc.bridge_helpers__bridge_aerobraking_integrator_name
label: _bridge_aerobraking_integrator_name
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _bridge_aerobraking_integrator_name
  lines:
  - 135
  - 135
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  type: String
  units: n/a
  description: Return value of `_bridge_aerobraking_integrator_name`.
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

# _bridge_aerobraking_integrator_name

## Purpose
Returns the name of the numerical integrator the legacy aerobraking bridge should request, defaulting to `"Julia"`.

## Design & Implementation
Reads `_bridge_optional_field(args, :integrator, "Julia")` and normalises the result through `String(...)`, so a `Symbol` such as `:Julia` is accepted. The annotated return type is `::String`. Unlike most sibling accessors this one takes no `mission` argument because the integrator choice is purely a runtime setting.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_bridge_aerobraking_integrator_name`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.bridge_helpers__make_aerobraking_runtime_settings|_make_aerobraking_runtime_settings]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:159-159`

**Downstream**

- `callees` → [[gnc.bridge_helpers__bridge_optional_field|_bridge_optional_field]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:136-136`
<!-- vulcan:connections:end -->

## Limitations
The value is not validated against the integrator names the solver policy recognises, so an unknown name is only rejected downstream. `String(...)` throws `MethodError` for numeric or `nothing` values. The `"Julia"` default is a hard-coded literal duplicated from the solver configuration layer rather than shared as a constant.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 135.
