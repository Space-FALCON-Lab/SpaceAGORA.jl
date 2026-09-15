---
id: gnc.bridge_helpers__bridge_get_solution
label: _bridge_get_solution
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _bridge_get_solution
  lines:
  - 221
  - 221
inputs:
- id: args
  type: Any
  units: n/a
  required: false
  description: Positional argument `args` (default `nothing`).
- id: solution
  type: Any
  units: n/a
  required: false
  description: Keyword argument `solution` (default `nothing`).
- id: cnf
  type: Any
  units: n/a
  required: false
  description: Keyword argument `cnf` (default `nothing`).
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
  description: Return value of `_bridge_get_solution`. Returns `solution` or `getproperty(args,
    :solution)` or `getproperty(cnf, :solution)`.
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

# _bridge_get_solution

## Purpose
Retrieves the trajectory `solution` object the aerobraking bridge should analyse, checking the keyword, the typed `args.solution` field, and finally the control state `cnf.solution`.

## Design & Implementation
Signature `_bridge_get_solution(args=nothing; solution=nothing, cnf=nothing)`. Returns `solution` when non-`nothing`; else `args.solution` when present; else `cnf.solution` when `cnf` has that property. Otherwise throws `ArgumentError("Solution state `solution` not found. ...")`. The three-tier priority lets callers that already hold `cnf` avoid threading the solution separately.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | no | Positional argument `args` (default `nothing`). |
| in | `solution` | Any | n/a | no | Keyword argument `solution` (default `nothing`). |
| in | `cnf` | Any | n/a | no | Keyword argument `cnf` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_bridge_get_solution`. Returns `solution` or `getproperty(args, :solution)` or `getproperty(cnf, :solution)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:12-12`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:10-10`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:13-13`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only `hasproperty`/`getproperty` access is supported, so dictionary-backed `args` are not searched. Field values that are present but `nothing` are returned as-is, deferring the failure to the caller. There is no verification that the returned object is actually an ODE solution type.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 221.
