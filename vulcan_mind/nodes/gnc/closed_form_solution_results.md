---
id: gnc.closed_form_solution_results
label: results
kind: function
source:
  file: src/gnc/guidance/aerobraking/common/closed_form_solution.jl
  symbol: results
  lines:
  - 403
  - 403
inputs:
- id: solution
  type: Solution
  units: n/a
  required: true
  description: Positional argument `solution`.
- id: t_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `t_cf`.
- id: h_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `h_cf`.
- id: gamma_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `γ_cf`.
- id: v_cf
  type: Any
  units: n/a
  required: true
  description: Positional argument `v_cf`.
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
  description: Return value of `results`.
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

# results

## Purpose
Accumulates one segment of closed-form aerobraking trajectory output into the campaign-wide solution record.

## Design & Implementation
Appends the time, altitude, flight-path angle and velocity arrays onto the corresponding `t_cf`, `h_cf`, `γ_cf` and `v_cf` vectors of `solution.closed_form`. Using `append!` rather than assignment means successive passes extend one continuous history rather than replacing it, which is what lets the campaign be reconstructed end to end afterwards.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `solution` | Solution | n/a | yes | Positional argument `solution`. |
| in | `t_cf` | Any | n/a | yes | Positional argument `t_cf`. |
| in | `h_cf` | Any | n/a | yes | Positional argument `h_cf`. |
| in | `gamma_cf` | Any | n/a | yes | Positional argument `γ_cf`. |
| in | `v_cf` | Any | n/a | yes | Positional argument `v_cf`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `results`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.closed_form_solution_closed_form|closed_form]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:14-14`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The four arrays are appended independently with no length cross-check, so a caller supplying mismatched segment lengths silently desynchronises the history.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/common/closed_form_solution.jl` line 403.
