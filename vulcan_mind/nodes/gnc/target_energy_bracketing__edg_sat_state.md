---
id: gnc.target_energy_bracketing__edg_sat_state
label: _edg_sat_state
kind: function
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: _edg_sat_state
  lines:
  - 180
  - 180
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: i
  type: Int
  units: n/a
  required: true
  description: Positional argument `i`.
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
  description: 'Return value of `_edg_sat_state`. Returns `hasproperty(u, :sc) ? u.sc[i]
    : u`.'
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

# _edg_sat_state

## Purpose
`_edg_sat_state` extracts the state sub-object for spacecraft `i` from the integrator state `u`, supporting both the multi-spacecraft layout (a container with an `sc` property indexed by spacecraft) and a single flat state vector. It lets the bracketing logic be written once for both cases.

## Design & Implementation
Declared `@inline` with signature `(u, i::Int)`. It returns `u.sc[i]` when `hasproperty(u, :sc)` is true and otherwise returns `u` itself. The result is passed to `_edg_pos_vel_mass` and to the control module's heat-load lookup. No copying is performed; the returned object aliases the integrator state.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_sat_state`. Returns `hasproperty(u, :sc) ? u.sc[i] : u`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.target_energy_bracketing__edg_run_target_energy_bracketing_bang|_edg_run_target_energy_bracketing!]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:281-281`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
In the single-state fallback the index `i` is ignored, so calling with `i > 1` on a flat state silently returns the first (only) spacecraft. `hasproperty` is evaluated at runtime and is not type-stable across the two layouts, which can introduce dynamic dispatch in the guidance hot path. No check confirms that `u.sc` has at least `i` elements.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl` line 180.
