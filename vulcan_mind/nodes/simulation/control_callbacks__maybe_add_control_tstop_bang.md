---
id: simulation.control_callbacks__maybe_add_control_tstop_bang
label: _maybe_add_control_tstop!
kind: function
source:
  file: src/simulation/callbacks/control_callbacks.jl
  symbol: _maybe_add_control_tstop!
  lines:
  - 1
  - 1
inputs:
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
- id: tstop
  type: Float64
  units: n/a
  required: true
  description: Positional argument `tstop`.
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
  type: Nothing
  units: n/a
  description: Return value of `_maybe_add_control_tstop!`; mutates `integrator` in
    place. Returns `nothing`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _maybe_add_control_tstop!

## Purpose
Registers a burn boundary time as a solver stop point, but only when doing so is both meaningful and legal for the current integrator. Mutates the integrator's tstop set and returns `nothing`.

## Design & Implementation
Marked `@inline`. It first rejects any `tstop` that is not finite or is not strictly greater than `integrator.t`, since the solver cannot stop in the past. It then queries `DiffEqBase.get_tstops_max(integrator)` inside a `try`/`catch` that falls back to `Inf` when the integrator type does not implement that query. The stop is added through `DiffEqBase.add_tstop!` only when `tstop <= current_max + 1e-9`, the tolerance absorbing floating-point drift at the span endpoint, and only when `applicable(DiffEqBase.add_tstop!, integrator, tstop)` confirms the method exists for this integrator type.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `tstop` | Float64 | n/a | yes | Positional argument `tstop`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_maybe_add_control_tstop!`; mutates `integrator` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/control_callbacks.jl`
- [[simulation.control_callbacks__register_control_tstops_bang|_register_control_tstops!]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:20-20`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The `catch` that returns `Inf` is unconditional, so a genuine integrator failure during the tstops query is silently reinterpreted as an unbounded span and the stop is added anyway. The `1e-9` second slack is absolute, so for very large mission times it is smaller than the representable spacing and effectively becomes an exact comparison. Rejections are silent: a burn time beyond the integration span produces no warning, and the manoeuvre boundary is simply not resolved.

## Provenance
Mapped from `src/simulation/callbacks/control_callbacks.jl` line 1.
