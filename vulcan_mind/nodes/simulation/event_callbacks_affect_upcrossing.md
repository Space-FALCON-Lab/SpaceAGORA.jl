---
id: simulation.event_callbacks_affect_upcrossing_bang
label: affect_upcrossing!
kind: function
source:
  file: src/simulation/callbacks/event_callbacks.jl
  symbol: affect_upcrossing!
  lines:
  - 145
  - 145
inputs:
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
- id: idx
  type: Int64
  units: n/a
  required: true
  description: Positional argument `idx`.
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
  description: Return value of `affect_upcrossing!`; mutates `integrator` in place.
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

# affect_upcrossing!

## Purpose
Upcrossing handler of the drag-state callback (line 145). When satellite `idx` climbs above the entry interface it flags the satellite as out of the atmosphere, relaxes the integrator's maximum step and tolerances for vacuum propagation, invalidates the vacuum-predicted GRAM cache, and schedules any event-driven thruster controls that should fire at atmospheric exit.

## Design & Implementation
Signature `affect_upcrossing!(integrator, idx::Int64)`. It writes `p.shared_buffers.in_atmosphere[idx] = false` and `in_atmosphere_sample_t[idx] = Float64(integrator.t)`. If `idx <= length(p.shared_buffers.vacuum_gram_caches)` and the cache is not `nothing`, it sets `cache.valid = false` so the next entry rebuilds the cache from fresh state instead of interpolating stale data. `integrator.opts.dtmax` becomes `p.args.integration_tolerances.dt_max_orbit`, and `reltol`/`abstol` are replaced by the pair returned from `_callback_tolerances_for_phase(reltol, abstol, p.args, false)`. Finally `schedule_event_driven_thruster_controls!(integrator, idx)` is invoked. Verbose output is gated by `callback_verbose`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `idx` | Int64 | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `affect_upcrossing!`; mutates `integrator` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/event_callbacks.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:151-151`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:148-148`
- `callees` → [[simulation.assembly__callback_tolerances_for_phase|_callback_tolerances_for_phase]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:161-161`
- `callees` → [[simulation.control_callbacks_schedule_event_driven_thruster_controls_bang|schedule_event_driven_thruster_controls!]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:169-169`
- `callees` → [[simulation.event_callbacks_affect_downcrossing_bang|affect_downcrossing!]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:172-172`
- `callees` → [[simulation.registry_callback_verbose|callback_verbose]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:147-147`
<!-- vulcan:connections:end -->

## Limitations
Because `integrator.opts` is a single object shared across satellites, relaxing tolerances for one exiting satellite loosens them for any satellite still inside the atmosphere. The tolerance adjustment starts from the current `integrator.opts` values, so repeated entry/exit cycles depend on `_callback_tolerances_for_phase` being idempotent; nothing here checks that. No guard exists for `dt_max_orbit <= 0`. Cache invalidation silently no-ops when the cache vector is shorter than `idx`.

## Provenance
Mapped from `src/simulation/callbacks/event_callbacks.jl` line 145.
