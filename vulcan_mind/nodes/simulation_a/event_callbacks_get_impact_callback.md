---
id: simulation_a.event_callbacks_get_impact_callback
label: get_impact_callback
kind: function
source:
  file: src/simulation/callbacks/event_callbacks.jl
  symbol: get_impact_callback
  lines:
  - 3
  - 35
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Number of spacecraft in the run; sets the length of the vector-valued
    root function and of the active-flag array consulted on a crossing.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: impact_callback
  type: VectorContinuousCallback
  units: n/a
  description: Root-finding callback that deactivates a spacecraft on surface impact
    and terminates the solve once every spacecraft is inactive.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation-a
origin: agent
---
# get_impact_callback

## Purpose
`get_impact_callback` builds the terminal event that detects surface impact for each spacecraft. It is the one callback `get_callbacks` installs unconditionally, including in gravity-backbone split mode, because a solve must never continue integrating a vehicle that has already reached the ground.

## Theory & Math
The event function for spacecraft $i$ is the signed altitude margin above the impact shell,

$$g_i(u,t) = \lVert r_{ii}^{(i)} \rVert - R_{p,e} - h_{\mathrm{impact}}$$

with $h_{\mathrm{impact}} = 50{,}000$ m held in `IMPACT_ALTITUDE_M` and $R_{p,e}$ the planet equatorial radius. Impact is the downcrossing $g_i \to 0^-$, which the solver isolates by root-finding within the step rather than by post-hoc detection, so the reported impact time is accurate to the root-finder tolerance rather than to the step size.

## Model & Assumptions
Impact is defined against a spherical shell at the planet equatorial radius plus a fifty-kilometre margin rather than against a terrain model, so the event is conservative for oblate bodies and for landing sites above the reference datum. The margin is a module-level constant, which means all spacecraft in a run share one impact altitude. Only spacecraft still flagged active in `p.is_active` react to a crossing, so a vehicle that has already impacted cannot retrigger the event.

## Design & Implementation
The condition writes one residual per spacecraft into the solver-provided output vector, making this a `VectorContinuousCallback` with `num_sats` roots. Only the downcrossing affect is supplied; the upcrossing slot is `nothing`. On a crossing the spacecraft is marked inactive, and in gravity-backbone split mode its velocity component is zeroed in place so the frozen vehicle contributes nothing further to the backbone state. When every flag has gone false the callback prints an explicit `termination_cause=impact` line — the `Terminated` retcode is shared with the orbit-count stop, so the cause is always stated — and calls `terminate!`. Verbose per-spacecraft reporting is gated on the run's verbosity setting.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `num_sats` | Int | n/a | yes | Number of spacecraft in the run; sets the length of the vector-valued root function and of the active-flag array consulted on a crossing. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `impact_callback` | VectorContinuousCallback | n/a | — | Root-finding callback that deactivates a spacecraft on surface impact and terminates the solve once every spacecraft is inactive. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:152-152`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:16-16`
- `callees` → [[simulation.event_callbacks_affect_downcrossing_bang|affect_downcrossing!]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:12-12`
- `callees` → [[simulation.event_callbacks_condition_bang|condition!]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:4-4`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:8-8`
- `callees` → [[simulation.registry_callback_verbose|callback_verbose]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:15-15`
<!-- vulcan:connections:end -->

## Limitations
The fixed fifty-kilometre shell will report impact well above ground for high-altitude terrain and well below the surface for a body whose polar radius differs sharply from the equatorial one. Deactivation is a mutation of shared run state, so it is not safe to evaluate the affect concurrently. Terminating on the last active spacecraft means a run configured with an already-inactive vehicle set can stop on the very first crossing.

## Provenance
Mapped from `src/simulation/callbacks/event_callbacks.jl:2-35`.
