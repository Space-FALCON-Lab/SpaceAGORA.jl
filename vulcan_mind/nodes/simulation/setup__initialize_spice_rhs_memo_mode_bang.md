---
id: simulation.setup__initialize_spice_rhs_memo_mode_bang
label: _initialize_spice_rhs_memo_mode!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_spice_rhs_memo_mode!
  lines:
  - 1417
  - 1417
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  description: Return value of `_initialize_spice_rhs_memo_mode!`; mutates `p` in
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

# _initialize_spice_rhs_memo_mode!

## Purpose
Captures whether the single-entry SPICE RHS memo is enabled for this run, so the ephemeris samplers read a `Ref` rather than the environment on every call.

## Design & Implementation
Assigns `_spice_rhs_memo_enabled()`, which parses the corresponding `SPACEAGORA_*` variable, into `shared_buffers.spice_rhs_memo_enabled[]`. Returns `nothing`. The solar and third-body samplers consult this flag before touching the memo's lock.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_spice_rhs_memo_mode!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:201-201`

**Downstream**

- `callees` → [[simulation.setup__spice_rhs_memo_enabled|_spice_rhs_memo_enabled]] · `callers` · call · `src/simulation/engine/setup.jl:1418-1418`
<!-- vulcan:connections:end -->

## Limitations
Captured once at setup; toggling the environment variable during a run has no effect, which is deliberate but differs from the older read-at-use behaviour some scripts assumed.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1417.
