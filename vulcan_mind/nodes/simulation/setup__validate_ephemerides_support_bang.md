---
id: simulation.setup__validate_ephemerides_support_bang
label: _validate_ephemerides_support!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _validate_ephemerides_support!
  lines:
  - 108
  - 108
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
  type: Nothing
  units: n/a
  description: Return value of `_validate_ephemerides_support!`; mutates `args` in
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

# _validate_ephemerides_support!

## Purpose
Rejects scenario configurations that pair the analytic `SimpleEphemeridesModel` with effectors needing real planetary positions (active N-body gravity or solar radiation pressure), steering users to `SpiceEphemeridesModel`.

## Design & Implementation
Reads `args.environment_model.ephemerides_model`; when it `isa SimulationModel.SimpleEphemeridesModel` it checks `_has_active_nbody_effector(dynamic_effectors)` and then `_has_active_srp_effector(dynamic_effectors)`, throwing a distinct `ArgumentError` for each with a message recommending `SpiceEphemeridesModel()`. Any other ephemerides model passes without inspection. Returns `nothing`; no mutation despite the `!` suffix.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_validate_ephemerides_support!`; mutates `args` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1687-1687`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:166-166`

**Downstream**

- `callees` → [[environment.simple_ephemerides_spiceephemeridesmodel|SpiceEphemeridesModel]] · `callers` · call · `src/simulation/engine/setup.jl:114-114`
- `callees` → [[simulation.setup__has_active_nbody_effector|_has_active_nbody_effector]] · `callers` · call · `src/simulation/engine/setup.jl:111-111`
- `callees` → [[simulation.setup__has_active_srp_effector|_has_active_srp_effector]] · `callers` · call · `src/simulation/engine/setup.jl:117-117`
<!-- vulcan:connections:end -->

## Limitations
Only the two named effector classes are checked; other effectors that quietly depend on ephemerides (for example custom tidal models) are not caught. The N-body check runs first, so a configuration with both problems reports only the N-body one.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 108.
