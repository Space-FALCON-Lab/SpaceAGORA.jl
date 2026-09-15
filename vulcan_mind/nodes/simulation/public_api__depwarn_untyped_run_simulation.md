---
id: simulation.public_api__depwarn_untyped_run_simulation
label: _depwarn_untyped_run_simulation
kind: function
source:
  file: src/simulation/engine/public_api.jl
  symbol: _depwarn_untyped_run_simulation
  lines:
  - 10
  - 10
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
  description: Return value of `_depwarn_untyped_run_simulation`. Returns `nothing`.
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

# _depwarn_untyped_run_simulation

## Purpose
Emits the deprecation warning for callers that still pass an untyped argument to `run_simulation`, telling them to construct a `SimulationConfiguration` first. It exists so the warning text and its cost live outside the hot typed path.

## Design & Implementation
Marked `@noinline` precisely to keep the string interpolation and warning machinery out of the inlined fast path. It calls `Base.depwarn` with the module-level constant `_RUN_SIMULATION_TYPED_BOUNDARY_DEPRECATION` plus `" Got $(typeof(args))."`, tags the warning with the symbol `:run_simulation`, and passes `force=true` so the message appears even when Julia is started without `--depwarn=yes`. It returns `nothing`; the caller is expected to follow it with `_require_simulation_configuration`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_depwarn_untyped_run_simulation`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/public_api.jl:24-24`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because `force=true` bypasses the user's depwarn setting, the message cannot be silenced through the normal Julia flag and will repeat on every untyped call, which is noisy in loops or batch sweeps. The warning is advisory only and does not itself reject the argument, so removing the companion type check would let untyped input through.

## Provenance
Mapped from `src/simulation/engine/public_api.jl` line 10.
