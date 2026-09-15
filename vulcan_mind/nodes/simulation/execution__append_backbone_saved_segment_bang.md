---
id: simulation.execution__append_backbone_saved_segment_bang
label: _append_backbone_saved_segment!
kind: function
source:
  file: src/simulation/engine/execution.jl
  symbol: _append_backbone_saved_segment!
  lines:
  - 67
  - 67
inputs:
- id: times_acc
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `times_acc`.
- id: data_acc
  type: Vector{SimulationModel.SaveData}
  units: n/a
  required: true
  description: Positional argument `data_acc`.
- id: sol
  type: Any
  units: n/a
  required: true
  description: Positional argument `sol`.
- id: save_fields
  type: Any
  units: n/a
  required: true
  description: Positional argument `save_fields`.
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
  description: Return value of `_append_backbone_saved_segment!`; mutates `times_acc`
    in place. Returns `nothing`.
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

# _append_backbone_saved_segment!

## Purpose
Reconstructs saved-output snapshots from a solved segment when the `:gravity_backbone_split` mode is active. In that mode the `SavingCallback` cannot be relied upon, so this routine walks the solution's own time grid and appends one `SaveData` snapshot per accepted step into the accumulators that later feed the results DataFrame.

## Design & Implementation
Arguments are the accumulators `times_acc::Vector{Float64}` and `data_acc::Vector{SimulationModel.SaveData}`, the segment solution `sol`, the resolved `save_fields`, and the parameter object `p`. If `length(sol.t) <= 1` it returns immediately, so a segment that produced only its initial point contributes nothing. Otherwise it forms `integrator_view = (p=p,)`, a minimal NamedTuple mimicking the integrator interface, and for `idx in 2:length(sol.t)` pushes `Float64(sol.t[idx])` and `SimulationCallbacks._save_snapshot(save_fields, sol.u[idx], t_sample, integrator_view)`. The loop is `@inbounds` and the function is `@inline`; it mutates both accumulator vectors in place and returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `times_acc` | Vector{Float64} | n/a | yes | Positional argument `times_acc`. |
| in | `data_acc` | Vector{SimulationModel.SaveData} | n/a | yes | Positional argument `data_acc`. |
| in | `sol` | Any | n/a | yes | Positional argument `sol`. |
| in | `save_fields` | Any | n/a | yes | Positional argument `save_fields`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_append_backbone_saved_segment!`; mutates `times_acc` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:359-359`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/execution.jl:77-77`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/engine/execution.jl:78-78`
- `callees` → [[simulation.save_fields__save_snapshot|_save_snapshot]] · `callers` · call · `src/simulation/engine/execution.jl:79-79`
<!-- vulcan:connections:end -->

## Limitations
Index 1 is always skipped on the assumption that it duplicates the previous segment's last point, so the very first segment's initial condition is never saved. No duplicate check is performed against `times_acc[end]`, unlike `_append_checkpoint_saved_segment!`. The `(p=p,)` view only satisfies snapshot functions that access `integrator.p`; a save field that reads `integrator.t` or `integrator.u` will throw a `FieldError`.

## Provenance
Mapped from `src/simulation/engine/execution.jl` line 67.
