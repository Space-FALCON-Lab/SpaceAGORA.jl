---
id: simulation.rhs_calibration__run_rhs_sweep_bang
label: _run_rhs_sweep!
kind: function
source:
  file: src/simulation/engine/rhs_calibration.jl
  symbol: _run_rhs_sweep!
  lines:
  - 252
  - 252
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: u0
  type: Any
  units: n/a
  required: true
  description: Positional argument `u0`.
- id: dynamic_effectors
  type: Any
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
- id: verbose
  type: Bool
  units: n/a
  required: true
  description: Positional argument `verbose`.
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
  description: Return value of `_run_rhs_sweep!`; mutates `p` in place. Returns `best_plan,
    best_elapsed`.
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

# _run_rhs_sweep!

## Purpose
Times each candidate execution plan by repeatedly evaluating the full RHS `spacecraft_dynamics!` and returns the fastest plan with its mean wall-clock cost.

## Design & Implementation
Signature `_run_rhs_sweep!(p, u0, dynamic_effectors, verbose::Bool)`. Fetches `n_warmup`, `n_timed`, and the candidate list; if only the baseline candidate exists it returns `(nothing, 0.0)` immediately. It allocates `du = zero(u0)` once, then for each candidate sets `p.shared_buffers.rhs_plan_override[] = candidate`, runs `n_warmup` untimed calls, brackets `n_timed` calls with `time_ns()`, and computes `elapsed_mean` in nanoseconds per call. A `try`/`catch`/`finally` block warns and skips a candidate that throws, and always resets the override to `nothing`. The minimum `elapsed_mean` wins; with `verbose` it prints per-candidate and best-plan lines in ms/call. Returns `(best_plan, best_elapsed)` where `best_plan` may be `nothing` if every candidate failed.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `u0` | Any | n/a | yes | Positional argument `u0`. |
| in | `dynamic_effectors` | Any | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `verbose` | Bool | n/a | yes | Positional argument `verbose`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_run_rhs_sweep!`; mutates `p` in place. Returns `best_plan, best_elapsed`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_rhs_calibration_calibrate_rhs_plan_if_needed__calibrate_rhs_plan_if_needed_bang|_calibrate_rhs_plan_if_needed!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:339-339`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:282-282`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:261-261`
- `callees` → [[simulation.rhs_calibration__rhs_calibrate_n_timed|_rhs_calibrate_n_timed]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:254-254`
- `callees` → [[simulation.rhs_calibration__rhs_calibrate_n_warmup|_rhs_calibrate_n_warmup]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:253-253`
- `callees` → [[simulation.rhs_calibration__rhs_plan_candidates|_rhs_plan_candidates]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:255-255`
- `callees` → [[simx.engine_dynamics_rhs_spacecraft_dynamics_bang|spacecraft_dynamics!]] · `callers` · call · `src/simulation/engine/rhs_calibration.jl:276-276`
<!-- vulcan:connections:end -->

## Limitations
Every candidate evaluates the RHS at `t = 0.0` with the initial state `u0`, so plans are ranked on a single point that may not represent later, denser-atmosphere phases. Timing uses the mean, not median, so one GC pause can flip the winner. The sweep mutates `p.shared_buffers.rhs_plan_override` and any RHS-side memo caches as a side effect, and `du` is discarded. Wall-clock timing is sensitive to other load on the machine.

## Provenance
Mapped from `src/simulation/engine/rhs_calibration.jl` line 252.
