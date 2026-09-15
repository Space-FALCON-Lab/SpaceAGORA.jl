---
id: gnc.lqmpc_rpo_ref_preview
label: rpo_ref_preview
kind: function
source:
  file: src/gnc/control/rpo_mpc/lqmpc.jl
  symbol: rpo_ref_preview
  lines:
  - 125
  - 125
inputs:
- id: plan
  type: RPOPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: t_elapsed_s
  type: Real
  units: n/a
  required: true
  description: Positional argument `t_elapsed_s`.
- id: dt
  type: Real
  units: n/a
  required: true
  description: Positional argument `dt`.
- id: horizon
  type: Int
  units: n/a
  required: true
  description: Positional argument `horizon`.
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
  description: Return value of `rpo_ref_preview`. Returns `out`.
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

# rpo_ref_preview

## Purpose
Extracts the window of the RPO reference trajectory the MPC should track over its horizon, starting from the current elapsed time.

## Design & Implementation
Allocates a six-by-`horizon + 1` matrix and returns it zero when the plan has no reference columns. Otherwise it converts `t_elapsed_s / dt` to a one-based start column, guarding `dt` below by 1e-9 to avoid division by zero and clamping the index into the plan's range. For each preview column it copies position from `plan.r_ref_rtn` and velocity from `plan.v_ref_rtn` at `min(start_idx + j, n_ref)`, so once the plan runs out the last reference point is repeated as a hold.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RPOPlan | n/a | yes | Positional argument `plan`. |
| in | `t_elapsed_s` | Real | n/a | yes | Positional argument `t_elapsed_s`. |
| in | `dt` | Real | n/a | yes | Positional argument `dt`. |
| in | `horizon` | Int | n/a | yes | Positional argument `horizon`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_ref_preview`. Returns `out`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:16-16`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:130-130`
<!-- vulcan:connections:end -->

## Limitations
The reference is assumed to be sampled at exactly `dt`, so a plan generated at a different rate is read at the wrong speed; the end-of-plan hold means the controller keeps driving toward the final point indefinitely with no signal that the plan has been exhausted.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/lqmpc.jl` line 125.
