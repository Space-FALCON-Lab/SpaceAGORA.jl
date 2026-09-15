---
id: simulation.setup__effector_satellite_share_budget
label: _effector_satellite_share_budget
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_satellite_share_budget
  lines:
  - 527
  - 527
inputs:
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
- id: budget
  type: Int
  units: n/a
  required: true
  description: Positional argument `budget`.
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
  type: Int
  units: n/a
  description: Return value of `_effector_satellite_share_budget`.
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

# _effector_satellite_share_budget

## Purpose
Divides the inner thread budget among concurrently evaluated satellites to produce the number of threads each satellite's effector loop may claim without oversubscribing the pool.

## Theory & Math
With $n$ satellites and thread budget $b$, concurrency is $c = \max(1, \min(n, b))$ and the per-satellite share is $s = \max\!\left(1, \left\lfloor b / c \right\rfloor\right)$.

## Design & Implementation
`_effector_satellite_share_budget(num_sats::Int, budget::Int)::Int` computes `sat_concurrency = max(1, min(max(1, num_sats), max(1, budget)))`, the number of satellites that can run at once, then returns `max(1, fld(max(1, budget), sat_concurrency))`. Every intermediate is clamped to at least 1 so zero or negative inputs cannot produce a zero divisor or zero share. Called from `_dynamic_effector_thread_decision` before the outer-parallel halving.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `budget` | Int | n/a | yes | Positional argument `budget`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_effector_satellite_share_budget`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__dynamic_effector_thread_decision|_dynamic_effector_thread_decision]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:671-671`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Assumes satellites are evaluated with full concurrency up to the budget, which is only true under `satellite_batch` execution; in serial-over-satellites modes the share is unnecessarily conservative. Integer floor division discards remainder threads rather than distributing them.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 527.
