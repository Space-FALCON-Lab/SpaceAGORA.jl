---
id: parcore.types_policytelemetry
label: PolicyTelemetry
kind: struct
source:
  file: src/parallel/policy/types.jl
  symbol: PolicyTelemetry
  lines:
  - 11
  - 66
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Base threading primitives; this file is included first in ParallelPolicy
    and depends on nothing else in the package.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: telemetry
  type: PolicyTelemetry
  units: n/a
  description: Mutable counter block holding decision totals, per-source breakdowns,
    last-decision fields, elapsed-time sums, quantum accounting and persistent hint
    statistics.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parcore
origin: agent
---

# PolicyTelemetry

## Purpose
`PolicyTelemetry` is the counter block that records what the inner-threading policy did. Every decision, every observation and every hint lookup writes into an instance of it, and the whole adaptive system is observable through this one record.

## Model & Assumptions
Fields fall into five groups. Aggregate totals count decisions, threading-enabled decisions, proposed versus dispatched threading and route discards. Per-source pairs count decisions and threading-enabled decisions separately for density, control, multibody and other sources. A block of last-decision fields mirrors the full argument set of the most recent decision, including mode, threshold, item count, budget, the outer-active and heavy-work flags and the resulting choice. Timing fields accumulate total, threaded and serial elapsed nanoseconds. The remaining fields carry quantum accounting and the persistent hint signature, allotment, confidence and regret.

## Design & Implementation
The record is built with `Base.@kwdef` and every field has a default, so a fresh `PolicyTelemetry()` is a zeroed counter block. It is mutable because it is updated in place under `_policy_telemetry_lock` rather than replaced. The same file declares the other policy state types — `AdaptiveControllerState` with its per-source desire and window accumulators, `AdaptiveChoiceStats` for the bandit, `_HintLayerStatsAccumulator` for reporting, and `PolicyContext`, which owns a telemetry instance plus the per-source adaptive state, decision signature and decision allotment dictionaries. Below the type declarations the file creates the module-level singletons: the telemetry lock, the task-local storage key, the global context reference, the persistent-foreach lock and pool registries, the hint lock and state reference, and the spin-barrier pool type with its own registry and lock.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Base threading primitives; this file is included first in ParallelPolicy and depends on nothing else in the package. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `telemetry` | PolicyTelemetry | n/a | — | Mutable counter block holding decision totals, per-source breakdowns, last-decision fields, elapsed-time sums, quantum accounting and persistent hint statistics. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/types.jl`
- [[parallel.policy_telemetry_reset_policy_telemetry_bang|reset_policy_telemetry!]] · `callees` → `callers` · call · `src/parallel/policy/policy_telemetry.jl:79-79`
- [[parallel.types_policycontext|PolicyContext]] · `callees` → `callers` · call · `src/parallel/policy/types.jl:90-90`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The counters are `Int64` and never reset themselves, so a very long campaign must reset explicitly to keep the elapsed sums interpretable. The last-decision fields describe only the most recent call, so under concurrent decisions from several sources they interleave and cannot be attributed. Adding a counter means touching this record, the recorder and the snapshot function together, since none of the three is derived from the others.

## Provenance
Mapped from `src/parallel/policy/types.jl:11-66`.
