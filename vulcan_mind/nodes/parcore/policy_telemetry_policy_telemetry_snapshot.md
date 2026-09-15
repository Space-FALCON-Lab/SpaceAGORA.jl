---
id: parcore.policy_telemetry_policy_telemetry_snapshot
label: policy_telemetry_snapshot
kind: function
source:
  file: src/parallel/policy/policy_telemetry.jl
  symbol: policy_telemetry_snapshot
  lines:
  - 87
  - 151
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ParallelPolicy namespace supplying the telemetry lock and the active
    policy context accessor.
- id: telemetry
  type: PolicyTelemetry
  units: n/a
  required: true
  description: Mutable telemetry record held by the active policy context, read field
    by field under the telemetry lock.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: snapshot
  type: NamedTuple
  units: n/a
  description: Immutable copy of the policy telemetry counters, safe to inspect, log
    or compare after the measured region has finished.
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

# policy_telemetry_snapshot

## Purpose
`policy_telemetry_snapshot` returns a consistent copy of the inner-threading telemetry for the active policy context. Tests and performance reports read it to assert how many decisions were made, how many actually threaded, how much time went to threaded versus serial regions, and what the persistent hint layer contributed.

## Model & Assumptions
The snapshot is taken under `_policy_telemetry_lock`, so it reflects one coherent instant rather than a torn read across a live counter set. It is scoped to the active policy context, which means a snapshot taken inside `with_policy_context` describes only that scope's work; outside a scope it describes the global context.

## Design & Implementation
The file holds the write and read ends of telemetry. `_record_policy_decision!` takes the full decision tuple — source, mode, threshold, item count, budget, adaptive flag, desire, allotment, the outer-active, allow-with-outer, heavy-only and heavy-work flags, the resulting threading choice, and the hint signature, allotment, confidence and regret — and updates both the aggregate counters and the per-source breakdown for density, control, multibody and other sources. It separates the proposed threading count from the dispatched count so that a route override is visible instead of appearing as a policy that declined to thread. `reset_policy_telemetry!` restores a fresh record, and `policy_telemetry_snapshot` copies the fields out.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ParallelPolicy namespace supplying the telemetry lock and the active policy context accessor. |
| in | `telemetry` | PolicyTelemetry | n/a | yes | Mutable telemetry record held by the active policy context, read field by field under the telemetry lock. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `snapshot` | NamedTuple | n/a | — | Immutable copy of the policy telemetry counters, safe to inspect, log or compare after the measured region has finished. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:452-452`

**Downstream**

- `callees` → [[parallel.context__active_policy_context|_active_policy_context]] · `callers` · call · `src/parallel/policy/policy_telemetry.jl:89-89`
<!-- vulcan:connections:end -->

## Limitations
Counters are absolute since the last reset, so comparing two regions requires snapshotting before and after and differencing. The per-source breakdown recognises a fixed set of source categories and folds everything else into a single other bucket, so a new call site is invisible until it is added. Because the snapshot is a copy taken under a lock, sampling it at high frequency contends with the decision path it is measuring.

## Provenance
Mapped from `src/parallel/policy/policy_telemetry.jl:87-151`.
