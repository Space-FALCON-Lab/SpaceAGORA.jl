---
id: parallel.policy_telemetry_reset_policy_telemetry_bang
label: reset_policy_telemetry!
kind: function
source:
  file: src/parallel/policy/policy_telemetry.jl
  symbol: reset_policy_telemetry!
  lines:
  - 76
  - 76
inputs:
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
  description: Return value of `reset_policy_telemetry!`. Returns `nothing`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# reset_policy_telemetry!

## Purpose
Clears the threading policy's accumulated telemetry and learned state so a benchmark or test starts from a known baseline.

## Design & Implementation
Under `_policy_telemetry_lock` it replaces `ctx.telemetry` with a freshly constructed `PolicyTelemetry()` and empties three dictionaries on the active context: `adaptive_state`, which holds per-source adaptive controllers, and `decision_signature` and `decision_allotment`, which memoise the signature and chosen allotment of each decision site. Emptying rather than reassigning keeps any existing references to those dictionaries valid. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `reset_policy_telemetry!`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/policy_telemetry.jl`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:169-169`

**Downstream**

- `callees` → [[parallel.context__active_policy_context|_active_policy_context]] · `callers` · call · `src/parallel/policy/policy_telemetry.jl:78-78`
- `callees` → [[parcore.types_policytelemetry|PolicyTelemetry]] · `callers` · call · `src/parallel/policy/policy_telemetry.jl:79-79`
<!-- vulcan:connections:end -->

## Limitations
It resets only the context that is active at the moment of the call, so telemetry accumulated under a different policy context survives; the persistent hint file on disk is untouched, and a subsequent run reloads the learned allotments this reset appeared to discard.

## Provenance
Mapped from `src/parallel/policy/policy_telemetry.jl` line 76.
