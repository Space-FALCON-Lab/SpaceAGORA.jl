---
id: flow.parallel
label: Parallel routing & thread policy
kind: group
inputs:
- id: policy_queries
  type: decision requests
  units: n/a
  description: Should this loop run threaded, and with what allotment.
outputs:
- id: thread_policy
  type: PolicyDecision
  units: n/a
  description: Use-threads verdict, allotment and mode, plus telemetry.
tags:
- master-flow
charts:
- master
origin: agent
opens: parallel
---

# Parallel routing & thread policy

## Purpose
Decides, for every parallelisable loop — density callbacks, effector sums, campaign samples — whether to thread it and how many workers to use, learning from measured timings and persisting what it learned across runs.

## Design & Implementation
`src/parallel/` provides the profile definitions (R0 through R5), the outer-route bandit that picks serial, threaded or process backends for campaigns, the inner thread-policy decision with adaptive allotments, persistent hints keyed by workload signature, per-context telemetry, and the Distributed worker pool bootstrapped with SpaceAGORA, GRAMSuite and kernels.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `policy_queries` | decision requests | n/a | — | Should this loop run threaded, and with what allotment. |
| out | `thread_policy` | PolicyDecision | n/a | — | Use-threads verdict, allotment and mode, plus telemetry. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.campaigns|Campaigns]] · `policy_queries` → `policy_queries` · feedback · `src/simulation/campaigns/adaptive_routing.jl`

**Downstream**

- `thread_policy` → [[flow.callbacks|Integration callbacks]] · `thread_policy` · dataflow · `src/simulation/callbacks/density_callbacks/config.jl`
- `thread_policy` → [[flow.campaigns|Campaigns]] · `thread_policy` · dataflow · `src/simulation/campaigns/adaptive_routing.jl`
- `thread_policy` → [[flow.rhs|Dynamics right-hand side]] · `thread_policy` · dataflow · `src/simulation/engine/setup.jl`
<!-- vulcan:connections:end -->

## Limitations
Learned allotments are keyed by a coarse workload signature that includes the machine label only if set, so runs on different machines pollute one history; the persistent hint file is saved only at exit.
