---
id: parcore.outer_route_metrics_record_outer_route_feedback_bang
label: record_outer_route_feedback!
kind: function
source:
  file: src/parallel/routing/outer_route_metrics.jl
  symbol: record_outer_route_feedback!
  lines:
  - 7
  - 58
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ParallelProfiles namespace supplying OuterRouteFeatures, OuterRouteTuning
    and the signature hierarchy helper.
- id: state
  type: OuterRouteState
  units: n/a
  required: true
  description: Mutable adaptive routing history, updated in place under its own lock.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: feedback
  type: Nothing
  units: n/a
  description: 'Side effect only: per-signature, per-route sample, success, failure
    and elapsed-time statistics are accumulated for future route selection.'
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

# record_outer_route_feedback!

## Purpose
`record_outer_route_feedback!` records what actually happened when a chosen outer route ran, so that later adaptive route selection has empirical runtime and reliability statistics instead of only static thresholds. It is the feedback half of the outer routing loop, mirroring what observation tracking does for inner threading.

## Theory & Math
Each route accumulates sufficient statistics for a mean and variance of elapsed time: sample count, sum and sum of squares. Failures are folded in as a fixed penalty $p$ from the tuning record, so a route that fails contributes $p$ seconds and $p^2$ to the squared sum. The squared sum is floored at the Cauchy-Schwarz bound

$$\sum t^2 \ge \frac{\left(\sum t\right)^2}{n},$$

which keeps the derived variance non-negative even when a caller supplies only the first moment.

## Model & Assumptions
Only the routes `:none`, `:threads` and `:process` are accepted; anything else returns without recording. Success and failure counts are clamped non-negative and their sum must be positive or the call is a no-op. When the caller omits the sum of squares it arrives as `NaN`, and the function substitutes the approximation from the first moment; when a value is supplied, the larger of the supplied and approximated values is kept so an understated input cannot produce a negative variance.

## Design & Implementation
Statistics are written into a hierarchy of signatures rather than a single key: `_outer_route_signature_hierarchy` derives progressively coarser descriptions of the workload from the `OuterRouteFeatures` record, and the same observation is credited to every level. That is what lets an unseen specific workload still fall back on evidence from its general class. The whole update runs under the `OuterRouteState` lock, using `get!` to create the per-signature bucket and the per-route `OuterRouteStats` on first use.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ParallelProfiles namespace supplying OuterRouteFeatures, OuterRouteTuning and the signature hierarchy helper. |
| in | `state` | OuterRouteState | n/a | yes | Mutable adaptive routing history, updated in place under its own lock. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `feedback` | Nothing | n/a | — | Side effect only: per-signature, per-route sample, success, failure and elapsed-time statistics are accumulated for future route selection. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.adaptive_routing__record_campaign_route_feedback_bang|_record_campaign_route_feedback!]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:219-219`

**Downstream**

- `callees` → [[parallel.outer_route_selection__outer_route_signature_hierarchy|_outer_route_signature_hierarchy]] · `callers` · call · `src/parallel/routing/outer_route_metrics.jl:41-41`
- `callees` → [[parallel.outer_route_state_outerroutestats|OuterRouteStats]] · `callers` · call · `src/parallel/routing/outer_route_metrics.jl:48-48`
- `callees` → [[parallel.outer_route_state_outerroutetuning|OuterRouteTuning]] · `callers` · call · `src/parallel/routing/outer_route_metrics.jl:15-15`
<!-- vulcan:connections:end -->

## Limitations
The failure penalty is a single tuning constant, so a fast failure and a timeout are indistinguishable in the statistics. Crediting every level of the signature hierarchy means coarse levels are dominated by whichever specific workloads ran most often, which can bias a first decision for a genuinely different workload in the same class. Elapsed times are wall clock and are never aged out, so history from an unrepresentative machine persists until the state is reset.

## Provenance
Mapped from `src/parallel/routing/outer_route_metrics.jl:7-58`.
