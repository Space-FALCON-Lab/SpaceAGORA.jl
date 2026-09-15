---
id: parallel.select_outer_route_select_outer_route_bang
label: select_outer_route!
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: select_outer_route!
  lines:
  - 535
  - 610
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: route_policy
  type: NamedTuple
  units: n/a
  required: true
  description: Normalized profile configuration and adaptive route policy.
- id: default_route
  type: Symbol
  units: n/a
  required: true
  description: Fallback route returned by default_outer_route.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: route
  type: Symbol
  units: n/a
  description: Selected outer execution route recorded in route state and passed to
    campaign setup.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
- routing
charts:
- parallel
origin: agent
---

# select_outer_route!

## Purpose
`select_outer_route!` chooses and records the route used by an adaptive campaign. It combines normalized profile policy, route candidates, current feedback, and a fallback route, then mutates the route state so later campaign calls can observe the decision and its outcomes.

## Theory & Math
Selection is a discrete optimization step over candidate routes. Feature values such as recent duration, success rate, and worker pressure are compared against tuning thresholds; the chosen route is the candidate with the best policy score subject to availability. The state update is a feedback transition, not a change to the spacecraft equations.

## Model & Assumptions
Candidate features must be comparable and recorded with the same units across routes. The route-state object is assumed to be owned by the current campaign context, so concurrent mutation requires the caller’s synchronization policy. A missing or unusable feature vector should fall back to `default_route`.

## Design & Implementation
The function resolves candidates, evaluates the policy, chooses a route, and writes the selected symbol into `OuterRouteState`. `record_outer_route_feedback!` later updates metrics used by subsequent calls. The selected route is consumed by adaptive campaign code and may determine whether `ensure_process_workers!` grows a process pool before sampling begins.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `route_policy` | NamedTuple | n/a | yes | Normalized profile configuration and adaptive route policy. |
| in | `default_route` | Symbol | n/a | yes | Fallback route returned by default_outer_route. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `route` | Symbol | n/a | — | Selected outer execution route recorded in route state and passed to campaign setup. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.outer_route_selection__best_candidate_confidence|_best_candidate_confidence]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:530-530`
- [[simulation.adaptive_routing__campaign_route_plan|_campaign_route_plan]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:150-150`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:603-603`
- `callees` → [[parallel.default_outer_route|default_outer_route]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:543-543`
- `callees` → [[parallel.outer_route_selection__best_candidate_confidence|_best_candidate_confidence]] · `callers` · feedback · `src/parallel/routing/outer_route_selection.jl:588-588`
- `callees` → [[parallel.outer_route_selection__outer_route_signature_hierarchy|_outer_route_signature_hierarchy]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:566-566`
- `callees` → [[parallel.outer_route_selection__outer_route_stats_snapshot_internal|_outer_route_stats_snapshot_internal]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:571-571`
- `callees` → [[parallel.outer_route_selection__under_sampled_candidate|_under_sampled_candidate]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:583-583`
- `callees` → [[parallel.outer_route_selection_outer_route_candidates|outer_route_candidates]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:554-554`
- `callees` → [[parallel.outer_route_state_outerroutetuning|OuterRouteTuning]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:538-538`
<!-- vulcan:connections:end -->

## Limitations
Adaptive selection is only as good as the feedback window and tuning values. A route can be selected from stale metrics, and the mutation is not a transaction across process-pool creation. The function does not guarantee that a selected route will remain available for the entire campaign.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl:528-610`.
