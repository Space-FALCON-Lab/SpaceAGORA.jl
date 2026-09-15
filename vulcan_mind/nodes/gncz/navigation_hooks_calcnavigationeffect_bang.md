---
id: gncz.navigation_hooks_calcnavigationeffect_bang
label: calcNavigationEffect!
kind: function
source:
  file: src/gnc/navigation/navigation_hooks.jl
  symbol: calcNavigationEffect!
  lines:
  - 17
  - 19
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: NavigationHooks module namespace that exports the hook and includes
    the rendezvous geometry and distance files.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: navigation_effect
  type: Nothing
  units: n/a
  description: No return value on a successful extension; the fallback method raises
    a method error naming the unsupported navigation model.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncz
origin: agent
---

# calcNavigationEffect!

## Purpose
`calcNavigationEffect!` is the single dispatch point through which typed periodic callbacks invoke navigation and sensor estimator models. Declaring it in one module gives every navigation model a common signature of model, state, parameters, time, and spacecraft index, and lets the callback machinery call navigation without knowing which estimator a scenario configured.

## Model & Assumptions
The method defined here is deliberately a fallback that throws a method error for its own argument tuple. That converts a missing extension into an immediate, self-describing failure at the first callback invocation rather than a silent no-operation that would leave the estimator state untouched for an entire run. Implementations are expected to mutate state in place and return nothing, matching the exclamation-mark naming convention used across the guidance and control hooks.

## Design & Implementation
The module also re-exports the rendezvous and proximity operations geometry and distance interface, then includes the station geometry, chaser geometry, combined reference geometry, mesh distance, clearance, surface frame, and standoff query files in dependency order. That ordering matters because the geometry types must exist before the distance routines that dispatch on them are parsed. Grouping the hook with the geometry it serves keeps the navigation surface of the GNC module in one namespace.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | NavigationHooks module namespace that exports the hook and includes the rendezvous geometry and distance files. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `navigation_effect` | Nothing | n/a | — | No return value on a successful extension; the fallback method raises a method error naming the unsupported navigation model. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.navigation_hooks_navigationhooks|NavigationHooks]] · `callees` → `callers` · call · `src/gnc/navigation/navigation_hooks.jl:12-12`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/navigation_guidance_callbacks.jl:10-10`
- [[simulation_a.navigation_guidance_callbacks_get_navigation_callbacks|get_navigation_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/navigation_guidance_callbacks.jl:10-10`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The hook carries no return contract beyond mutation, so a model that forgets to write its estimate cannot be detected here. Time is fixed to a double precision scalar and the spacecraft index to a machine integer, which prevents dispatching the hook under automatic differentiation without a widened signature.

## Provenance
Mapped from `src/gnc/navigation/navigation_hooks.jl:1-28`.
