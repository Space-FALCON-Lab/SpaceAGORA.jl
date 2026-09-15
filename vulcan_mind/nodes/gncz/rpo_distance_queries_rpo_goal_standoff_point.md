---
id: gncz.rpo_distance_queries_rpo_goal_standoff_point
label: rpo_goal_standoff_point
kind: function
source:
  file: src/gnc/navigation/rpo_nav/distances/rpo_distance_queries.jl
  symbol: rpo_goal_standoff_point
  lines:
  - 2
  - 7
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: NavigationHooks namespace providing the static vector and norm operations
    used to build the standoff goal.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: goal_point
  type: SVector{3, Float64}
  units: m
  description: Body-frame goal position offset from the given surface point along
    the unit surface normal by the requested standoff distance.
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

# rpo_goal_standoff_point

## Purpose
`rpo_goal_standoff_point` converts a chosen point on the target station surface into the position a chaser should actually fly to. Planners need a goal that sits clear of the structure, and this routine defines that goal as a fixed offset along the local surface normal, which keeps the approach corridor perpendicular to the surface.

## Theory & Math
The goal is $g = p + d\,\hat{n}$ with $\hat{n} = n/\lVert n \rVert$, where $p$ is the surface point, $n$ the supplied surface normal, and $d$ the standoff distance. Normalising inside the routine means the caller may pass an unnormalised normal, such as the raw difference vector produced by the point-cloud normal estimator, without introducing a scale error into the offset.

## Model & Assumptions
The surface point, the normal, and the resulting goal are all expressed in the station body frame, and the standoff distance is interpreted along the normal only, so no lateral bias is applied. A normal whose magnitude is at or below machine epsilon is rejected with an argument error rather than being replaced by a default direction, because a degenerate normal indicates that the upstream geometry query itself failed.

## Design & Implementation
Inputs are converted to static three-vectors so the computation is allocation free and can be called inside a planner scoring loop. The standoff distance is widened to double precision, allowing integer or rational arguments from configuration without a separate conversion at the call site.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | NavigationHooks namespace providing the static vector and norm operations used to build the standoff goal. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `goal_point` | SVector{3, Float64} | m | — | Body-frame goal position offset from the given surface point along the unit surface normal by the requested standoff distance. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/navigation/rpo_nav/distances/rpo_distance_queries.jl:6-6`
<!-- vulcan:connections:end -->

## Limitations
The routine does not check that the resulting goal is itself clear of the station, so a normal pointing into a concave pocket can place the goal inside the keepout volume. It also ignores chaser dimensions, so the standoff must already include any body half extent the caller cares about.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/distances/rpo_distance_queries.jl:1-8`.
