---
id: gnc.planner_comparison_rpo_comparison_artifact_slug
label: rpo_comparison_artifact_slug
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_artifact_slug
  lines:
  - 1126
  - 1126
inputs:
- id: value
  type: Any
  units: n/a
  required: true
  description: Positional argument `value`.
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
  description: 'Return value of `rpo_comparison_artifact_slug`. Returns `isempty(token)
    ? "case" : token`.'
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

# rpo_comparison_artifact_slug

## Purpose
Turns an arbitrary label into a filesystem-safe lowercase token for naming generated HTML and CSV artifacts, replacing any run of disallowed characters with a single underscore.

## Design & Implementation
`rpo_comparison_artifact_slug(value)` lowercases `String(value)`, applies `replace(..., r"[^a-z0-9_=-]+" => "_")`, strips leading and trailing underscores with `strip(token, '_')`, and returns `"case"` if the result is empty. Allowed characters are ASCII letters, digits, underscore, equals, and hyphen.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `value` | Any | n/a | yes | Positional argument `value`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_comparison_artifact_slug`. Returns `isempty(token) ? "case" : token`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Non-ASCII letters are replaced rather than transliterated, so labels in other scripts collapse to underscores and may collide as "case". Distinct labels differing only in punctuation (for example "a.b" and "a-b") can map to different or identical slugs unpredictably ("a_b" versus "a-b"). It is not currently called by the output writers in this file, which use the planner symbol directly in file names.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 1126.
