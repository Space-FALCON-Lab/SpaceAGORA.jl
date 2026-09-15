---
id: gncy.guidance_hooks_guidancehooks
label: GuidanceHooks
kind: struct
source:
  file: src/gnc/guidance/guidance_hooks.jl
  symbol: GuidanceHooks
  lines:
  - 1
  - 94
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: gnc_namespace
  type: Module
  units: n/a
  required: true
  description: Parent GNC namespace that supplies configuration, abstract types, navigation
    geometry, and reference-system utilities.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: guidance_api
  type: Module
  units: n/a
  description: Exported guidance surface covering RPO planning, replanning, planner
    comparison, and aerobraking strategy dispatch.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncy
origin: agent
---

# GuidanceHooks

## Purpose
`GuidanceHooks` is the module that assembles the entire guidance layer. It gathers the RPO path-planning stack, the planner-comparison harness, the replanning logic, and the aerobraking strategy dispatch into one namespace and exports the names the rest of the simulator calls.

## Model & Assumptions
The module treats guidance as a set of files included in dependency order rather than as separate submodules, so every included file shares one namespace and can call the others without qualification. Control-side entry points are reached lazily through `_control_module`, which fetches the sibling `ControlHooks` module from the cached parent module reference at call time. That indirection exists because guidance and control are mutually referential and cannot both be loaded first.

## Design & Implementation
Line 23 aliases `config` to `Structure` and line 24 caches `_PARENT` from `parentmodule`. Three inlined forwarders wrap `asim_ctrl`, `asim_ctrl_rf`, and `control_solarpanels_heatrate` from the control module. The export list spans `calcGuidanceEffect!`, the full family of PSO settings structs, the RRT-Connect and RRT-Star settings, the CHOMP and STOMP planner entry points, the planner-comparison configuration and batch runner, the replanning configuration and decision functions, and the aerobraking strategy tags with their input and output records. Includes run from line 63 through line 93, ordered so that PSO parameters load before sampling, sampling before costs, costs before planning, and the aerobraking strategy files after the solvers they call.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `gnc_namespace` | Module | n/a | yes | Parent GNC namespace that supplies configuration, abstract types, navigation geometry, and reference-system utilities. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `guidance_api` | Module | n/a | — | Exported guidance surface covering RPO planning, replanning, planner comparison, and aerobraking strategy dispatch. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[gnc.guidance_hooks__control_asim_ctrl|_control_asim_ctrl]] · `callers` · call · `src/gnc/guidance/guidance_hooks.jl:29-29`
- `callees` → [[gnc.guidance_hooks__control_asim_ctrl_rf|_control_asim_ctrl_rf]] · `callers` · call · `src/gnc/guidance/guidance_hooks.jl:32-32`
- `callees` → [[gnc.guidance_hooks__control_module|_control_module]] · `callers` · call · `src/gnc/guidance/guidance_hooks.jl:26-26`
- `callees` → [[gnc.guidance_hooks__control_solarpanels_heatrate|_control_solarpanels_heatrate]] · `callers` · call · `src/gnc/guidance/guidance_hooks.jl:35-35`
<!-- vulcan:connections:end -->

## Limitations
Because everything shares one namespace, a name collision between two included files is only caught as a method redefinition at load time. The lazy control-module lookup defers a genuine wiring error until the first guidance call rather than at load. Include order is load-bearing and is not enforced by any declaration, so reordering the include block breaks compilation in ways that point at the wrong file.

## Provenance
Mapped from guidance_hooks.jl lines 1-94.
