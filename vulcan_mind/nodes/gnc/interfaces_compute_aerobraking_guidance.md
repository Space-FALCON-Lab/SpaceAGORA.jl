---
id: gnc.interfaces_compute_aerobraking_guidance
label: compute_aerobraking_guidance
kind: function
source:
  file: src/gnc/guidance/aerobraking/interfaces.jl
  symbol: compute_aerobraking_guidance
  lines:
  - 27
  - 27
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
  type: Any
  units: n/a
  description: Return value of `compute_aerobraking_guidance`.
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

# compute_aerobraking_guidance

## Purpose
Generic entry point for aerobraking guidance and the fallback method of the dispatch surface. Concrete strategies add methods narrowing the first argument; this definition catches any strategy that has not supplied one and converts the omission into an immediate, explicit error.

## Design & Implementation
The signature `compute_aerobraking_guidance(::AbstractAerobrakingStrategy, ::AerobrakingGuidanceInput)` ignores both arguments and unconditionally executes `throw(MethodError(compute_aerobraking_guidance, (AbstractAerobrakingStrategy, AerobrakingGuidanceInput)))`. Because the declared parameter is the abstract supertype, Julia's dispatch prefers any more specific method, so this body runs only when a subtype without its own method reaches the call. The hand-constructed `MethodError` names the function and the abstract argument tuple, so the message points at the missing method rather than at an argument conversion failure.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `compute_aerobraking_guidance`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.dispatcher_dispatch_aerobraking_guidance|dispatch_aerobraking_guidance]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/dispatcher.jl:11-11`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The thrown `MethodError` reports the abstract types rather than the concrete strategy that was actually passed, which hides which subtype is unimplemented. Defining a catch-all method also means callers never get a compile-time or load-time signal that a strategy lacks an implementation; the failure appears only when guidance is first invoked, potentially well into a long propagation.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/interfaces.jl` line 27.
