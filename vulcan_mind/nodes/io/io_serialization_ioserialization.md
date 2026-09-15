---
id: io.io_serialization_ioserialization
label: IOSerialization
kind: module
source:
  file: src/io/serialization/io_serialization.jl
  symbol: IOSerialization
  lines:
  - 1
  - 1
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
  description: Value produced by this symbol.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- io
charts:
- io
origin: agent
---

# IOSerialization

## Purpose
The module holding checkpoint persistence, giving a long simulation a crash-safe way to write its state and resume from it.

## Design & Implementation
Imports SHA for digests, TOML for the human-readable manifest, `Serialization` for the binary payload and `Dates` for the creation stamp, and takes checkpoint paths from `IOConfig`. It exports five entry points: the atomic file writer, the digest helper, and write, load and clear operations over a checkpoint pair. Splitting each checkpoint into a serialized data file and a TOML manifest carrying its size and SHA-256 means a reader can detect truncation or corruption before attempting to deserialize.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.io|IOConfig]] · `api` → `module_api` · call · `src/io/serialization/io_serialization.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Loading does not verify the payload against the digest the manifest records, so the integrity information is written but not enforced on the resume path.

## Provenance
Mapped from `src/io/serialization/io_serialization.jl` line 1.
