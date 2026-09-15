---
id: output.checkpoint
label: Checkpoint (.jls + manifest)
kind: external
inputs:
- id: checkpoint
  type: serialized state + TOML
  units: n/a
  description: Written at each segment boundary.
outputs: []
tags:
- master-flow
charts:
- master
origin: agent
---

# Checkpoint (.jls + manifest)

## Purpose
A resumable snapshot of a running simulation: the integrator state, time and solver mode serialized to disk, with a TOML manifest recording its size and digest, written after every checkpoint segment.

## Design & Implementation
Produced by `_write_checkpoint!` when `checkpoint_enabled` is set with a positive interval; consumed by `_load_checkpoint` when `resume_from_checkpoint` is set, and cleared by `_clear_checkpoint!` on a clean finish. A `gravity_backbone_split` checkpoint can only be resumed in that mode.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `checkpoint` | serialized state + TOML | n/a | — | Written at each segment boundary. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.write_results|Write results & checkpoints]] · `checkpoint` → `checkpoint` · dataflow · `src/io/serialization/io_serialization.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Loading trusts the payload after a key check without validating the digest; the format is Julia-version specific.
