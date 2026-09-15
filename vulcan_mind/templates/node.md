---
id: REPLACE.node_id
label: REPLACE_symbol_name
kind: function
source:
  file: REPLACE/path/to/source.ext
  symbol: REPLACE_symbol_name
  lines: [1, 1]
inputs:
  - id: REPLACE_input
    type: REPLACE_Type
    units: "REPLACE units or n/a"
    required: true
    description: REPLACE — what this input is, in the caller's terms.
outputs:
  - id: REPLACE_output
    type: REPLACE_Type
    units: "REPLACE units or n/a"
    description: REPLACE — what this output is.
tags: []
charts: [master]
origin: agent
---

# REPLACE_symbol_name

## Purpose
One or two sentences: what this unit is responsible for, in the system's terms.

## Theory & Math
The governing equations, in LaTeX. Define every symbol. Delete this section only
if the unit genuinely involves no mathematics.

$$
y = f(x)
$$

where $x$ is ... and $y$ is ...

## Model & Assumptions
- Every assumption that would change results if violated.
- State the regime in which each holds.

## Design & Implementation
How it actually works: algorithm, solver, data structures, control flow.
Reference real identifiers and real line ranges from the source.

## Interface (ICD)
<!-- vulcan:icd:begin -->
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
<!-- vulcan:connections:end -->

## Limitations
Known failure modes, valid ranges, numerical caveats. Be specific: state the
tolerance, the divergence condition, the regime where the model stops holding.

## Provenance
Mapped from `REPLACE/path/to/source.ext`.
