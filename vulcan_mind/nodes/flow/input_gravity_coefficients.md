---
id: input.gravity_coefficients
label: Gravity harmonics coefficients
kind: external
inputs: []
outputs:
- id: coefficient_file
  type: CSV
  units: n/a
  description: Spherical-harmonics C and S coefficients, normalisation convention,
    reference radius.
tags:
- master-flow
charts:
- master
origin: agent
---

# Gravity harmonics coefficients

## Purpose
A spherical-harmonics coefficient file — for example the LP165P lunar field or a Mars or Earth model — giving the degree, order, cosine and sine coefficients the `GravitationalHarmonicsModel` evaluates.

## Design & Implementation
Read once when the model is constructed in `src/dynamics/coupled/perturbations.jl`, converted to fully normalised form whatever convention the file used, and memoised by a key built from the path, degree, order, normalisation, J2 source and planet so repeated construction reuses the matrices. The reference radius is inferred from the filename or planet when not given.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `coefficient_file` | CSV | n/a | — | Spherical-harmonics C and S coefficients, normalisation convention, reference radius. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `coefficient_file` → [[flow.environment|Environment sampling]] · `coefficient_file` · dataflow · `src/dynamics/coupled/perturbations.jl`
<!-- vulcan:connections:end -->

## Limitations
The file format is the project's own CSV layout, not a standard such as ICGEM; and the reference-radius inference is a filename heuristic that silently uses the equatorial radius for an unrecognised file.
