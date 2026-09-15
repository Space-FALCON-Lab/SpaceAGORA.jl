---
id: gnc.path_retiming_rpo_curvature_from_samples
label: rpo_curvature_from_samples
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_retiming.jl
  symbol: rpo_curvature_from_samples
  lines:
  - 72
  - 72
inputs:
- id: samples
  type: Any
  units: n/a
  required: true
  description: Positional argument `samples`.
- id: s_vals
  type: Any
  units: n/a
  required: true
  description: Positional argument `s_vals`.
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
  description: Return value of `rpo_curvature_from_samples`. Returns `κ`.
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

# rpo_curvature_from_samples

## Purpose
Estimates the unsigned geometric curvature κ (1/m) at each sample of an arc-length-parameterised RPO path, using central finite differences over neighbouring samples. `rpo_retime_path` uses κ to cap speed via `sqrt(a_max / κ)`.

## Theory & Math
With $\mathbf{r}' \approx \dfrac{\mathbf{p}_{j+1} - \mathbf{p}_{j-1}}{s_{j+1} - s_{j-1}}$ and $\mathbf{r}'' \approx \dfrac{2}{s_{j+1}-s_{j-1}}\left(\dfrac{\mathbf{p}_{j+1}-\mathbf{p}_j}{s_{j+1}-s_j} - \dfrac{\mathbf{p}_j-\mathbf{p}_{j-1}}{s_j-s_{j-1}}\right)$, curvature is $\kappa_j = \dfrac{\lVert \mathbf{r}' \times \mathbf{r}'' \rVert}{\lVert \mathbf{r}' \rVert^3}$.

## Design & Implementation
After converting inputs, `κ = zeros(n)` is returned unchanged when `n < 3`. For interior `j in 2:n-1` the forward, backward, and central spacings `ds2`, `ds1`, `ds` are computed; if any is `<= eps(Float64)` the curvature is set to `0.0` and the sample skipped. Otherwise the first derivative `r′` is the central difference over `ds`, and `r″` is twice the difference of forward and backward secant slopes divided by `ds`. Curvature is `norm(cross(r′, r″)) / norm(r′)^3`, guarded by `norm(r′) > eps(Float64)`. Endpoints are filled by copying their neighbours: `κ[1] = κ[2]`, `κ[end] = κ[end-1]`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `samples` | Any | n/a | yes | Positional argument `samples`. |
| in | `s_vals` | Any | n/a | yes | Positional argument `s_vals`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_curvature_from_samples`. Returns `κ`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.path_retiming_rpo_retime_path|rpo_retime_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:141-141`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The second-derivative stencil is not the standard non-uniform-grid formula (it divides by `ds` rather than the mean spacing), so on strongly non-uniform sampling κ is biased. Curvature at a sharp corner scales with sample spacing and can be enormous, driving the retimer to near-zero speed. Endpoint extrapolation by copying can be wrong where the path straightens abruptly. `cross` assumes 3-row input; 2-D paths will throw. No check that `length(s_vals) == size(samples, 2)`.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_retiming.jl` line 72.
