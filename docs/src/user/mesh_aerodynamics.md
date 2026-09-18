# Mesh Aerodynamics

The aerodynamic models in `SimulationModel` describe a spacecraft as boxes
(`AerodynamicCoefficientfM`, the Hart closed forms for a rectangular prism
summed over links). When a CAD model of the vehicle exists, the mesh
aerodynamics path derives free-molecular force and moment coefficients from
that geometry instead, then fits an analytic surrogate the integrator
evaluates at run time without touching the mesh.

```julia
using SpaceAGORA

panels = mesh_aero_panels("path/to/spacecraft.glb";
    scale=1.0, rotation_deg=(0, 0, 90))              # the viewer's model transform: long axis along body x
surrogate = fit_mesh_aero_surrogate(panels; degree=8, speed_ratios=(3, 5, 8, 12, 20), verbose=true)
write_mesh_aero_surrogate("output/spacecraft_aero.json", surrogate)

effector = AerodynamicCoefficientMeshSurrogate(surrogate; wall_temperature_k=300.0)
# ... dynamic_effectors=(GravitationalHarmonicsModel(...), effector) in place of AerodynamicCoefficientfM()
```

No GRAM, SPICE or additional package is needed: the readers, the panel method
and the fit use only what the package already carries.

## Panel method

`mesh_aero_panels` reads STL, OBJ or glTF/GLB (embedded buffers, no Draco)
through the same readers the viewer uses (`load_model_triangles`), applies
the same scale, XYZ Euler rotation, centering and `articulations` (parts of
the model posed first, see `articulate_triangles`), and keeps per-facet
centroids, outward normals (from the vertex winding, or `outward_normals=true`
to flip normals toward the outside for convex meshes with mixed winding) and
areas. Coefficients are normalized by the dynamic pressure of the
wind-relative airspeed and `reference_area_m2` (default: total surface area
over four, the mean projected area of a convex body); moments by
`reference_area_m2 * reference_length_m` (default length: the bounding-box
diagonal), about `moment_reference_m` (default: the origin).

`panel_aero_coefficients(panels, vhat, s; sigma_n, sigma_t, tw_ratio)` sums
the Schaaf and Chambré pressure and shear over the facets for one airspeed
direction `vhat` in the mesh frame (the direction the body moves relative to
the gas) and speed ratio `s = V / sqrt(2 R T)`, with normal and tangential
momentum accommodation `sigma_n`, `sigma_t` (1 is fully diffuse) and the
wall-to-freestream temperature ratio `tw_ratio`. Facets behind other facets
along the flow are found with a depth buffer (`panel_shadow_mask`) and
contribute nothing, so a dish in front of a bus shields it and a panel
behind another is dead. The Hart box formulas are the same integrand with the
back faces dropped, and a box mesh reproduces them away from grazing incidence. At exactly
grazing faces the panel model retains tangential shear that the box model
omits.

## Surrogate

The wall temperature enters the free-molecular flux only through the diffuse
re-emission term, linearly in `sqrt(Tw/T)`, so each of the six coefficients
splits exactly as `A(vhat, s) + sqrt(Tw/T) B(vhat, s)`.
`fit_mesh_aero_surrogate` tabulates `A` and `B` over `n_directions`
Fibonacci-spread directions and the listed speed ratios and fits real
spherical harmonics of degree `degree` in the direction times a polynomial of
degree `poly_degree` in `1/s` by least squares. The result is a
`MeshAeroSurrogate`: two `6 x n` coefficient matrices, the normalization, the
moment reference point, the accommodation coefficients, the fitted
speed-ratio range (evaluation clamps to it) and metadata with the residuals,
both over the training samples and over fresh holdout directions. Read the
holdout error against the coefficient scale before choosing a degree. Fit
quality and construction time depend on geometry, shadowing resolution and
the sampled directions and speed ratios. Validate the fit for the vehicle
and flight conditions you will use.

`mesh_aero_coefficients(surrogate, vhat, s; tw_ratio)` evaluates the fitted
coefficients without loading or resampling the mesh. `write_mesh_aero_surrogate` and
`read_mesh_aero_surrogate` store it as JSON.

### Limits and validation

`degree` is limited to `MESH_AERO_MAX_DEGREE` (20): the evaluator keeps its
direction basis in a fixed buffer of that size. The coefficient count uses
checked integer arithmetic to reject overflowing degree combinations. The
limit, `poly_degree >= 0`,
the coefficient shapes and finiteness, the positive normalization, the finite
reference point, accommodation coefficients in `[0, 1]` and
`0 < speed_ratio_min <= speed_ratio_max` are checked by the one
`MeshAeroSurrogate` constructor, so there is no construction route that
skips them. `fit_mesh_aero_surrogate` checks the degrees before it tabulates
anything. Fitting requires at least `(degree + 1)^2` directions and
`poly_degree + 1` distinct speed ratios; a larger total sample count cannot
compensate for too few directions or speeds. The holdout checks fresh
directions at a training speed, so it does not establish interpolation
accuracy between sampled speeds. `read_mesh_aero_surrogate` checks every field of the file
(schema, integer degrees, six numeric coefficient rows of equal length,
numeric normalization and range, a three-component reference point) and
raises an `ArgumentError` naming the file for anything else. A file that
declares a degree above the limit is rejected before its coefficient arrays
are read.

## Runtime effector

`AerodynamicCoefficientMeshSurrogate(Dict(link_index => surrogate);
wall_temperature_k)` carries one surrogate per link, keyed by the link's
integer position in `spacecraft.links`; the single-surrogate form puts a
whole-vehicle fit on the actual root, regardless of its position in that list. Links without a surrogate carry no
aerodynamic load, so a rigid vehicle uses one whole-vehicle mesh, and a
vehicle whose panels articulate under aerobraking guidance uses one mesh per
link, composed through the link poses the way the box model is (with no
shadowing between links).

Per link and per RHS call the effector forms the airspeed direction in that
link's frame, the speed ratio and `wall_temperature_k / T`, evaluates the
surrogate, and scales by the dynamic pressure and the reference area (and
length for the moment). Attitude follows the other aerodynamic effectors:
with `orientation_sim` the propagated root attitude and the stored child
attitudes place the links; without it the root is held in the
velocity-aligned frame (x along the airspeed, z toward nadir) and each
link's configured quaternion selects its incidence. It uses the `wrench`
interface described in [Adding a Force or Torque](custom_effector.md); the
older `calcForceTorque` hook is not provided.

### Reference points and torque

Torque is returned about the root centre of mass in the root body frame, and
only when the attitude is propagated. A link's fitted moment is about its
surrogate's `moment_reference_m`, a point `P` in the mesh (link) frame. The
effector moves it to the link origin `O` with

```text
M_O = M_P + (P - O) x F
```

then rotates it into the root frame through the link attitude and adds the
lever arm of the link's root-frame offset `body.r`. The torque about the root
therefore does not depend on where the fit put its reference point; two
surrogates of the same mesh fitted about different points give the same
motion. What the engine does assume, as for every other effector, is that the
link origin is the link's centre of mass. Position the mesh accordingly:
`center=true` puts the bounding-box centre at the origin, which is only an
approximation of the mass centre, and a model whose origin is elsewhere needs
to be translated before fitting.

### Diagnostics in radial flight

The drag, lift and cross save fields are projections of the total force.
Drag (along the airspeed) is defined whenever there is airspeed. Lift and
cross are built from the orbit normal `r x v`, which vanishes in radial
flight (a vertical descent or ascent, and to within a relative tolerance of
`1e-9` of it), and are then written as zero. The force applied to the
spacecraft is computed before that split and is the same either way: a
vertical descent decelerates under drag exactly as a horizontal pass does.

### What is not modelled

Shadowing between links; per-link atmosphere sampling (one density and
temperature sample per spacecraft, as for the box model's default); thermal
incidence (the heating callback keeps its stored link angles); torque in the
fixed-attitude mode. The effector is not on the automatic effector-level
threading allowlist.

## Validation

`test/unit/dynamics/mesh_aero_tests.jl` checks a sphere against the
closed-form free-molecular drag, a box against the Hart formulas the box
model uses, a plate behind another plate, the harmonic basis, the fit and
its JSON round trip, the moment-reference translation (exact in the panel
method, invariant across fits and hand-composed for an offset, rotated child
link), radial flight, malformed constructor and file inputs, a propagated
attitude run whose rate build-up matches the plate torque, a radial descent
whose energy loss matches the recorded drag work, and a fixed-attitude run of
a cube against the same cube under `AerodynamicCoefficientfM`.
`scripts/dev/aero/fit_mesh_aero_surrogate.jl` fits a model file from the
command line. See the [Mesh Aerodynamics API](../generated/mesh_aerodynamics_api.md)
for the supported constructors and fitting functions.
