# Mesh Aerodynamics

The built-in aerodynamic models describe a spacecraft as boxes
(`AerodynamicCoefficientfM` sums the Hart closed forms for a rectangular prism
over the links). When a CAD model of the vehicle exists, the mesh aerodynamics
path derives free-molecular force and moment coefficients from that geometry
instead, then fits an analytic surrogate that the integrator evaluates at run
time without touching the mesh.

With a model file of your own in place of the path, the whole workflow is:

```julia
using SpaceAGORA

panels = mesh_aero_panels("path/to/spacecraft.glb";
    scale=1.0, rotation_deg=(0, 0, 90))              # the viewer's model transform: long axis along body x
surrogate = fit_mesh_aero_surrogate(panels; degree=8, speed_ratios=(3, 5, 8, 12, 20), verbose=true)
write_mesh_aero_surrogate("output/spacecraft_aero.json", surrogate)

effector = AerodynamicCoefficientMeshSurrogate(surrogate; wall_temperature_k=300.0)
# ... dynamic_effectors=(GravitationalHarmonicsModel(...), effector) in place of AerodynamicCoefficientfM()
```

It needs no GRAM, no SPICE and no extra package: the model readers, the panel
method and the fit are all part of SpaceAGORA.

## Panel method

`mesh_aero_panels` reads STL, OBJ or glTF/GLB (with embedded buffers; Draco
compression is not supported) through the same readers the viewer uses
(`load_model_triangles`). It applies the viewer's model transform: the scale,
the XYZ Euler rotation in degrees, centering (`center=true` by default) and
`articulations`, the parts of the model posed first (see
`articulate_triangles`). It keeps each facet's centroid, outward normal and
area, and drops facets smaller than `min_area_m2` (default 0). The normals come
from the vertex winding; `outward_normals=true` instead flips every normal that
points toward the mesh's centroid, which repairs mixed winding on a convex mesh
but is wrong on concave parts such as the inside of a dish.

Coefficients are normalized by the dynamic pressure of the wind-relative
airspeed and by `reference_area_m2`, which defaults to the total surface area
over four (the mean projected area of a convex body). Moments are normalized
by `reference_area_m2 * reference_length_m`, where the length defaults to the
bounding-box diagonal, and are taken about `moment_reference_m`, which
defaults to the origin.

`panel_aero_coefficients(panels, vhat, s; sigma_n, sigma_t, tw_ratio)` sums
the Schaaf and Chambré pressure and shear over the facets. `vhat` is the
airspeed direction in the mesh frame (the direction the body moves relative to
the gas), and `s = V / sqrt(2 R T)` is the speed ratio. `sigma_n` and
`sigma_t` are the normal and tangential momentum accommodation coefficients
(1 is fully diffuse, the default), and `tw_ratio` is the wall-to-freestream
temperature ratio. Facets that lie behind other facets along the flow are
found with a depth buffer (`panel_shadow_mask`) and contribute nothing, so a
dish in front of a bus shields it and a panel behind another adds no load. The
Hart box formulas are the same integrand with the back faces dropped, and a
box mesh reproduces them away from grazing incidence. At exactly grazing faces
the panel model keeps tangential shear that the box model omits.

## Surrogate

The wall temperature enters the free-molecular flux only through the diffuse
re-emission term, linearly in `sqrt(Tw/T)`, so each of the six coefficients
splits exactly as `A(vhat, s) + sqrt(Tw/T) B(vhat, s)`.
`fit_mesh_aero_surrogate` tabulates `A` and `B` over `n_directions` (default
1200) Fibonacci-spread directions and the listed speed ratios (default
`(3, 5, 8, 12, 20)`). It then fits, by least squares, real spherical harmonics
of degree `degree` (default 8) in the direction times a polynomial of degree
`poly_degree` (default 2) in `1/s`.

The result is a `MeshAeroSurrogate`: two `6 x n` coefficient matrices, the
normalization, the moment reference point, the accommodation coefficients,
the fitted speed-ratio range (evaluation clamps to it) and metadata with the
fit residuals, both over the training samples and over `holdout_directions`
(default 200) fresh directions. Read the holdout error against the coefficient
scale before choosing a degree. Fit quality and the time the fit takes depend
on the geometry, the shadowing resolution and the sampled directions and speed
ratios, so validate the fit for the vehicle and flight conditions you will
use.

`mesh_aero_coefficients(surrogate, vhat, s; tw_ratio)` evaluates the fitted
coefficients without loading or resampling the mesh.
`write_mesh_aero_surrogate` and `read_mesh_aero_surrogate` store a surrogate
as JSON.

### Limits and validation

`degree` is limited to `MESH_AERO_MAX_DEGREE` (20), because the evaluator
keeps its direction basis in a fixed buffer of that size. The coefficient
count uses checked integer arithmetic, so degree combinations that would
overflow are rejected. The single `MeshAeroSurrogate` constructor checks the
degree limit, `poly_degree >= 0`, the coefficient shapes and finiteness, the
positive normalization, the finite reference point, accommodation coefficients
in `[0, 1]` and `0 < speed_ratio_min <= speed_ratio_max`, so no way of
building a surrogate skips these checks.

`fit_mesh_aero_surrogate` checks the degrees before it tabulates anything.
Fitting requires at least `(degree + 1)^2` directions and `poly_degree + 1`
distinct speed ratios; a larger total sample count cannot make up for too few
directions or speeds. The holdout check uses fresh directions at the middle
training speed ratio, so it does not establish accuracy between the sampled
speeds.

`read_mesh_aero_surrogate` checks every field of the file (schema, integer
degrees, six numeric coefficient rows of equal length, numeric normalization
and range, a three-component reference point) and raises an `ArgumentError`
naming the file for anything else. A file that declares a degree above the
limit is rejected before its coefficient arrays are read.

## Runtime effector

`AerodynamicCoefficientMeshSurrogate(Dict(link_index => surrogate);
wall_temperature_k=300.0)` carries one surrogate per link, keyed by the link's
integer position in `spacecraft.links`. The single-surrogate form puts a
whole-vehicle fit on the actual root link, wherever it sits in that list.
Links without a surrogate carry no aerodynamic load. A rigid vehicle therefore
uses one whole-vehicle mesh, and a vehicle whose panels articulate under
aerobraking guidance uses one mesh per link, composed through the link poses
the way the box model is, with no shadowing between links.

At each right-hand-side call, for each link, the effector forms the airspeed
direction in that link's frame, the speed ratio and
`wall_temperature_k / T`, evaluates the surrogate, and scales by the dynamic
pressure and the reference area (and the reference length for the moment).
Attitude follows the other aerodynamic effectors. With `orientation_sim`, the
propagated root attitude and the stored child attitudes place the links.
Without it, the root is held in the velocity-aligned frame (x along the
airspeed, z toward nadir) and each link's configured quaternion selects its
incidence. The effector uses the `wrench` interface described in
[Adding a Force or Torque of Your Own](custom_effector.md).

### Reference points and torque

Torque is returned about the root center of mass in the root body frame, and
only when the attitude is propagated. A link's fitted moment is about its
surrogate's `moment_reference_m`, a point `P` in the mesh (link) frame. The
effector moves it to the link origin `O` with

```text
M_O = M_P + (P - O) x F
```

then rotates it into the root frame through the link attitude and adds the
lever arm of the link's root-frame offset `body.r`. The torque about the root
therefore does not depend on where the fit put its reference point: two
surrogates of the same mesh fitted about different points give the same
motion. What the engine does assume, as for every other effector, is that the
link origin is the link's center of mass, so position the mesh accordingly.
`center=true` puts the bounding-box center at the origin, which only
approximates the center of mass; a model whose center of mass is elsewhere
needs to be translated before fitting.

### Diagnostics in radial flight

The drag, lift and cross save fields are projections of the total force. Drag
(along the airspeed) is defined whenever there is airspeed. Lift and cross are
built from the orbit normal `r x v`, which vanishes in radial flight (a
vertical descent or ascent, or within a relative tolerance of `1e-9` of one),
and they are then written as zero. The force applied to the spacecraft is
computed before that split and is the same either way: a vertical descent
decelerates under drag as a horizontal pass does.

### What is not modelled

Shadowing between links; per-link atmosphere sampling (one density and
temperature sample per spacecraft, as for the box model's default); thermal
incidence (the heating callback keeps its stored link angles); torque in the
fixed-attitude mode. For the supported threaded routes and their
restrictions, see [Parallel execution](#Parallel-execution).

## Validation

`test/unit/dynamics/mesh_aero_tests.jl` checks:

- a sphere against the closed-form free-molecular drag, and a box against the
  Hart formulas the box model uses;
- a plate behind another plate (shadowing);
- the harmonic basis, the fit and its JSON round trip;
- the moment-reference translation: exact in the panel method, the same across
  fits, and composed by hand for an offset, rotated child link;
- radial flight, and malformed constructor and file inputs;
- a run with propagated attitude whose rate build-up matches the plate torque,
  and a radial descent whose energy loss matches the recorded drag work;
- a fixed-attitude run of a cube against the same cube under
  `AerodynamicCoefficientfM`.

`test/unit/dynamics/mesh_aero_threading_tests.jl` checks the threaded routes.
`scripts/dev/aero/fit_mesh_aero_surrogate.jl` fits a model file from the
command line. See the
[Mesh Aerodynamics API](../generated/mesh_aerodynamics_api.md) for the
supported constructors and fitting functions.

## Parallel execution

The mesh surrogate supports the engine's automatic and explicitly selected
threaded routes. Evaluation reads the fitted coefficients and the spacecraft
geometry, and its temporary basis and force buffers belong to each call. Treat
the surrogate matrices, the link map and the geometry as read-only while a run
is active.

A tuple that contains more than one mesh or box aerodynamic effector uses the
serial route, because those effectors share the drag, lift and cross
diagnostic slots. Use one mesh effector with its per-link surrogate map. Other
model and thread-budget restrictions still apply. Being eligible for threading
does not by itself make a run faster; the execution policy chooses the route
for the workload.
