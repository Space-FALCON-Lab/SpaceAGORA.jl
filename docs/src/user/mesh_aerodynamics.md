# Mesh Aerodynamics

The aerodynamic models in `SimulationModel` describe a spacecraft as boxes
(`AerodynamicCoefficientfM`, the Hart closed forms for a rectangular prism
summed over links). When a CAD model of the vehicle exists, the mesh
aerodynamics path derives free-molecular force and moment coefficients from
that geometry instead, then fits an analytic surrogate the integrator
evaluates at run time without touching the mesh.

```julia
using SpaceAGORA

panels = mesh_aero_panels("data/models/magellan_nasa_3d_resources.glb";
    scale=1.0, rotation_deg=(0, 0, 90))              # the viewer's model transform: long axis along body x
surrogate = fit_mesh_aero_surrogate(panels; degree=8, speed_ratios=(3, 5, 8, 12, 20), verbose=true)
write_mesh_aero_surrogate("output/magellan_aero.json", surrogate)

effector = AerodynamicCoefficientMeshSurrogate(surrogate; wall_temperature_k=300.0)
# ... dynamic_effectors=(GravitationalHarmonicsModel(...), effector) in place of AerodynamicCoefficientfM()
```

## Panel method

`mesh_aero_panels` reads STL, OBJ or glTF/GLB (embedded buffers, no Draco)
through the same readers the viewer uses, applies the same scale, XYZ Euler
rotation and centring, and keeps per-facet centroids, outward normals (from
the vertex winding, or `outward_normals=true` to flip normals toward the
outside for convex meshes with mixed winding) and areas. Coefficients are
normalised by the dynamic pressure of the wind-relative airspeed,
`reference_area_m2` (default: total surface area over four, the mean projected
area of a convex body); `articulations` pose parts of the model first (the
Magellan wings turned broadside, see `articulate_triangles`) and `reference_length_m` for moments (default: the
bounding-box diagonal), about `moment_reference_m` (default: the origin,
which is the link centre of mass once the model is centred).

`panel_aero_coefficients(panels, vhat, s; sigma_n, sigma_t, tw_ratio)` sums
the Schaaf and Chambré pressure and shear over the facets for one airspeed
direction `vhat` in the mesh frame (the direction the body moves relative to
the gas) and speed ratio `s = V / sqrt(2 R T)`, with normal and tangential
momentum accommodation `sigma_n`, `sigma_t` (1 is fully diffuse) and the
wall-to-freestream temperature ratio `tw_ratio`. Facets behind other facets
along the flow are found with a depth buffer (`panel_shadow_mask`) and
contribute nothing, so a dish in front of a bus shields it and a panel
behind another is dead. The Hart box formulas are the same integrand with the
back faces dropped, and a box mesh reproduces them to round-off.

## Surrogate

The wall temperature enters the free-molecular flux only through the diffuse
re-emission term, linearly in `sqrt(Tw/T)`, so each of the six coefficients
splits exactly as `A(vhat, s) + sqrt(Tw/T) B(vhat, s)`.
`fit_mesh_aero_surrogate` tabulates `A` and `B` over `n_directions`
Fibonacci-spread directions and the listed speed ratios and fits real
spherical harmonics of degree `degree` in the direction times a polynomial of
degree `poly_degree` in `1/s` by least squares. The result is a
`MeshAeroSurrogate`: two `6 x n` coefficient matrices, the normalisation, the
accommodation coefficients, the fitted speed-ratio range (evaluation clamps
to it) and metadata with the residuals, both over the training samples and
over fresh holdout directions. Read the holdout error against the coefficient
scale before trusting a degree: a smooth body fits at degree 6, a body with
thin panels and edges wants 10 or more. The fit takes one shadow buffer per
direction and speed ratio; a few thousand facets take seconds, a
half-million-facet model minutes.

`mesh_aero_coefficients(surrogate, vhat, s; tw_ratio)` evaluates it in a few
hundred multiply-adds. `write_mesh_aero_surrogate` and
`read_mesh_aero_surrogate` store it as JSON.

## Runtime effector

`AerodynamicCoefficientMeshSurrogate(Dict(link_index => surrogate);
wall_temperature_k)` carries one surrogate per link, keyed by the link's
position in `spacecraft.links` (1 is the root); the single-surrogate form
puts a whole-vehicle fit on the root. Links without a surrogate carry no
aerodynamic load, so a rigid vehicle uses one whole-vehicle mesh, and a
vehicle whose panels articulate under aerobraking guidance uses one mesh
per link, composed through the link poses the way the box model is (with no
shadowing between links, as before).

Per link and per RHS call the effector forms the airspeed direction in that
link's frame, the speed ratio and `wall_temperature_k / T`, evaluates the
surrogate, and scales by the dynamic pressure and the reference area (and
length for the moment). Attitude follows the other aerodynamic effectors:
with `orientation_sim` the propagated root attitude and the stored child
attitudes place the links; without it the root is held in the
velocity-aligned frame (x along the airspeed, z toward nadir, the attitude
the viewer draws) and each link's configured quaternion selects its
incidence. Torque about the root centre of mass is returned only when the
attitude is propagated. The drag, lift and cross save fields are filled from
the total force as for the box model.

## Validation

`test/unit/dynamics/mesh_aero_tests.jl` checks a sphere against the
closed-form free-molecular drag, a box against the Hart formulas the box
model uses, a plate behind another plate, the harmonic basis, the fit and
its JSON round trip, and a short simulation of a cube with the surrogate
against the same cube under `AerodynamicCoefficientfM`.
`scripts/dev/aero/fit_mesh_aero_surrogate.jl` fits a model file from the
command line, and the Magellan and Cassini viewer demos switch to the mesh
surrogate with `SPACEAGORA_DEMO_MESH_AERO=1`, where the SPICE ghost gives a
direct comparison of the two aerodynamic models against the flown trajectory.
