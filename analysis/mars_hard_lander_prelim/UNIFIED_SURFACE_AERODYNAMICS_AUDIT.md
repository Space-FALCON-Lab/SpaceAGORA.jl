# Unified Surface Aerodynamics Audit

## Implemented architecture

One regime-neutral surface representation now owns the geometric inputs used
by both endpoint laws:

```text
AerodynamicGeometry
  surface.positions_body_m
  surface.inward_normals_body
  surface.area_weights_m2
  reference_area_m2
  reference_length_m
```

The sphere-cone constructor supplies a spherical cap, conical frustum, and rear
base disk. `aerodynamic_plate_surface` supplies an exact rectangular plate
primitive, and `combine_aerodynamic_surfaces` assembles multiple primitives
under one force/moment normalization. A thin plate exposed on both sides must
be represented by two coincident surfaces with opposite inward normals.

The old `Newtonian*` geometry names remain aliases, and
`sphere_cone_newtonian_geometry` still returns the forebody without the base.
This preserves the Grant-Braun results. The new
`sphere_cone_aerodynamic_geometry` returns a closed surface by default.

## Free-molecular endpoint

The implemented local law is the Schaaf-Chambre pressure and tangential-shear
model reproduced on page 133 of NASA CR-182076, *Aerodynamic
Preliminary Analysis System II, Part I - Theory*:

<https://ntrs.nasa.gov/api/citations/19910013767/downloads/19910013767.pdf>

For each surface element,

```text
S   = sqrt(gamma/2) M
sn  = S dot(Vinf, n_inward)
Cp  = complete pressure expression on NASA CR-182076 page 133
Ct  = complete tangential expression on NASA CR-182076 page 133
dCF = (Cp n_inward + Ct_vector) dA / Aref
dCM = cross(r, dCF) / lref
```

All signed incidences are evaluated; there is no Newtonian leeward clipping in
the free-molecular endpoint. The implementation uses separate normal and
tangential momentum accommodation coefficients, wall temperature, freestream
temperature, and molecular speed ratio. The accommodation convention is the
standard Schaaf-Chambre convention documented by Hart et al.: zero is specular
and one is fully accommodated/diffuse.

Hart et al. source:
<https://repository.gatech.edu/bitstreams/aceec4be-2b83-4a6a-9841-8fa4dd36afc9/download>

NASA CR-182076's legend reverses the normal-coefficient endpoint labels even
though its printed equation has the standard Schaaf-Chambre algebra. The code
follows the equation and the Hart definition, not that conflicting legend.
This choice is locked by specular and diffuse limiting-case tests.

NASA TM-101600 independently supports the surface-element architecture: it
describes a vehicle as many small flat plates, computes pressure and shear
force and all moments for each element, and integrates them over the body:

<https://ntrs.nasa.gov/api/citations/19890015753/downloads/19890015753.pdf>

## Knudsen number and transitional bridge

The automatic model uses

```text
lambda = mu/rho sqrt(pi/(2 R T))
Kn     = lambda/Lref
```

This is the viscosity-based mean-free-path definition obtained by rewriting
NASA TM-2006-214397 Eq. (21) with `p=rho R T`:

<https://ntrs.nasa.gov/api/citations/20060049127/downloads/20060049127.pdf>

The transition model blends the complete force and moment vectors before
deriving body- and wind-axis coefficients:

```text
C = (1-Pb) Ccontinuum + Pb Cfree-molecular
Pb = sin(phi)^2
phi = (pi/2) (log10(Kn)-log10(Kn_cont))
             / (log10(Kn_fm)-log10(Kn_cont))
```

`Pb` is clamped exactly to zero below `Kn_cont` and one above `Kn_fm`. The
defaults are `Kn_cont=0.001` and `Kn_fm=10`, giving the classic logarithmic
sine-squared bridge. The governing form comes from Wilmoth, Blanchard, and
Moss, *Rarefied Transitional Bridging of Blunt Body Aerodynamics*:

<https://ntrs.nasa.gov/api/citations/20040095951/downloads/20040095951.pdf>

The bounds are configurable because the paper shows that bridge parameters
depend on vehicle, attitude, and coefficient. It recommends DSMC information
inside the transition region when accuracy beyond an engineering bridge is
required.

## Operational behavior

`AerodynamicSurfaceModel` is the regime-neutral sampled effector. The historic
`AerodynamicCoefficientNoBallisticFlight` name is an alias. Supported selectors
are:

- `:continuum`: regular or modified Newtonian only;
- `:free_molecular`: Schaaf-Chambre only;
- `:transitional`: force/moment bridge at a supplied or computed Knudsen number;
- `:automatic`: continuum, bridge, or free molecular based on Knudsen number.

The sampled path uses wind-relative speed consistently for Mach, speed ratio,
dynamic pressure, and direction. It computes Knudsen number from the sampled
density and temperature plus dynamic viscosity and a characteristic length.
Both may be overridden. Any nonzero free-molecular contribution requires an
explicit positive wall temperature. Positive representable densities are not
discarded through an arbitrary machine-epsilon cutoff; an infinite Knudsen
number cleanly selects the free-molecular endpoint.

The old `AerodynamicCoefficientfM` rectangular-prism model remains available
only for compatibility. It is not used by `AerodynamicSurfaceModel`.

## Verification

Tests cover:

- unchanged Grant-Braun regular and modified Newtonian reference values;
- closed-sphere-cone base area and unchanged forward continuum coefficients;
- direct NASA free-molecular equation values;
- high-speed specular and diffuse limits;
- the closed-form Schaaf-Chambre diffuse-sphere drag coefficient;
- exact plate force and moment about an offset reference;
- finite three-axis sphere-cone forces and moments;
- bridge endpoints, midpoint, monotonicity, and vector blending;
- mean free path and Knudsen number;
- the operational automatic-Knudsen wrench path; and
- invalid regimes, accommodation values, temperatures, and bridge bounds.

## Validity limits

This is a complete surface integration of the selected endpoint laws, not a
general high-fidelity flow solver.

- The continuum endpoint is still pressure-only Newtonian theory. It does not
  add continuum skin friction, base pressure, separated wake flow, real-gas
  chemistry, ablation, or shock-layer interaction.
- The free-molecular endpoint assumes collisionless freestream molecules and
  spatially uniform wall temperature and accommodation coefficients for one
  evaluation. These quantities can be varied between calls, not per element.
- The surface container does not perform ray-traced inter-component
  visibility. It is exact for externally exposed convex surfaces. For concave
  or overlapping component assemblies, callers must supply only the exposed
  gas-facing surfaces or add a separate visibility/occlusion model.
- The sine-squared transition is an endpoint interpolation, not a Boltzmann,
  DSMC, or Navier-Stokes solution. It can miss coefficient-specific behavior,
  wake rarefaction, and static-stability changes in the transition regime.
- Planet viscosity values in the current environment model are constants.
  Automatic Knudsen selection therefore uses an engineering estimate unless a
  temperature-dependent viscosity is supplied through the model override.
- The characteristic length defaults to the geometry reference length. A
  different rarefaction length scale must be configured when appropriate.
