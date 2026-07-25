# Grant-Braun Sphere-Cone Cn/Ca Audit

## Scope and source

This audit covers the full pressure-only Newtonian aerodynamic implementation
used by the Mars hard-lander preliminary analysis and package effector. The
primary source is Michael J. Grant and Robert D. Braun, *Analytic Hypersonic
Aerodynamics for Conceptual Design of Entry Vehicles*, AIAA 2010-1212,
doi:10.2514/6.2010-1212. The authors use regular Newtonian theory for the
published analytic results and state that modified Newtonian results are
obtained with the appropriate pressure multiplier.

Primary-source PDF:
<https://engineering.purdue.edu/~mjgrant/48th-aiaa-aerospace-science.pdf>

Grant and Braun specify that modified Newtonian results use the appropriate
stagnation-pressure multiplier, but do not print its equation. The implemented
perfect-gas normal-shock expression was therefore checked against Eqs. (6) and
(7) of NASA TN D-2231, *Force-Coefficient and Moment-Coefficient Correlations
and Air-Helium Simulation for Spherically Blunted Cones*:
<https://ntrs.nasa.gov/api/citations/19650000987/downloads/19650000987.pdf>

## Source completeness

The paper publishes the governing surface integrals, coordinate conventions,
shadow rule, reference normalization, and basic-shape parametrizations. It does
not publish the generated Matlab closed-form routines used for its validation
plots. Therefore, an exact transcription of those omitted symbolic routines is
not possible from the paper alone.

The implementation evaluates the published continuous-surface integrals
directly with deterministic Gaussian quadrature over the exact parametrized
surfaces. It is not a flat-panel approximation: surface locations, inward
normals, and differential-area Jacobians come from the continuous sphere and
cone equations.

## Audited equations and conventions

The body-axis freestream direction is

```text
Vinf = [-cos(alpha)cos(beta), -sin(beta), -sin(alpha)cos(beta)]
```

For every point on a continuously parametrized surface,

```text
incidence = dot(Vinf, inward_unit_normal)
Cp = pressure_scale * incidence^2    when incidence > 0
Cp = 0                               when incidence <= 0
```

`pressure_scale` is 2 for regular Newtonian theory and `Cp_max` for modified
Newtonian theory. Force and moment coefficients follow the paper's surface
integrals:

```text
Cforce  = integral(Cp * inward_unit_normal * dA) / Aref
Cmoment = integral(Cp * cross(r, inward_unit_normal) * dA)
          / (Aref * lref)
```

The sphere-cone is the superposition of the paper's spherical segment and
conical frustum, normalized by the common base area `pi*rb^2`. Coordinates use
the paper's aircraft body axes: `+x` forward, `+y` right, and `+z` down.

## Findings

1. The preliminary compact expression is not the complete Grant-Braun
   sphere-cone result. It omits part of the spherical cap's angle-dependent
   axial and lateral/normal force. It agrees at zero incidence but diverges as
   total incidence increases.
2. The preliminary configuration exposed `:modified_newtonian`, but the old
   evaluator ignored the selector and always returned the regular Newtonian
   result. The nominal modified-Newtonian configuration was therefore not
   actually modified Newtonian.
3. The legacy package wrapper was incorrect. It evaluated `cos(delta^2)` in
   `CN` instead of `cos(delta)^2` and subtracted an undocumented `0.15` from
   `CD`. Its unused `Cp_max = 2` assignment did not implement modified
   Newtonian theory.
4. Neither path implemented sideslip, surface shadowing, side force, moments,
   stability derivatives, or operational six-degree-of-freedom force mapping.

The port now has one package-owned continuous-surface model with explicit
regular and modified pressure laws, arbitrary angle of attack and sideslip,
pointwise shadow clipping, all three force and moment coefficients, wind-axis
coefficients, stability derivatives, continuous basic-surface superposition, and
input/geometry validation. The sampled aerodynamic wrench applies the result
using stage-consistent atmosphere, wind-relative Mach, dynamic pressure, and
attitude. The preliminary analysis and legacy coefficient wrapper delegate to
this implementation.

For the paper's Table 1 sphere-cone (`rn=0.3 ft`, `rb=1.25 ft`, `delta=70 deg`),
regular Newtonian evaluation gives:

```text
alpha=10 deg, beta=0 deg:  CA=1.7170715356, CN=0.0398739686
alpha=0 deg, beta=20 deg:  CA=1.5737901746, CS=-0.0749385482, CN=0
```

The second result follows the sphere-cone curves in the paper's Figure 16.

As a numerical convergence check, the default 12-by-48 sphere/cone quadrature
was compared with a 30-by-240 evaluation of the Table 1 geometry over
`alpha=-90:10:90 deg` and `beta=-40:10:40 deg`. The largest absolute differences
were `5.3e-6` among the force coefficients and `8.4e-6` among the moment
coefficients. Callers can increase both quadrature orders when tighter accuracy
is required near a moving shadow boundary.

## Validity limits

This Newtonian endpoint is appropriate as a conceptual-design pressure model
when all of the following assumptions are acceptable. The separate
free-molecular implementation and Knudsen bridge are audited in
`UNIFIED_SURFACE_AERODYNAMICS_AUDIT.md`.

- hypersonic continuum flow;
- a tangent sphere-cone geometry;
- base-area normalization and the Grant-Braun axis convention; and
- a calorically perfect gas with a supplied constant `gamma` for the modified
  pressure correction.

High-angle shadowing and nonzero sideslip are now evaluated, subject to the
paper's convex-surface assumption. The Newtonian endpoint should still not be treated as
correct for mutual shadowing between interacting components, concave geometry,
separated or wake flow, base drag, skin friction, rarefied/transitional flow,
real-gas chemistry, ablation, or shock-layer interaction. Modified Newtonian
theory improves the stagnation pressure level; it does not add those missing
flow physics.

For numerical continuity, `modified_newtonian_cp_max` supplies an isentropic
continuation at and below Mach 1. That continuation prevents a trajectory from
becoming singular, but Newtonian surface coefficients are not physically
validated in that regime and should be replaced or blended with an appropriate
subsonic/supersonic aerodynamic model.
