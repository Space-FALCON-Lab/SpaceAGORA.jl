# Geometry helper ownership and compatibility

These boundaries describe the maintained calculations. They do not change the
reference surfaces, scientific settings, saved-output meanings or native policy.

## Bowring altitude

`SimulationModel.Geodesy` owns the common scalar Bowring latitude and final
ellipsoidal-height expressions. `geodetic_altitude` first forms Cartesian
coordinates from its Float64 radius and direction. `FrameTransforms.rtolatlong`
accepts a Float64 Cartesian `SVector`, retains its Float64 flattening conversion,
and returns `[altitude_m, geodetic_latitude_rad, east_longitude_rad]`. The wrappers
retain their original arithmetic, prime-vertical radius calculation and return
conversions; neither normalizes a caller's position merely to reuse another API.

The direction-based altitude/radius inversion remains shared by initialization
and guidance. Directions are expected to be unit vectors; callers own validation.
Negative inversion targets still return `NaN` internally, while the planet-aware
`InitialCondition` constructor still raises `ArgumentError`. The constructor's
altitudes refer to the reference ellipsoid along its transformed apsis directions.

The optional topography branch remains in `rtolatlong`: it subtracts the supplied
surface radius from `norm(position)` and passes the includer module's `args` when
present (otherwise `nothing`) to the existing six-argument callback. It still
computes geodetic latitude/longitude before that callback. Callback failures are
not intercepted. The overload accepting an ephemerides model still delegates to
the default ellipsoidal path; `rtolatlongrad` remains a separate geocentric-radius
and latitude calculation. No sphere/ellipsoid/topography policy is substituted.

Production imports the one Geodesy module. The internal `reference_system.jl`
leaf expects its enclosing aggregator to supply `Geodesy` and `OrbitalGeometry`;
it does not orchestrate source loading. Standalone test sandboxes raw-include
`test/helpers/reference_system_sandbox.jl`, which loads the same canonical numeric
sources and then the frame leaf into the includer. It imports neither SpaceAGORA
nor RuntimeServices, preserving local `args` and the early missing-lock probe.
The frame file contains no duplicate formulas or test-specific loading policy.

Focused tests use kernel-free CSPICE `recgeo`, analytic spherical cases and actual
topography callbacks. Bowring latitude is approximate at finite height: the
specified Earth/Mars/Venus grid from 0 to 500 km uses a 1e-8 rad latitude bound
and retains the 1 mm height criterion. The pinned original and current results
match exactly on that grid; no improved or arbitrary-height accuracy is claimed.
Existing constructor/apsis/frame probes remain applicable.
This is bounded geometry compatibility, not a new scientific campaign or a claim
that spherical osculating periapsis and instantaneous geodetic altitude coincide.

## Apsis and orbital-element rotations

`SimulationModel.OrbitalGeometry` owns the internal Float64 matrix
`_perifocal_from_inertial(i, ω, Ω)`, with inclination, argument of periapsis and
ascending node in radians. It maps inertial components to perifocal components;
its transpose maps perifocal vectors into inertial axes. The original nine
component expressions and their evaluation order are preserved. The owner retains
the constructor's literal `@SMatrix` construction and inline annotation: replacing
that construction with an array conversion changes some final normalized direction
bits even when the matrix entries match. Bounded exact-base differential checks
therefore cover both actual callers, including signed zero and nonfinite inputs.

`orbitalelemtorv` retains its original angular-momentum, position and velocity
scaling, multiplication by the transpose, and returned mutable vectors.
`_initial_condition_apsis_direction_ii` still constructs its own anomaly direction,
multiplies by that transpose and normalizes the result. It does not substitute a
unit orbital radius or inherit an orbit's speed. Initial-condition degree/radian
conversion, frame precedence and altitude validation remain caller-owned.

The matrix owner is loaded before the spacecraft model; FrameTransforms imports
the same module. The standalone test helper loads that owner's exact source in
its sandbox before the frame leaf. No new root export or public API is added.
The legacy plotting matrix and example-specific element interpretations remain
separate: their expression order and data contracts have not been declared
interchangeable. Tests discriminate transpose direction with elementary rotations
and check conserved two-body state invariants without propagation.

## Legacy velocity-frame compatibility

`FrameTransforms.latlongtoOE(LATLONGH, planet, γ, α, v)` keeps its established
inertial-speed behavior. `LATLONGH` supplies geodetic latitude, east longitude
(radians), and height above the reference ellipsoid (metres). Speed `v` is in m/s;
`γ` is measured above the local geodetic tangent plane, and azimuth `α` is measured
from north toward east. Those local directions describe how an **inertial**
velocity vector is resolved at the current location. They do not make the speed
ground-relative or air-relative. At nonzero flattening the geodetic up direction
need not equal the radial direction.

The caller must supply `planet.L_PI` for its chosen epoch. The helper has no time
argument and does not update that orientation. It transforms the constructed
position and rotates the velocity components with `L_PI'`, without adding
`ω × position`. It returns `[a, e, i, Ω, ω, ν]` in metres/dimensionless/radians.
The existing positional signature, printed diagnostics, validation behavior and
calculation are unchanged. The old test's vis-viva and speed-preserving roundtrip
already distinguish this contract; new fixtures additionally cover rotating and
zero-spin equatorial east/west cases and non-equatorial tilted frames.

A planet-relative Cartesian velocity instead requires the complete existing
`r_pintor_i(position, velocity, planet)` transform, including its transport term.
The supported Cartesian initial-condition path accepts an explicitly prepared
inertial state; it does not call this legacy geodetic helper. Wind is a third
quantity and needs its own explicitly established convention. No new relative-
speed overload is introduced, and the separate descent-kit initialization is
not reinterpreted.

The tracked current-base search found this helper's definition, export and tests,
not a production caller. That bounded inventory does not establish that external
or historical consumers never used it. No existing scientific product is
reclassified, regenerated or declared unaffected by this contract clarification.
