# Geometry helper ownership

Which internal module owns the shared geodetic and orbital-rotation expressions,
what each caller keeps for itself, and which tests check them. None of these
helpers is exported from `SpaceAGORA`. Sharing an expression does not change a
reference surface, a saved-output meaning or any scientific setting.

## Geodetic latitude and altitude

`SimulationModel.Geodesy` (`src/foundation/core/numerics/geodesy.jl`) owns the
scalar Bowring expressions for geodetic latitude (`_bowring_latitude`) and
ellipsoidal height (`_bowring_altitude`). Two callers use them:

- `FrameTransforms.rtolatlong(r_p, planet, spherical_harmonic_topography=false)`
  takes a planet-fixed `SVector{3, Float64}` in metres and returns
  `[altitude_m, geodetic_latitude_rad, east_longitude_rad]`. It computes the
  flattening in Float64 (`_planet_flattening`) even when the planet's radii are
  `Float32` or `BigFloat`. A position that is not a Float64 `SVector` is a
  `MethodError`.
- `Geodesy.geodetic_altitude(radius, u_pp, planet)` builds the Cartesian position
  from a Float64 radius and a planet-fixed direction, then applies the same
  expressions.

Each caller keeps its own input construction, prime-vertical radius and return
conversion. Neither normalizes the caller's position in order to reuse the other.

`Geodesy.radius_for_geodetic_altitude(target_altitude_m, u_pp, planet)` inverts
`geodetic_altitude` along a direction by bisection. The planet-aware
`InitialCondition` constructor in `SpacecraftModels` and the thruster guidance
functions (`src/models/gnc/guidance/maneuvers/thruster_guidance_functions.jl`)
share it. Directions must be unit vectors, and the callers validate them. A
negative target returns `NaN`; the `InitialCondition` constructor instead throws
`ArgumentError("Oblate InitialCondition altitudes must be nonnegative; …")`. Its
altitudes are measured from the reference ellipsoid along the constructor's
transformed apsis directions.

With `spherical_harmonic_topography=true`, `rtolatlong` still computes geodetic
latitude and longitude first. It then returns `norm(r_p)` minus the surface radius
from `planet.topography_function(args, Clm_topo, Slm_topo, lat, lon, A_topo)`,
where `args` is the including module's `args` global if it has one and `nothing`
otherwise. Errors from the callback propagate. The method that takes an
ephemerides model, `rtolatlong(r_p, planet, ephemerides_model)`, uses the
ellipsoid path. `rtolatlongrad` is a separate calculation that returns radius,
geocentric latitude and longitude. No helper chooses between sphere, ellipsoid
and topography on the caller's behalf.

Production code loads one `Geodesy` and one `OrbitalGeometry` module, in
`src/orchestration/simulation/composition/simulation_model.jl`, and
`FrameTransforms` imports both. The frame file
`src/foundation/core/interfaces/reference_system.jl` expects its enclosing module
to supply them and loads nothing itself. Standalone test sandboxes include
`test/helpers/reference_system_sandbox.jl` instead: it includes the same
`geodesy.jl` and `orbital_geometry.jl` and then the frame file into the
including module, without loading SpaceAGORA or `RuntimeServices`. That keeps
the includer's own `args` visible and leaves the missing-lock probe in
`test/probes/reference_system_probes.jl` reachable.

`test/unit/geometry_altitude_contract_tests.jl` compares `rtolatlong` and
`geodetic_altitude` with CSPICE `recgeo` (which needs no kernels) for Earth, Mars
and Venus, at latitudes from −90° to 90° and heights of 0, 150 and 500 km. It
requires altitude within 1 mm (`atol=1e-3`), latitude within 1e-8 rad and
longitude within 1e-14 rad. It also checks spherical cases analytically (heights
of −500 m, 0 and 1e9 m, within 2e-7 m), that the sandbox and package copies return
identical results, the validation of negative targets, and the topography branch
with real callbacks. One-pass Bowring latitude is approximate at finite height:
the 1e-8 rad bound is established on that grid only, and nothing here claims
accuracy at other heights. Nor does it claim that the spherical osculating
periapsis and the instantaneous geodetic altitude coincide.

## Orbital-element rotation

`SimulationModel.OrbitalGeometry` (`src/foundation/core/numerics/orbital_geometry.jl`)
owns `_perifocal_from_inertial(i, ω, Ω)`, a Float64 `SMatrix` built from
inclination, argument of periapsis and ascending node in radians. It maps inertial
components to perifocal components, and its transpose maps perifocal vectors into
inertial axes. Keep the literal `@SMatrix` construction and its expression order:
building the matrix from an array instead changes the last bits of some
normalized directions even when the entries agree.

Two callers use it:

- `FrameTransforms.orbitalelemtorv(oe, planet)` keeps its own angular-momentum,
  position and velocity scaling, multiplies by the transpose, and returns position
  and velocity as mutable `Vector{Float64}`.
- `SpacecraftModels._initial_condition_apsis_direction_ii(i, ω, Ω, ν)` builds its
  own anomaly direction, multiplies by the transpose and normalizes the result. It
  does not use a unit orbit radius or take an orbit's speed. Degree-to-radian
  conversion, frame precedence and altitude validation stay with the
  `InitialCondition` constructor.

The element-to-inertial matrix in `scripts/plotting/plot_data.jl` and the
element handling in individual examples are separate code. Nobody has shown that
their expression order or data contracts match this one, so they do not use it.

`test/unit/apsis_rotation_contract_tests.jl` checks the transpose direction with
elementary rotations, orthogonality and the determinant, both callers against
independently composed active rotations, and the two-body invariants (angular
momentum and vis-viva energy) without propagating.

## Legacy geodetic velocity: `latlongtoOE`

`FrameTransforms.latlongtoOE(LATLONGH, planet, γ, α, v)` converts a geodetic
location and a speed to orbital elements:

- `LATLONGH` is geodetic latitude and east longitude in radians, then height
  above the reference ellipsoid in metres.
- `v` is the inertial speed in m/s. `γ` is measured above the local geodetic
  tangent plane, and the azimuth `α` from north toward east, both in radians.
- `planet.L_PI` must hold the inertial-to-planet-fixed rotation for the caller's
  epoch. The function has no time argument and does not update it.
- It returns `[a, e, i, Ω, ω, ν]`, in metres, dimensionless and radians, and
  prints its inputs and intermediate values.

The local angles say how an inertial velocity vector is resolved at that
location; they do not make the speed ground-relative or air-relative. At nonzero
flattening, geodetic up is not the radial direction. The position goes through
`r_pintor_i` with zero velocity, and the velocity components are only rotated
with `planet.L_PI'`; no `ω × r` term is added.

For a planet-relative Cartesian velocity, use `r_pintor_i(position, velocity, planet)`,
which includes the transport term. `CartesianInitialCondition` takes an
explicitly prepared inertial state and does not call `latlongtoOE`. Wind is a
third quantity and needs its own stated convention.

`test/unit/legacy_velocity_frame_contract_tests.jl` covers equatorial east and
west cases with and without planet spin (speed, vis-viva and the difference from
the transport-term reading) and tilted `L_PI` at non-equatorial latitudes against
CSPICE `georec`. `test/probes/reference_system_probes.jl` checks the vis-viva
semi-major axis and the apsis eccentricity for one case.

No code in this repository calls `latlongtoOE` apart from its tests; it is still
exported from `FrameTransforms`. That says nothing about scripts outside the
repository.
