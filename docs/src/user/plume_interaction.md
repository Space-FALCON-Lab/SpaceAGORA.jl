# Plume interaction

[`PlumeSurfaceInteractionModel`](@ref) adds a modeled ground-effect force to
an engine whose held thrust is available as `control.actuators.thrust_n`, a
finite, nonnegative vector in simulation spacecraft order. The effector does
not apply the engine's primary thrust; the control model must already supply
it. A run may have at most one plume effector, which can cover several
spacecraft and writes the seven `sc{i}_plume_*` diagnostic columns (see
[the plume columns](output_columns.md#Plume-interaction)).

Use this model when you already have a landing configuration and want to look
at its plume diagnostics. The [lunar landing guide](lunar_landing.md) gives a
complete scenario; this page explains how to configure the plume model and
read what it reports.

The closure represents a Gaussian pressure footprint, wall shear, erosion by
excess shear, a bounded ejecta speed and a small axial augmentation of the
thrust. The defaults are illustrative phenomenological parameters, not an
Apollo flight calibration, and the model does not solve gas dynamics, regolith
transport or a changing surface. [`PlumeSurfaceConfig`](@ref) holds the
parameters. Its `expansion_ratio`, `chamber_pressure_pa`, `exit_mach`,
`bulk_density_kg_m3` and `cohesion_pa` are recorded metadata that do not enter
the closure.

## Set the terrain and engine geometry

The terrain uses planetocentric latitude and longitude in degrees and a
reference sphere in meters. Pass `reference_radius_m` when the mission uses an
explicit datum, and use the same radius for guidance and touchdown. An explicit
radius must match a DEM terrain's declared radius. Without it, DEM terrain uses
its own radius and other terrain uses the planet's equatorial radius. For
example, `PlumeSurfaceInteractionModel(control, terrain;
reference_radius_m=landing_radius)` keeps a landing and its plume on the same
reference sphere.

At each state, the local radial tangent plane at the terrain sample below the
spacecraft gives the radial clearance, and dividing that by the cosine between
the exhaust axis and the downward vertical gives the slant height. The plume
quantities use the height of the nozzle exit plane, the slant height less
`nozzle_offset_m`. This approximation does not ray-trace the terrain, move the
impact point of a tilted plume, or account for local slopes and planetary
curvature. Horizontal and upward exhaust do not interact with the surface.

With attitude, the exhaust points along body +z and the augmentation acts
along body -z; without attitude, the exhaust points radially inward. The
scalar `nozzle_offset_m` lies along this exhaust line from the integrated
vehicle reference point, so the force line passes through that point and the
modeled torque is zero. An off-axis engine needs an additional geometry and
torque model.

## Read the diagnostics across a run

Evaluating the force does not change the diagnostics. An accepted-step callback
samples the trajectory and integrates the erosion rate by trapezoidal
quadrature, using the thrust held across each interval; a throttle change at
the end of an interval applies to the next one. A new run resets the integral,
and checkpoint segments inside one run keep it. A separately resumed run starts
a new diagnostic integral at its resume time. Its results table holds the
whole run (see [Outputs](outputs.md)), so the cumulative mass column drops to
zero at the resume time and counts only the resumed interval from there.

The cumulative mass is a diagnostic, not an ODE state controlled by the solver
tolerances. Show convergence by reducing `dt_max_orbit` and, for an
atmospheric phase, `dt_max_atmosphere`. The default save fields recompute the
instantaneous quantities from the saved state and the thrust held over the
interval, and the saved mass is interpolated along the trapezoid of the
accepted interval. A save must therefore fall within the current accepted
interval; these fields do not reconstruct earlier history from the effector's
current state.

## Choose compatible stopping and solver settings

The model checks its configuration when it is built: values must be finite,
lengths and material scales positive, coefficients nonnegative,
`expansion_ratio` at least 1, `plume_half_angle_deg` strictly between 0 and
90, the ejecta speed bounds in order and `ground_effect_max_fraction` within
[0, 1]. The thrust vector and the plume state must match the run's spacecraft
count. A full propagation still needs its own stopping or touchdown policy
suited to the mission; this diagnostic effector does not end a run at impact.

The `gravity_backbone_split` solver mode refuses this nonconservative,
planet-frame-dependent effector.
