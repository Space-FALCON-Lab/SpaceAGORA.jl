# Plume interaction

[`PlumeSurfaceInteractionModel`](@ref) adds a modeled ground-effect force to an
engine whose held thrust is available as `control.actuators.thrust_n`. Supply
that vector in simulation spacecraft order. The primary engine thrust must
already be supplied by the control model. One plume effector publishes the seven
`sc{i}_plume_*` diagnostic columns; it may cover multiple spacecraft.

The closure represents a Gaussian pressure footprint, wall shear, excess-shear
erosion, bounded ejecta speed and a small axial augmentation. Defaults are
illustrative phenomenological parameters. They have not been established as an
Apollo flight calibration, and the model does not solve gas dynamics, regolith
transport or a changing surface. Engine and soil metadata fields that do not
enter the closure are identified in [`PlumeSurfaceConfig`](@ref).

The supplied terrain uses planetocentric latitude/longitude in degrees and a
reference sphere in meters. For DEM terrain its declared radius is used. For
other terrain the planet's equatorial radius supplies that reference. At each
state the local radial tangent plane at the subspacecraft terrain sample gives
the radial clearance; dividing by the downward exhaust cosine gives slant
height. This approximation does not ray trace the terrain, relocate the tilted
impact point or account for local slopes and planetary curvature. Horizontal
and upward exhaust do not interact with the surface.

The attitude quaternion places exhaust along body +z and augmentation along
body -z. When attitude is absent, exhaust points radially inward. The scalar
`nozzle_offset_m` lies along this exhaust line from the integrated vehicle
reference point. The force line therefore crosses that point and the modeled
torque is zero. An off-axis engine needs an additional geometry/torque model.

Pure wrench evaluations do not alter diagnostics. An accepted-step callback
samples the trajectory and integrates erosion rate by trapezoidal quadrature,
using the thrust held across that interval. Endpoint throttle changes apply to
the following interval. A new run resets the integral; checkpoint segments
inside one run retain it. A separately resumed run starts a new diagnostic
integral at its resume time.

The cumulative mass is a diagnostic, not an extra ODE state controlled by solver
tolerances. Demonstrate convergence by reducing `dt_max_orbit` and, for an
atmospheric phase, `dt_max_atmosphere`. Default save fields recompute instantaneous
quantities from the supplied saved state and held interval thrust. Saved mass
uses the accepted-interval trapezoidal interpolation. Save queries must occur in
the current accepted interval; these fields do not reconstruct arbitrary past
history from the current mutable effector state.

The model validates finite configuration values, positive length/material scales,
nonnegative coefficients and matching spacecraft counts. A full propagation
still needs an independent stopping/touchdown policy appropriate to the mission;
this diagnostic effector does not terminate an impact.

The existing `gravity_backbone_split` policy rejects this nonconservative,
planet-frame-dependent effector. The plume callback uses the engine's canonical
state accessors, but that does not make the force eligible for that solver mode.
