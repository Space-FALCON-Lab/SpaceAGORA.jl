# Orbital scenes with MuJoCo

SpaceAGORAMuJoCo connects SpaceAGORA's gravitational models to MuJoCo's rigid-body,
joint and contact dynamics. It is useful when the local mechanical scene matters:
a spacecraft with moving links, two vehicles approaching each other, or bodies
that may come into contact. You can advance a scene on its own, or give it
responsibility for selected spacecraft during a SpaceAGORA run.

The implementation lives in the separate
[SpaceAGORAMuJoCo.jl repository](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl).
Its [user guide](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl/blob/988bd957de44799f8aacb0ef87f6bedf4fab16c6/docs/user_guide.md)
contains a complete two-spacecraft example, its extension into a SpaceAGORA
configuration, and explanations of resetting, forcing and reading scene states.
The companion owns these examples and the native implementation, so the
instructions can evolve with the package that executes them.

## Install the compatible pair

Follow the [optional companion installation](installation_environment.md#Optional-MuJoCo-companion)
to create a separate environment, then load `SpaceAGORA` and
`SpaceAGORAMuJoCo` explicitly. The current native binding requires Julia 1.12
and Linux x86_64, using the pinned MuJoCo 3.11.0 artifact. Importing SpaceAGORA
alone does not load the native library.

Use the source revisions recorded by the companion's installer and
[release preparation notes](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl/blob/988bd957de44799f8aacb0ef87f6bedf4fab16c6/docs/releases/v0.1.0.md).
The core and companion are qualified as a pair; a matching version-family
constraint is insufficient to establish compatibility. The connection follows
the shared [adapter standard](../maintainer/adapter_contract.md) and version 1.0.0
of [ExternalPropagationServices](../maintainer/external_propagation_services.md).

## Understand the coordinates and ownership

The scene uses a moving origin whose axes remain parallel to the inertial frame.
This keeps native positions near the mechanical scene while a separate chief
state carries the common orbital motion. It is a translating frame, not a
rotating radial-transverse-normal frame. The public body states are absolute
inertial positions and velocities, in metres and metres per second.

The scene quaternion uses scalar-first order `(w, x, y, z)`, while SpaceAGORA's
spacecraft quaternion uses scalar-last order `(x, y, z, w)`. Both describe the
active rotation from body to inertial axes. The companion's conversion helpers
reorder the components; they do not reverse the rotation. Body angular velocity
is expressed in body axes, in radians per second.

In a coupled run, each assigned spacecraft has one propagation owner. MuJoCo
advances the scene-owned bodies, including their native joints and contacts.
The core does not add its aerodynamic, solar-pressure or control loads to those
bodies. In particular, attaching the scene to a configuration does not transfer
the configuration's control commands into MuJoCo actuators. These restrictions
should guide the choice of a demonstration or research case.

## Interpret the output

The native scene advances in fixed steps, while the core solver can take
adaptive steps. Between scene boundaries, the spacecraft state supplied to the
core uses a short extrapolation. Ordinary saved columns can also contain the
core solver's interpolation. The companion's `scene_body_pose_save_field`
instead reads the current last completed native pose, with a scalar-last
quaternion. It does not reconstruct a pose at an earlier saved-row timestamp.
For a time-matched comparison, accepted solver steps must be aligned with the
requested native boundaries; setting the output period alone is insufficient.
The companion guide and coupling tests explain this timing requirement.

The overall scene method is first order in its fixed step, because gravitational
wrenches are held during each native step. Repeat a representative case with
smaller steps and compare the positions, attitudes and contact response that
matter to the study. Tightening core solver tolerances does not reduce this
native stepping error. The companion's
[architecture and validation notes](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl/blob/main/docs/architecture.md)
describe the measured orbital agreement and the remaining attitude-comparison
limits. They support those tested cases, rather than establishing an error bound
for every mechanical scene.

Use a fixed mission duration for the supported coupling. External owners do not
support core checkpoint/resume or process-pool campaigns. Atmosphere-interface
events are unsupported in coupled runs, while core impact detection skips the
owned entries. The full guide
explains these limits alongside the examples, including which solver modes and
parallel execution routes are supported.
