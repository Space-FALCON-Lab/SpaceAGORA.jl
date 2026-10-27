# Orbital scenes with MuJoCo

!!! note "Optional package"
    Requires the optional `SpaceAGORAMuJoCo` companion; see
    [Install the compatible pair](#Install-the-compatible-pair).

SpaceAGORAMuJoCo combines SpaceAGORA's gravity models with MuJoCo's rigid-body,
joint and contact dynamics. Use it when the local mechanics matter: a
spacecraft with moving links, two vehicles approaching each other, or bodies
that may touch. You can step a scene on its own, or let it propagate selected
spacecraft during a SpaceAGORA run.

The implementation lives in the separate
[SpaceAGORAMuJoCo.jl repository](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl).
Its [user guide at the recorded source revision](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl/blob/7b5a00318f6883d7cdc7c7e7c0c18830f96ce57f/docs/user_guide.md)
has a complete two-spacecraft example, the same example carried into a
SpaceAGORA configuration, and explanations of resetting, forcing and reading
scene states. The examples are kept in the companion, next to the native code
that runs them, so they change together with it.

## Install the compatible pair

Follow the [optional companion installation](installation_environment.md#Optional-MuJoCo-companion)
to create a separate environment, then load `SpaceAGORA` and
`SpaceAGORAMuJoCo` explicitly. The native binding requires Julia 1.12 and Linux
x86_64 and uses the pinned MuJoCo 3.11.0 artifact. Loading SpaceAGORA alone
does not load the native library.

For SpaceAGORA 0.7.1, use the companion source that the core checkout records,
revision `7b5a00318f6883d7cdc7c7e7c0c18830f96ce57f`. Its package version is
still 0.5.0, but the existing v0.5.0 release tag predates compatibility with
this core. From the core checkout, prepare a separate environment:

```text
julia --startup-file=no scripts/optional_packages.jl setup mujoco /path/to/mujoco-env
```

If you have already initialized the companion checkout, keep any local work
and update it to the revision core records before running setup; see
[Set up the companion sources recorded by the updated core](upgrading.md#Set-up-the-companion-sources-recorded-by-the-updated-core).
Core and companion are tested as a pair: a matching version range alone does
not show that two revisions work together. The connection follows the shared
[adapter standard](../reference/adapter_contract.md) and version 2.0.0 of
[ExternalPropagationServices](../reference/external_propagation_services.md).

## Understand the coordinates and ownership

The scene uses a moving origin whose axes stay parallel to the inertial frame.
This keeps the native positions small, near the mechanical scene, while a
separate chief state carries the common orbital motion. It is a translating
frame, not a rotating radial-transverse-normal frame. The body states the
public interface takes and returns are absolute inertial positions and
velocities, in meters and meters per second.

The scene quaternion is scalar-first, `(w, x, y, z)`, and SpaceAGORA's
spacecraft quaternion is scalar-last, `(x, y, z, w)`. Both describe the same
attitude; for SpaceAGORA's quaternion, `v_body = rot(q) * v_inertial`, with
`rot` defined in the
[external propagation reference](../reference/external_propagation_services.md).
The companion's conversion helpers reorder the components; they do not invert
the rotation. Body angular velocity is in body axes, in radians per second.

In a coupled run, each assigned spacecraft has exactly one propagation owner.
MuJoCo advances the bodies the scene owns, including their native joints and
contacts, and the core adds none of its aerodynamic, solar-pressure or control
loads to them. Core refuses a control effector that acts on a scene-owned
spacecraft, so the configuration's control commands never reach MuJoCo
actuators. Keep these limits in mind when you choose a demonstration or
research case.

## Interpret the output

The native scene advances in fixed steps, while the core solver can take
adaptive steps. Between scene steps, the spacecraft state the core sees is a
short extrapolation, and ordinary saved columns can also contain the core
solver's interpolation. The companion's `scene_body_pose_save_field` instead
reads the last completed native pose, with a scalar-last quaternion. It does
not reconstruct the pose at an earlier saved-row time. For a time-matched
comparison, the accepted solver steps must line up with the native step
boundaries; setting the output period alone does not do that. The companion
guide and coupling tests explain how.

The scene method as a whole is first order in its fixed step, because the
gravitational wrenches are held constant during each native step. Repeat a
representative case with smaller steps and compare the positions, attitudes
and contact response that matter to your study. Tightening the core solver
tolerances does not reduce this stepping error. The companion's
[architecture and validation notes](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl/blob/7b5a00318f6883d7cdc7c7e7c0c18830f96ce57f/docs/architecture.md)
describe the measured orbital agreement and the remaining limits of the
attitude comparisons. They support the cases tested there; they are not an
error bound for every mechanical scene.

Core checks a coupled run before it starts. Use a fixed mission duration:
orbit-count termination, atmosphere-interface and entry events, and touchdown
events on owned spacecraft are refused, and core impact detection skips the
owned spacecraft. Coupled runs support the solver modes `:tsit5`,
`:auto_stiff`, `:rodas5p` and `:dp8`. Guidance, navigation and control periods
must be whole multiples of the scene step. An owned spacecraft cannot also use
core articulated joints, compliant attachments or a robot arm. External owners
do not support checkpointing or resume, `run_constellation_ensemble` refuses
them, and the companion guide lists process-pool campaigns as unsupported. The
guide explains these limits next to its examples, including which parallel
routes work.
