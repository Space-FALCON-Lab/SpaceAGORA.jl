# Quaternion helper contracts

The canonical, robotics and compliant-multibody quaternion helpers all use
scalar-last quaternions `[x, y, z, w]`, but they do not share every mathematical
contract. `OrientationQuaternionMath`, in
`src/foundation/core/numerics/orientation_quaternion_utils.jl`, owns the Float64
orientation operations that `Robotics`, `ClothMultibody` and
`ClothRobotArmDynamics` use. (`CompliantMultibody` and `CompliantRobotArmDynamics`
are the public names of the last two modules.) Those modules import the functions
under their existing private names; nothing new is exported.

`test/unit/robotics/quaternion_contract_tests.jl` pins the behavior below. It keeps
orientation operations, which project to unit length, separate from the raw
products that derivatives need.

## Current behavior

| Operation | Canonical `QuaternionMath` | `Robotics` | `ClothMultibody` and `ClothRobotArmDynamics` |
| --- | --- | --- | --- |
| Multiplication | `quat_mult` is raw; input scale is kept. | `_quat_mul` is raw; forward kinematics normalizes the resulting link quaternion. | `_quat_mul` normalizes both operands and the result. `_quat_raw_mul` keeps scale for body-rate derivatives. |
| Unit projection | `project_unit_quaternion` requires a finite squared norm strictly greater than `eps(Float64)`. | `_unit_quat` requires a finite norm strictly greater than `eps(Float64)`. | Same norm cutoff as `Robotics`. |
| Rejected projection input | Identity for zero, nonfinite or below-cutoff input. | Identity. | Identity. |
| Rotation | `rot` is passive: a positive quarter turn about z maps x to −y. It does not normalize, so doubling q scales the matrix by four. | `_rot` projects its input and is active: the same quaternion maps x to +y. | Same as `Robotics`. |
| Zero quaternion | `rot` returns the zero matrix. | `_rot` returns the identity after projection. | Identity after projection. |
| Quaternion derivative | `DynamicsRotational.quaternion_derivative` is raw in both q and the body rate. | Forward kinematics only consumes orientations. | All three RHS paths project the state quaternion, then form a raw product with the body rate. |

The two cutoffs give different answers near zero. For `q = [1e-10, 0, 0, 0]`,
canonical projection returns the identity while the local projectors return
`[1, 0, 0, 0]`. The test covers each cutoff exactly, its neighbouring
floating-point values, and a NaN or infinity in each component. These are the
behaviors callers rely on today, not a recommendation.

## Shared orientation helpers

`OrientationQuaternionMath` holds `_unit_quat`, `_rot`, `_quat_conj`, `_quat_mul`,
`_quat_raw_mul` and `_axis_angle_error`.

- `_unit_quat` converts with `SVector{4, Float64}(q)` before taking `norm`, so
  tuples and other numeric element types are accepted and converted.
- `_quat_raw_mul` accepts only two `SVector{4, Float64}` arguments. `Robotics`
  imports it under its `_quat_mul` name.
- The cloth modules import the normalized `_quat_mul`, which projects the first
  input, then the second, applies the raw product, and projects the result.
- `_axis_angle_error` projects its input, negates it only when the scalar part is
  strictly negative, and returns `2 * vector_part` when `norm(vector_part) <= 1e-12`.

Raw products are never redirected to the canonical `quat_mult`.

## Where each consumer uses them

- Canonical multiplication and rotation are in
  `src/foundation/core/numerics/quaternion_utils.jl`; canonical rotational
  kinematics is in `src/models/dynamics/rotational/attitude_kinematics.jl`.
- `Robotics.cloth_fk`, in `src/models/vehicle/robotics/robotics.jl`, normalizes
  each raw composition before using its active rotation for link positions.
- `ClothMultibody.compliant_joint_loads` uses normalized products for rest and
  error attitudes. `compliant_multibody_dynamics` uses the raw product for the
  quaternion derivative after `compliant_state_parts` projects the state. The
  compiled `compliant_body_derivatives!` path uses the same projector and raw
  product; its mount and in-place load calculations keep the active rotation and
  normalized composition.
- `ClothRobotArmDynamics` uses normalized products for rest and error attitudes.
  `assign_coupled_cloth_robot_arm_rhs!` uses the raw product after
  `_coupled_body_state` projects the state.

On all three cloth RHS paths, a zero body rate gives a zero quaternion
derivative, doubling the rate doubles the derivative, and a body rate that is not
parallel to the rotation axis fixes the multiplication order. Replacing the raw
product with a normalized one would change all three results.

## What the test covers

The test calls the package helpers themselves, forward kinematics and the three
RHS paths. It also checks that the shared functions are the same objects in every
consumer, the input conversions, the raw-product dispatch restriction, shortest
error vectors, and the separate articulated-body helpers. A one-body model
without joints and a stationary one-link arm isolate quaternion behavior without a
trajectory solve, external data, SPICE kernels or native GRAM, and the test checks
that evaluating a derivative leaves the input state unchanged.

`test/probes/shared_math_ownership_probes.jl` covers canonical ownership and frame
conventions. Neither file shows that the orientation helpers could replace the
canonical or articulated ones, and neither validates trajectories.

## Helpers that stay separate

| Owner | Operations and why they stay separate |
| --- | --- |
| `QuaternionMath` and `DynamicsRotational` | Canonical projection uses the squared-norm cutoff; rotation is passive and raw in scale. Canonical multiplication accepts general vectors, and its expression order differs from the typed orientation product. The body-rate derivative keeps its cross/dot form and is raw in the state quaternion. |
| `ArticulatedBody` | `_qmul`, `_qconj`, `_qnormalize`, `_rotmat`, `_axis_quat`, `_axis_rotmat` and `_axis_angle` keep generic element types, including dual numbers and `BigFloat`. Normalization has no identity fallback, rotation divides by the raw squared norm, and the error vector uses a different threshold and formula. |
| `Robotics` | `_quat_from_axis_angle` is the only owner of its norm-checked Float64 axis-angle conversion. Forward kinematics projects the raw composition. |
| `ControlHooks` | Robot-arm conjugation and axis-angle projection use canonical projection and an axis-projected angle, not the shared projector or error vector. |
| `FrameTransforms` | LVLH-to-quaternion construction and `rotate_vector_by_quaternion` keep their frame-specific conversion and vector formula; the vector formula does not project its input. |
| `SceneVisualization` | Short-arc interpolation, sign continuity and velocity-aligned frames work on time histories or frame definitions, not single orientations. |
| `SimulationCallbacks` | The accepted-step projection keeps its tolerance and activity rules; ball-joint projection skips an invalid norm instead of substituting the identity. |

The cloth `_rest_child_parent_quat` functions keep their own names but delegate
conjugation and composition to `OrientationQuaternionMath`. Landing guidance
converts direction-cosine matrices with the canonical `dcm_to_quaternion`; its
wrapper adds no arithmetic of its own.
