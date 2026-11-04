# Quaternion helper contracts

The canonical, robotics and cloth helpers use scalar-last quaternions, but they do not share every mathematical contract. `OrientationQuaternionMath`, in `src/foundation/core/numerics/orientation_quaternion_utils.jl`, now owns the equivalent Float64 operations used by Robotics, ClothMultibody and ClothRobotArmDynamics. The existing private bindings import these functions; the module adds no root exports. The tests in `test/unit/robotics/quaternion_contract_tests.jl` characterize the helpers and their direct consumers while preserving the distinction between orientation operations and raw products used for derivatives.

## Current behavior

| Operation | Canonical QuaternionMath | Robotics | ClothMultibody and ClothRobotArmDynamics |
| --- | --- | --- | --- |
| Multiplication | `quat_mult` is raw; input scale is retained. | `_quat_mul` is raw; forward kinematics explicitly normalizes the resulting link quaternion. | `_quat_mul` normalizes both operands and the result. `_quat_raw_mul` retains scale for body-rate derivatives. |
| Unit projection | `project_unit_quaternion` requires a finite squared norm strictly greater than `eps(Float64)`. | `_unit_quat` requires a finite norm strictly greater than `eps(Float64)`. | Same norm-based cutoff as Robotics. |
| Rejected projection input | Identity quaternion for zero, nonfinite or below-cutoff input. | Identity quaternion. | Identity quaternion. |
| Rotation | `rot` is passive for the tested unit inputs: a positive z quarter-turn maps x to negative y. It does not normalize input; scaling q by two scales the matrix by four. | `_rot` projects input and maps x to positive y for the same quaternion. | Same active direction and projection behavior as Robotics. |
| Zero quaternion rotation | Zero matrix, because `rot` uses the raw quaternion. | Identity matrix after projection. | Identity matrix after projection. |
| Quaternion derivative | `DynamicsRotational.quaternion_derivative` is raw in both q and body angular rate. | Forward kinematics is an orientation consumer. | All three RHS paths project the state quaternion, then use a raw product with the body angular rate. |

The cutoff difference is observable. For `q = [1e-10, 0, 0, 0]`, canonical projection returns identity while the local projectors return `[1, 0, 0, 0]`. The fixture covers each exact cutoff and its adjacent representable values; it also checks nonfinite inputs in each component. These fallback behaviors are existing contracts, not recommendations to change them.

## Shared orientation ownership

The shared owner contains `_unit_quat`, `_rot`, `_quat_conj`, `_quat_mul`, `_quat_raw_mul` and `_axis_angle_error`. Its projector first converts with `SVector{4, Float64}(q)`, then evaluates `norm`; tuples and convertible non-Float64 inputs keep that conversion behavior. The raw product still accepts only two `SVector{4, Float64}` inputs. Robotics imports this raw product under its existing `_quat_mul` binding. The cloth modules import the separately named normalized product, which projects the first input, then the second input, evaluates the unchanged raw expression, then projects the result. Raw products are not redirected to canonical multiplication.

The shared error-vector operation projects the input, negates it only when the scalar part is strictly negative, and retains the `norm(vector_part) <= 1e-12` branch and its `2 * vector_part` result. The identity fallback has the same Float64 components in all consumers. No consumer call order, state projection, derivative, or public export changes with this ownership move.

## Consumer boundaries

- Canonical multiplication and rotation belong to `src/foundation/core/numerics/quaternion_utils.jl`. Canonical rotational kinematics belongs to `src/models/dynamics/rotational/attitude_kinematics.jl`.
- `Robotics.cloth_fk`, in `src/models/vehicle/robotics/robotics.jl`, normalizes a raw composition before using its active rotation for link positions.
- `ClothMultibody.compliant_joint_loads` uses normalized orientation products for rest/error attitudes. `compliant_multibody_dynamics` uses a raw product for the quaternion derivative after `compliant_state_parts` projects the state. The compiled `compliant_body_derivatives!` path uses the same projector and raw product; mount and in-place load calculations retain active rotation and normalized orientation composition.
- `ClothRobotArmDynamics` uses normalized products for rest/error attitudes. `assign_coupled_cloth_robot_arm_rhs!` uses a raw product after `_coupled_body_state` projects the state.

For all three cloth RHS paths, zero angular rate must produce zero quaternion derivative. Doubling the body rate doubles the derivative, and a nonparallel body rate pins multiplication order. Replacing the raw product with a normalized orientation product would change these observable results.

## Fixture scope

The fixture calls the actual package helpers, forward kinematics and all three direct RHS paths. It also checks the shared function identities, preserved input conversions and raw dispatch restrictions, shortest error vectors, and the deliberately separate articulated numeric types. A one-body model without joints and a stationary one-link reference isolate quaternion behavior without a trajectory solve, external data, SPICE kernels or native GRAM execution. It checks that derivative evaluation leaves the input state unchanged.

These tests supplement `test/probes/shared_math_ownership_probes.jl`, which already covers canonical ownership and frame conventions. They do not establish that orientation helpers can be substituted for canonical or generic articulated helpers, or validate trajectories or performance. The following distinct owners remain intentional.


## Operations that retain distinct ownership

| Owner | Retained operations and reason |
| --- | --- |
| `QuaternionMath` and `DynamicsRotational` | Canonical projection uses the squared-norm cutoff; rotation is passive and raw in scale. Canonical multiplication accepts general vectors, and its expression evaluation differs from the typed orientation product. The body-rate derivative uses its existing cross/dot expression and remains raw in the state quaternion. |
| `ArticulatedBody` | `_qmul`, `_qconj`, `_qnormalize`, `_rotmat`, `_axis_quat`, `_axis_rotmat` and `_axis_angle` preserve generic numeric types, including dual numbers and BigFloat. Normalization has no identity fallback; rotation divides by the raw squared norm, and the error-vector branch uses a different threshold and formula. These are not Float64 orientation helpers. |
| `Robotics` | `_quat_from_axis_angle` remains the sole owner of its norm-checked Float64 axis-angle conversion. Forward kinematics explicitly projects the raw composition result. |
| `ControlHooks` | Robot-arm conjugation and axis-angle projection use canonical projection and an axis-projected angle, not the shared norm projector or error vector. |
| `FrameTransforms` | LVLH-to-quaternion construction and `rotate_vector_by_quaternion` keep their frame-specific conversion and vector formula. The vector formula does not project its input. |
| `SceneVisualization` | Short-arc interpolation, sign continuity, and velocity-aligned frame construction operate on time histories or frame definitions. They remain separate from orientation projection. |
| `SimulationCallbacks` | Accepted-step projection keeps its integration tolerance and activity policy; ball-joint projection skips invalid norms rather than substituting identity. |

The cloth `_rest_child_parent_quat` adapters retain their named consumer convention while delegating conjugation and composition to the shared owner. Landing guidance already delegates passive DCM conversion to canonical `dcm_to_quaternion`; its named wrapper is not a second arithmetic implementation.
