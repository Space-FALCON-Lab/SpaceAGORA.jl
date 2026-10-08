# Quaternion helper contracts

The canonical, robotics and cloth helpers use scalar-last quaternions, but they do not share every mathematical contract. The tests in `test/unit/robotics/quaternion_contract_tests.jl` characterize the existing implementations and direct consumers. They preserve the distinction between orientation operations and raw products used for derivatives.

## Current behavior

| Operation | Canonical QuaternionMath | Robotics | ClothMultibody and ClothRobotArmDynamics |
| --- | --- | --- | --- |
| Multiplication | `quat_mult` is raw; input scale is retained. | `_quat_mul` is raw; forward kinematics explicitly normalizes the resulting link quaternion. | `_quat_mul` normalizes both operands and the result. `_quat_raw_mul` retains scale for body-rate derivatives. |
| Unit projection | `project_unit_quaternion` requires a finite squared norm strictly greater than `eps(Float64)`. | `_unit_quat` requires a finite norm strictly greater than `eps(Float64)`. | Same norm-based cutoff as Robotics. |
| Rejected projection input | Identity quaternion for zero, nonfinite or below-cutoff input. | Identity quaternion. | Identity quaternion. |
| Rotation | `rot` is passive for the tested unit inputs: a positive z quarter-turn maps x to negative y. It does not normalize input; scaling q by two scales the matrix by four. | `_rot` projects input and maps x to positive y for the same quaternion. | Same active direction and projection behavior as Robotics. |
| Zero quaternion rotation | Zero matrix, because `rot` uses the raw quaternion. | Identity matrix after projection. | Identity matrix after projection. |
| Quaternion derivative | `DynamicsRotational.quaternion_derivative` is raw in both q and body angular rate. | Forward kinematics is an orientation consumer. | Both actual RHS callers project the state quaternion, then use a raw product with the body angular rate. |

The cutoff difference is observable. For `q = [1e-10, 0, 0, 0]`, canonical projection returns identity while the local projectors return `[1, 0, 0, 0]`. The fixture covers each exact cutoff and its adjacent representable values; it also checks nonfinite inputs in each component. These fallback behaviors are existing contracts, not recommendations to change them.

## Consumer boundaries

- Canonical multiplication and rotation belong to `src/core/numerics/quaternion_utils.jl`. Canonical rotational kinematics belongs to `src/dynamics/rotational/attitude_kinematics.jl`.
- `Robotics.cloth_fk`, in `src/vehicle/robotics/robotics.jl`, normalizes a raw composition before using its active rotation for link positions.
- `ClothMultibody.compliant_joint_loads` uses normalized orientation products for rest/error attitudes. `compliant_multibody_dynamics` uses a raw product for the quaternion derivative after `compliant_state_parts` projects the state.
- `ClothRobotArmDynamics` uses normalized products for rest/error attitudes. `assign_coupled_cloth_robot_arm_rhs!` uses a raw product after `_coupled_body_state` projects the state.

For both cloth RHS callers, zero angular rate must produce zero quaternion derivative. Doubling the body rate doubles the derivative, and a nonparallel body rate pins multiplication order. Replacing the raw product with a normalized orientation product would change these observable results.

## Fixture scope

The fixture calls the actual package helpers, forward kinematics and both direct RHS paths. A one-body model without joints and a stationary one-link reference isolate quaternion behavior without a trajectory solve, external data, SPICE kernels or native GRAM execution. It checks that derivative evaluation leaves the input state unchanged.

These tests supplement `test/probes/shared_math_ownership_probes.jl`, which already covers canonical ownership and frame conventions. They do not establish that local helpers can be substituted for canonical helpers, validate trajectories, or settle mathematical ownership. A future consolidation must preserve each consumer's normalization, direction and threshold requirements and account for concurrent dynamics changes.
