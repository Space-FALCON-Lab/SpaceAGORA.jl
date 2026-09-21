"""
    quaternion_derivative(omega_body, q)

Scalar-last quaternion kinematics for a BODY-frame angular velocity:
`qdot = 0.5 * q ⊗ [omega_body, 0]`, whose vector part is
`qs*omega + qv × omega = qs*omega - omega × qv`.

The pre-2026-07 implementation used `+ omega × qv` — the INERTIAL-rate
composition — while [`angular_acceleration`](@ref) integrates Euler's
equations with body-frame torques and gyroscopic term on the same state.
That frame mixing broke rigid-body invariants: a torque-free asymmetric
tumble (I = diag(1.5, 1.0, 2.0), |omega| ~ 0.04 rad/s) drifted its inertial
angular momentum by 68% over 2000 s; with the body-rate form the same test
conserves to integrator accuracy. The probe suite pins this invariant.
"""
@inline function quaternion_derivative(
    omega_body::SVector{3, Float64},
    q::AbstractVector{<:Real}
)::SVector{4, Float64}
    qv = SVector{3, Float64}(Float64(q[1]), Float64(q[2]), Float64(q[3]))
    qs = Float64(q[4])

    vector_part = qs * omega_body - cross(omega_body, qv)
    scalar_part = -dot(omega_body, qv)

    return 0.5 * SVector{4, Float64}(
        vector_part[1],
        vector_part[2],
        vector_part[3],
        scalar_part,
    )
end
