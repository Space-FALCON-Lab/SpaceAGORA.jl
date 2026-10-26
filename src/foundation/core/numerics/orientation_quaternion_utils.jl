"""Internal Float64 quaternion policy shared by robotics and compliant multibody dynamics.

Scalar-last, active rotations use a norm-based identity fallback. Orientation
products project their inputs and output; the typed raw product retains scale
for forward-kinematics composition and body-rate derivatives. This policy is
intentionally distinct from QuaternionMath and generic articulated dynamics.
"""
module OrientationQuaternionMath

using LinearAlgebra
using StaticArrays

const _Q_IDENTITY = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)

"""Project with the existing norm cutoff after conversion to Float64."""
@inline function _unit_quat(q)::SVector{4, Float64}
    qv = SVector{4, Float64}(q)
    nq = norm(qv)
    return isfinite(nq) && nq > eps(Float64) ? qv / nq : _Q_IDENTITY
end

"""Conjugate the projected orientation, including its identity fallback."""
@inline function _quat_conj(q)::SVector{4, Float64}
    qv = _unit_quat(q)
    return SVector{4, Float64}(-qv[1], -qv[2], -qv[3], qv[4])
end

"""Raw scalar-last product of two Float64 static quaternions; never projects."""
@inline function _quat_raw_mul(a::SVector{4, Float64}, b::SVector{4, Float64})::SVector{4, Float64}
    ax, ay, az, aw = a
    bx, by, bz, bw = b
    return SVector{4, Float64}(
        aw * bx + bw * ax + ay * bz - az * by,
        aw * by + bw * ay + az * bx - ax * bz,
        aw * bz + bw * az + ax * by - ay * bx,
        aw * bw - ax * bx - ay * by - az * bz,
    )
end

"""Project both operands, multiply in order, then project the orientation result."""
@inline function _quat_mul(a, b)::SVector{4, Float64}
    aq = _unit_quat(a)
    bq = _unit_quat(b)
    return _unit_quat(_quat_raw_mul(aq, bq))
end

"""Active body-to-world rotation after norm-based projection."""
@inline function _rot(q_raw)::SMatrix{3, 3, Float64}
    q = _unit_quat(q_raw)
    x, y, z, w = q
    return SMatrix{3, 3, Float64}(
        1 - 2y^2 - 2z^2, 2x * y + 2z * w, 2x * z - 2y * w,
        2x * y - 2z * w, 1 - 2x^2 - 2z^2, 2y * z + 2x * w,
        2x * z + 2y * w, 2y * z - 2x * w, 1 - 2x^2 - 2y^2,
    )
end

"""Shortest orientation error, retaining the existing 1e-12 small-vector branch."""
@inline function _axis_angle_error(q_err_raw)::SVector{3, Float64}
    q_err = _unit_quat(q_err_raw)
    q_err = q_err[4] < 0.0 ? -q_err : q_err
    v = SVector{3, Float64}(q_err[1], q_err[2], q_err[3])
    nv = norm(v)
    if nv <= 1.0e-12
        return 2.0 * v
    end
    θ = 2.0 * atan(nv, q_err[4])
    return (θ / nv) * v
end

end # module OrientationQuaternionMath
