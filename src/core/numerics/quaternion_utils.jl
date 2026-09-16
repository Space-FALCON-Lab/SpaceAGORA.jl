using LinearAlgebra
using StaticArrays

const IDENTITY_QUATERNION = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)

### Quaternion Function ###
# Quaternion multiplication (scalar-last convention: q = [qx, qy, qz, qw]).
# Returns SVector to avoid heap allocation.
@inline function quat_mult(q::AbstractVector, p::AbstractVector)
    q0 = q[4]; q1 = q[1]; q2 = q[2]; q3 = q[3]
    p0 = p[4]; p1 = p[1]; p2 = p[2]; p3 = p[3]
    s  = q0*p0 - q1*p1 - q2*p2 - q3*p3
    v1 = q0*p1 + p0*q1 + q2*p3 - q3*p2
    v2 = q0*p2 + p0*q2 + q3*p1 - q1*p3
    v3 = q0*p3 + p0*q3 + q1*p2 - q2*p1
    return SVector{4, Float64}(v1, v2, v3, s)
end

@inline function project_unit_quaternion(q::AbstractVector{<:Real})::SVector{4, Float64}
    q_unit = SVector{4, Float64}(Float64(q[1]), Float64(q[2]), Float64(q[3]), Float64(q[4]))
    qnorm2 = dot(q_unit, q_unit)
    if !(isfinite(qnorm2) && qnorm2 > eps(Float64))
        return IDENTITY_QUATERNION
    end
    return q_unit / sqrt(qnorm2)
end

## skew matrix
@inline function hat(ω::AbstractVector{<:Float64})::SMatrix{3, 3, Float64}
    """
    skew matrix, The matrix equivalent of cross product
    params:
    ω -> vector ∈ R³
    return:
    [ω]x
    """
    w1 = ω[1]; w2 = ω[2]; w3 = ω[3]
    # Column-major: col1=[0,w3,-w2], col2=[-w3,0,w1], col3=[w2,-w1,0]
    return SMatrix{3, 3, Float64}(0, w3, -w2, -w3, 0, w1, w2, -w1, 0)
end

# Rotation matrix from quaternion (scalar-last: q = [qx,qy,qz,qw]).
# Uses column-major SMatrix constructor to avoid a temporary matrix allocation.
@inline function rot(q::AbstractVector{<:Float64})::SMatrix{3, 3, Float64}
    """
    Convert quaternion to rotation matrix
    """
    q1, q2, q3, q4 = q[1], q[2], q[3], q[4]
    q1q1 = q1*q1; q2q2 = q2*q2; q3q3 = q3*q3; q4q4 = q4*q4
    q1q2 = q1*q2; q1q3 = q1*q3; q1q4 = q1*q4
    q2q3 = q2*q3; q2q4 = q2*q4; q3q4 = q3*q4
    # Column-major entries (col1, col2, col3):
    return SMatrix{3, 3, Float64}(
        q1q1-q2q2-q3q3+q4q4, 2*(q1q2-q3q4),          2*(q1q3+q2q4),
        2*(q1q2+q3q4),        -q1q1+q2q2-q3q3+q4q4,   2*(q2q3-q1q4),
        2*(q1q3-q2q4),         2*(q2q3+q1q4),          -q1q1-q2q2+q3q3+q4q4
    )
end

# Builds the attitude-error quaternion (scalar-last convention) without slice allocations.
@inline function error_quaternion(current::SVector{4, Float64}, target::SVector{4, Float64})
    a, b, c, s = current[1], current[2], current[3], current[4]
    # q_matrix = [s*I(3) - hat([a,b,c]), [a,b,c]'; -[a,b,c], s] (column-major)
    q_matrix = SMatrix{4, 4, Float64}(
         s, -c,  b, -a,
         c,  s, -a, -b,
        -b,  a,  s, -c,
         a,  b,  c,  s
    )
    qd = SVector{4, Float64}(-target[1], -target[2], -target[3], target[4])
    result = q_matrix * qd
    return result / norm(result)
end

function qToEulerAngles(q)
    # // this implementation assumes normalized quaternion
    # // converts to Euler angles in 3-2-1 sequence
    qx, qy, qz, qw = normalize(q)
    # roll (x-axis rotation)
    sinr_cosp = 2 * (qw * qx + qy * qz);
    cosr_cosp = 1 - 2 * (qx * qx + qy * qy);
    roll = atan(cosr_cosp,sinr_cosp);# - π/2;

    # pitch (y-axis rotation)
    sinp = sqrt(1 + 2 * (qw * qy - qx * qz));
    cosp = sqrt(1 - 2 * (qw * qy - qx * qz));
    pitch = 2 * atan(cosp,sinp) - π/2;

    # yaw (z-axis rotation)
    siny_cosp = 2 * (qw * qz + qx * qy);
    cosy_cosp = 1 - 2 * (qy * qy + qz * qz);
    yaw = atan(cosy_cosp,siny_cosp);# - π/2;

    return [roll, pitch, yaw]
end

function dcm_to_quaternion(dcm::SMatrix{3, 3, Float64})
    """
    Convert a Direction Cosine Matrix (DCM) to a quaternion.
    
    # Arguments
    - `dcm`: A 3x3 Direction Cosine Matrix (DCM) represented as a StaticVector.
    
    # Returns
    - A quaternion represented as a StaticVector{4, Float64}.
    """
    tr = dcm[1, 1] + dcm[2, 2] + dcm[3, 3]
    i = argmax([dcm[1, 1], dcm[2, 2], dcm[3, 3], tr])
    if i == 4
        q = SVector{4, Float64}([dcm[2, 3] - dcm[3, 2], dcm[3, 1] - dcm[1, 3], dcm[1, 2] - dcm[2, 1], 1+tr])
        return q / norm(q)
    else
        q = zeros(4)
        q[i] = 1+2*dcm[i, i] - tr
        for j in 1:4
            if j != i
                if j != 4
                    q[j] = dcm[i, j] + dcm[j, i]
                else
                    idx = i+1 > 3 ? (i+1)%3 : i+1
                    second_idx = i+2 > 3 ? (i+2)%3 : i+2
                    
                    q[j] = dcm[idx, second_idx] - dcm[second_idx, idx]
                end
            end
        end
        return q / norm(q)
    end
end
