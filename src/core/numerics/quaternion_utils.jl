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

@inline function quaternion_norm_error(q::AbstractVector{<:Real})::Float64
    q_unit = SVector{4, Float64}(Float64(q[1]), Float64(q[2]), Float64(q[3]), Float64(q[4]))
    qnorm2 = dot(q_unit, q_unit)
    if !isfinite(qnorm2)
        return Inf
    end
    return abs(qnorm2 - 1.0)
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

## Left — scalar-first convention: Q[1] is scalar, Q[2:4] are vector components.
# Replaces sparse(lIs,lJs,Vs) with a stack-allocated SMatrix (no GC pressure).
@inline function L(Q)
    """
    Left Multiply quaternion
    q2 ⨂ q1 -> L(q2)q1
    params:
    Q -> unit quaternion (scalar-first)
    returns:
    L(Q) -> 4x4 SMatrix
    """
    s = Q[1]; a = Q[2]; b = Q[3]; c = Q[4]
    # Column-major layout: col1=[s,a,b,c], col2=[-a,s,c,-b], col3=[-b,-c,s,a], col4=[-c,b,-a,s]
    return SMatrix{4,4,Float64}(
         s, a, b, c,
        -a, s, c,-b,
        -b,-c, s, a,
        -c, b,-a, s
    )
end

function L!(Lmat, Q)
    Lmat .= L(Q)
end

## Right — scalar-first convention: Q[1] is scalar, Q[2:4] are vector components.
@inline function R(Q)
    """
    Right Multiply quaternion
    q2 ⨂ q1 -> R(q1)q2
    params:
    Q -> unit quaternion (scalar-first)
    returns:
    R(Q) -> 4x4 SMatrix
    """
    s = Q[1]; a = Q[2]; b = Q[3]; c = Q[4]
    # Column-major layout: col1=[s,a,b,c], col2=[-a,s,-c,b], col3=[-b,c,s,-a], col4=[-c,-b,a,s]
    return SMatrix{4,4,Float64}(
         s, a, b, c,
        -a, s,-c, b,
        -b, c, s,-a,
        -c,-b, a, s
    )
end

function R!(Rmat, Q)
    Rmat .= R(Q)
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

## take angle components (4, 3)
H = SMatrix{4,3}([zeros(1, 3); I]);
S = SMatrix{4,3}([[1 0 0]; zeros(3,3)])
## Inverse
T = SMatrix{4,4}(Diagonal([1.0; -1; -1; -1]))

## Attitude Jacobian
function G(Q)
    return L(Q) * H
end

function Q̄(q1, q2)
    """
    Average quaternion
    params:
    q1 -> unit quaternion
    q2 -> unit quaternion
    returns: 
    q̄ -> unit quaternion
    """
    q̂ = L(q1)' * q2
    normv = norm(q̂[1:3])
    if normv == 0.0
        return q1
    else
        r = q̂[1:3] / normv
        θ = (2 * atan(normv, q̂[4]))
    end
    # r = q̂[2:4] / (norm(q̂[2:4]) + 1e-10)
    # θ = 2 * atan(norm(q̂[2:4]), q̂[1])
    q̄ = L(q1) * [cos(θ / 4); sin(θ / 4) * r]
    # values, vectors = eigen(0.5.*([q1 q2] * [q1 q2]'))
    q̄ ./= norm(q̄)

    return q̄

end


@inline function phi_from_q(q)
    """
    axis angle from quaternion (scalar last)
    params:
    q -> unit quaternion
    returns:
    ϕ -> axis (r) * angle (θ)
    """
    v1, v2, v3, s = q[1], q[2], q[3], q[4]
    normv = sqrt(v1*v1 + v2*v2 + v3*v3)
    inv_n = 1.0 / (normv + eps())
    θ = 2 * atan(normv, s)
    return SVector{3, Float64}(v1*inv_n*θ, v2*inv_n*θ, v3*inv_n*θ)
end

@inline function q_from_phi(ϕ)
    """
    quaternion from axis angle (scalar last)
    params:
    ϕ -> axis (r) * angle (θ)
    returns:
    q -> unit quaternion
    """
    θ = norm(ϕ)
    if θ == 0.0
        return SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
    else
        inv_θ = 1.0 / θ
        sθ2 = sin(θ * 0.5)
        return SVector{4, Float64}(ϕ[1]*inv_θ*sθ2, ϕ[2]*inv_θ*sθ2, ϕ[3]*inv_θ*sθ2, cos(θ * 0.5))
    end
end

function dphi_from_q(q)
    v = @views q[1:3]
    # v = SVector(q[2], q[3], q[4])
    s = q[4]
    normv = norm(v)
    return [-2*v/(s^2*(1.0 + normv^2/s^2)) (2 * atan(normv, s))/(normv+1e-10)*Matrix(I, 3,3) .- v*v'*(2 * atan(normv, s))/(normv+1e-10)^3 .+ 2*v*v'/((normv+1e-10)*s*(normv+1e-10)*(1.0 + normv^2/s^2))] 
end

function cayley_map(q) 
    return q[1:3]/q[4]
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


function EulerAnglesToq(θ)
    roll, pitch, yaw = θ[1],θ[2],θ[3]
    cr = cos(roll * 0.5)
    sr = sin(roll * 0.5)
    cp = cos(pitch * 0.5)
    sp = sin(pitch * 0.5)
    cy = cos(yaw * 0.5)
    sy = sin(yaw * 0.5)

    qw = cr * cp * cy + sr * sp * sy
    qx = sr * cp * cy - cr * sp * sy
    qy = cr * sp * cy + sr * cp * sy
    qz = cr * cp * sy - sr * sp * cy

    # check for norm close to 1 (allowing a small tolerance)
    norm = sqrt(qw^2 + qx^2 + qy^2 + qz^2)
    if abs(norm - 1) > 1e-6
        println("Warning: Quaternion is not normalized. Norm is: ", norm)
    end
    return [qx, qy, qz, qw]
end

function rotation_matrix(link)
    q = link.q  # Assuming link.q stores the quaternion
    # Quaternion components
    qx, qy, qz, qw = q[1], q[2], q[3], q[4]
    
    # Rotation matrix from quaternion
    R = [1 - 2*qy^2 - 2*qz^2    2*qx*qy - 2*qz*qw    2*qx*qz + 2*qy*qw;
         2*qx*qy + 2*qz*qw    1 - 2*qx^2 - 2*qz^2    2*qy*qz - 2*qx*qw;
         2*qx*qz - 2*qy*qw    2*qy*qz + 2*qx*qw    1 - 2*qx^2 - 2*qy^2]
    return R
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

@inline function Ψ(q::SVector{4, Float64})
    """
    Attitude kinematics matrix Ψ = [s*I(3) - [v]×; -vᵀ] (scalar-last: q=[v,s]).

    # Arguments
    - `q`: unit quaternion [qx,qy,qz,qw].

    # Returns
    - SMatrix{4, 3, Float64}
    """
    a, b, c, s = q[1], q[2], q[3], q[4]
    # Column-major: col1=[s,-c,b,-a], col2=[c,s,-a,-b], col3=[-b,a,s,-c]
    return SMatrix{4, 3, Float64}(
         s, -c,  b, -a,
         c,  s, -a, -b,
        -b,  a,  s, -c
    )
end

@inline function Ξ(q::SVector{4, Float64})
    """
    Attitude kinematics matrix Ξ = [s*I(3) + [v]×; -vᵀ] (scalar-last: q=[v,s]).

    # Arguments
    - `q`: unit quaternion [qx,qy,qz,qw].

    # Returns
    - SMatrix{4, 3, Float64}
    """
    a, b, c, s = q[1], q[2], q[3], q[4]
    # Column-major: col1=[s,c,-b,-a], col2=[-c,s,a,-b], col3=[b,-a,s,-c]
    return SMatrix{4, 3, Float64}(
         s,  c, -b, -a,
        -c,  s,  a, -b,
         b, -a,  s, -c
    )
end
